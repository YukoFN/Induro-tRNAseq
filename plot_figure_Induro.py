import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns
import colorsys
import matplotlib.colors as mc
from math import ceil
import glob
from scipy.stats import pearsonr
from scipy import stats



#####################################################################################

# Induro-tRNAseq vs mim-tRNAseq
directory = ""
sampledata = "sample_data_Induro.txt"
organism = 'Homo_sapiens'
sample = ["Induro_42C_1h","Induro_42C_2h","Induro_42C_16h","Induro_55C_1h","Induro_55C_2h","Induro_55C_16h"]
folder = directory.split("/")[-2]

#####################################################################################



# [0] Mapping rates
def plot_mapping_rate(directory, sampledata):
    # sample data
    df = pd.read_csv('{}{}'.format(directory,sampledata),sep="\t",header=None,names=['file','sample'])
    # print(df.head())
    for i in df[df['file'].str.contains('.fastq.gz')].index:
        df.loc[i,'file_unpaired_uniq'] = df.loc[i,'file'][:-9]+'.unpaired_uniq'
    for i in df[df['file'].str.contains('.fq.gz')].index:
        df.loc[i,'file_unpaired_uniq'] = df.loc[i,'file']+'.unpaired_uniq'

    with open('{}align/mapping_stats.txt'.format(directory)) as f:
        n = False
        for line in f:
            if '*' in line:
                n = True
            elif n == True:
                if '.fq' in line or '.fastq' in line:
                    file = line.rstrip('\n') # Remove "\n"
                elif "Uniquely" in line:
                    unique = line.rstrip("\n").split(": ")[1]
                    uread = int(unique.split(" (")[0])
                elif "Multi-mapping" in line:
                    multi = line.rstrip("\n").split(": ")[1]
                    mread = int(multi.split(" (")[0])
                elif "Unmapped" in line:
                    Unmap = line.rstrip("\n").split(": ")[1]
                    Uread = int(Unmap.split(" (")[0])
                elif 'Total' in line:
                    total = int(line.rstrip("\n").split(": ")[1])
                    if len(df[df['file'] == file]) > 0:
                        d = df[df['file'] == file].index[0]
                        df.loc[d,'Uniquely_mapped_reads'] = uread
                        df.loc[d,'Uniquely_mapped_rates'] = round(uread/total*100,1)
                        df.loc[d,'Multi-mapped_reads'] = mread
                        df.loc[d,'Multi-mapped_rates'] = round(mread/total*100,1)
                        df.loc[d,'Unmapped_reads'] = Uread
                        df.loc[d,'Unmapped_rates'] = round(Uread/total*100,1)
                        df.loc[d,'Total_reads'] = total
                    n = False

    # preparation for plot
    color = ['#d7eaf3','#199FB1','#0D5C75']
    sample = sorted(set(df['sample']))
    x = ['Uniquely mapped reads','Multi-mapped reads','Unmapped reads']

    # plot
    fig = plt.figure(figsize=(5,4*len(sample)))
    sub = 0
    # for cumulative bar graph
    umr = [] # uniquely mapped rates
    mmr = [] # multi-mapped rates
    ur = [] # unmapped rates
    tnr = [] # Total number of reads
    for s in sample:
        dfs = df[df['sample']==s]
        umr.append(sum(dfs["Uniquely_mapped_rates"])/len(dfs))
        mmr.append(sum(dfs["Multi-mapped_rates"])/len(dfs))
        ur.append(sum(dfs["Unmapped_rates"])/len(dfs))
        tnr.append(sum(dfs["Total_reads"])/len(dfs))
        rn = 0 # replicate number
        for ss in dfs.index:
            rn += 1
            df.loc[ss,"replicate_num"] = rn
            y = [df.loc[ss,"{}_reads".format(col)] for col in ['Uniquely_mapped','Multi-mapped','Unmapped']]
            # print("sample:",s,"bam:",df.loc[ss,"file"])
            sub += 1
            ax = fig.add_subplot(len(df),1,sub)
            ax.barh([0,1,2],y,tick_label=x,color=color,edgecolor="black")
            for p,r in zip([0,1,2],y):
                plt.text(r+50000,p-0.15,"{}%".format(round(r/df.loc[ss,"Total_reads"]*100)))
            plt.xlim(0,int(max(df['Total_reads'].astype("int"))))
            plt.xlabel('Number of reads',fontsize=10)
            plt.yticks(fontsize=10)
            plt.title('{}, Total '.format(s)+"{0:,}".format(df.loc[ss,"Total_reads"]),fontsize=10)
            plt.tight_layout()
    plt.savefig('{}/figures/mapping_rate_each.png'.format(directory))

    # plot cumulative bar graph
    fig = plt.figure(figsize=(5,0.6*len(sample)))
    left = [0 for n in ur]
    for y, c, l in zip([ur,mmr,umr],reversed(color),reversed(x)):
        plt.barh(range(len(y)),y,left=left,edgecolor="black",color=c,tick_label=sample,label=l)
        left = [l + n for l, n in zip(left,y)]
        print("left:",left)
        if l == "Uniquely mapped reads":
            for p,r in zip(range(len(y)),y):
                plt.text(r/2,p-0.15,"{}%".format(round(r,1)))
    plt.xlabel("Percentage of reads",fontsize=10)
    plt.legend(loc="lower center",bbox_to_anchor=(.5,1))
    plt.tight_layout()
    plt.savefig('{}/figures/mapping_rate_average.png'.format(directory))
    df["replicate_num"] = df["replicate_num"].astype("int")
    return df.loc[:,['file', 'sample', "replicate_num", 'file_unpaired_uniq']]

# [1] 3'-CCA end analysis
def CCAanalysis(directory,sampledic):
    df = pd.read_csv("{}CCAanalysis/CCAcounts.csv".format(directory),sep="\t")
    df = df[df['gene'].str.contains(organism)].reset_index(drop=True)
    df['gene'] = df["gene"].str.replace("{}_".format(organism),"")
    df['isoacceptor'] = df["gene"].str.replace("mito_tRNA-".format(organism),"mito-")
    df['isoacceptor'] = df["isoacceptor"].str.replace("tRNA-".format(organism),"")
    df['isoacceptor'] = df['isoacceptor'].str.replace("/","").str.replace(r"[0-9]+","").str[:-1]

    print(df.head())
    # plot
    sample = sorted(set(sampledic["sample"]))
    fig = plt.figure(figsize=(2*len(sample),10))
    sub = 0
    for s in sample:
        sub += 1
        print("sample:",s)
        # For isoacceptor
        genei = []
        cci = []
        ci = []
        absenti = []
        ccai = []
        for a in sorted(set(df["isoacceptor"]),reverse=True):
            dfa = df[df['isoacceptor']==a]
            # For isoacceptor
            cc = []
            c = []
            absent = []
            for g in sorted(set(dfa['gene'])):
                dfg = dfa[dfa['gene']==g]
                for f in set(sampledic[sampledic['sample']==s]['file']):
                    #print("file:",f)
                    dfs = dfg[dfg['sample'].str.contains(f.split(".")[0])].reset_index(drop=True)
                    #print(len(dfs),"rows")
                    if len(dfs) == 4:
                        sumg = sum(dfs['count'])
                        absent_ = dfs[dfs['end']=='Absent'].reset_index(drop=True).loc[0,'count']/sumg*100
                        absent.append(absent_)
                        c_ = absent_ + dfs[dfs['end']=='C'].reset_index(drop=True).loc[0,'count']/sumg*100
                        c.append(c_)
                        cc_ = c_ + dfs[dfs['end']=='CC'].reset_index(drop=True).loc[0,'count']/sumg*100
                        cc.append(cc_)
                    else:
                        print("error!")
                        print(dfs)
            # For isoacceptor
            genei.append(a)
            cci.append(sum(cc)/len(cc))
            ci.append(sum(c)/len(c))
            absenti.append(sum(absent)/len(absent))
            allin = [sum(cc)/len(cc),sum(c)/len(c),sum(absent)/len(absent)]
            ccai.append(100-sum(allin))

        y = [absenti,ci,cci,ccai]
        ax = fig.add_subplot(1,len(sample),sub)
        bottom = [0 for g in range(len(genei))]
        print("len(bottom) :",len(bottom))
        for y_, l, c in zip(y,["3'-N","3'-NC","3'-NCC","3'-NCCA"],['#012E54','#0F4C81','#92B0CA','#E9F4FF']):
            ax.barh(range(1,len(genei)+1),y_,left=bottom,label=l,color=c,align="center")
            bottom = [x + y for (x, y) in zip(bottom, y_)]
            print(l,"ave.",round(sum(y_)/len(y_),2),"%")
            if sub == 1:
                ax.set_yticklabels(genei,ha='right')
            else:
                ax.set_yticklabels([])
        ax.set_xlim([0,100])
        ax.set_yticks(range(1,len(genei)+1))
        ax.set_ylim(0,len(bottom)+1)
        ax.set_xlabel('Percentage(%)',fontsize=10)
        ax.set_xticks([0,20,40,60,80,100])
        ax.tick_params(axis='x',labelsize=10)
        ax.tick_params(axis='y',labelsize=10)
        ax.set_title(s,fontsize=10)
    plt.legend(loc="upper right",bbox_to_anchor=(2,1))
    plt.tight_layout()
    plt.savefig("{}figures/CCAanalysis_isoacceptors.png".format(directory))

# [2] ercentage of full length reads for each isodecoder
def full_length(directory,sample):
    out = pd.DataFrame()
    o = -1
    for s in sample:
        print("sample:",s)
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        df["cov"] = df["cov"].astype("float")
        for i in set(df["isodecoder"]):
            dfi = df[df["isodecoder"]==i]
            for r in set(dfi["replicate_num"]):
                o += 1
                dfr = dfi[dfi["replicate_num"]==r].reset_index(drop=True)
                cov = dfr.loc[0,"cov"]
                if np.isnan(cov):
                    print("NaN in",i)
                    break
                if cov != min(dfr["cov"]):
                    print("min(coverage) does not equal coverage at position 1 in",i)
                    print("min cov",cov)
                    print(dfr.head())
                out.loc[o,"condition"] = s
                out.loc[o,"cyto_mito"] = dfr.loc[0,"cyto_mito"]
                out.loc[o,"full_length(%)"] = cov/max(dfr["cov"])*100
    print("done")
    out.to_csv("{}full_length_reads.csv".format(directory),index=False)

    # plot
    sns.set()
    sns.set_style("ticks")
    color = ["#D6FA8C","#A5D721","#82B300","#FDB777","#FD9346","#FF6200"]
    fig = plt.figure(figsize=(4,3*0.6*len(sample)))
    a = 0
    for cm in ["all","cyto","mito"]:
        a += 1
        ax = fig.add_subplot(3, 1, a)
        if cm == "all":
            data = out
        else:
            data = out[out["cyto_mito"]==cm].reset_index(drop=True)
        flierprops = dict(marker='o', markersize=1, alpha=0.5)
        sns.boxplot(x="condition",y="full_length(%)",data=data,ax=ax,palette=color,order=sample,flierprops=flierprops)
        for l in range(len(sample)-1):
            signal = test(data[data["condition"]==sample[l]]["full_length(%)"],data[data["condition"]==sample[l+1]]["full_length(%)"],sample[l],sample[l+1])
            if signal != "":
                plt.text(l+0.5,102,signal,ha="center",fontsize=20)
                plt.plot([l+0.1,l+0.9],[107,107],c="black")
        ax.set_title(f"{cm} tRNAs")
        plt.ylim(-10,120)
        plt.xticks(rotation=45)
        plt.tight_layout()
    plt.savefig("{}figures/Full-length_reads.png".format(directory))

# [3] Percentage of misincorporation proportion for each isodecoder
def ave_mismatch(directory,sample):
    out = pd.DataFrame()
    o = -1
    for s in sample:
        print("sample:",s)
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        df = df[df["prop_type"]=="mismatch"].reset_index(drop=True)
        print(df.head())
        for i in set(df["isodecoder"]):
            dfi = df[df["isodecoder"]==i]
            for r in set(dfi["replicate_num"]):
                o += 1
                dfr = dfi[dfi["replicate_num"]==r].reset_index(drop=True)
                out.loc[o,"cyto_mito"] = dfr.loc[0,"cyto_mito"]
                out.loc[o,"condition"] = s
                out.loc[o,"ave_mismatch"] = sum(dfr["proportion"])/len(dfr)
    print("done")
    print(out.head())
    out.to_csv("{}ave_misincorporation.csv".format(directory),index=False)

    # plot
    sns.set()
    sns.set_style("ticks")
    color = ["#D6FA8C","#A5D721","#82B300","#FDB777","#FD9346","#FF6200"]
    fig = plt.figure(figsize=(4,3*0.6*len(sample)))
    a = 0
    for cm in ["all","cyto","mito"]:
        a += 1
        ax = fig.add_subplot(3, 1, a)
        if cm == "all":
            data = out
        else:
            data = out[out["cyto_mito"]==cm].reset_index(drop=True)
        flierprops = dict(marker='o', markersize=1, alpha=0.5)
        sns.boxplot(x="condition",y="ave_mismatch",data=data,ax=ax,palette=color,order=sample,flierprops=flierprops)
        for l in range(len(sample)-1):
            signal = test(data[data["condition"]==sample[l]]["ave_mismatch"],data[data["condition"]==sample[l+1]]["ave_mismatch"],sample[l],sample[l+1])
            if signal != "":
                plt.text(l+0.5,round(max(data["ave_mismatch"]),0)+1,signal,ha="center",fontsize=20)
                plt.plot([l+0.1,l+0.9],[round(max(data["ave_mismatch"]),0)+1,round(max(data["ave_mismatch"]),0)+1],c="black")
        ax.set_title(f"{cm} tRNAs")
        plt.ylim(-0.01,0.08)
        plt.xticks(rotation=45)
        plt.tight_layout()
    plt.savefig("{}figures/ave_misincorporation.png".format(directory))

# [a] Summary table of RTstop, misincorporation, & readthrough proportions at individual position
def RTstop_mismatch(directory,sampledic,organism):
    sample = sorted(set(sampledic['sample']))
    print("sample:",sample)
    folder = directory.split("/")[-2].replace("3rd_","")

    ind = -1
    for s in sample:
        print("sample:",s)
        out = pd.DataFrame()
        for file in ['mismatch','RTstop','readthrough']:
            print("file:",file)
            df = pd.read_csv('{}mods/{}Table.csv'.format(directory,file),sep="\t")
            df = df[df['isodecoder'].str.contains(organism)].reset_index(drop=True)
            df['isodecoder'] = df['isodecoder'].str.replace("{}_".format(organism),"")
            df["isoacceptor"] = df['isodecoder'].str.replace("/","").str.replace(r"[0-9]+","").str[:-1]
            print(df.head())

            for rep in sampledic[sampledic["sample"]==s].index:
                repnum = sampledic.loc[rep,"replicate_num"]
                bam = sampledic.loc[rep,"file"]
                if "DN" in s:
                    bam = bam.split(".")[0]
                dfr = df[df["bam"] == f"{folder}/{bam}.unpaired_uniq.bam"]
                print("replicate:",repnum,bam)
                print("length:",len(dfr))
                print(dfr.head())
                for i in sorted(set(df["isodecoder"])):
                    # print("isodecoder:",i)
                    dfi = dfr[dfr["isodecoder"]==i]
                    # print(round(sum(dfi["proportion"]),0))
                    # print(len(dfi),"rows")
                    for p in sorted(set(dfi["pos"])):
                        # print("pos:",p)
                        ind += 1
                        dfp = dfi[dfi["pos"]==p].reset_index(drop=True)
                        # print(dfp)
                        # print("len:",len(dfp))
                        if file == "mismatch":
                            if len(dfp) != 4:
                                print("Warning: >4 bases at one position.")
                                # print(dfp.head())
                        elif file == "RTstop":
                            if len(dfp) != 1:
                                print("Warning: >1 bases at one position.")
                        out.loc[ind,"isodecoder"] = dfp.loc[0,"isodecoder"]
                        out.loc[ind,"isoacceptor"] = dfp.loc[0,"isoacceptor"]
                        if "mito" in dfp.loc[0,"isodecoder"]:
                            out.loc[ind,"cyto_mito"] = "mito"
                        elif "mito" not in dfp.loc[0,"isodecoder"]:
                            out.loc[ind,"cyto_mito"] = "cyto"
                        out.loc[ind,"pos"] = dfp.loc[0,"pos"]
                        out.loc[ind,"canon_pos"] = dfp.loc[0,"canon_pos"]
                        out.loc[ind,"cov"] = dfp.loc[0,"cov"]
                        out.loc[ind,"condition"] = dfp.loc[0,"condition"]
                        out.loc[ind,"replicate_num"] = repnum
                        out.loc[ind,"prop_type"] = file
                        if file == "readthrough":
                            if p == 0 or p == 1:
                                out.loc[ind,"proportion"] = 1
                            else:
                                rtin = dfi[dfi["pos"]==p-1].index
                                out.loc[ind,"proportion"] = float(dfi[dfi["pos"]==p-1]["cov"]/dfp.loc[0,"cov"])
                            ind += 1
                            if "mito" in dfp.loc[0,"isodecoder"]:
                                out.loc[ind,"cyto_mito"] = "mito"
                            elif "mito" not in dfp.loc[0,"isodecoder"]:
                                out.loc[ind,"cyto_mito"] = "cyto"
                            out.loc[ind,"isodecoder"] = dfp.loc[0,"isodecoder"]
                            out.loc[ind,"isoacceptor"] = dfp.loc[0,"isoacceptor"]
                            out.loc[ind,"pos"] = dfp.loc[0,"pos"]
                            out.loc[ind,"canon_pos"] = dfp.loc[0,"canon_pos"]
                            out.loc[ind,"cov"] = dfp.loc[0,"cov"]
                            out.loc[ind,"condition"] = dfp.loc[0,"condition"]
                            out.loc[ind,"replicate_num"] = repnum
                            out.loc[ind,"prop_type"] = "RTstop_eachpos"
                            out.loc[ind,"proportion"] = 1-out.loc[ind-1,"proportion"]
                        else:
                            out.loc[ind,"proportion"] = sum(dfp["proportion"])

        print('done')
        out.to_csv('{}RTstop_mismatch_summarized_{}.csv'.format(directory,s),index=False)

# [b] Add info of annotated tRNA modificartions in Nedialkova's paper
def annotation_mod(directory,sample):
    # mod info
    dfm = pd.read_csv("modContext.csv")
    dfm["canon_pos"] = dfm["canon_pos"].astype(str).str.replace("20.3","20a").str.replace("20.6","20b")#.str.split(".")[0]
    dfp = dfm["canon_pos"].astype(str).str.split(".",expand=True).rename(columns={0:"canon_pos"})
    dfm["canon_pos"] = dfp["canon_pos"]
    modifications = ["I","Y","acp3U","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
    print("modifications:",set(dfm["Annotated modification"]))
    print(set(dfm["canon_pos"]))
    print(dfm.head())
    for s in sample:
        print("sample;",s)
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        df = df[(df["prop_type"]=="mismatch") | (df["prop_type"]=="RTstop_eachpos")].reset_index(drop=True)
        for m in dfm.index:
            # print(dfm.loc[m,:])
            cpos = dfm.loc[m,"canon_pos"]
            iso = dfm.loc[m,"isodecoder"]
            mod = dfm.loc[m,"Annotated modification"]
            df_ = df[(df["canon_pos"]==cpos) & (df["isodecoder"]==iso)]
            for din in df_.index:
                # print(df.loc[din,:])
                df.loc[din,"Annotated modification"] = mod
        df = df.fillna(0)
        print(df.head())
        df.to_csv("{}RTstop_mismatch_summarized_{}_annotated_modification.csv".format(directory,s),index=False)

# [4] RTstop & misincorporation pattern at modification sites
def detectable_modifications(directory,sample,cal):
    for cm in ["cyto","mito","all"]:
        print(cm,"tRNAs")
        modifications = ["acp3U","I","Y","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
        # plot 1
        sns.set()
        sns.set_style("ticks")
        fig = plt.figure(figsize=(2.5*len(sample),4*len(modifications)))
        color = ["#3CA2C8","#DB4C77"]
        a = 0
        out = pd.DataFrame()
        o = -1
        for m in modifications:
            print("modification:",m)
            for s in sample:
                a += 1 # for plot
                print("sample:",s)
                df = pd.read_csv("{}RTstop_mismatch_summarized_{}_annotated_modification.csv".format(directory,s))
                df["mod_pos"] = ["a" for n in df.index]
                dfm = df[df["Annotated modification"]==m].loc[:,["isodecoder","pos","canon_pos","Annotated modification","cyto_mito"]]
                if cm != "all":
                    dfm = dfm[dfm["cyto_mito"]==cm]
                modlist_ = []
                if len(dfm) > 0:
                    for i in dfm.index:
                        mod = dfm.loc[i,"Annotated modification"]
                        iso = dfm.loc[i,"isodecoder"]
                        pos = dfm.loc[i,"pos"]
                        cpos = dfm.loc[i,"canon_pos"]
                        modlist = list(set(df.loc[i-2:i+2,"Annotated modification"]))
                        modlist.remove("0")
                        modlist.remove(mod)
                        modlist_.extend(modlist)
                        for n in range(-1,2):
                            df.loc[i+n,"mod_pos"] = n
                        modlist1 = list(set(df.loc[i-1:i+1,"Annotated modification"]))
                        modlist1.remove("0")
                        modlist1.remove(mod)
                        o += 1
                        pt = df.loc[i,"prop_type"]
                        # print(pt)
                        out.loc[o,"isodecoder"] = iso
                        out.loc[o,"canon_pos"] = cpos
                        out.loc[o,"Annotated modification"] = mod
                        out.loc[o,"condition"] = df.loc[i,"condition"]
                        out.loc[o,"prop_type"] = pt
                        if pt == "mismatch":
                            # print(df.loc[i,:])
                            out.loc[o,"proportion"] = df.loc[i,"proportion"]
                        else:
                            if cal == "ave":
                                out.loc[o,"proportion"] = sum(df.loc[i-1:i+1,"proportion"])/3
                            elif cal == "max":
                                out.loc[o,"proportion"] = max(df.loc[i-1:i+1,"proportion"])
                        if len(modlist1) > 0:
                            newmod = set(modlist1) & set(modifications)
                            if len(newmod) > 0:
                                out.loc[o,"near_mod"] = modlist1[0]
                            else:
                                out.loc[o,"near_mod"] = "0"
                        else:
                            out.loc[o,"near_mod"] = "0"
                    # plot 1
                    df = df.fillna("a")
                    data = df[(df["mod_pos"] != "a")]
                    ax = fig.add_subplot(len(modifications),len(sample), a)
                    if len(data) > 0:
                        data1 = data[data["prop_type"] == "RTstop_eachpos"]
                        flierprops = dict(marker='o', markersize=1, alpha=0.5)
                        sns.boxplot(x="mod_pos",y="proportion",data=data1,ax=ax,color=color[1],flierprops=flierprops)
                        plt.plot([-0.5,2.5],[0.1,0.1],linestyle="dashed",color="gray",linewidth=0.5)
                        data2 = data[data["prop_type"] == "mismatch"]
                        data2["proportion"] = -data2["proportion"]
                        sns.boxplot(x="mod_pos",y="proportion",data=data2,ax=ax,color=color[0],flierprops=flierprops)
                        plt.plot([-0.5,2.5],[-0.1,-0.1],linestyle="dashed",color="gray",linewidth=0.5)
                    ax.set_title("{}x{} with {}".format(int(len(data)/10),m,set(modlist_)))
                    plt.ylim(-1.1,1.4)
                    plt.legend(["RTstop","Misincorporation"],loc="lower center",bbox_to_anchor=(0.5,1.2))
        plt.tight_layout()
        plt.savefig("{}figures/Detectable_modifications_{}-tRNAs.png".format(directory,cm))
    out.to_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,cal),index=False)

    # plot 2
    modifications = ["acp3U","I","Y","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
    out = pd.read_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,cal))
    print(out.head())
    for near in ["all","near"]:
        fig2 = plt.figure(figsize=(8,3*len(modifications)))
        color2 = ["#D6FA8C","#A5D721","#82B300","#FDB777","#FD9346","#FF6200"]
        b = 0
        for m in modifications:
            print("modification:",m)
            data2 = out[out["Annotated modification"] == m]
            # print(data2.head())
            print(len(data2)/len(set(data2["condition"]))/len(set(data2["isodecoder"])),"modifications")
            print(set(data2["canon_pos"]))
            if near == "near":
                data2 = data2[data2["near_mod"]==0]
                print("near")
                print(len(data2),"modifications")
            if len(data2) > 0:
                for pt in ["mismatch","RTstop_eachpos"]:
                    print("prop_type:",pt)
                    data3 = data2[data2['prop_type']==pt]
                    if len(data3) > 0:
                        b += 1
                        ax2 = fig2.add_subplot(len(modifications),2, b)
                        flierprops = dict(marker='o', markersize=1, alpha=0.5)
                        sns.boxplot(x="condition",y="proportion",data=data3,ax=ax2,palette=color2,order=sample, flierprops=flierprops)
                        for l in range(len(sample)-1):
                            signal = test(data3[data3["condition"]==sample[l]]["proportion"],data3[data3["condition"]==sample[l+1]]["proportion"],sample[l],sample[l+1])
                            if signal != "":
                                plt.text(l+0.5,round(max(data3["proportion"]),1)+0.1,signal,ha="center",fontsize=12)
                                plt.plot([l+0.1,l+0.9],[round(max(data3["proportion"]),1)+0.1,round(max(data3["proportion"]),1)+0.1],c="black")
                        ax2.set_title("{} at {}".format(pt,m))
                        plt.ylim(-0.1,1.4)
        plt.tight_layout()
        plt.savefig("{}figures/Detectable_modifications_all_conditions_{}_mods_{}.png".format(directory,near,cal))






# Statistic analysis
def test(a,b,a_,b_):
    # print(a_,"vs",b_)
    A_var = np.var(a, ddof=1)
    B_var = np.var(b, ddof=1)
    A_df = len(a) - 1
    B_df = len(b) - 1
    if B_var == 0:
        return "error"
    else:
        f = A_var / B_var
        one_sided_pval1 = stats.f.cdf(f, A_df, B_df)
        one_sided_pval2 = stats.f.sf(f, A_df, B_df)
        two_sided_pval = min(one_sided_pval1, one_sided_pval2) * 2
        signal = ""
        if round(two_sided_pval, 3) < 0.05: # Welch's t-test
            t, pval = stats.ttest_ind(a,b,equal_var=False)
            if pval < 0.01:
                signal = "**"
                print("Welch's t-test p-value: ",pval,a_,"vs",b_)
            elif pval < 0.05:
                signal = "*"
                print("Welch's t-test p-value: ",pval,a_,"vs",b_)
        else: # Student's t-test
            t, pval = stats.ttest_ind(a,b,equal_var=True)
            if pval < 0.01:
                signal = "**"
                print("Student's t-test p-value: ",pval,a_,"vs",b_)
            elif pval < 0.05:
                signal = "*"
                print("Student's t-test p-value: ",pval,a_,"vs",b_)
        return signal






# [0] Mapping rates
sampledic = plot_mapping_rate(directory, sampledata)
# print("sampledic:",sampledic)
print("sample:",sample)

# [1] 3'-CCA end analysis
CCAanalysis(directory,sampledic)

# [2] Percentage of full length reads for each isodecoder
full_length(directory,sample)

# [3] Percentage of misincorporation proportion for each isodecoder
ave_mismatch(directory,sample)    

# [a] Summary table of RTstop, misincorporation, & readthrough proportions at individual position
RTstop_mismatch(directory,sampledic,organism)

# [b] Add info of annotated tRNA modificartions in Nedialkova's paper
annotation_mod(directory,sample)

# [4] RTstop & misincorporation pattern at modification sites
detectable_modifications(directory,sample,"max")

