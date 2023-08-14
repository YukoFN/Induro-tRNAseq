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
directory = '/Volumes/YMH_tRNAseq/3rd_tRNA-seq/3rd_HEK293T_vs_K562_020223/'
sampledata = 'sample_data_HEK293T-K562.txt'
organism = 'Homo_sapiens'
sample2 = ["K562_totalRNA","K562_tRNA_DN"]

#####################################################################################



# [0] Mapping rates
def plot_mapping_rate(directory, sampledata):
    # sample data
    df = pd.read_csv('{}{}'.format(directory,sampledata),sep="\t",header=None,names=['file','sample'])
    # print(df.head())
    # df['file_unpaired_uniq'] = [0 for i in df.index]
    for i in df[df['file'].str.contains('.fastq.gz')].index:
        df.loc[i,'file_unpaired_uniq'] = df.loc[i,'file'][:-9]+'.unpaired_uniq'
    for i in df[df['file'].str.contains('.fq.gz')].index:
        df.loc[i,'file_unpaired_uniq'] = df.loc[i,'file']+'.unpaired_uniq'
    # print(df[df['file_unpaired_uniq']==0])

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
    # for g in ["tRNA-",r"[0-9]+",r"-$",r"-/$",r"-//$"]: #$をつけると、末尾で一致するものになる。
    #     df['isoacceptor'] = df['isoacceptor'].str.replace(g,"")

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
                        #cca_ = 100#dfs[dfs['end']=='CA'].reset_index(drop=True).loc[0,'count']/sumg*100
                        #ccad.append(cca_)
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
                # ax.get_yaxis().set_tick_params(pad=45)
            else:
                ax.set_yticklabels([])
        ax.set_xlim([0,100])
        ax.set_yticks(range(1,len(genei)+1))
        ax.set_ylim(0,len(bottom)+1)
        #ax.set_ylabel('All human tRNA isodecoders',fontsize=10)
        ax.set_xlabel('Percentage(%)',fontsize=10)
        ax.set_xticks([0,20,40,60,80,100])
        ax.tick_params(axis='x',labelsize=10)
        ax.tick_params(axis='y',labelsize=10)
        ax.set_title(s,fontsize=10)
    plt.legend(loc="upper right",bbox_to_anchor=(2,1))
    plt.tight_layout()
    plt.savefig("{}figures/CCAanalysis_isoacceptors.png".format(directory))

# [2] Differential abundance analysis
def differential_expression(directory,sampledic,log10,organism):
    # tRNA: "tRNA-Xxx-XXX-X", or "tRNA-Xxx-XXX", or "mito_tRNA-Xxx-XXX"

    for file in ['cyto_Isodecoder','cyto_Anticodon','mito_Isodecoder']:
        print("file:",file)
        if file == 'cyto_Anticodon':path = "{}DESeq2/cyto/anticodon/".format(directory)
        elif file == 'cyto_Isodecoder':path = "{}DESeq2/cyto/isodecoder/".format(directory)
        elif file == 'mito_Isodecoder':path = "{}DESeq2/organelle/isodecoder/".format(directory)
        files = glob.glob("{}*results.csv".format(path))

        fig = plt.figure(figsize=(2.9,3.1*(len(files)-1)))
            
        n = 0

        for f in files:
            f_ = f.split(path)[1]
            ff_ = [f_.split('vs')[0],f_.split("vs")[1].split('_diff')[0]]
            ff = sorted(ff_)
            print('sample combination',ff)
            f1 = ff[0]
            f2 = ff[1]
            col1 = list(sampledic[sampledic['sample'] == f1]['file_unpaired_uniq'])
            col2 = list(sampledic[sampledic['sample'] == f2]['file_unpaired_uniq'])
            df = pd.read_csv(f,usecols=['Gene','log2FoldChange','padj']+col1+col2)
            # print(df.head())
            print("len(file):",len(df))
            df = df[df['Gene'].str.contains(organism)].reset_index(drop=True)
            df['{}_ave'.format(f1)] = df.loc[:,col1].sum(axis=1)/len(col1)
            df['{}_ave_log10'.format(f1)] = np.log10(df['{}_ave'.format(f1)])
            df['{}_ave'.format(f2)] = df.loc[:,col2].sum(axis=1)/len(col2)
            df['{}_ave_log10'.format(f2)] = np.log10(df['{}_ave'.format(f2)])
            df.replace([np.inf,-np.inf],0,inplace=True)
            # print(df.head())
            # print(df.columns)
            if file == "mito_Isodecoder":
                print("Ser in {}".format(f1),df[df["Gene"].str.contains("Ser")]['{}_ave'.format(f1)])
                print("Ser in {}".format(f2),df[df["Gene"].str.contains("Ser")]['{}_ave'.format(f2)])

            #plot
            n += 1
            ax = fig.add_subplot(len(files),1,n) #all tRNA
                
            if log10 == 'log10':
                xcol = '{}_ave_log10'.format(f1)
                ycol = '{}_ave_log10'.format(f2)
            else:
                xcol = '{}_ave'.format(f1)
                ycol = '{}_ave'.format(f2)
            x = df[xcol]
            y = df[ycol]
            # threshold
            maxval = ceil(max([max(x),max(y)])/1)*1
            ax.plot([0,6],[0,6],linewidth=0.5,c='gray')
            r, p_value = pearsonr(df['{}_ave'.format(f1)],df['{}_ave'.format(f2)])
            ax.text(maxval-2,0.5,"R = {}".format(round(r,2)))

            dfo = df
            alpha = 0.3
            xx = dfo[xcol]
            yy = dfo[ycol]
            ax.scatter(xx,yy,c='darkgray',alpha=alpha,edgecolor="darkgray",label="all")#c='lightgray'                        
            ax.set_ylim(0,6)#maxval)
            ax.set_xlim(0,6)#maxval)
            ax.tick_params(axis='x',labelsize=10)
            ax.tick_params(axis='y',labelsize=10)
            ax.set_xlabel(xcol,fontsize=10)
            ax.set_ylabel(ycol,fontsize=10)
            ax.set_title("{}, All {} tRNAs".format(file,len(df)),fontsize=10)
            plt.subplots_adjust(hspace=0.5,wspace=0.5)
            plt.tight_layout()
        plt.savefig('{}figures/differential_expression_{}.png'.format(directory,file))

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

# [b] Table for heatmap of misincorporation and RTStop
def RTstop_mismatch_for_heatmap(directory,sample,isotype):
    position = [f"{n}" for n in range(21)]+["20a","20b"]+[f"{n}" for n in range(21,46)]+[f"e{n}" for n in range(1,20)]+[f"{n}" for n in range(46,76)]
    position.remove("e7")
    position.remove("e13")
    for s in sample:
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        print("sample:",s)
        for proptype in ["mismatch","RTstop_eachpos"]:
            print("prop_type:",proptype)
            dft = df[df["prop_type"]==proptype].reset_index(drop=True)
            print(dft.head())
            # print(dfp.head())
            out = pd.DataFrame(columns=position, dtype=object)
            # print(out.head())
            if isotype == "isoacceptor":
                iso_ = sorted(set(df["isoacceptor"]))
            elif isotype == "isodecoder":
                iso_ = sorted(set(df["isodecoder"]))
            for iso in iso_:
                if isotype == "isoacceptor":
                    dfi = dft[dft["isoacceptor"]==iso].reset_index(drop=True)
                elif isotype == "isodecoder":
                    dfi = dft[dft["isodecoder"]==iso].reset_index(drop=True)
                # print(iso)
                # print(dfi.head())
                if len(dfi) > 0:
                    for pos in position:
                        dfp = dfi[dfi["canon_pos"]==pos]
                        if len(dfp) > 0:
                            out.loc[iso,pos] = sum(dfp["proportion"])/len(dfp)
                        else:
                            out.loc[iso,pos] = -1
            out.to_csv("{}heatmap_{}_{}_{}.csv".format(directory,proptype,isotype,s))

# [3] heatmap of misincorporation and RTStop
def heatmap_RTsop_mismatch(directory,sample):
    for s in sample:
        print("sample:",s)
        fig = plt.figure(figsize=(20,25))
        n = 0
        for proptype in ["mismatch","RTstop_eachpos"]:
            out = pd.read_csv("{}heatmap_{}_{}.csv".format(directory,proptype,s),index_col=0)
            out_mask = (out == -1)
            iso_ = list(out.index)
            n += 1
            ax = fig.add_subplot(2,1,n)
            if proptype == "mismatch":
                cmap = "YlGnBu"
            else:
                cmap = "YlOrRd"
            sns.heatmap(out,square=True,ax=ax,cmap='Greys',cbar=False)
            sns.heatmap(out,square=True,ax=ax,cmap=cmap,vmin=0,vmax=1,mask=out_mask,linewidths=0.2)#,linecolor="gray")
            ax.set_xlabel("Position")
            ax.set_ylabel("tRNA isoacceptors")
            ax.set_title("{} in {}".format(proptype,s))
            ax.set_xticks([i+0.5 for i in range(len(out.columns))])
            ax.set_xticklabels(list(out.columns))
            ax.set_yticks([i+0.5 for i in range(len(iso_))])
            ax.set_yticklabels(iso_)
            plt.xticks(fontsize=10,rotation=90)
            plt.yticks(fontsize=10)
        plt.tight_layout()
        fig.savefig("{}figures/heatmap_{}".format(directory,s))

# [c] Add info of annotated tRNA modificartions in Nedialkova's paper
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

# [d] Summary table of RTstop & misincorporation proportions at modification sites
def detectable_modifications(directory,sample,cal):
    for cm in ["cyto","mito","all"]:
        print(cm,"tRNAs")
        modifications = ["acp3U","I","Y","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
        out = pd.DataFrame()
        o = -1
        for m in modifications:
            print("modification:",m)
            for s in sample:
                print("sample:",s)
                df = pd.read_csv("{}RTstop_mismatch_summarized_{}_annotated_modification.csv".format(directory,s))
                df["mod_pos"] = ["a" for n in df.index]
                # print(df.head())
                # print(set(df["Annotated modification"]))
                dfm = df[df["Annotated modification"]==m].loc[:,["isodecoder","pos","canon_pos","Annotated modification","cyto_mito"]]#.drop_duplicates()
                if cm != "all":
                    dfm = dfm[dfm["cyto_mito"]==cm]
                # print(dfm.head())
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
                        # print(mod,iso,pos)
                        # for "out"
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

    out.to_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,cal),index=False)

# [4] Comparison of readthrough proportion at annotated modification sites between two samples
def compare_readthrough(directory,sample):
    am = "max"
    # plot 1
    modifications = ["m1A","m1G","m2,2G","m3C","I","m1I","Y","acp3U","ms2i6A/ms2t6A"]
    out = pd.read_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,am))
    out = out[(out["condition"]==sample[0]) | (out["condition"]==sample[1])]
    print(out.head())
    near = "all"
    fig2 = plt.figure(figsize=(6,10))
    color = ["#B564E3","#D4D4D4"]
    b = 0
    for pt in ["mismatch","RTstop_eachpos","readthrough"]:
        print("prop_type:",pt)
        if pt == "readthrough":
            data3 = out[out['prop_type']=="RTstop_eachpos"].reset_index(drop=True)
            data3[pt] = 3-data3["proportion"]
            y ="readthrough"
        else:
            data3 = out[out['prop_type']==pt].reset_index(drop=True)
            y = "proportion"
        if len(data3) > 0:
            b += 1
            ax2 = fig2.add_subplot(3,1, b)
            flierprops = dict(marker='o', markersize=1, alpha=0.5)
            sns.boxplot(x="Annotated modification",y=y,data=data3,ax=ax2,palette=color,hue="condition",hue_order=sample,order=modifications,linewidth=1, flierprops=flierprops)
            signal = test(data3[data3["condition"]==sample[0]][y],data3[data3["condition"]==sample[1]][y],sample[0],sample[1])
            if signal != "":
                plt.text(l+0.5,round(max(data3["proportion"]),1)+0.1,signal,ha="center",fontsize=12)
                plt.plot([l+0.1,l+0.9],[round(max(data3[y]),1)+0.1,round(max(data3[y]),1)+0.1],c="black")
            ax2.set_title(pt,fontsize=10)
            plt.legend(loc="lower center",bbox_to_anchor=(0.5,1.2))
            if pt == "mismatch":
                plt.ylim(-0.1,1.1)
            else:
                if max(data3[y]) < 1:
                    plt.ylim(-0.1,1.1)
                else:
                    plt.ylim(-0.1,round(max(data3[y]),1)+0.5)
    plt.tight_layout()
    plt.savefig("{}figures/Detectable_modifications_all_conditions_{}_mods_readthrough_{}.png".format(directory,near,am))

    # plot 2
    data = pd.DataFrame()
    d  = -1
    for m in modifications:
        print("modification",m)
        for s in sample:
            d += 1
            data3 = out[(out["proportion"] > 0.1) & (out["condition"]==s) & (out["Annotated modification"]==m)].loc[:,["isodecoder","canon_pos","Annotated modification"]]
            # print(data3.head())
            print((~data3.duplicated()).sum(),"modifications")
            data.loc[d,"condition"] = s
            data.loc[d,"Annotated modification"] = m
            data.loc[d,"num_detected_mod"] = (~data3.duplicated()).sum()
    # print(data.head())
    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot(111)
    sns.barplot(x="Annotated modification",y="num_detected_mod",data=data,ax=ax,palette=color,hue="condition",order=modifications,linewidth=1,edgecolor="black")
    for m, x in zip(modifications,range(0,len(modifications))):
        y1 = data[(data["Annotated modification"]==m) & (data["condition"]==sample[0])].reset_index(drop=True).loc[0,"num_detected_mod"]
        y2 = data[(data["Annotated modification"]==m) & (data["condition"]==sample[1])].reset_index(drop=True).loc[0,"num_detected_mod"]
        plt.text(x-0.2,y1+1,int(y1),ha="center")
        plt.text(x+0.2,y2+1,int(y2),ha="center")
    plt.ylim(0,60)
    plt.legend(loc="lower center",bbox_to_anchor=(0.5,1.2))
    plt.tight_layout()
    plt.savefig("{}figures/Detectable_modifications_all_conditions_{}_mods_detected_mod_num_{}.png".format(directory,near,am))

# [e] Proportion of read identity (A,T,G,C) at annotated modifications with upstream and downstream base info
def annotation_mod_misincorporated_base(directory,sample):
    # mod info
    dfm = pd.read_csv("modContext.csv")
    dfm["canon_pos"] = dfm["canon_pos"].astype(str).str.replace("20.3","20a").str.replace("20.6","20b")
    dfp = dfm["canon_pos"].astype(str).str.split(".",expand=True).rename(columns={0:"canon_pos"})
    dfm["canon_pos"] = dfp["canon_pos"]
    modifications = ["I","Y","acp3U","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
    dfm = dfm[dfm["Annotated modification"].isin(modifications)].reset_index(drop=True)
    print("modifications:",set(dfm["Annotated modification"]))
    print(set(dfm["canon_pos"]))
    print(dfm.head())
    for s in sample:
        print("sample;",s)
        # output file
        out = pd.DataFrame()
        o = -1
        # mismatch proportion
        df = pd.read_csv("{}mods/mismatchTable.csv".format(directory,s),"\t")
        df = df[df['isodecoder'].str.contains(organism)].reset_index(drop=True)
        df['isodecoder'] = df['isodecoder'].str.replace("{}_".format(organism),"")
        df["isoacceptor"] = df['isodecoder'].str.replace("/","").str.replace(r"[0-9]+","").str[:-1]
        df = df[df["condition"]==s]
        print(df.head())
        replicate = set(df["bam"])
        print(len(replicate),"replicates")
        for rep, r in zip(replicate,range(1,len(replicate)+1)):
            print(rep,r)
            dfr = df[df["bam"]==rep].reset_index(drop=True)
            for m in dfm.index:
                cpos = dfm.loc[m,"canon_pos"]
                iso = dfm.loc[m,"isodecoder"]
                mod = dfm.loc[m,"Annotated modification"]
                identity = dfm.loc[m,"identity"]
                upstream = dfm.loc[m,"upstream"]
                downstream = dfm.loc[m,"downstream"]
                dfrm = dfr[(dfr["isodecoder"]==iso) & (dfr["canon_pos"]==cpos)].reset_index(drop=True)
                for di in dfrm.index:
                    o += 1
                    out.loc[o,"isodecoder"] = iso
                    out.loc[o,"canon_pos"] = cpos
                    out.loc[o,"Annotated modification"] = mod
                    out.loc[o,"original"] = identity
                    out.loc[o,"upstream"] = upstream
                    out.loc[o,"downstream"] = downstream
                    out.loc[o,"replicate"] = r
                    out.loc[o,"type"] = dfrm.loc[di,"type"]
                    if identity == dfrm.loc[di,"type"]:
                        out.loc[o,"proportion"] = 1-sum(dfrm["proportion"])
                    else:
                        out.loc[o,"proportion"] = dfrm.loc[di,"proportion"]
        out["type"] = out["type"].str.replace("T","U")
        print("done.")
        out.to_csv("{}RTstop_mismatch_summarized_{}_annotated_modification_misincorporated_base.csv".format(directory,s),index=False)

# [5] Context-dependency of read identity
def plot_base_preference(directory,sample,type_):
    modifications = ["m1A","m1G","m2,2G","m3C","I","m1I","Y","acp3U","ms2i6A/ms2t6A"]
    identity = ["A","U","G","C"]
    for s in sample:
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}_annotated_modification_misincorporated_base.csv".format(directory,s))
        print(df.head())
        sns.set()
        sns.set_style("ticks")
        sns.set_palette("Set3")
        fig = plt.figure(figsize=(2*len(modifications),10))
        m = 0
        for mod in modifications:
            m += 1
            print("modification:",mod)
            dfm = df[df["Annotated modification"]==mod]
            p = -1
            for pos in sorted(set(dfm[type_])):
                p += 1
                dfp = dfm[dfm[type_]==pos]
                ax = fig.add_subplot(4,len(modifications), m+(p*len(modifications)))
                flierprops = dict(marker='o', markersize=1, alpha=0.5)
                sns.boxplot(x="type",y="proportion",data=dfp,ax=ax,order=identity,flierprops=flierprops)
                for i in range(0,len(identity)-1):
                    signal = test(dfp[dfp["type"]==identity[i]]["proportion"],dfp[dfp["type"]==identity[i+1]]["proportion"],identity[i],identity[i+1])
                    if signal == "error":
                        plt.text(i+0.5,1.1,signal,ha="center",fontsize=8)
                        plt.plot([i+0.1,i+0.9],[1.05,1.05],c="black")
                    elif signal != "":
                        plt.text(i+0.5,1.05,signal,ha="center",fontsize=12)
                        plt.plot([i+0.1,i+0.9],[1.05,1.05],c="black")
                if type_ == "upstream":
                    ax.set_title(f"{pos}-{mod}")
                elif type_ == "downstream":
                    ax.set_title(f"{mod}-{pos}")
                elif type_ == "canon_pos":
                    ax.set_title(f"{mod}{pos}")
                ax.set_ylim(-0.1,1.2)
        plt.tight_layout()
        plt.savefig("{}figures/{}_preference_{}.png".format(directory,type_,s))

# [6] RTstop and misincorporation in individual tRNA isoacceptor
def RTstop_mismatch_individual_isoacceptor(directory,sample,tRNA,canon_pos,mod):
    fig = plt.figure(figsize=(5,2.7*len(sample)))
    n = 0
    for s in sample:
        print("sample:",s)
        n += 1
        ax = fig.add_subplot(len(sample),1,n)
        for label, proptype, color in zip(["Misincorporation","RTstop"],["mismatch","RTstop_eachpos"],["#3CA2C8","#DB4C77"]):
            df = pd.read_csv("{}heatmap_{}_isoacceptor_{}.csv".format(directory,proptype,s),index_col=0)
            df = df[df.index == tRNA].reset_index(drop=True)
            df = df.replace(-1,0)
            df_ = df.transpose()
            df_ = df_.rename(columns={0:"proportion"})
            df_["x"] = [n for n in range(len(df_))]
            print(df_.head())
            if proptype == "mismatch" and canon_pos != False:
                x_ = df_.loc[canon_pos,"x"]
                y_ = df_.loc[canon_pos,"proportion"]
                if y_ > 0.9:
                    plt.text(x_,y_+0.05,"{}:{}".format(mod,round(y_,2)),ha="center",fontsize=8,color=color)
                else:
                    plt.text(x_,y_+0.07,"{}{}{}".format(mod,"\n",round(y_,2)),ha="center",fontsize=8,color=color)
            ax.plot(df_["x"],df_["proportion"],color="black",linewidth=0.8)
            ax.fill_between(df_["x"],df_["proportion"],[0 for n in range(len(df_))],color=color,alpha=0.5,label=label)
            ax.set_ylim(0,1)
        ax.plot([0,len(df_)],[0.1,0.1],linestyle="dashed",color="gray",linewidth=0.5)
        plt.text(len(df_)-5,0.12,"0.1",color="gray",fontsize=8)
        ax.set_ylim(0,1.2)
        ax.set_xlim(0,len(df_))
        ax.set_xticks([9,36,39,77])
        ax.set_xticklabels(["9","34","37","58"])
        ax.set_xlabel("Canonical position",fontsize=10)
        ax.set_ylabel("Proportion",fontsize=10)
        ax.set_title(f"{tRNA} in {s}",fontsize=10)
        plt.legend(loc="lower center",bbox_to_anchor=(0.5,1.2))
    plt.tight_layout()
    fig.savefig("{}figures/RTstop_mismatch_in_{}".format(directory,tRNA))

# [7] Read identity at specific site (ex. position 37 in tRNAPro(AGG))
def search_misincorporated_base2(directory,sampledic,tRNA,pos,identity,mod,isotype,organism,sample):
    dfm = pd.read_csv('{}mods/mismatchTable.csv'.format(directory),sep="\t")
    dfm = dfm[dfm['isodecoder'].str.contains(organism)==True].reset_index(drop=True)
    dfm['isodecoder'] = dfm['isodecoder'].str.replace("{}_".format(organism),"")
    dfm = dfm.replace({'canon_pos': {'20a':'20.3', '20b':'20.6'}})
    dfm = dfm[(dfm['canon_pos'].str.contains("-") == False) & (dfm['canon_pos'].str.contains('[a-z]+') == False)].reset_index(drop=True)
    dfm['canon_pos'] = dfm['canon_pos'].astype("float").astype("str")
    dfm['type'] = dfm['type'].str.replace("T","U")
    dfm = dfm[(dfm["isodecoder"].str.contains(tRNA)) & (dfm["canon_pos"]=="{}".format(pos))]
    dfm = dfm[dfm["condition"]==sample]
    dfmt = dfm[~(dfm['isodecoder'].str.contains("mito"))].reset_index(drop=True)
    
    #plot each sample
    print(dfmt.columns)
    print("dfmt")
    print(dfmt)

    sns.set()
    sns.set_style('ticks')
    sns.set_palette('Set3')
    if isotype == 'isodecoder':
        fig = plt.figure(figsize=(4*len(set(dfmt['isodecoder'])),2*len(sample)))
    elif isotype == 'isoacceptor':
        for t in set(dfmt['isodecoder']):
            for b in set(dfmt["bam"]):
                inde = dfmt[(dfmt["isodecoder"]==t)&(dfmt["bam"]==b)&(dfmt["type"]==identity)].index
                dfmt.loc[inde,'proportion'] = 1-sum(dfmt[(dfmt["bam"]==b)&(dfmt["isodecoder"]==t)]["proportion"])
            print("dfmt")
            print(dfmt)
        fig = plt.figure(figsize=(2.7,3))

    if isotype == 'isodecoder':
        for t in set(dfmt['isodecoder']):
            n += 1
            print(t)
            dfmstpt = dfmstp[dfmstp['isodecoder'] == t].reset_index(drop=True)
            boxplot(fig,dfmstpt,t,mod,len(sample),len(set(dfmt['isodecoder'])),s)

    elif isotype == 'isoacceptor':
        ax = fig.add_subplot(111)
        sns.barplot(x='type', y='proportion', data=dfmt, ax=ax, edgecolor="black",linewidth=1, order=["A","U","G","C"])
        sns.stripplot(x='type', y='proportion', data=dfmt, jitter=True, color='black', ax=ax, s=2, order=["A","U","G","C"])
        plt.ylim(-0.1,1.1)
        plt.xlabel('Base',fontsize=12)
        plt.ylabel('Proportion at {}{}'.format(mod,int(pos)),fontsize=12)
        plt.subplots_adjust(wspace=0.8,hspace=0.6)

    plt.tight_layout()
    plt.savefig("{}figures/misincorporation_preference_{}_at_{}_{}.png".format(directory,tRNA,mod,sample))

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
sample = sorted(set(sampledic['sample']))
print("sample:",sample)

# [1] 3'-CCA end analysis
CCAanalysis(directory,sampledic)

# [2] Differential abundance analysis
differential_expression(directory,sampledic,'log10',organism)

# [a] Summary table of RTstop, misincorporation, & readthrough proportions at individual position
RTstop_mismatch(directory,sampledic,organism)

# [b] Table for heatmap of misincorporation and RTStop
RTstop_mismatch_for_heatmap(directory,sample,"isoacceptor")

# [3] heatmap of misincorporation and RTStop
heatmap_RTsop_mismatch(directory,sample)

# [c] Add info of annotated tRNA modificartions in Nedialkova's paper
annotation_mod(directory,sample)

# [d] Summary table of RTstop & misincorporation proportions at modification sites
detectable_modifications(directory,sample,"max")

# [4] Comparison of readthrough proportion at annotated modification sites between two samples
compare_readthrough(directory,sample2)

# [e] Proportion of read identity (A,T,G,C) at annotated modifications with upstream and downstream base info
annotation_mod_misincorporated_base(directory,sample)

# [5] Context-dependency of read identity
for type_ in ["upstream", "downstream", "canon_pos"]:
    plot_base_preference(directory,sample2,type_)

# [6] RTstop and misincorporation in individual tRNA isoacceptor
RTstop_mismatch_individual_isoacceptor(directory,sample,"tRNA-Pro-AGG","34","I")
RTstop_mismatch_individual_isoacceptor(directory,sample,"tRNA-Pro-TGG","37","m1G")
RTstop_mismatch_individual_isoacceptor(directory,sample,"tRNA-Pro-CGG","37","m1G")

# [7] Read identity at specific site (ex. position 37 in tRNAPro(AGG))
search_misincorporated_base2(directory,sampledic,"tRNA-Pro-AGG",34.0,"A","I","isoacceptor",organism,"K562_totalRNA")

