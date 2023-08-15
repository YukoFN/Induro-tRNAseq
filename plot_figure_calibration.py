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
directory = "/Volumes/YMH_tRNAseq/5th_tRNA-seq/Calibration/"
sampledata = "sample_data_calibration.txt"
organism = "Homo_sapiens"
folder = directory.split("/")[-2]
sample = ["G100_m1G0","G75_m1G25","G25_m1G75","G0_m1G100"]

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

# [a] Summary table of RTstop, misincorporation, & readthrough proportions at individual position
def RTstop_mismatch(directory,sampledic,organism):
    sample = sorted(set(sampledic['sample']))
    print("sample:",sample)

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

# [1] Misincorporation proportion as a function of the G9:m1G9 ratio of the transcript of human mt-Leu(TAA)
def calibration(directory,sample,calib):
    x = []
    y = []
    for s, c in zip(sample,calib):
        # print("sample:",s)
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        df = df[df["isoacceptor"] == "mito_tRNA-Leu-TAA"]
        df = df[df["canon_pos"]=="9"]
        df = df[df["prop_type"]=="mismatch"].reset_index()
        # print(len(df),"rows")
        # print(df.head())
        x += [c for l in df.index]
        y += [df.loc[l,"proportion"] for l in df.index]
    print("x:",x)
    print("y:",y)
    # plot
    grid = plt.GridSpec(1,4)
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(grid[0,:3])
    ax.scatter(x,y,c='blue',alpha=0.3,edgecolor="black",s=25)
    b, a = np.polyfit(x, y, deg=1)
    xseq = np.linspace(0, 1, num=100)
    ax.plot(xseq, a + b * xseq, color="k", lw=1)
    r, p_value = pearsonr(x,y)
    ax.text(0.1,0.8,"R = {}".format(round(r,2)))
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(-0.1,1.1)
    ax.tick_params(axis='x',labelsize=10)
    ax.tick_params(axis='y',labelsize=10)
    ax.set_xlabel("m1G9 level",fontsize=10)
    ax.set_ylabel("Misincorporation",fontsize=10)

    # Native tRNA
    x_ = []
    y_ = []
    file = ["/Volumes/YMH_tRNAseq/4th_tRNA-seq/Induro/RTstop_mismatch_summarized_Induro_42C_16h.csv","/Volumes/YMH_tRNAseq/3rd_tRNA-seq/3rd_HEK293T_vs_K562_020223/RTstop_mismatch_summarized_K562_totalRNA.csv"]
    for f in file:
        dfn = pd.read_csv(f)
        dfn = dfn[(dfn["isoacceptor"] == "mito_tRNA-Leu-TAA") & (dfn["canon_pos"]=="9") & (dfn["prop_type"]=="mismatch")].reset_index()
        print(dfn.head())
        x_ += [0 for l in dfn.index]
        y_ += [dfn.loc[l,"proportion"] for l in dfn.index]
    print("x_:",x_)
    print("y_:",y_)
    ax2 = fig.add_subplot(grid[0,3:])
    ax2.scatter(x_,y_,c='blue',alpha=0.3,edgecolor="black",s=25)
    ax2.set_ylim(-0.1,1.1)
    ax2.set_xlim(-0.5,0.5)
    ax2.set_xticks([0])
    ax2.axes.yaxis.set_ticklabels([])
    ax2.tick_params(axis='x',labelsize=10)
    ax2.tick_params(axis='y',labelsize=10)
    plt.tight_layout()
    plt.savefig('{}figures/Calibration_m1G9_mito_tRNA-Leu-TAA.png'.format(directory))




# [0] Mapping rates
sampledic = plot_mapping_rate(directory, sampledata)
# print("sampledic:",sampledic)
print("sample:",sample)

# [a] Summary table of RTstop, misincorporation, & readthrough proportions at individual position
RTstop_mismatch(directory,sampledic,organism)

# [1] Misincorporation proportion as a function of the G9:m1G9 ratio of the transcript of human mt-Leu(TAA)
calibration(directory,sample,[0,0.25,0.75,1])

