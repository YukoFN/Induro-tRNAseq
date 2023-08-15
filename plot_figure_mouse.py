#!/usr/bin/env python
# coding: utf-8

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
directory = '/Volumes/YMH_tRNAseq/3rd_tRNA-seq/3rd_mouse_wt_vs_mut/'
sampledata = 'sample_data_mouse.txt'
organism = 'Mus_musculus'

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
    fig = plt.figure(figsize=(5,2*len(sample)))
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
    fig = plt.figure(figsize=(5,1*len(sample)))
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
    fig = plt.figure(figsize=(3.5*len(sample),10))
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

# [2] Differential abundance analysis
def differential_expression(directory,sampledic,log10,organism):
    # tRNA: "tRNA-Xxx-XXX-X", or "tRNA-Xxx-XXX", or "mito_tRNA-Xxx-XXX"

    for file in ['cyto_Isodecoder','cyto_Anticodon','mito_Isodecoder']:
        print("file:",file)
        if file == 'cyto_Anticodon':path = "{}DESeq2/cyto/anticodon/".format(directory)
        elif file == 'cyto_Isodecoder':path = "{}DESeq2/cyto/isodecoder/".format(directory)
        elif file == 'mito_Isodecoder':path = "{}DESeq2/organelle/isodecoder/".format(directory)

        f1 = sampledic.loc[0,"sample"]
        f2 = sampledic.loc[1,"sample"]
        ff = [f1,f2]
        ff_ = ff
        col1 = [sampledic.loc[0,"file_unpaired_uniq"]]
        col2 = [sampledic.loc[1,"file_unpaired_uniq"]]
        df = pd.read_csv(f"{path}singleRep-normCounts.csv",usecols=["Gene"]+col1+col2).rename(columns={col1[0]:f1,col2[0]:f2})
        df = df[df['Gene'].str.contains(organism)].reset_index(drop=True)
        df['{}_log10'.format(f1)] = np.log10(df[f1])
        df['{}_log10'.format(f2)] = np.log10(df[f2])
        df.replace([np.inf,-np.inf],0,inplace=True)
        print(df.head())
        print("len(file):",len(df))

        # plot
        fig = plt.figure(figsize=(2.9,3.1))
        ax = fig.add_subplot(111)

        if log10 == 'log10':
            xcol = '{}_log10'.format(f1)
            ycol = '{}_log10'.format(f2)
        else:
            xcol = f1
            ycol = f2

        x = df[xcol]
        y = df[ycol]
        ax.scatter(x,y,c='darkgray',alpha=0.3,edgecolor="darkgray",s=10)
        ax.plot([0,6],[0,6],linewidth=0.5,c='gray')
        r, p_value = pearsonr(df[f1],df[f2])
        ax.text(4,0.5,"R = {}".format(round(r,2)))
        ax.set_ylim(0,6)#maxval)
        ax.set_xlim(0,6)#maxval)
        ax.tick_params(axis='x',labelsize=10)
        ax.tick_params(axis='y',labelsize=10)
        ax.set_xlabel(xcol,fontsize=10)
        ax.set_ylabel(ycol,fontsize=10)
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
                    # print(len(dfi),"rows")
                    if i == "tRNA-Arg-TCT-4":
                        print("tRNA: ",i)
                    for p in sorted(set(dfi["pos"])):
                        ind += 1
                        dfp = dfi[dfi["pos"]==p].reset_index(drop=True)
                        # print(dfp)
                        # print("len:",len(dfp))
                        if file == "mismatch":
                            if i == "tRNA-Arg-TCT-4" and dfp.loc[0,"canon_pos"] == "50":
                                print(i,dfp.loc[0,"canon_pos"])
                                print(dfp.loc[0,:])
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
                            if i == "tRNA-Arg-TCT-4" and dfp.loc[0,"canon_pos"] == "50" and file == "mismatch":
                                out.loc[ind,"proportion"] = 1 - sum(dfp["proportion"])
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
            out = pd.read_csv("{}heatmap_{}_isoacceptor_{}.csv".format(directory,proptype,s),index_col=0)
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
            ax.spines[:].set_visible(True)
            plt.xticks(fontsize=10,rotation=90)
            plt.yticks(fontsize=10)
        plt.tight_layout()
        fig.savefig("{}figures/heatmap_isoacceptor_{}".format(directory,s))

# [4] RTstop and misincorporation in individual tRNA isodecoder
def RTstop_mismatch_individual_isodecoder(directory,sample,tRNA,canon_pos,mod):
    fig = plt.figure(figsize=(5,2.7*len(sample)))
    n = 0
    for s in sample:
        print("sample:",s)
        n += 1
        ax = fig.add_subplot(len(sample),1,n)
        for label, proptype, color in zip(["Misincorporation","RTstop"],["mismatch","RTstop_eachpos"],["#3CA2C8","#DB4C77"]):
            df = pd.read_csv("{}heatmap_{}_isodecoder_{}.csv".format(directory,proptype,s),index_col=0)
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
                elif y_ < 0.1:
                    plt.text(x_,y_+0.15,"{}:{}".format(mod,round(y_,2)),ha="center",fontsize=8,color=color)
                else:
                    plt.text(x_,y_+0.07,"{}{}{}".format(mod,"\n",round(y_,2)),ha="center",fontsize=8,color=color)
            ax.plot(df_["x"],df_["proportion"],color="black",linewidth=0.8)
            ax.fill_between(df_["x"],df_["proportion"],[0 for n in range(len(df_))],color=color,alpha=0.5,label=label)
            ax.set_ylim(0,1)
        ax.plot([0,len(df_)],[0.1,0.1],linestyle="dashed",color="gray",linewidth=0.5)
        plt.text(len(df_)-5,0.12,"0.1",color="gray",fontsize=8)
        ax.set_ylim(0,1.2)
        ax.set_xlim(0,len(df_))
        ax.set_xticks([9,28,34,69,77])
        ax.set_xticklabels(["9","26","32","50","58"])
        ax.set_xlabel("Canonical position",fontsize=10)
        ax.set_ylabel("Proportion",fontsize=10)
        ax.set_title(f"{tRNA} in {s}",fontsize=10)
        plt.legend(loc="lower center",bbox_to_anchor=(0.5,1.2))
    plt.tight_layout()
    fig.savefig("{}figures/RTstop_mismatch_in_{}".format(directory,tRNA))

# [5] RTstop and misincorporation in all tRNA isodecoders
def RTstop_mismatch_all_isoacceptor(directory,sample):
    fig = plt.figure(figsize=(5,2.7*len(sample)))
    n = 0
    for s in sample:
        print("sample:",s)
        n += 1
        ax = fig.add_subplot(len(sample),1,n)
        out = pd.DataFrame()
        for label, proptype, color in zip(["Misincorporation","RTstop"],["mismatch","RTstop_eachpos"],["#3CA2C8","#DB4C77"]):
            df = pd.read_csv("{}heatmap_{}_isoacceptor_{}.csv".format(directory,proptype,s),index_col=0)
            df = df.replace(-1,0)
            print(df.head())
            print(df.columns)
            print(df.describe().loc["mean"])
            y = list(df.describe().loc["mean"])
            print(len(df.columns)==len(y))
            x = [n for n in range(len(y))]
            ax.plot(x,y,color="black",linewidth=0.8)
            ax.fill_between(x,y,[0 for n in range(len(y))],color=color,alpha=0.5,label=label)
            ax.set_ylim(0,1)
        ax.plot([0,len(y)],[0.1,0.1],linestyle="dashed",color="gray",linewidth=0.5)
        ax.set_ylim(0,1.1)
        ax.set_xlim(0,len(y))
        ax.set_xticks([9,28,34,36,39,77])
        ax.set_xticklabels(["9","26","32","34","37","58"])
        ax.set_xlabel("Canonical position",fontsize=10)
        ax.set_ylabel("Proportion",fontsize=10)
        ax.set_title(f"All tRNA in {s}",fontsize=10)
        plt.legend(loc="lower center",bbox_to_anchor=(0.5,1.2))
    plt.tight_layout()
    fig.savefig("{}figures/RTstop_mismatch_in_all_tRNAs.png".format(directory))




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
RTstop_mismatch_for_heatmap(directory,sample,"isodecoder")

# [3] heatmap of misincorporation and RTStop
heatmap_RTsop_mismatch(directory,sample)

# [4] RTstop and misincorporation in individual tRNA isodecoder
RTstop_mismatch_individual_isodecoder(directory,sample,"tRNA-Arg-TCT-4","50","U")

# [5] RTstop and misincorporation in all tRNA isodecoders
RTstop_mismatch_all_isoacceptor(directory,sample)



