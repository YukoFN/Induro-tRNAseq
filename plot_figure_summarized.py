import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns
import colorsys
import matplotlib.colors as mc
from math import ceil #小数点以下を切り上げする
import glob
from scipy.stats import pearsonr
from scipy import stats





# [0] plot mapping rates
def plot_mapping_rate(directory, sampledata):
    # sample data
    df = pd.read_csv('{}{}'.format(directory,sampledata),sep="\t",header=None,names=['file','sample'])
    print(df.head())
    # df['file_unpaired_uniq'] = [0 for i in df.index]
    for i in df[df['file'].str.contains('.fastq.gz')].index:
        df.loc[i,'file_unpaired_uniq'] = df.loc[i,'file'][:-9]+'.unpaired_uniq'
    for i in df[df['file'].str.contains('.fq.gz')].index:
        df.loc[i,'file_unpaired_uniq'] = df.loc[i,'file']+'.unpaired_uniq'
    print(df[df['file_unpaired_uniq']==0])

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
    # for 積み立てグラフ
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

    # plot 積み立てグラフ
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

# [1] plot CCAcount
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
        for y_, l, c in zip(y,['Absent','C','CC','CCA'],['#012E54','#0F4C81','#92B0CA','#E9F4FF']):
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

# [2] countをscatter plotする
def differential_expression(directory,sampledic,tRNA,log10,organism,significance):
    # tRNA: "tRNA-Xxx-XXX-X", or "tRNA-Xxx-XXX", or "mito_tRNA-Xxx-XXX"
    
    if len(tRNA)>0:
        if 'mito' in tRNA[0]:searchtype = 'mito_Isodecoder'
        else:
            if tRNA[0].count('-') == 2:searchtype = 'cyto_Anticodon'
            elif tRNA[0].count('-') == 3:searchtype = 'cyto_Isodecoder'
    else:
        print("searchtype = False")
        searchtype = False
    
    
    for file in ['cyto_Isodecoder']:#,'cyto_Anticodon','mito_Isodecoder']:

        print("file:",file)
        if file == 'cyto_Anticodon':path = "{}DESeq2/cyto/anticodon/".format(directory)
        elif file == 'cyto_Isodecoder':path = "{}DESeq2/cyto/isodecoder/".format(directory)
        elif file == 'mito_Isodecoder':path = "{}DESeq2/organelle/isodecoder/".format(directory)
        files = glob.glob("{}*results.csv".format(path))

        if searchtype == file:
            fig = plt.figure(figsize=(2.9+len(tRNA),3.1*(len(files))))
            grid = plt.GridSpec(len(files),5+len(tRNA))
        else:
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
            print(df.head())
            print(df.columns)
            # print("total # of detected all tRNAs in {}:".format(f1),len(df[(df['{}_ave'.format(f1)]>0)]))
            # print("total # of detected cyto-tRNAs in {}:".format(f1),len(df[(df['{}_ave'.format(f1)]>0) & ~(df["Gene"].str.contains("mito"))]))
            # print("total # of detected mito-tRNAs in {}:".format(f1),len(df[(df['{}_ave'.format(f1)]>0) & (df["Gene"].str.contains("mito"))]))
            # print("total # of detected all tRNAs in {}:".format(f2),len(df[(df['{}_ave'.format(f2)]>0)]))
            # print("total # of detected cyto-tRNAs in {}:".format(f2),len(df[(df['{}_ave'.format(f2)]>0) & ~(df["Gene"].str.contains("mito"))]))
            # print("total # of detected mito-tRNAs in {}:".format(f2),len(df[(df['{}_ave'.format(f2)]>0) & (df["Gene"].str.contains("mito"))]))
            if file == "mito_Isodecoder":
                print("Ser in {}".format(f1),df[df["Gene"].str.contains("Ser")]['{}_ave'.format(f1)])
                print("Ser in {}".format(f2),df[df["Gene"].str.contains("Ser")]['{}_ave'.format(f2)])

            #plot
            n += 1
            if len(tRNA) > 0:
                if len(df[df['Gene'] == '{}_{}'.format(organism,tRNA[0])])>0:
                    ax = fig.add_subplot(grid[n-1,0:5]) #all tRNA
                    axt  = fig.add_subplot(grid[n-1,5:]) #specific tRNA
            else:
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
            maxval = ceil(max([max(x),max(y)])/1)*1 #ceil: 小数点の切り上げ
            #ax.plot([0,maxval],[0,maxval],linewidth=0.5,c='gray')
            ax.plot([0,6],[0,6],linewidth=0.5,c='gray')
            r, p_value = pearsonr(df['{}_ave'.format(f1)],df['{}_ave'.format(f2)])
            #ax.text(maxval-1.5,0.5,"r2 = {}".format(round(r*r,2)))
            ax.text(maxval-2,0.5,"R = {}".format(round(r,2)))

            if significance == True:
                padj = 0.001
                alpha = 0.3
                # Others
                dfo = df[(df['padj']>padj)]
            else:
                dfo = df
                alpha = 0.3
            xx = dfo[xcol]
            yy = dfo[ycol]
            print("len(xx)",len(xx))
            print("type(xx)",type(xx))
            # print(xx)
            print("len(yy)",len(yy))
            print("type(xx)",type(yy))
            # print(yy)
            ax.scatter(xx,yy,c='darkgray',alpha=alpha,edgecolor="darkgray",label="all")#c='lightgray'

            if significance == True:
                # Up tRNAs
                dfu = df[(df['log2FoldChange']>0) & (df['padj']<=padj)]
                xx = dfu[xcol]
                yy = dfu[ycol]
                # 薄い色の組み合わせ。pink and green, up: #FF91A4, down : #83C760
                # 濃い色の組み合わせ。pink and green, up: #FB607F, down : #66B032
                # 濃い色の組み合わせ。pink and blue, up: #FB607F, down : #2694BF
                
                if ff != ff_:
                    ax.scatter(xx,yy,label='{} up tRNAs'.format(len(dfu)),c='#FB607F',alpha=alpha)#,marker='^',edgecolor='firebrick',linewidth=0.5,c='lightgray')
                else:
                    ax.scatter(xx,yy,label='{} down tRNAs'.format(len(dfu)),c='#2694BF',alpha=alpha)#,marker='v',edgecolor='darkblue',linewidth=0.5,c='lightgray')

                # Down tRNAs
                dfd = df[(df['log2FoldChange']<0) & (df['padj']<=padj)]
                xx = dfd[xcol]
                yy = dfd[ycol]
                if ff != ff_:
                    ax.scatter(xx,yy,label='{} down tRNAs'.format(len(dfd)),c='#2694BF',alpha=alpha)#,marker='v',edgecolor='darkblue',linewidth=0.5,c='lightgray')
                else:
                    ax.scatter(xx,yy,label='{} up tRNAs'.format(len(dfd)),c='#FB607F',alpha=alpha)#,marker='^',edgecolor='firebrick',linewidth=0.5,c='lightgray')

            # specific tRNA
            if len(tRNA) > 0:
                if len(df[df['Gene'] == '{}_{}'.format(organism,tRNA[0])])>0:
                    tn = 0
                    colort = ['#FDEA79','#ABD6DFFF']
                    mint = []
                    fold = []
                    tnx = []
                    for t in tRNA:
                        tn += 1
                        dft = df[df['Gene'] == '{}_{}'.format(organism,t)].reset_index(drop=True)
                        #print(dft)
                        if len(dft) > 0:
                            ttype = False
                            if (dft.loc[0,'log2FoldChange'] < 0) & (dft.loc[0,'padj']<=padj):
                                if ff != ff_:ttype='down'
                                else:ttype='up'
                            elif (dft.loc[0,'log2FoldChange'] > 0) & (dft.loc[0,'padj']<=padj):
                                if ff != ff_:ttype='up'
                                else:ttype='down'
                            if ttype == 'up':
                                marker='^'
                                c='#FB607F'
                            elif ttype == 'down':
                                marker='v'
                                c='#66B032'
                            else:
                                marker='o'
                                c='lightgray'
                            ax.scatter(dft[xcol],dft[ycol],label=t.replace("T","U"),marker=marker,c=c,edgecolor='black',linewidth=0.5)
                            # add scatter plot
                            tnxx = tn/(len(tRNA)+1)
                            tnx += [tnxx]
                            scat1 = axt.scatter([tnxx for i in range(len(col1))],np.log10(dft[col1])
                                        ,c=colort[0],edgecolor='black',linewidth=0.5) #label=f1
                            axt.plot([(tn-0.4)/(len(tRNA)+1),(tn+0.4)/(len(tRNA)+1)],
                                     [dft['{}_ave_log10'.format(f1)],dft['{}_ave_log10'.format(f1)]],c=colort[0])
                            scat2 = axt.scatter([tnxx for i in range(len(col2))],np.log10(dft[col2])
                                        ,c=colort[1],edgecolor='black',linewidth=0.5) #label=f2
                            axt.plot([(tn-0.4)/(len(tRNA)+1),(tn+0.4)/(len(tRNA)+1)],
                                     [dft['{}_ave_log10'.format(f2)],dft['{}_ave_log10'.format(f2)]],c=colort[1])
                            mint += [min(dft.loc[0,['{}_ave_log10'.format(f) for f in [f1,f2]]])]
                            fold += ['x{}'.format(round(dft.loc[0,'{}_ave'.format(f2)]/dft.loc[0,'{}_ave'.format(f1)],2))]

                    mint_ = [min(mint),min(mint)]
                    for x, m, f in zip(tnx, mint_,fold):
                        axt.text(x,m-1,f,horizontalalignment="center",verticalalignment="center",size=8)

                    axt.legend([scat1,scat2],[f1,f2],loc='lower center', borderaxespad=0,bbox_to_anchor=(0.5,1.02))
                    axt.set_xlim(0,1),
                    axt.set_ylim(min(mint)-0.7)
                    axt.set_xticks([tn/(len(tRNA)+1) for tn in range(1,len(tRNA)+1)])
                    axt.set_xticklabels([t.replace("tRNA-",'').replace("T","U") for t in tRNA],rotation=45)
                    axt.set_ylabel('Normalized count',fontsize=10)

                        
            ax.set_ylim(0,6)#maxval)
            ax.set_xlim(0,6)#maxval)
            ax.tick_params(axis='x',labelsize=10)
            ax.tick_params(axis='y',labelsize=10)
            ax.set_xlabel(xcol,fontsize=10)
            ax.set_ylabel(ycol,fontsize=10)
            ax.set_title("{}, All {} tRNAs".format(file,len(df)),fontsize=10)
            if significance == True:
                ax.legend()
            plt.subplots_adjust(hspace=0.5,wspace=0.5)
            plt.tight_layout()
        #plt.show()
        if significance == True:
            plt.savefig('{}figures/differential_expression_{}_padj{}_alpha{}.png'.format(directory,file,padj,alpha))
        else:
            plt.savefig('{}figures/differential_expression_{}.png'.format(directory,file))
def differential_expression_mouse(directory,sampledic,log10,organism,significance):
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

# [4] Isodecoder毎に、modificationの判定をする
# 各isodecoderごとに、各位置のRTstopとmismatchのproportionをまとめた表
def RTstop_mismatch(directory,folder,sampledic,organism):#10%のthreshold判定をなくしてみた
    sample = sorted(set(sampledic['sample']))
    print("sample:",sample)
    
    #output files
    # ol = ['{}'.format(n) for n in range(20)]
    # ol.extend(['20','20.3','20.6'])
    # ol.extend(['{}'.format(n) for n in range(21,76)])
    # out = pd.DataFrame()

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
            # print(set(df['isoacceptor']))

            for rep in sampledic[sampledic["sample"]==s].index:
                repnum = sampledic.loc[rep,"replicate_num"]
                bam = sampledic.loc[rep,"file"]#.split(".")[0]
                # print(sampledic.loc[rep,"file"])
                # print(bam)
                dfr = df[df["bam"] == f"{folder}{bam}.unpaired_uniq.bam"]
                # print(set(dfr["bam"]))
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
                                # print(dfi[dfi["pos"]==p-1])
                                # print(dfp["cov"])
                                # print(dfi[dfi["pos"]==p-1]["cov"]/dfp.loc[0,"cov"])
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
    # return out
# [5] % Full length reads
def full_length(directory,sample):
    out = pd.DataFrame()
    o = -1
    for s in sample:
        print("sample:",s)
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        df["cov"] = df["cov"].astype("float")
        # print(df.head())
        # print(df[df["canon_pos"]=="0"])
        for i in set(df["isodecoder"]):
            # print("isodecoder:",i)
            dfi = df[df["isodecoder"]==i]
            # print(dfi.head())
            for r in set(dfi["replicate_num"]):
                o += 1
                dfr = dfi[dfi["replicate_num"]==r].reset_index(drop=True)
                # print("replicate:",r)
                cov = dfr.loc[0,"cov"]
                if np.isnan(cov):
                    print("NaN in",i)
                    break
                if cov != min(dfr["cov"]):
                    print("min(coverage) does not equal coverage at position 1 in",i)
                    print("min cov",cov)
                    print(dfr.head())
                    # print(dfr.tail())
                # print(cov)
                # print(max(dfr["cov"]))
                # print(cov/max(dfr["cov"])*100)
                out.loc[o,"condition"] = s
                out.loc[o,"cyto_mito"] = dfr.loc[0,"cyto_mito"]
                out.loc[o,"full_length(%)"] = cov/max(dfr["cov"])*100
    print("done")
    # print(out.head())
    # print(out.tail())
    out.to_csv("{}full_length_reads.csv".format(directory),index=False)

    # plot
    sns.set()
    sns.set_style("ticks")
    # color = ["#d0dfed","#7fabdb","#2376c8","#d9c0d3","#d97aa6","#d3095f"]
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
        # sns.boxplot(x="condition",y="full_length(%)",data=data,ax=ax,palette=color,order=sample,flierprops=flierprops)
        sns.violinplot(x="condition",y="full_length(%)",data=data,ax=ax,palette=color,order=sample,flierprops=flierprops)
        for l in range(len(sample)-1):
            # print(sample[l],"vs",sample[l+1])
            signal = test(data[data["condition"]==sample[l]]["full_length(%)"],data[data["condition"]==sample[l+1]]["full_length(%)"],sample[l],sample[l+1])
            if signal != "":
                plt.text(l+0.5,102,signal,ha="center",fontsize=20)
                plt.plot([l+0.1,l+0.9],[107,107],c="black")
        ax.set_title(f"{cm} tRNAs")
        plt.ylim(-10,120)
        plt.xticks(rotation=45)
        plt.tight_layout()
    plt.savefig("{}figures/Full-length_reads_version3.png".format(directory))

# [6] Sum of misincorporation proportions
def ave_mismatch(directory,sample):
    out = pd.DataFrame()
    o = -1
    for s in sample:
        print("sample:",s)
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        df = df[df["prop_type"]=="mismatch"].reset_index(drop=True)
        print(df.head())
        # print(df[df["canon_pos"]=="0"])
        for i in set(df["isodecoder"]):
            # print("isodecoder:",i)
            dfi = df[df["isodecoder"]==i]
            # print(dfi.head())
            for r in set(dfi["replicate_num"]):
                o += 1
                dfr = dfi[dfi["replicate_num"]==r].reset_index(drop=True)
                # print("len:",len(dfr))
                # print(dfr.head())
                # out.loc[o,"isodecoder"] = i
                out.loc[o,"cyto_mito"] = dfr.loc[0,"cyto_mito"]
                out.loc[o,"condition"] = s
                out.loc[o,"ave_mismatch"] = sum(dfr["proportion"])/len(dfr)
    print("done")
    print(out.head())
    # print(out.tail())
    out.to_csv("{}ave_misincorporation.csv".format(directory),index=False)

    # plot
    sns.set()
    sns.set_style("ticks")
    # color = ["#d0dfed","#7fabdb","#2376c8","#d9c0d3","#d97aa6","#d3095f"]
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
            # print(sample[l],"vs",sample[l+1])
            signal = test(data[data["condition"]==sample[l]]["ave_mismatch"],data[data["condition"]==sample[l+1]]["ave_mismatch"],sample[l],sample[l+1])
            if signal != "":
                plt.text(l+0.5,round(max(data["ave_mismatch"]),0)+1,signal,ha="center",fontsize=20)
                plt.plot([l+0.1,l+0.9],[round(max(data["ave_mismatch"]),0)+1,round(max(data["ave_mismatch"]),0)+1],c="black")
        ax.set_title(f"{cm} tRNAs")
        plt.ylim(-0.01,0.08)
        plt.xticks(rotation=45)
        plt.tight_layout()
    plt.savefig("{}figures/ave_misincorporation_version2.png".format(directory))

# [7] Nedialkova's paperのannotated tRNA modificartions を、自分の結果にannotateする。
def annotation_mod(directory,sample):
    # mod info
    dfm = pd.read_csv("/Volumes/YMH_tRNAseq/modContext_updated.csv")
    dfm["canon_pos"] = dfm["canon_pos"].astype(str).str.replace("20.3","20a").str.replace("20.6","20b")#.str.split(".")[0]
    dfp = dfm["canon_pos"].astype(str).str.split(".",expand=True).rename(columns={0:"canon_pos"})
    dfm["canon_pos"] = dfp["canon_pos"]
    modifications = ["I","Y","acp3U","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
    # dfm = dfm[dfm["Annotated modification"].isin(modifications)].reset_index(drop=True)
    print("modifications:",set(dfm["Annotated modification"]))
    print(set(dfm["canon_pos"]))
    print(dfm.head())
    for s in sample:
        print("sample;",s)
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}.csv".format(directory,s))
        # df["canon_pos"] = df["canon_pos"].str.replace("20a","20.3").str.replace("20b","20.6")
        df = df[(df["prop_type"]=="mismatch") | (df["prop_type"]=="RTstop_eachpos")].reset_index(drop=True)
        # print(set(df["canon_pos"]))
        # print(df.head())
        # print(set(df["prop_type"]))
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

# [8] RTstop & misincorporation pattern at modification sites
def detectable_modifications(directory,sample,cal,plot2):
    
    if plot2 == False:
        for cm in ["cyto","mito","all"]:
            print(cm,"tRNAs")
            # modifications = ['ms2i6A/ms2t6A', 'm1A', 'D', 't6A', 'm1G', 'm3C', 'm2,2G', 'τm5s2U', 'Q', 'U', 'm5U', 'Y', 'acp3U', 'f5C', 'i6A', 'C', 'Ψ', 'm5C', 'G', 'm1I', 'm2G', 'I']
            modifications = ["acp3U","I","Y","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
            # plot 1
            sns.set()
            sns.set_style("ticks")
            fig = plt.figure(figsize=(2.5*len(sample),4*len(modifications)))
            # color = ["#00b5a1","#ff416d"]
            color = ["#3CA2C8","#DB4C77"]#,"#F9FFF8","lightgray"] #label = ['mismatch_here',"RTstop_here",'Other','RTstop_before']
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
                                # if len(modlist) > 0:
                                #     for ml in modlist:
                                        # if ml not in modifications:
                                        #     print(modlist)
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
                            # sns.boxplot(x="mod_pos",y="proportion",data=data,ax=ax,hue="prop_type",palette=color)
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
            plt.savefig("{}figures/Detectable_modifications_{}-tRNAs_version2.png".format(directory,cm))
        out.to_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,cal),index=False)

    if plot2 == True:
        # plot 2
        modifications = ["acp3U","I","Y","m1A","m1G","m1I","m2,2G","m3C","ms2i6A/ms2t6A"]
        out = pd.read_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,cal))
        print(out.head())
        for near in ["all","near"]:
            fig2 = plt.figure(figsize=(8,3*len(modifications)))
            # color2 = ["#d0dfed","#7fabdb","#2376c8","#d9c0d3","#d97aa6","#d3095f"]
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
                            # sns.stripplot(x="condition",y="proportion",data=data3,ax=ax2,palette=color2,order=sample,s=3)
                            for l in range(len(sample)-1):
                                # print(sample[l],"vs",sample[l+1])
                                signal = test(data3[data3["condition"]==sample[l]]["proportion"],data3[data3["condition"]==sample[l+1]]["proportion"],sample[l],sample[l+1])
                                if signal != "":
                                    plt.text(l+0.5,round(max(data3["proportion"]),1)+0.1,signal,ha="center",fontsize=12)
                                    plt.plot([l+0.1,l+0.9],[round(max(data3["proportion"]),1)+0.1,round(max(data3["proportion"]),1)+0.1],c="black")
                            ax2.set_title("{} at {}".format(pt,m))
                            plt.ylim(-0.1,1.4)
                            # ax2.get_legend().remove()
                            # if pt == "mismatch":
                            #     plt.ylim(-0.1,1.4)
                            # else:
                            #     if max(data3["proportion"]) < 1:
                            #         plt.ylim(-0.1,1.1)
                            #     else:
                            #         plt.ylim(-0.1,round(max(data3["proportion"]),1)+0.5)
            plt.tight_layout()
            plt.savefig("{}figures/Detectable_modifications_all_conditions_{}_mods_{}_version2.png".format(directory,near,cal))

# [9] Readthrough comparison between two samples
def compare_readthrough(directory,sample):
    for am in ["ave","max"]:
        # plot 1
        modifications = ["m1A","m1G","m2,2G","m3C","I","m1I","Y","acp3U","ms2i6A/ms2t6A"]
        out = pd.read_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,am))
        out = out[(out["condition"]==sample[0]) | (out["condition"]==sample[1])]
        print(out.head())
        for near in ["all","near"]:
            fig2 = plt.figure(figsize=(6,10))
            color = ["#B564E3","#D4D4D4"]#["#f7feae","#00718b"]
            b = 0
            if near == "near":
                data2 = out[out["near_mod"]=="0"]
            else:
                data2 = out
            for pt in ["mismatch","RTstop_eachpos","readthrough"]:
                print("prop_type:",pt)
                if pt == "readthrough":
                    data3 = data2[data2['prop_type']=="RTstop_eachpos"].reset_index(drop=True)
                    data3[pt] = 3-data3["proportion"]
                    y ="readthrough"
                else:
                    data3 = data2[data2['prop_type']==pt].reset_index(drop=True)
                    y = "proportion"
                if len(data3) > 0:
                    b += 1
                    ax2 = fig2.add_subplot(3,1, b)
                    flierprops = dict(marker='o', markersize=1, alpha=0.5)
                    sns.boxplot(x="Annotated modification",y=y,data=data3,ax=ax2,palette=color,hue="condition",order=modifications,linewidth=1, flierprops=flierprops)
                    signal = test(data3[data3["condition"]==sample[0]][y],data3[data3["condition"]==sample[1]][y],sample[0],sample[1])
                    if signal != "":
                        plt.text(l+0.5,round(max(data3["proportion"]),1)+0.1,signal,ha="center",fontsize=12)
                        plt.plot([l+0.1,l+0.9],[round(max(data3[y]),1)+0.1,round(max(data3[y]),1)+0.1],c="black")
                    ax2.set_title(pt,fontsize=10)
                    # ax2.get_legend().remove()
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
                    data3 = data2[(data2["proportion"] > 0.1) & (data2["condition"]==s) & (data2["Annotated modification"]==m)].loc[:,["isodecoder","canon_pos","Annotated modification"]]
                    # print(data3.head())
                    # 重複を除いた行のカウント
                    print((~data3.duplicated()).sum(),"modifications")
                    data.loc[d,"condition"] = s
                    data.loc[d,"Annotated modification"] = m
                    data.loc[d,"num_detected_mod"] = (~data3.duplicated()).sum()
            # print(data.head())
            fig = plt.figure(figsize=(6,3))
            # color = ["#f7feae","#00718b"]
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

# [10] heatmap of misincorporation and RTStop
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
            for iso in iso_: #["mito_tRNA-Lys-TTT"]:
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
                    # print(out.loc[iso, out.isnull().any()])
            out.to_csv("{}heatmap_{}_{}_{}.csv".format(directory,proptype,isotype,s))
def heatmap_RTsop_mismatch(directory,sample):
    for s in sample:
        print("sample:",s)
        fig = plt.figure(figsize=(20,25))
        n = 0
        for proptype in ["mismatch","RTstop_eachpos"]:
            out = pd.read_csv("{}heatmap_{}_{}.csv".format(directory,proptype,s),index_col=0)
            out_mask = (out == -1) # proportion == -1部分だけ抽出 -> 別にプロット
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
            #ax.set_xticks(list(range(len(pos))))
            #ax.set_xticklabels(pos)
            ax.set_xticks([i+0.5 for i in range(len(out.columns))])
            ax.set_xticklabels(list(out.columns))
            ax.set_yticks([i+0.5 for i in range(len(iso_))])
            ax.set_yticklabels(iso_)
            plt.xticks(fontsize=10,rotation=90)
            plt.yticks(fontsize=10)
        plt.tight_layout()
        fig.savefig("{}figures/heatmap_{}".format(directory,s))

# [11] Annotation of known modifications & misincorporated bases
def annotation_mod_misincorporated_base(directory,sample):
    # mod info
    dfm = pd.read_csv("/Volumes/YMH_tRNAseq/modContext_updated.csv")
    dfm["canon_pos"] = dfm["canon_pos"].astype(str).str.replace("20.3","20a").str.replace("20.6","20b")#.str.split(".")[0]
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
        # print(set(df["canon_pos"]))
        # print(set(df["prop_type"]))
        for rep, r in zip(replicate,range(1,len(replicate)+1)):
            print(rep,r)
            dfr = df[df["bam"]==rep].reset_index(drop=True)
            for m in dfm.index:
                # print(dfm.loc[m,:])
                cpos = dfm.loc[m,"canon_pos"]
                iso = dfm.loc[m,"isodecoder"]
                mod = dfm.loc[m,"Annotated modification"]
                identity = dfm.loc[m,"identity"]
                upstream = dfm.loc[m,"upstream"]
                downstream = dfm.loc[m,"downstream"]
                dfrm = dfr[(dfr["isodecoder"]==iso) & (dfr["canon_pos"]==cpos)].reset_index(drop=True)
                # print(mod,"in",iso,"at",cpos)
                # print(dfrm.head())
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
def plot_upstream_position_preference(directory,sample,samplenum,modifications):
    # modifications = ["m1A","m1G","m2,2G","m3C","I","m1I","Y","acp3U","ms2i6A/ms2t6A"]
    identity = ["A","U","G","C"]
    for s, snum in zip(sample, samplenum):
        df = pd.read_csv("{}RTstop_mismatch_summarized_{}_annotated_modification_misincorporated_base.csv".format(directory,s))
        print(df.head())
        for mod in modifications:
            print("modification:",mod)
            dfm = df[df["Annotated modification"]==mod]
            canonpos = sorted(set(dfm["canon_pos"]))
            sns.set()
            sns.set_style("ticks")
            sns.set_palette("Set3")
            fig = plt.figure(figsize=(2*len(canonpos),10))
            p = 0
            for pos in canonpos:
                p += 1
                dfp = dfm[dfm["canon_pos"] == pos]
                u = -1
                for up in sorted(set(dfm["upstream"])):
                    u += 1
                    dfu = dfp[dfp["upstream"]==up]
                    print("upstream:",up,"at",pos)
                    print(dfu)
                    ax = fig.add_subplot(4,len(canonpos), p+(u*len(canonpos)))
                    flierprops = dict(marker='o', markersize=1, alpha=0.5)
                    sns.boxplot(x="type",y="proportion",data=dfu,ax=ax,order=identity,flierprops=flierprops)
                    for i in range(0,len(identity)-1):
                        signal = test(dfu[dfu["type"]==identity[i]]["proportion"],dfu[dfu["type"]==identity[i+1]]["proportion"],identity[i],identity[i+1])
                        if signal == "error":
                            plt.text(i+0.5,1.1,signal,ha="center",fontsize=8)
                            plt.plot([i+0.1,i+0.9],[1.05,1.05],c="black")
                        elif signal != "":
                            plt.text(i+0.5,1.05,signal,ha="center",fontsize=12)
                            plt.plot([i+0.1,i+0.9],[1.05,1.05],c="black")
                    ax.set_title(f"{mod}{pos}, x{len(dfu)/4/snum}")
                    ax.set_ylabel(f"Upstream {up}")
                    # ax.get_legend().remove()
            plt.tight_layout()
            plt.savefig("{}figures/{}_upstream_preference_{}.png".format(directory,mod,s))

# [12] plot individual tRNA isoacceptor
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
                # print(df_[df_.index==canon_pos])
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

# 有意差検定
def test(a,b,a_,b_):
    # print(a_,"vs",b_)
    A_var = np.var(a, ddof=1)  # Aの不偏分散
    B_var = np.var(b, ddof=1)  # Bの不偏分散
    A_df = len(a) - 1  # Aの自由度
    B_df = len(b) - 1  # Bの自由度
    if B_var == 0:
        return "error"
    else:
        f = A_var / B_var  # F比の値
        one_sided_pval1 = stats.f.cdf(f, A_df, B_df)  # 片側検定のp値 1
        one_sided_pval2 = stats.f.sf(f, A_df, B_df)   # 片側検定のp値 2
        two_sided_pval = min(one_sided_pval1, one_sided_pval2) * 2  # 両側検定のp値
        # print('F:       ', round(f, 3))
        # print('p-value: ', round(two_sided_pval, 3)) #p<0.05なら、不等分散 -> Welchのt検定
        signal = ""
        # Welchのt-test
        if round(two_sided_pval, 3) < 0.05:
            t, pval = stats.ttest_ind(a,b,equal_var=False)
            if pval < 0.01:
                signal = "**"
                print("Welch's t-test p-value: ",pval,a_,"vs",b_)
            elif pval < 0.05:
                signal = "*"
                print("Welch's t-test p-value: ",pval,a_,"vs",b_)
        else:
            t, pval = stats.ttest_ind(a,b,equal_var=True)
            if pval < 0.01:
                signal = "**"
                print("Student's t-test p-value: ",pval,a_,"vs",b_)
            elif pval < 0.05:
                signal = "*"
                print("Student's t-test p-value: ",pval,a_,"vs",b_)
        return signal


# 3rd run
# HEK293T vs K562 including Nedialkova's data

# directory = '/Volumes/YMH_tRNAseq/3rd_tRNA-seq/3rd_HEK293T_vs_K562_020223/'
# sampledata = 'sample_data_HEK293T-K562.txt'
# organism = 'Homo_sapiens'
# sample = ["HEK293T_totalRNA","HEK293T_tRNA","HEK293T_tRNA_DN","K562_totalRNA","K562_tRNA","K562_tRNA_DN"]
# sample = ["K562_totalRNA","K562_tRNA_DN"]

# directory = "/Volumes/YMH_tRNAseq/4th_tRNA-seq/Induro/"
# sampledata = "sample_data_Induro_for_analysis.txt"
# organism = 'Homo_sapiens'
# sample = ["Induro_42C_1h","Induro_42C_2h","Induro_42C_16h","Induro_55C_1h","Induro_55C_2h","Induro_55C_16h"]

# directory = '/Volumes/YMH_tRNAseq/3rd_tRNA-seq/3rd_mouse_wt_vs_mut/'
# sampledata = 'sample_data_mouse.txt'
# organism = 'Mus_musculus'
# sample = ["mouse_WT","mouse_mut"]

directory = "/Volumes/YMH_tRNAseq/5th_tRNA-seq/Calibration/"
sampledata = "sample_data_calibration.txt"
organism = "Homo_sapiens"
folder = "Calibration/"
sample = ["G100_m1G0","G75_m1G25","G25_m1G75","G0_m1G100"]

# [0] plot mapping rates for all samples
sampledic = plot_mapping_rate(directory, sampledata)
print("sampledic:",sampledic)

# [1] plot CCAcount for all samples
# CCAanalysis(directory,sampledic)

# [2] countをscatter plotする
# differential_expression(directory,sampledic,[],'log10',organism,False)
# for mouse
# differential_expression_mouse(directory,sampledic,'log10',organism,False)

# [4] Isodecoder毎に、modificationの判定をする
# RTstop_mismatch(directory,folder,sampledic,organism)

# [5] % Full length reads
# full_length(directory,sample)
# out = pd.read_csv("{}full_length_reads.csv".format(directory))
# print(out.head())

# [6] Sum of misincorporation proportions
# ave_mismatch(directory,sample)    

# [7] Nedialkova's paperのannotated tRNA modificartions を、自分の結果にannotateする。
# annotation_mod(directory,sample)

# [8] RTstop & misincorporation pattern at modification sites
# detectable_modifications(directory,sample,"max",False)
# detectable_modifications(directory,sample,"ave",True)


# [8.5] look into m1A and m1G
def specific_modification(directory,sample):
    for cal in ["ave","max"]:
        df = pd.read_csv("{}Detectable_modifications_mismtach_{}RTstop.csv".format(directory,cal))
        # print(df.head())
        for m in ["m1A","m1G"]:
            print("modification:",m)
            dfm = df[df["Annotated modification"] == m].reset_index(drop=True)
            # print(dfm.head())
            fig = plt.figure(figsize=(7,2.7*len(set(dfm["canon_pos"]))))
            color = ["#D6FA8C","#A5D721","#82B300","#FDB777","#FD9346","#FF6200"]
            flierprops = dict(marker='o', markersize=1, alpha=0.5)
            b = 0
            for p in sorted(set(dfm["canon_pos"])):
                print("canonical position:",p)
                dfp = dfm[dfm["canon_pos"] == p].reset_index(drop=True)
                print(dfp.head())
                for t in sorted(set(dfp["prop_type"])):
                    dft = dfp[dfp["prop_type"] == t]
                    b += 1
                    ax = fig.add_subplot(len(set(dfm["canon_pos"])),2, b)
                    sns.boxplot(x="condition",y="proportion",data=dft,ax=ax,palette=color,order=sample, flierprops=flierprops)
                    # sns.stripplot(x="condition",y="proportion",data=dft,ax=ax,palette=color,order=sample,s=3)
                    # plt.plot([-0.5,5.5],[0.1,0.1],linestyle="dashed",color="gray",linewidth=0.5)
                    for l in range(len(sample)-1):
                        # print(sample[l],"vs",sample[l+1])
                        signal = test(dft[dft["condition"]==sample[l]]["proportion"],dft[dft["condition"]==sample[l+1]]["proportion"],sample[l],sample[l+1])
                        if signal != "":
                            plt.text(l+0.5,round(max(dft["proportion"]),1)+0.1,signal,ha="center",fontsize=12)
                            plt.plot([l+0.1,l+0.9],[round(max(dft["proportion"]),1)+0.1,round(max(dft["proportion"]),1)+0.1],c="black")
                    ax.set_ylabel(t)
                    ax.set_title("{}{}x{},{}".format(m,p,len(set(dft["isodecoder"])),t))
                    ax.set_ylim(-0.1,1.4)
                    # plt.ylim(-0.1,1.4)
                plt.tight_layout()
            plt.savefig("{}figures/{}_position_dependent_RTstop-mismatch_{}.png".format(directory,m,cal))
                            
# specific_modification(directory,sample)

# [9] Readthrough comparison between two samples
# compare_readthrough(directory,["K562_totalRNA","K562_tRNA_DN"])

# [10] heatmap of misincorporation and RTStop
# RTstop_mismatch_for_heatmap(directory,sample,"isoacceptor")
# RTstop_mismatch_for_heatmap(directory,sample,"isodecoder")
# heatmap_RTsop_mismatch(directory,["K562_totalRNA"])
# heatmap_RTsop_mismatch(directory,["mouse_mut"])

# [11] Annotation of known modifications & misincorporated bases
# annotation_mod_misincorporated_base(directory,sample)
# plot_base_preference(directory,sample,"canon_pos")
# plot_upstream_position_preference(directory,sample,[4,2],["m1A","m1G"])
# plot_base_preference(directory,sample,"downstream")

# [12] plot individual tRNA
# RTstop_mismatch_individual_isoacceptor(directory,sample,"mito_tRNA-Leu-TAA","9","m1G")
# RTstop_mismatch_individual_isoacceptor(directory,sample,"tRNA-Pro-TGG","37","m1G")
# RTstop_mismatch_individual_isoacceptor(directory,sample,"tRNA-Pro-CGG","37","m1G")
# RTstop_mismatch_individual_isoacceptor(directory,sample,"tRNA-Pro-AGG","34","I")

# これは、使っているheatmap用のデータがisoacceptor毎にまとまってしまっているから、isodecoder用には別の関数を準備する必要がある。
# [12] plot individual tRNA isodecoder
def RTstop_mismatch_individual_isodecoder(directory,sample,tRNA,canon_pos,mod,special):
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
            if special == True and proptype == "mismatch":
                # print(df_.head())
                # print(df_.index)
                df_.loc[canon_pos,"proportion"] = 1 - df_.loc[canon_pos,"proportion"]
            df_["x"] = [n for n in range(len(df_))]
            print(df_.head())
            if proptype == "mismatch" and canon_pos != False:
                x_ = df_.loc[canon_pos,"x"]
                y_ = df_.loc[canon_pos,"proportion"]
                # print(df_[df_.index==canon_pos])
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
        # ax.set_xticks([9,36,39,69,77])
        # ax.set_xticklabels(["9","34","37","50","58"])
        ax.set_xticks([9,28,34,69,77])
        ax.set_xticklabels(["9","26","32","50","58"])
        ax.set_xlabel("Canonical position",fontsize=10)
        ax.set_ylabel("Proportion",fontsize=10)
        ax.set_title(f"{tRNA} in {s}",fontsize=10)
        plt.legend(loc="lower center",bbox_to_anchor=(0.5,1.2))
    plt.tight_layout()
    fig.savefig("{}figures/RTstop_mismatch_in_{}".format(directory,tRNA))
# RTstop_mismatch_individual_isodecoder(directory,sample,"tRNA-Arg-TCT-4","50","U",True)


# For calibration
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
    # Fit linear regression via least squares with numpy.polyfit
    # It returns an slope (b) and intercept (a)
    # deg=1 means linear fit (i.e. polynomial of degree 1)
    b, a = np.polyfit(x, y, deg=1)
    # Create sequence of 100 numbers from 0 to 100 
    xseq = np.linspace(0, 1, num=100)
    # Plot regression line
    ax.plot(xseq, a + b * xseq, color="k", lw=1)
    # ax.plot([0,6],[0,6],linewidth=0.5,c='gray')
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

calibration(directory,sample,[0,0.25,0.75,1])
