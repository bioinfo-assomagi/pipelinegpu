import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Process, Lock
import argparse
import csv
import datetime
import glob
import os
import datetime
import subprocess
import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.nonparametric.kde import KDEUnivariate
from statsmodels.nonparametric.kernel_density import KDEMultivariate
from statsmodels.stats.diagnostic import kstest_normal
from statsmodels.nonparametric.kernel_regression import KernelReg
import scipy.stats as st
from scipy.interpolate import interp1d
from os import listdir, system
from os.path import isfile, join
import shutil
import time
# from cStringIO import StringIO
from io import StringIO
pd.options.mode.chained_assignment = None
import time, warnings
from sklearn_pandas import DataFrameMapper
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from Pipes.Pipe import Pipe
from Pipes.ParallelPipe import ParallelPipe

import config
import dir_tree


warnings.filterwarnings("ignore",category =RuntimeWarning)
warnings.filterwarnings("ignore",category =UserWarning)

np.errstate(invalid='ignore')

geno37 = join('/home/magi/','dataset/GENOME/37/all_chr37.fa')
geno38 = join('/home/magi/','dataset/GENOME/38/all_chr38.fa')

HGMD_37 = join('/home/magi/','dataset/HGMD/37/HGMD_pro_2014_2.bed')
HGMD_38 = join('/home/magi/','dataset/HGMD/38/LiftOverHGMD_pro_2014_2.bed')



class VariantCallPipeIndel(ParallelPipe):

    BUCHIARTIFICIALI = join('/home/magi/', 'PROJECT/diagnosys/bin/BUCHIARTIFICIALI.txt')

    def __init__(self):
        super().__init__()

    def process(self, **kwargs):
        self.thread_id = os.getpid()

        self.sample = kwargs.pop("sample")
        self.panel = kwargs.pop("panel", None)
        self.genome_type = kwargs.pop("genome", "geno38")
        self.dest = kwargs.pop("dest")

        self.thread_print("TESTING VARIANT Calling INDEL PIPE: {}".format(self.sample.name))

        self.run()

        kwargs.update({"sample": self.sample, "genome": self.genome_type, "panel": self.panel, "dest": self.dest})
        return kwargs


      

    def ClusterJustControl(self, param,sample_x,PCA_folder,CONTROL_folder,PCA_NEW,bedfile, folder_bam):
        #CREATE CONTROLS DATA
        # -mode StartWithBam
        # -inputDir /home/magi/PROJECT/diagnosys/RESULT/21_Jul_2022_CANCER/bam/
        # -useSampleAsControl
        # -controlsDir /home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/PCA
        # -outputDir /home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/NEWENTRY
        # -bed /home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/BED/bed_CANCER_TWIST91957849.bed
        bashCommand0=(' '.join(['CoNVaDING.pl -mode StartWithBam','-inputDir',join(folder_bam,'inanalysis', sample_x),
                        '-useSampleAsControl',
                        '-controlsDir',PCA_folder,
                        '-outputDir',PCA_NEW,
                        '-bed',bedfile]))
        process = subprocess.Popen(bashCommand0.split(), stdout=subprocess.PIPE)
        output, error = process.communicate() #DEBUG_ONLY
        return True

    def ClusterSetting(self, param,sample_x,PCA_folder,CONTROL_folder,folder_indel,bedfile):
        print ('PCA ANALYSIS!!!\n')

        files=glob.glob(join(PCA_folder,'*_final.normalized.coverage.txt'))
        DATA=pd.DataFrame()
        #READ NORMALIZED CONTROLS
        for file in files:
            t_=pd.read_csv(file, sep='\t')
            t_['sample']=file.split('/')[-1].split('.n')[0]
            t_['coordinates']=t_['CHR']+':'+t_['START'].map(str)+'-'+t_['STOP'].map(str)
            t_.set_index('coordinates',inplace=True)
            t_=t_[['NORMALIZED_TOTAL']]
            t_=t_.T
            t_.index = [file.split('/')[-1].split('.n')[0]]
            DATA=pd.concat([DATA,t_])

        path=join(join(folder_indel,sample_x),sample_x+'_final.normalized.coverage.txt')

        t_=pd.read_csv(path, sep='\t')
        t_['sample']=sample_x		#path.split('/')[-1].split('.n')[0]
        t_['coordinates']=t_['CHR']+':'+t_['START'].map(str)+'-'+t_['STOP'].map(str)
        t_.set_index('coordinates',inplace=True)
        t_=t_[['NORMALIZED_TOTAL']]
        t_=t_.T
        t_.index = [path.split('/')[-1].split('.n')[0]]
        DATA=pd.concat([DATA,t_])
        #################SCALING DATA###########
        mapper = DataFrameMapper([(DATA.columns, StandardScaler())])
        scaled_features = mapper.fit_transform(DATA.copy())
        scaled_features_df = pd.DataFrame(scaled_features, index=DATA.index, columns=DATA.columns)
        #################PCA####################
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(scaled_features_df)
        principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'], index=scaled_features_df.index)
        cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
        scaled_features_df['cluster']=cluster.fit_predict(principalDf)
        principalDf['cluster']=list(scaled_features_df['cluster'])
        clust_n=scaled_features_df.loc[sample_x+'_final','cluster']
        #GET CONTROLS CORRESPONDING TO SAMPLE CLUSTERS
        files.append(path)
        principalDf['paths']=files

        try: controls_to_load=list(principalDf[principalDf['cluster']==clust_n]['paths'])
        except: controls_to_load=list(principalDf[principalDf['cluster']==clust_n[0]]['paths'])

        del controls_to_load[-1]
        for file in controls_to_load:
            shutil.copy(file,CONTROL_folder)

        return
    

    def ConvadingCore(self, param,sample_x,CONTROL_folder,folder_bam,folder_indel,bedfile):
        #NORMALIZE BAM COUNTS
        print ('NORMALIZE BAM COUNTS!!!\n')
        bashCommand1=(' '.join(['perl','/home/magi/PROJECT/diagnosys/bin/CoNVaDING.pl -mode StartWithBam','-inputDir',join(folder_bam,'inanalysis', sample_x),
                            '-controlsDir',CONTROL_folder,'-outputDir',join(folder_indel,sample_x),'-bed',bedfile]))
        process = subprocess.Popen(bashCommand1.split(), stdout=subprocess.PIPE)
        output, error = process.communicate() #DEBUG_ONLY

        #CHOOSE CONTROLS
        print( 'CHOOSE CONTROLS!!!\n')
        bashCommand2=(' '.join(['perl','/home/magi/PROJECT/diagnosys/bin/CoNVaDING.pl -mode StartWithMatchScore ','-inputDir ',join(folder_indel,sample_x),
                            '-sexChr','-controlSamples 30','-controlsDir',CONTROL_folder,'-outputDir',join(folder_indel,sample_x)]))
        process = subprocess.Popen(bashCommand2.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        #ANALYZE CNV
        print ('ANALYZE CNV!!!\n')
        # NOTE: if the control group is not good, the cutoff of 3.5 and -2.6 are used
        # 
        # bashCommand3=(' '.join(['CoNVaDING.pl -mode StartWithBestScore ','-inputDir',join(folder_indel,sample_x),'-controlsDir ',CONTROL_folder,
        # 						'-sexChr','-zScoreCutOffLow=-3.2','-zScoreCutOffHigh 3.2 ','-outputDir',join(folder_indel,sample_x)]))
        bashCommand3=(' '.join(['perl','/home/magi/PROJECT/diagnosys/bin/CoNVaDING.pl -mode StartWithBestScore ','-inputDir',join(folder_indel,sample_x),'-controlsDir ',CONTROL_folder,
                            '-sexChr','-zScoreCutOffLow=-2.6','-zScoreCutOffHigh=3.5','-outputDir',join(folder_indel,sample_x)]))

        process = subprocess.Popen(bashCommand3.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return
    

    def IntegrateResult(self, param,sample_x,TOTALbedfile,lastFOLDER,phenotype, folder_indel):
        result_totalist = pd.DataFrame(columns=['CHR','START','STOP','GENE','NUMBER_OF_TARGETS','NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST','ABBERATION'])
        calls = pd.DataFrame(columns=['CHR','START','STOP','GENE','NUMBER_OF_TARGETS','NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST','ABBERATION'])
        phenotype_ = phenotype[phenotype['sample'].astype(str) == str(sample_x)]
        a = phenotype_[['malattia','gene']].drop_duplicates()
        x1 = pd.DataFrame(a['gene'])
        #x2= pd.DataFrame({'gene':pd.Series(['CTD-3074O7.11','TM4SF2'])})
        x = x1 #.append(x2)

        if self.dest == 'r':
            dest = 'rovereto'
            path_django55 = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/INDEL/'
            path_django = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/INDEL/'
            path_download55 = '/home/magi/VIRTUAL/MAGIS/DOWNLOADS/NGSINFO/Indel'
        elif self.dest == 'b':
            dest = 'bolzano'
            path_django55 = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/INDEL/'
            path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/INDEL/'
            path_download55 = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/Indel'
        elif self.dest == 's':
            dest = 'sanfelice'
            path_django55 = '/home/magi/VIRTUAL/SANFELICE/NGS_RESULT/INDEL/'
            path_download55 = '/home/magi/VIRTUAL/SANFELICE/DOWNLOADS/NGSINFO/Indel'
        elif self.dest == 'z':
            dest = 'ricerca'
            path_django55 = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/INDEL/'
            path_django = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/INDEL/'
            path_download55 = '/home/magi/VIRTUAL/RICERCA/DOWNLOADS/NGSINFO/Indel'
        elif self.dest == 'p':
            dest = 'privato'
            path_django55 = '/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/'
            path_download55 = '/home/magi/VIRTUAL38/apimagi_prod/DOWNLOADS/NGSINFO/Indel'
            path_django133 = 'bioinfo@192.168.1.133:/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/'
            path_download133 = 'bioinfo@192.168.1.133:/home/magi/VIRTUAL38/apimagi_prod/DOWNLOADS/NGSINFO/Indel'
            path_django = 'bioinfo@192.168.1.120:/home/magi/VIRTUAL/SANFELICE/NGS_RESULT/INDEL/'
            path_download = 'bioinfo@192.168.1.120:/home/magi/VIRTUAL/SANFELICE/DOWNLOADS/NGSINFO/Indel'
        if not os.path.exists(path_download55):
            os.makedirs(path_download55)
        if not os.path.exists(path_django55):
            os.makedirs(path_django55)

        info_all_sample = join(folder_indel,'INDELs_info.csv')
        info_all_sample_download = join(path_download55,'INDELs_%s') % (str(lastFOLDER))

        #INTEGRATE RESULTS WITH INFO FROM NORM COVERAGE AND BED FILE
        full_bed= TOTALbedfile    #pd.read_csv(sys.argv[3], sep='\t')
        coverage=pd.read_csv(join(folder_indel,sample_x,sample_x+'_final.normalized.coverage.txt'), sep='\t',header=0)
        callsshort=pd.read_csv(join(folder_indel,sample_x,sample_x+'_final.best.score.shortlist.txt'), sep='\t',header=0)
        callslong=pd.read_csv(join(folder_indel,sample_x,sample_x+'_final.best.score.longlist.txt'), sep='\t',header=0)
        calls2=pd.read_csv(join(folder_indel,sample_x,sample_x+'_final.best.score.totallist.txt'), sep='\t',header=0)
        calls_full={}

        calls2b=calls2[calls2['ABBERATION']!='.']
        clong = callslong[callslong['NUMBER_OF_TARGETS']>=3]
        clong = pd.merge(clong[['CHR','GENE','NUMBER_OF_TARGETS', 'NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST']], calls2b[['CHR','START','STOP','GENE','AUTO_RATIO','AUTO_ZSCORE','ABBERATION']], on=['CHR','GENE'], how='left')

        if len(clong) == 0:
            calls1 = callsshort
            calls1['AUTO_RATIO'] = -999
            calls1['AUTO_ZSCORE'] = -999
        else:
            calls1=pd.concat([callsshort,clong])

        calls1.drop_duplicates(inplace=True)
        result_homdel = calls2[calls2['ABBERATION']=='HOM_DEL']
        result_homdup = calls2[calls2['ABBERATION']=='HOM_DUP']
        if ((len(result_homdel)==0) & (len(result_homdup)!=0)): result_totalist = result_homdup
        elif ((len(result_homdel)!=0) & (len(result_homdup)==0)): result_totalist = result_homdel
        elif ((len(result_homdel)==0) & (len(result_homdup)==0)): result_totalist = pd.DataFrame(columns=['CHR','START','STOP','GENE',
                                                'NUMBER_OF_TARGETS','NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST','ABBERATION','AUTO_RATIO','AUTO_ZSCORE'])
        else: result_totalist = pd.concat([result_homdel,result_homdup])
        result_totalist['NUMBER_OF_TARGETS'] = 1
        result_totalist['NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST'] = 1
        result_totalist = result_totalist[['CHR','START','STOP','GENE','NUMBER_OF_TARGETS','NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST','ABBERATION','AUTO_RATIO','AUTO_ZSCORE']]
        if ((len(result_totalist)==0) & (len(calls1)!=0)): calls = calls1
        elif ((len(result_totalist)!=0) & (len(calls1)==0)): calls = result_totalist
        elif ((len(result_totalist)==0) & (len(calls1)==0)): calls = pd.DataFrame(columns=['CHR','START','STOP','GENE','NUMBER_OF_TARGETS','NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST','ABBERATION','AUTO_RATIO','AUTO_ZSCORE'])
        else: calls = pd.concat([result_totalist,calls1])

        if len(calls) == 0:
            result4=pd.DataFrame(columns=['sample_id','hgvs','chrom','GENE','exone','strand','refseq', 'ABBERATION',
                        'NUMBER_OF_TARGETS','REGION_COV','AVG_TOTAL_COV','AVG_GENE_COV',
                        'NORMALIZED_TOTAL','NORMALIZED_GENE','AUTO_RATIO','AUTO_ZSCORE'])
            result4['sample_id'] = sample_x
            result4['hgvs'] = -999
        else:
            j=0
            for index, row in calls.iterrows():
                strand=full_bed[full_bed['GENE']==row.GENE]['strand'].unique()
                if strand[0]==1:
                    temp=full_bed[(full_bed['START']>=row.START)&(full_bed['#CHROM']==(row.CHR))]
                    temp=temp[(temp['END']<=row.STOP)&(full_bed['#CHROM']==row.CHR)]
                else:
                    temp=full_bed[(full_bed['START']<=row.START)&(full_bed['#CHROM']==(row.CHR))]
                    temp=temp[(temp['END']>=row.STOP)&(full_bed['#CHROM']==(row.CHR))]

                for i,r in temp.iterrows():
                    final_r=pd.concat([row.drop(['START','GENE']),r])
                    calls_full[i]=final_r
                    j=j+1
            result=pd.DataFrame(calls_full).T
            if len(result) != 0:
                result['STOP']=result['END']
                result.drop(['END'],axis=1,inplace=True)
                result['sample_id']=sample_x
                result['chrom'] = result['#CHROM']
                result['hgvs']=result['#CHROM'].astype(str)+':'+result['START'].astype(str)+'-'+result['STOP'].astype(str)
                result['START'] = result['START'].astype(int)
                result['STOP'] = result['STOP'].astype(int)
                result=result[['sample_id','hgvs','chrom','START','STOP','GENE','exone','strand','refseq', 'ABBERATION', 'NUMBER_OF_TARGETS','NUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST']]
                result=result.reset_index(drop=True)
                result2 = pd.merge(result, coverage, on=['START','STOP', 'GENE'], how='left')
                result2['NUMBER_OF_TARGETS'] = result2['NUMBER_OF_TARGETS'].astype('int')
                result2['CHR'] = result2['chrom']
                result3b = pd.merge(result2,calls2[['CHR','START','GENE','AUTO_RATIO','AUTO_ZSCORE']],on = ['CHR','START','GENE'], how='left')
                result3bnew=result3b[abs(result3b['AUTO_ZSCORE'])>=2.0]
                result3=result3b[['sample_id','hgvs','chrom','GENE','exone','strand','refseq', 'ABBERATION',
                            'NUMBER_OF_TARGETS','REGION_COV','AVG_TOTAL_COV','AVG_GENE_COV',
                            'NORMALIZED_TOTAL','NORMALIZED_GENE','AUTO_RATIO','AUTO_ZSCORE']][((result3b['NUMBER_OF_TARGETS'].astype(int)>1) | (result3b['ABBERATION'] == 'HOM_DUP') | (result3b['ABBERATION'] == 'HOM_DEL'))]
                ##################################################
                result4B = result3[result3['GENE'].isin(x['gene'])]
                result4C = result3[result3['ABBERATION']=='HOM_DEL']
                result4D = result4B[(result4B['chrom']!='chrX')&(result4B['ABBERATION']!='HOM_DEL')]
                result4 = result4D[result4D['GENE'].isin(x['gene'])]
            #######################################################
            else:
                result4=pd.DataFrame(columns=['sample_id','hgvs','chrom','GENE','exone','strand','refseq', 'ABBERATION',
                            'NUMBER_OF_TARGETS','REGION_COV','AVG_TOTAL_COV','AVG_GENE_COV',
                            'NORMALIZED_TOTAL','NORMALIZED_GENE','AUTO_RATIO','AUTO_ZSCORE'])
                result4['sample_id'] = sample_x
                result4['hgvs'] = -999
        sizegene = pd.DataFrame(result4.groupby(by=['GENE']).size().reset_index())
        sizegene['size'] = sizegene[0]
        sizegene.drop(0,axis=1,inplace=True)
        sizegene = sizegene[sizegene['size']>=3]
        result5A = result4[result4['GENE'].isin(sizegene['GENE'])]
        try: result5=pd.concat([result5A,result4C]).drop_duplicates()
        except:result5=result5A.drop_duplicates()
        result5 = result5[result5['GENE'].isin(x['gene'])]
        #result5.to_csv(join(path_django,sample_x+'_prefinal_indel.csv'),sep='\t',index=None)
        result5.to_csv(join(folder_indel,sample_x,sample_x+'_prefinal_indel.csv'),sep='\t',index=None)
        #system(' '.join(['scp',join(folder_indel,sample_x,sample_x+'_final_indel.csv'), "bioinfo@192.168.1.133:/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/"]))
        #system(' '.join(['cp',join(folder_indel,sample_x,sample_x+'_final_indel.csv'), "/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/"]))
        #print ('---->'+str(sample_x)+' ('+str(len(result5))+')'+'<----')
        cl = open(info_all_sample,'a')
        cl.write(sample_x)
        cl.write('\t')
        cl.write(str(len(result5)))
        cl.write('\n')
        if len(result5) >= 1: cl.write(str(result5[['sample_id','GENE','exone','ABBERATION','NUMBER_OF_TARGETS']]))
        cl.write('\n--------------------------------\n')
        cl.close()
        cl_download = open(info_all_sample_download,'a')
        cl_download.write(sample_x)
        cl_download.write('\t')
        cl_download.write(str(len(result5)))
        cl_download.write('\n')
        if len(result5) >= 1: cl_download.write(str(result5[['sample_id','GENE','exone','ABBERATION','NUMBER_OF_TARGETS']]))
        cl_download.write('\n--------------------------------\n')
        cl_download.close()
        return result5
    

    def cov_graph(self, param,result,sample_x,folder_indel,CONTROL_folder,PCA_folder):
	    #print( 'GRAPHING!!!!! ',sample_x)
        if len(result) == 0:
            pass
        else:
            try:
                result['START'] = result['hgvs'].str.split(':').str.get(1).str.split('-').str.get(0)
                result = result.reset_index()
                control_fin=pd.DataFrame()
                controls_files=glob.glob(CONTROL_folder+'/*.txt')
                if len(controls_files) == 0:
                    controls_files=glob.glob(PCA_folder+'/*.txt')
                for file in controls_files:
                    temp=pd.read_csv(file,sep='\t')
                    control_fin=pd.concat([control_fin,temp])
                l=[]
                ticks=[]
                genes=[]
                fig1, ax1 = plt.subplots(figsize=(15,10))
                colors=['g','b','c','y','silver','orange','salmon','khaki','chocolate']
                k=0
                j=0
                for i in range(0,len(result)):
                    l.append(control_fin[control_fin['START'].astype(int)== int(result.loc[i,'START'])]['NORMALIZED_TOTAL'])
                    ticks.append(result.loc[i,'GENE']+'_'+str(result.loc[i,'exone']))
                    genes.append(result.loc[i,'GENE'])
                for i in range(0,len(result)):
                    ax1.boxplot(l, whis=1.5)
                    y = result.loc[i,'NORMALIZED_TOTAL']
                    x = i+1
                    if genes[i]==genes[k]:
                        plt.axvspan(i+0.5, i+1.5, facecolor=colors[j], alpha=0.5)
                    else:
                        k=i
                        if j < len(genes): j=j+1
                        else: j=0
                        plt.axvspan(i+0.5, i+1.5, facecolor=colors[j], alpha=0.5)
                    plt.plot(x, y, 'r.', alpha=1, ms=10, marker='x')
                plt.xticks(range(1,len(result)+1), ticks)

                folder_name = join(folder_indel,sample_x,'figure')
                if self.dest == 'r':
                    dest = 'rovereto'
                    folder_VIRTUAL = join('/home/magi/VIRTUAL/MAGIS/MAGIS/static/INDEL_image',sample_x)
                    if not os.path.exists(folder_VIRTUAL):
                        os.makedirs(folder_VIRTUAL)
                elif self.dest == 'b':
                    dest = 'bolzano'
                    folder_VIRTUAL = join('/home/magi/VIRTUAL/EUREGIO/EUREGIO/static/INDEL_image',sample_x)
                    if not os.path.exists(folder_VIRTUAL):
                        os.makedirs(folder_VIRTUAL)
                elif self.dest == 's':
                    dest = 'sanfelice'
                    folder_VIRTUAL = join('/home/magi/VIRTUAL/SANFELICE/SANFELICE/static/INDEL_image',sample_x)
                    if not os.path.exists(folder_VIRTUAL):
                        os.makedirs(folder_VIRTUAL)
                elif self.dest == 'z':
                    dest = 'ricerca'
                    folder_VIRTUAL = join('/home/magi/VIRTUAL/RICERCA/RICERCA/static/INDEL_image',sample_x)
                    if not os.path.exists(folder_VIRTUAL):
                        os.makedirs(folder_VIRTUAL)
                elif self.dest == 'p':
                    dest = 'privato'
                    folder_VIRTUAL  = join('/home/magi/VIRTUAL/SANFELICE/SANFELICE/static/INDEL_image',sample_x)
                    folder_VIRTUAL133 = join('bioinfo@192.168.1.133:/home/magi/VIRTUAL/SANFELICE/SANFELICE/static/INDEL_image',sample_x)
                    folder_VIRTUAL220 = join('bioinfo@192.168.1.120:/home/magi/VIRTUAL/SANFELICE/SANFELICE/static/INDEL_image',sample_x)
                    if not os.path.exists(folder_VIRTUAL):
                        os.makedirs(folder_VIRTUAL)
                if not os.path.exists(folder_name):
                    os.makedirs(folder_name)

                #name2 = name.replace(' ','')
                fig_name = sample_x+'_'+'cnv'+'.jpg'
                plt.savefig(join(folder_name,fig_name))
                #plt.savefig(main_folder+'results/'+sample_x+f'/{sample_x}_results.png')
                plt.savefig(join(folder_VIRTUAL,fig_name))
                plt.close(fig1)
                time.sleep(1)
            except: pass
        return
    

    def run(self):
        

        if self.panel == "trusightone":
            panel = self.panel
            project = self.panel

        else:

            panel = self.panel.upper()
            project = self.panel.upper()
            
            if 'OCULAR' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/OCULARE/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/OCULARE/BED/bed_OCULARE_TWIST91957849.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/OCULARE/BED/TOTAL_OCULARE_2022.bed',sep='\t',header=0)
            elif 'CANCER' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/BED/bed_CANCER_TWIST91957849.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/BED/TOTAL_CANCER_2022.bed',sep='\t',header=0)
            elif 'VASCULAR' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/VASCULAR/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/VASCULAR/BED/bed_VASCULAR.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/VASCULAR/BED/TOTAL_VASCULAR_2020.bed',sep='\t',header=0)
            elif 'NEUROLOGY' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/NEUROLOGY/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/NEUROLOGY/BED/bed_NEUROLOGY.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/NEUROLOGY/BED/TOTAL_NEUROLOGY_2020.bed',sep='\t',header=0)
            elif 'MIXED' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/MIXED1/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/MIXED1/BED/bed_MIXED1.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/MIXED1/BED/TOTAL_MIXED1_2020.bed',sep='\t',header=0)
            elif 'LYMPHOBESITY' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/LYMPHOBESITY/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/LYMPHOBESITY/BED/bed_LYMPHOBESITY.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/LYMPHOBESITY/BED/TOTAL_LYMPHOBESITY_2020.bed',sep='\t',header=0)
            elif 'INFERTILIT' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INFERTILITY/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INFERTILITY/BED/bed_INFERTILITY.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INFERTILITY/BED/TOTAL_INFERTILITY_2020.bed',sep='\t',header=0)
            elif 'INTEGRACARDIOSTANCHEZZA' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INTEGRACARDIOSTANCHEZZA/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INTEGRACARDIOSTANCHEZZA/BED/bed_INTEGRACARDIOSTANCHEZZA.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INTEGRACARDIOSTANCHEZZA/BED/TOTAL_INTEGRACARDIOSTANCHEZZA.bed',sep='\t',header=0)
            elif 'GENEOB' in self.panel:
                STAMPO = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/GENEOB/'
                bedfile = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/GENEOB/BED/bed_GENEOB.bed'
                TOTALbedfile = pd.read_csv('/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/GENEOB/BED/TOTAL_GENEOB.bed',sep='\t',header=0)
            else: pass


        folder_name = dir_tree.principal_directory.path
        folder_indel = dir_tree.principal_directory.indel.path
        print("Folder_indel = {}".format(folder_indel))
        folder_coverage = dir_tree.principal_directory.coverage.path
        folder_pheno = dir_tree.principal_directory.pheno.path
        folder_bam = dir_tree.principal_directory.bam.path
        print("Folder_bam = {}".format(folder_bam))
        input_phenotype = join(folder_name,'pheno/phenotype')
        lastFOLDER = folder_name.split('/')[-1]

        if self.genome_type != 'geno38':
            return

        hgmd = pd.read_csv(config.HGMD38,sep='\t',header=None,names=['#CHROM','START','END','hgmd'], encoding='ISO-8859-9')

        folder_COV = dir_tree.principal_directory.coverage.path
        sample_name = str(self.sample.name)
        coverage = pd.read_csv(os.path.join(folder_COV, sample_name, sample_name + "_all"), sep='\t', header=0) # TODO; read it from self.sample.coverage_all
        _stat_cov = join(folder_COV, sample_name, sample_name + '_stat_cov.csv')
        stat_cov = pd.read_csv(_stat_cov, sep='\t')

        sample_x = str(self.sample.name)

        if stat_cov[stat_cov['cutoff']=='depth_macro>=25'].iloc[0,3] >= 92:
        #if True:

            # Not an issue for parallelization, since each sample corresponds to one (his own) directory
            if not os.path.exists(join(folder_indel, sample_x)):
                os.makedirs(join(folder_indel, sample_x))
            else:
                print ('Folder already exists! skipping sample analysis for: ', sample_x)

            bam_path = self.sample.bam
            bai_path = self.sample.bai

            BAM = str(sample_x+'_final.bam')
            BAI = str(sample_x+'_final.bam.bai')

            CONTROL_folder = join(STAMPO,'CONTROLS')
            PCA_folder = join(STAMPO,'PCA')
            phenotype = pd.read_csv(input_phenotype, sep='\t', header=0,encoding='utf-8')
            phenotype['sample'] = phenotype['sample'].astype('str')
            phenotype['sample'].replace(r'\.202$','.2020',inplace=True,regex=True)
            phenotype['sample'] = phenotype['sample'].astype('str')
            self.thread_print(CONTROL_folder + ":" + bam_path)

            if not os.path.exists(join(folder_bam, 'inanalysis', sample_x)):
                os.makedirs(join(folder_bam, 'inanalysis', sample_x))

            try:
                shutil.move(bam_path, join(folder_bam,'inanalysis', sample_x))
                shutil.move(bai_path, join(folder_bam,'inanalysis', sample_x))
            except Exception as e: self.thread_print ('errore moving bam files to inanalysis: ' + str(e))

            self.thread_print('JUST REAL!!!\n')
            ##################################################################
            self.ConvadingCore(None, sample_x,PCA_folder,folder_bam,folder_indel,bedfile)
            # # # ##################################################################
            result = self.IntegrateResult(None, sample_x,TOTALbedfile,lastFOLDER,phenotype, folder_indel)
            # # # # # ##################################################################
            self.cov_graph(None,result,sample_x,folder_indel,CONTROL_folder,PCA_folder)

            if ('OCULARE' in project):
                #1################################################################
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/OCULARE/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/OCULARE/BED/bed_OCULARE_TWIST91957849.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, _CONTROL_folder_, PCA_NEW, _bedfile_, folder_bam)
                self.thread_print('JUST CONTROL OCULARE2!!!\n')
            if ('CANCER' in project): #|('OCULARE' in project)):
                ##1################################################################
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/CANCER/BED/bed_CANCER_TWIST91957849.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, _CONTROL_folder_, PCA_NEW, _bedfile_, folder_bam)
                self.thread_print('JUST CONTROL CANCER!!!\n')
                #2################################################################
                # _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/OCULARE/'
                # _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/OCULARE/BED/bed_OCULARE_TWIST91957849.bed'
                # _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                # PCA_folder = join(_STAMPO_,'PCA')
                # PCA_NEW = join(_STAMPO_,'NEWENTRY')
                # ClusterJustControl(args,sample_x,PCA_folder,_CONTROL_folder_,PCA_NEW,_bedfile_, folder_bam)
                # print ('JUST CONTROL OCULARE!!!\n')

            if ('INTEGRACARDIOSTANCHEZZA' in project):
                #1################################################################
                self.thread_print('JUST CONTROL INTEGRACARDIOSTANCHEZZA!!!\n')
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INTEGRACARDIOSTANCHEZZA/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INTEGRACARDIOSTANCHEZZA/BED/bed_INTEGRACARDIOSTANCHEZZA.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, CONTROL_folder, PCA_NEW, _bedfile_, folder_bam)
            elif ('GENEOB' in project):
                #1################################################################
                self.thread_print('JUST CONTROL GENEOB!!!\n')
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/GENEOB/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/GENEOB/BED/bed_GENEOB.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x,PCA_folder,CONTROL_folder,PCA_NEW,_bedfile_, folder_bam)
            elif (('MIXED1' in project) | ('INFERTILIT' in project)):
                #1################################################################
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/MIXED1/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/MIXED1/BED/bed_MIXED1.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, CONTROL_folder, PCA_NEW, _bedfile_, folder_bam)
                self.thread_print('JUST CONTROL MIXED1!!!\n')
                #2################################################################
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INFERTILITY/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/INFERTILITY/BED/bed_INFERTILITY.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, _CONTROL_folder_, PCA_NEW, _bedfile_, folder_bam)
                self.thread_print('JUST CONTROL INFERTILITY!!!\n')
            elif (('LYMPHOBESITY' in project)|('NEUROLOGY' in project)|('VASCULAR' in project)):
                #1################################################################
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/LYMPHOBESITY/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/LYMPHOBESITY/BED/bed_LYMPHOBESITY.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, _CONTROL_folder_, PCA_NEW, _bedfile_, folder_bam)
                self.thread_print('JUST CONTROL LYMPHOBESITY!!!\n')
                #2###################################################################
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/NEUROLOGY/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/NEUROLOGY/BED/bed_NEUROLOGY.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, _CONTROL_folder_, PCA_NEW, _bedfile_, folder_bam)
                self.thread_print('JUST CONTROL NEUROLOGY!!!\n')
                #3###################################################################
                _STAMPO_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/VASCULAR/'
                _bedfile_ = '/home/magi/PROJECT/diagnosys/bin/COVSTAMP/convading/VASCULAR/BED/bed_VASCULAR.bed'
                _CONTROL_folder_ = join(_STAMPO_,'CONTROLS')
                PCA_folder = join(_STAMPO_,'PCA')
                PCA_NEW = join(_STAMPO_,'NEWENTRY')
                self.ClusterJustControl(None, sample_x, PCA_folder, _CONTROL_folder_, PCA_NEW, _bedfile_, folder_bam)
                self.thread_print('JUST CONTROL VASCULAR!!!\n')
            try:
                #print (' '.join(['rm',(CONTROL_folder+'/*.txt')]))
                #system(' '.join(['rm',(CONTROL_folder+'/*.txt')])) TODO: move this to end pipe
                shutil.move(join(folder_bam,'inanalysis',sample_x, BAM), folder_bam)
                shutil.move(join(folder_bam,'inanalysis',sample_x, BAI), folder_bam)
            except Exception as e: self.thread_print(e)
            self.thread_print('---------------------------------------------')

        else:
            # pass
            # Not an issue for parallelization, since each sample corresponds to one (his own) directory
            if not os.path.exists(join(folder_indel,sample_x)):
                os.makedirs(join(folder_indel,sample_x))
            else:
                print ('Folder already exist! skipping sample analysis for: ',sample_x)

            if self.dest == 'r':
                dest = 'rovereto'
                path_django = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/INDEL/'
                path_download = '/home/magi/VIRTUAL/MAGIS/DOWNLOADS/NGSINFO/Indel'
            elif self.dest == 'b':
                dest = 'bolzano'
                path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/INDEL/'
                path_download = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/Indel'
            elif self.dest == 's':
                dest = 'sanfelice'
                path_django = '/home/magi/VIRTUAL/SANFELICE/NGS_RESULT/INDEL/'
                path_download = '/home/magi/VIRTUAL/SANFELICE/DOWNLOADS/NGSINFO/Indel'
            elif self.dest == 'z':
                dest = 'ricerca'
                path_django = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/INDEL/'
                path_download = '/home/magi/VIRTUAL/RICERCA/DOWNLOADS/NGSINFO/Indel'

            result = pd.DataFrame(columns=['sample_id','hgvs','chrom','GENE','exone','strand','refseq', 'ABBERATION',
                        'NUMBER_OF_TARGETS','REGION_COV','AVG_TOTAL_COV','AVG_GENE_COV',
                        'NORMALIZED_TOTAL','NORMALIZED_GENE','AUTO_RATIO','AUTO_ZSCORE'])
            
            result.to_csv(join(folder_indel,sample_x,sample_x+'_prefinal_indel.csv'),sep='\t',index=None)
            result.to_csv(join(path_django,sample_x+'_final_indel.csv'),sep='\t',index=None)
            #system(' '.join(['scp',join(folder_indel,sample_x,sample_x+'_final_indel.csv'), "bioinfo@192.168.1.133:/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/"]))
            #system(' '.join(['cp',join(folder_indel,sample_x,sample_x+'_final_indel.csv'), "/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/"]))

            self.sample.prefinal_indel = join(folder_indel,sample_x,sample_x+'_prefinal_indel.csv')
            self.sample.saveJSON()

            info_all_sample = join(folder_indel,'INDELs_info.csv')
            info_all_sample_download = join(path_download,'INDELs_%s') % (str(lastFOLDER))

            cl = open(info_all_sample,'a')
            cl.write(sample_x)
            cl.write('\t')
            cl.write('L\'analisi di screening per le CNV non e\' stata eseguita perche\' meno del 92% delle reads hanno un coverage di almeno 25X.')
            cl.write('\n')
            cl.write('\n--------------------------------\n')
            cl.close()
            cl_download = open(info_all_sample_download,'a')
            cl_download.write(sample_x)
            cl_download.write('\t')
            cl_download.write('L\'analisi di screening per le CNV in eterozigosi non e\' stata eseguita perche\' meno del 92% delle reads hanno un coverage di almeno 25X.')
            cl_download.write('\n')
            cl_download.write('\n--------------------------------\n')
            cl_download.close()