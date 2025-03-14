from Pipes.Pipe import Pipe
import utils
import os
import pandas as pd
import numpy as np
import csv
import glob
import config
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None

from Pipes.ParallelPipe import ParallelPipe

geno37 = os.path.join('/home/magi/','dataset/GENOME/37/all_chr37.fa')
geno38 = os.path.join('/home/magi/','dataset/GENOME/38/all_chr38.fa')

HGMD_37 = os.path.join('/home/magi/','dataset/HGMD/37/HGMD_pro_2014_2.bed')
HGMD_38 = os.path.join('/home/magi/','dataset/HGMD/38/LiftOverHGMD_pro_2014_2.bed')

import dir_tree

#cols = ["#CHROM", "POS", "C%", "G%", "T%", "A%", "ins%", "del%", "sum", "DEPTH", "GENE", "exone", "length", "strand", "refseq", "hgmd", "sample"]

class CoverageStatisticsPipe(ParallelPipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        self.thread_id = os.getpid()

        self.panel = kwargs.pop("panel")
        self.sample = kwargs.pop("sample")
        # self.sample.temp_dir = os.path.join(dir_tree.principal_directory.temp.path, str(self.sample.name))
        self.dest = kwargs.pop("dest")
        self.genome_type = kwargs.pop("genome", "geno38")
        self.lock = kwargs.pop("lock")


        self.run()

        kwargs.update({"sample": self.sample, "panel": self.panel, "dest": self.dest, "genome": self.genome_type, "lock": self.lock})


    def cov(self, lock, param_dummy, coverage_folder, name, gene_size, all_genes_exons):
        if self.panel == 'trusightone':
            panel = self.panel
            project = self.panel
        elif self.panel == 'VASCULAR' :
                panel = self.panel.upper()
                project = 'ANOMALIE_VASCOLARI'
        elif self.panel == 'LYMPHATIC' :
                panel = self.panel.upper()
                project ='MALFORMAZIONI_LINFATICHE'
        elif self.panel == 'LIMPHOBESITY' :
                panel = self.panel.upper()
                project ='MALFORMAZIONI_LINFATICHE'
        elif self.panel == 'LYMPHOBESITY' :
                panel = self.panel.upper()
                project ='MALFORMAZIONI_LINFATICHE'
        elif self.panel == 'MIXED1' :
                panel = self.panel.upper()
                project ='GENETEST1'
        elif self.panel == 'MIXED2' :
                panel = self.panel.upper()
                project ='GENETEST2'
        elif self.panel == 'NEUROLOGY' :
                panel = self.panel.upper()
                project ='SIFSR'
        elif self.panel == 'INTEGRACARDIOSTANCHEZZA':
                panel = self.panel.upper()
                project ='SIFSR'
        elif self.panel == "INFERTILITA":
                panel = self.panel.upper()
                project = self.panel.upper()
        elif self.panel == "INFERTILITY":
                panel = 'INFERTILITA'
                project = 'INFERTILITA'
        elif self.panel == "OBESITA":
                panel = self.panel.upper()
                project = self.panel.upper()
        else:
            panel = self.panel.upper()
            project = self.panel.upper()


        sample_x = str(self.sample.name)
        coverage_df = pd.read_csv((name),sep='\t',header=0)
        coverage_df['POS'] = coverage_df['POS'].astype(int)
        coverage_df['control'] = coverage_df['GENE']+'_'+coverage_df['exone'].astype(str)
        i=0
        coverage_df.sort_values(by=['#CHROM','POS'],inplace=True)
        coverage_df = coverage_df.reset_index()
        coverage_df['variante'] = np.nan

        if len(coverage_df) != 0:
            diff = []
            for index,row in coverage_df.iterrows():
                try:
                    diff = 1
                    coverage_df.loc[0,'diff'] = diff
                    if coverage_df['control'][i] == coverage_df['control'][i+1]:
                        diff = coverage_df['POS'][i+1]-coverage_df['POS'][i]
                        if diff > 1:
                            coverage_df.loc[index+1,'diff'] = diff
                    elif coverage_df['control'][i] != coverage_df['control'][i+1]:
                        diff = 1
                        coverage_df.loc[index+1,'diff'] = diff
                except KeyError:
                    self.thread_print('CONTROLLO TERMINATO!!!')
                i+=1

        if len(coverage_df) != 0:
            for index,row in coverage_df.iterrows():
                try:
                    if ((row['G%'] != 1.0) & (row['G%'] != 0.0)):
                        coverage_df.loc[index,'variante'] = ('POS='+str(row['POS'])+' '+'G='+str(row['G%']))
                    if ((row['C%'] != 1.0) & (row['C%'] != 0.0)):
                        coverage_df.loc[index,'variante'] = 'POS='+str(row['POS'])+' '+'C='+str(row['C%'])
                    if ((row['T%'] != 1.0) & (row['T%'] != 0.0)):
                        coverage_df.loc[index,'variante'] = 'POS='+str(row['POS'])+' '+'T='+str(row['T%'])
                    if ((row['A%'] != 1.0) & (row['A%'] != 0.0)):
                        coverage_df.loc[index,'variante'] = 'POS='+str(row['POS'])+' '+'A='+str(row['A%'])
                    if ((row['ins%'] != 1.0) & (row['ins%'] != 0.0)):
                        coverage_df.loc[index,'variante'] = 'POS='+str(row['POS'])+' '+'ins='+str(row['ins%'])
                    if ((row['del%'] != 1.0) & (row['del%'] != 0.0)):
                        coverage_df.loc[index,'variante'] = 'POS='+str(row['POS'])+' '+'del='+str(row['del%'])
                except KeyError:
                    self.thread_print('TROVA VARIANTE TERMINATO!!!')

            coverage_df['diff'].fillna(method='ffill',inplace=True)
            group_cov = coverage_df.groupby(['GENE','exone','refseq','diff'])
            size = pd.DataFrame(group_cov.size()).reset_index()
            size['size'] = size[0]
            size.drop(0,axis=1,inplace=True)
            df = pd.DataFrame()
            for index, group in group_cov:
                mean = group['DEPTH'].mean()
                group['mean'] = '{:,.1f}'.format(mean)
                group['variante'].fillna(method='ffill',inplace=True)
                group['variante'].fillna(method='bfill',inplace=True)
                df = pd.concat([df,group])
            try:
                df = df[['GENE','exone','mean','diff','variante']]
                df['mean'].fillna(100,inplace=True)
                maxmin = pd.DataFrame(group_cov['POS'].agg([np.max, np.min])).reset_index()
                limit = pd.merge(maxmin,df,on=['GENE','exone','diff'],how='left')
            except KeyError:
                limit = pd.DataFrame(group_cov['POS'].agg([np.max, np.min])).reset_index()

            print("LIMIT DF = {}".format(limit))
            limit['buco_exone'] = ((limit['max']-limit['min'])+1) # NOTE: changed from amin to min and amax to max
            limit_ = pd.merge(limit,size,on = ['GENE','exone','refseq','diff'],how='left')

            limit_['note'] = np.where(limit_['buco_exone']!=limit_['size'],'multiple','unique')
            limit_['buco_exone'] = limit_['size']
            limit_.drop('size',axis=1,inplace=True)

            limit2 = pd.merge(limit_,coverage_df,on=['GENE','exone','refseq','diff'],how='left')
            limit2.drop_duplicates(subset=['GENE','exone','refseq','diff'],inplace=True)
            limit2_ = pd.merge(limit2,all_genes_exons,on=['GENE','exone','refseq','length'],how='left')
            limit2_.drop_duplicates(subset=['GENE','refseq','exone','length','diff'],inplace=True)

            limit2_['len_exon'] = (abs((limit2_['START'])-(limit2_['END'])))+1
            limit2_['%for_exon'] = (100-((limit2_['buco_exone'].astype(float))/(limit2_['len_exon'].astype(float))*100)).map('{:,.2f}'.format)
    ###################################################################################
            group2 = limit2_.groupby(['GENE'])
            gene_sum = pd.DataFrame(group2['buco_exone'].agg([np.sum])).reset_index()
            gene_sum['buco_gene'] = gene_sum['sum']
            gene_sum.drop('sum',axis=1,inplace=True)
    ####################################################################################
            df4 = pd.merge(limit2_,gene_sum,on=['GENE'],how='left')
            df5 = pd.merge(df4,gene_size,on=['GENE'],how='outer')
            df5.fillna(0,inplace=True)
            df5['%for_gene'] = (100-((df5['buco_gene'].astype(float))/(df5['len_gene'].astype(float))*100)).map('{:,.2f}'.format)
            df5['start_buco'] = df5['min'].astype(int)
            df5['end_buco'] = df5['max'].astype(int)
            df5.drop(['max','min'],axis=1,inplace=True)

            df5['hgmd'] = df5['hgmd_x']
            df5['strand'] = df5['strand_x']

            df5.fillna(0,inplace=True)
            df5['%for_exon'].replace(0,100,inplace=True)
            df5['%for_gene'].replace(0,100,inplace=True)
            df5['START'].replace(0,-999,inplace=True)
            df5['END'].replace(0,-999,inplace=True)
            df5['exone'].replace(0,-999,inplace=True)

            df5['sample'] = str(sample_x)

            df5['Exon_Start'] = df5['START']
            df5['Exon_End'] = df5['END']
            df5.drop(['START','END'],axis=1,inplace=True)

            df5['CHROM'] = df5['#CHROM_x']
            df5.drop_duplicates(subset=['GENE','exone','Exon_Start','Exon_End','diff'],inplace=True)
            df5['variante'] = df5['variante_x']
            try:
                df5 = df5[['sample','GENE','strand','refseq','hgmd','exone','CHROM','Exon_Start','Exon_End','start_buco','end_buco',
                'len_gene','len_exon','buco_exone','note','buco_gene','%for_exon','%for_gene','mean','variante']]
            except KeyError:
                df5['mean'] = 10
                df5 = df5[['sample','GENE','strand','refseq','hgmd','exone','CHROM','Exon_Start','Exon_End','start_buco','end_buco',
                'len_gene','len_exon','buco_exone','note','buco_gene','%for_exon','%for_gene','mean','variante']]

            df5['strand'] = df5['strand'].astype(int)
            df5['exone'] = df5['exone'].astype(int)
            df5['Exon_Start'] = df5['Exon_Start'].astype(int)
            df5['Exon_End'] = df5['Exon_End'].astype(int)
            df5['start_buco'] = df5['start_buco'].astype(int)
            df5['end_buco'] = df5['end_buco'].astype(int)
            df5['len_exon'] = df5['len_exon'].astype(int)
            df5['len_gene'] = df5['len_gene'].astype(int)
            df5['buco_exone'] = df5['buco_exone'].astype(int)
            df5['buco_gene'] = df5['buco_gene'].astype(int)

            folder_stat = os.path.join(coverage_folder, str(sample_x))
            try:
                os.makedirs(folder_stat)
            except FileExistsError:
                pass
            # if not os.path.exists(folder_stat):
            #     os.makedirs(folder_stat)
            if self.dest == 'r':
                dest = 'rovereto'
                path_django = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/coverage/'
                path_igv = os.path.join('/var/www/html/IGV',dest,project)
            elif self.dest == 'b':
                dest = 'bolzano'
                path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/coverage/'
                path_igv = os.path.join('/var/www/html/IGV',dest,project)
            elif self.dest == 'z':
                dest = 'ricerca'
                path_django = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/coverage/'
                path_igv = os.path.join('/var/www/html/IGV',dest,project)


            df6 = df5[['CHROM','start_buco','end_buco','GENE','mean']]
            df7 = df6[df6['CHROM']!=0]

            if df7.empty:
                df7 = pd.DataFrame()
            csv_stat2 = os.path.join(folder_stat, sample_x+'_gene_cov.csv')
            df5.to_csv(os.path.join(path_django,sample_x+'_gene_cov.csv'),sep='\t',index=False,encoding='utf-8')
            df5.to_csv(csv_stat2,sep='\t',index=False,encoding='utf-8')
            #system(' '.join(['scp',csv_stat2,"bioinfo@192.168.1.133:/home/magi/VIRTUAL38/apimagi_dev/NGS_RESULT/coverage/"]))
            #system(' '.join(['cp',csv_stat2,"/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/coverage/"]))
            #print (path_igv,'AAAAAAAAAA')
            try:
                df7.to_csv(os.path.join(path_igv,sample_x+'_buchi.bed'),sep='\t',index=False,header=False,encoding='utf-8')
            except:
                path_igv = os.path.join('/var/www/html/IGV',dest,'ALLGENE')
                df7.to_csv(os.path.join(path_igv,sample_x+'_buchi.bed'),sep='\t',index=False,header=False,encoding='utf-8')
            return
        else:
            dfass2 = all_genes_exons[['GENE','strand','#CHROM']]
            dfass3 = dfass2.drop_duplicates()

            dfass3['sample'] = sample_x
            dfass3['refseq'] = 0
            dfass3['hgmd'] = 0
            dfass3['exone'] = 0
            dfass3['CHROM'] = dfass3['#CHROM']
            dfass3['Exon_Start'] = 0
            dfass3['Exon_End'] = 0
            dfass3['start_buco'] = 0
            dfass3['end_buco'] = 0
            dfass3['len_gene'] = 0
            dfass3['len_exon'] = 0
            dfass3['buco_exone'] = 0
            dfass3['note'] = 'unique'
            dfass3['buco_gene'] = 0
            dfass3['%for_exon'] = 100
            dfass3['%for_gene'] = 100
            dfass3['mean'] = 0
            dfass3['variante'] = 0

            dfass5 = dfass3[['sample','GENE','strand','refseq','hgmd','exone','CHROM','Exon_Start','Exon_End','start_buco','end_buco',
                    'len_gene','len_exon','buco_exone','note','buco_gene','%for_exon','%for_gene','mean','variante']]

            folder_stat = os.path.join(coverage_folder, str(sample_x))
            try:
                os.makedirs(folder_stat)
            except FileExistsError:
                pass
            if self.dest == 'r':
                dest = 'rovereto'
                path_django = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/coverage/'
                path_igv = os.path.join('/var/www/html/IGV',dest,project)
            elif self.dest == 'b':
                dest = 'bolzano'
                path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/coverage/'
                path_igv = os.path.join('/var/www/html/IGV',dest,project)
            elif self.dest == 'z':
                dest = 'ricerca'
                path_django = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/coverage/'
                path_igv = os.path.join('/var/www/html/IGV',dest,project)

            csv_stat2 = os.path.join(folder_stat, sample_x+'_gene_cov.csv')
            dfass5.to_csv(os.path.join(path_django,sample_x+'_gene_cov.csv'),sep='\t',index=False,encoding='utf-8')
            dfass5.to_csv(csv_stat2,sep='\t',index=False,encoding='utf-8')
            print ('COVERAGE ASSENTE!!! CONTROLLA I FILE O IL FENOTIPO!!!!')

    def cov_stat(self, param_dummy, sample_coverage_all, name_sex, folder, all_genes_exons):
        sexF = 0
        sexM = 0
        _sexpredict_ = 'Undefied'

        sample_x = str(self.sample.name)
        #print("SAMPLE COVERAGE ALL = {}".format(sample_coverage_all))
        coverage = pd.read_csv((sample_coverage_all),sep='\t',header=0)
        coverage['POS'] = coverage['POS'].astype(int)
        _all_macro = os.path.join(folder, str(sample_x), str(sample_x)+'_all_macro')
        all_macro = pd.read_csv(_all_macro,sep='\t',header=0)
        print("NAME SEX = {}".format(name_sex))
        sex = pd.read_csv((name_sex),sep='\t',header=0)

        sex['filt'] = 1

        for index, row in sex.iterrows():
            # print(row)
            if row['GENE'] == 'AMELX':
                # print('AMELX')
                if (int(row['POS']) < 11296910) | (int(row['POS']) > 11296920):
                    # print('in')
                    sex.loc[index,'filt'] = 0
            elif row['GENE'] == 'AMELY':
                if (int(row['POS']) < 6869900) | (int(row['POS']) > 6869920):
                    sex.loc[index,'filt'] = 0
            elif row['GENE'] == 'SRY':
                if (int(row['POS']) < 2787170) | (int(row['POS']) > 2787270):
                    sex.loc[index,'filt'] = 0

        sex = sex[sex['filt']==1]
        sex = sex.drop('filt', axis=1)
        sex.to_csv(name_sex,sep='\t')

        try: sexF = len(sex[sex['GENE'].isin(['AMELX'])])
        except: sexF=0

        try: sexM = len(sex[sex['GENE'].isin(['AMELY','SRY'])])
        except: sexM= 0

        if sexF == 11:
            if sexM == 122: _sexpredict_ = 'M'
            elif sexM == 0: _sexpredict_ = 'F'
            else: _sexpredict_ = 'Undefied'
        else: _sexpredict_ = 'Undefied'

        group_gene = coverage.groupby(['GENE'])
        gene_size = pd.DataFrame(group_gene.size()).reset_index()
        gene_size['len_gene'] = gene_size[0]
        gene_size.drop(0,axis=1,inplace=True)

        pheno = coverage['GENE'].drop_duplicates()
        cov_10 = len(coverage[coverage['DEPTH'] >= 10])
        cov_25 = len(coverage[coverage['DEPTH'] >= 25])
        cov_macro_25 = len(all_macro[all_macro['DEPTH'] >= 25])
        cov_40 = len(coverage[coverage['DEPTH'] >= 40])
        cov_50 = len(coverage[coverage['DEPTH'] >= 50])
        cov_100 = len(coverage[coverage['DEPTH'] >= 100])
        cov_200 = len(coverage[coverage['DEPTH'] >= 200])
        buchi_10 = len(coverage[coverage['DEPTH'] < 10])
        buchi_8 = len(coverage[coverage['DEPTH'] < 8])
        buchi_5 = len(coverage[coverage['DEPTH'] < 5])
        buchi_0 = len(coverage[coverage['DEPTH'] == 0])
        tot = cov_10+buchi_10
        tot_macro = all_macro.shape[0]

        try:
            cov_10_perc ='{:,.1f}'.format(float(cov_10)/float(tot)*100)
            cov_25_perc = '{:,.1f}'.format(float(cov_25)/float(tot)*100)
            cov_macro_25_perc = '{:,.1f}'.format(float(cov_macro_25)/float(tot_macro)*100)
            cov_40_perc = '{:,.1f}'.format(float(cov_40)/float(tot)*100)
            cov_50_perc = '{:,.1f}'.format(float(cov_50)/float(tot)*100)
            cov_100_perc = '{:,.1f}'.format(float(cov_100)/float(tot)*100)
            cov_200_perc = '{:,.1f}'.format(float(cov_200)/float(tot)*100)
            buchi_10_perc ='{:,.1f}'.format(float(buchi_10)/float(tot)*100)
            buchi_8_perc = '{:,.1f}'.format(float(buchi_8)/float(tot)*100)
            buchi_5_perc = '{:,.1f}'.format(float(buchi_5)/float(tot)*100)
            buchi_0_perc = '{:,.1f}'.format(float(buchi_0)/float(tot)*100)
            media = '{:,.1f}'.format(coverage['DEPTH'].mean())
        except ZeroDivisionError:
            print ('PROBLEME IN THIS SAMPLE:',sample_x,'!!!!!!!!!!!!')
            cov_10_perc = 0.0
            cov_25_perc = 0.0
            cov_40_perc = 0.0
            cov_50_perc = 0.0
            cov_100_perc = 0.0
            cov_200_perc = 0.0
            buchi_10_perc = 0.0
            buchi_8_perc = 0.0
            buchi_5_perc = 0.0
            buchi_0_perc = 0.0
            media = 0.0
            cov_macro_25_perc = 0.0

        index=['depth>=10','depth>=25','depth>=40','depth>=50','depth>=100','depth>=200','depth<10','depth<8','depth<5','depth=0','depth_macro>=25']
        all_cov = np.array([cov_10,cov_25,cov_40,cov_50,cov_100,cov_200,buchi_10,buchi_8,buchi_5,buchi_0,cov_macro_25])
        all_cov_perc = np.array([cov_10_perc,cov_25_perc,cov_40_perc,cov_50_perc,cov_100_perc,cov_200_perc,buchi_10_perc,
                    buchi_8_perc,buchi_5_perc,buchi_0_perc,cov_macro_25_perc])
        stat_ = pd.DataFrame({'cutoff':index,'count' : all_cov,'perc' : all_cov_perc,'cov_medio':media})
        stat_['sample'] = str(sample_x)
        stat_['sesso'] = str(_sexpredict_)
        stat_ = stat_[['sample','cutoff','count','perc','cov_medio','sesso']]

        folder_stat = os.path.join(folder, str(sample_x))
        try:
            os.makedirs(folder_stat)
        except FileExistsError:
            pass

        if self.dest == 'r':
            dest = 'rovereto'
            path_django = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/coverage/'
        elif self.dest == 'b':
            dest = 'bolzano'
            path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/coverage/'
        elif self.dest == 'z':
            dest = 'ricerca'
            path_django = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/coverage/'

        #stat_.to_csv(join(path_django,sample_x+'_stat_cov.csv'),sep='\t',index=False,encoding='utf-8')
        csv_stat = os.path.join(folder_stat,sample_x+'_stat_cov.csv')
        stat_.to_csv(csv_stat,sep='\t',index=False,encoding='utf-8')
        #system(' '.join(['scp',csv_stat,"bioinfo@192.168.1.133:/home/magi/VIRTUAL38/apimagi_dev/NGS_RESULT/coverage"]))
        #system(' '.join(['cp',csv_stat,"/home/magi/VIRTUAL38/apimagi_prod/NGS_RESULT/coverage"]))

        print ('--->'+str(sample_x)+'<---')

        try:
            print ('COVERAGE INFO PER:',str(sample_x))
            print ('Len tot:',len(coverage['DEPTH']),'Media: ','{:,.1f}'.format(coverage['DEPTH'].mean()))
            print ('Len COV  >=10: ',  cov_10,'-->', str(cov_10_perc)+'%')
            print ('Len COV  >=25: ',  cov_25,'-->', str(cov_25_perc)+'%')
            print ('Len COV  >=40: ',  cov_40,'-->', str(cov_40_perc)+'%')
            print ('Len COV  >=50: ',  cov_50,'-->', str(cov_50_perc)+'%')
            print ('Len COV >=100: ', cov_100,'-->', str(cov_100_perc)+'%')
            print ('Len COV >=200: ', cov_200,'-->', str(cov_200_perc)+'%')
            print ('Len BUCHI<10: ',buchi_10,'-->', str(buchi_10_perc)+'%')
            print ('Len BUCHI <8: ', buchi_8,'-->', str(buchi_8_perc)+'%')
            print ('Len BUCHI <5: ', buchi_5,'-->', str(buchi_5_perc)+'%')
            print ('Len BUCHI =0: ', buchi_0,'-->', str(buchi_0_perc)+'%')
            print ('Len Target >=25: ', buchi_0,'-->', str(cov_macro_25_perc)+'%')
            print ('Sesso Predetto: ', 'sesso','-->', str(_sexpredict_))

            if float(cov_10_perc) >= 95.0:
                COV_CONFORMITA = str(sample_x+'\t'+str(media)+'\t'+str(cov_10_perc)+'\t'+str(cov_25_perc)+'\t'+str(cov_macro_25_perc)+'\t'+str(_sexpredict_)+'\t'+''+'\n')
            else:
                COV_CONFORMITA = str(sample_x+'\t'+str(media)+'\t'+str(cov_10_perc)+'\t'+str(cov_25_perc)+'\t'+str(cov_macro_25_perc)+'\t'+str(_sexpredict_)+'\t'+'!'+'\n')
        except ZeroDivisionError: print (0)

        name = os.path.join(folder,sample_x,sample_x+'_buchi')
        self.cov(None, None, folder,name,gene_size,all_genes_exons)
        return COV_CONFORMITA


    def cov_graph(self, dummy, folder, name, exons, HGMD):
        
        sample_x = str(self.sample.name)
        print (' GET FIG OF: ',sample_x)

        try:
            HGMD['POS'] = HGMD['END']
            HGMD['present'] = 1
            HGMD['mutation'] = HGMD['hgmd'].str.split('|').str.get(2).str.split(';').str.get(0)
            HGMD.drop(['START','END','hgmd'],axis=1,inplace=True)
        except KeyError:
            HGMD = HGMD

        sample_x = str(self.sample.name)
        coverage = pd.read_csv((name),sep='\t',header=0)
        coverage.sort_values(by=['#CHROM','POS'],inplace=True)
        coverage = coverage.reset_index()
        coverage['POS'] = coverage['POS'].astype(int)
        new_coverage = pd.merge(coverage,exons, on=['#CHROM','GENE','exone','refseq','length','hgmd','strand'],how='left')
        self.thread_print("SAMPLE {} COVERAGE \n {}".format(self.sample.name, new_coverage))
        grouped = new_coverage.groupby(['GENE','exone','refseq'])

        for name, group in grouped:
            group_df = pd.DataFrame(group)

            group_df['genex'] = group_df['GENE']+' ex: '+group_df['exone'].astype(str)+'_'+group_df['refseq']
            group_df = group_df[['#CHROM','POS','DEPTH','genex','strand','refseq','START','END','sample']]

            group_df.reset_index(inplace=True)
            group_df.drop('index',axis=1,inplace=True)

            group_df['DEPTH_buchi'] = group_df['DEPTH']
            start = group_df['START'][0]
            end = group_df['END'][0]

            larger = np.arange(start,end)
            large = pd.DataFrame({'POS': larger})

            new_group_ = pd.merge(group_df, large, on='POS',how='outer')
            new_group = pd.merge(new_group_,HGMD,on=['#CHROM','POS'],how='left')
            new_group.drop_duplicates(inplace=True)
            new_group['DEPTH'].fillna(10,inplace=True)
            new_group['genex'].ffill(inplace=True)
            new_group['sample'].ffill(inplace=True)
            new_group['strand'].ffill(inplace=True)
            new_group['refseq'].ffill(inplace=True)
            new_group['#CHROM'].ffill(inplace=True)
            new_group['START'].ffill(inplace=True)
            new_group['END'].ffill(inplace=True)
            new_group['genex'].bfill(inplace=True)
            new_group['sample'].bfill(inplace=True)
            new_group['strand'].bfill(inplace=True)
            new_group['refseq'].bfill(inplace=True)
            new_group['#CHROM'].bfill(inplace=True)
            new_group['START'].bfill(inplace=True)
            new_group['END'].bfill(inplace=True)
            name = new_group['genex'][0]
            new_group['POS'] = new_group['POS'].astype(int)
            new_group['ten'] = 10

            #mask1 = new_group['POS'] < start+15
            #mask2 = new_group['POS'] > end-15

            mask1 = new_group['POS'] < start+5
            mask2 = new_group['POS'] > end-5

            new_group.loc[mask1,'ten'] = -1
            new_group.loc[mask2,'ten'] = -1
            new_group.sort_values(by=['POS'],inplace=True)
            new_group.reset_index(inplace=True)
            new_group.drop('index',axis=1,inplace=True)
            x_df = new_group[new_group['present']==1]
            print("X_DF just before plot: {}".format(x_df))
    #################################################################################################
    #################################################################################################
            fig = plt.figure(figsize=(10,8))
            ax1 = fig.add_subplot(111)
            plt.title(name,fontsize=30,loc='center',color='orange',weight='bold',style='italic')
            ax1.set_xticks(np.arange(0,len(new_group['POS']),1))
            ax1.set_xlim(0,(len(new_group['POS'])))
            ax1.set_ylim(-0.5,12)
            ax1.set_xticklabels(new_group['POS'].astype(int),rotation=75,size=9)
            #print (len(ax1.set_xticklabels()))
            for i in range(len(ax1.get_xticklabels())):
                #print (i)
                ax1.get_xticklabels()[i].set_visible(False)
                #print ('OOOOKKKK')
            if len(ax1.get_xticklabels()) <= 300:
                for k in (range(0,len(ax1.get_xticklabels()),10)):
                    ax1.get_xticklabels()[k].set_visible(True)
            elif (len(ax1.get_xticklabels())>300)&(len(ax1.get_xticklabels())<500):
                for k in (range(0,len(ax1.get_xticklabels()),20)):
                    ax1.get_xticklabels()[k].set_visible(True)
            elif (len(ax1.get_xticklabels())>500)&(len(ax1.get_xticklabels())<1000):
                for k in (range(0,len(ax1.get_xticklabels()),50)):
                    ax1.get_xticklabels()[k].set_visible(True)
            elif (len(ax1.get_xticklabels())>=1000):
                for k in (range(0,len(ax1.get_xticklabels()),100)):
                    ax1.get_xticklabels()[k].set_visible(True)
            ax1.set_xlabel('position',fontsize=9,fontweight='light')
            ax1.xaxis.set_label_coords(1.04, -0.0025)
            ax1.set_ylabel('coverage',fontsize=9,fontweight='light',rotation=90)
            ax1.yaxis.set_label_coords(-0.03, 0.5)
            label_CDS = 'CDS_15nt= '+str(int(new_group['START'][0]))+'-'+str(int(new_group['END'][0]))
            line1, = ax1.plot(new_group.DEPTH,'--',label=label_CDS)
            line2, = ax1.plot(new_group.ten,'g-',linewidth=5,label='CDS')
            line3, = ax1.plot(new_group.DEPTH_buchi,'r-',linewidth=5,label='Buchi')
            i=0
            for index,row in x_df.iterrows():
                if row['mutation'] == 'null':
                    hgmd = row['POS']
                else:
                    hgmd = 'POS:'+str(row['POS'])+'->'+row['mutation']

                y = np.arange(-1,11,1)
                xx = np.array([index])
                x = np.repeat(xx, len(y))
                if i <= 2: line4, = ax1.plot(x,y,'--',linewidth=2,alpha=0.5,color='grey',label=hgmd)
                else: line4, = ax1.plot(x,y,'--',linewidth=2,alpha=0.5,color='grey')
                i+=1
            ax1.fill_between(new_group.index,new_group.DEPTH_buchi,-1,facecolor='red', alpha=0.2)
            ax1.legend(title='Legend',ncol=2,loc='best',fontsize=9,bbox_to_anchor=(1,1))
    ######################################################################################################
    ######################################################################################################
            folder_name = os.path.join(folder,sample_x,'figure')
            if self.dest == 'r':
                dest = 'rovereto'
                folder_VIRTUAL = os.path.join('/home/magi/VIRTUAL/MAGIS/MAGIS/static/NGS_image',sample_x)
                try:
                    os.makedirs(folder_VIRTUAL)
                except FileExistsError:
                    pass
            elif self.dest == 'b':
                dest = 'bolzano'
                folder_VIRTUAL = os.path.join('/home/magi/VIRTUAL/EUREGIO/EUREGIO/static/NGS_image',sample_x)
                try:
                    os.makedirs(folder_VIRTUAL)
                except FileExistsError:
                    pass
            elif self.dest == 's':
                dest = 'sanfelice'
                folder_VIRTUAL = os.path.join('/home/magi/VIRTUAL/SANFELICE/SANFELICE/static/NGS_image',sample_x)
                try:
                    os.makedirs(folder_VIRTUAL)
                except FileExistsError:
                    pass
            elif self.dest == 'z':
                dest = 'ricerca'
                folder_VIRTUAL = os.path.join('/home/magi/VIRTUAL/RICERCA/RICERCA/static/NGS_image',sample_x)
                try:
                    os.makedirs(folder_VIRTUAL)
                except FileExistsError:
                    pass

            try:
                os.makedirs(folder_name)
            except FileExistsError:
                pass

            name2 = name.replace(' ','')
            fig_name = sample_x+'_'+name2+'.png'
            fig.savefig(os.path.join(folder_name,fig_name))
            fig.savefig(os.path.join(folder_VIRTUAL,fig_name))
            plt.close(fig)


    def run(self):
        folder_name = dir_tree.principal_directory.path
        coverage_folder = dir_tree.principal_directory.coverage.path
        coverage_folder_sample = os.path.join(coverage_folder, str(self.sample.name))
        sample_coverage_all = os.path.join(coverage_folder_sample, str(self.sample.name) + "_all")
        sample_buchi_all = os.path.join(coverage_folder_sample, str(self.sample.name) + "_buchi")
        sample_sex_final = os.path.join(coverage_folder_sample, str(self.sample.name) + "_final_sex")

        name_all = sample_coverage_all
        name_buchi = sample_buchi_all
        name_sex = sample_sex_final


        folder_control = dir_tree.principal_directory.control.path
        folder_final = dir_tree.principal_directory.final.path
        folder_pheno = dir_tree.principal_directory.pheno.path

        lastFOLDER = os.path.basename(os.path.normpath(folder_name))
        self.thread_print("FOLDERNAME FOR SAMPLE {} IS: {}".format(self.sample.name, folder_name))

        if self.dest == 'r':
            dest = 'rovereto'
            path_download = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/CONFORMITA_COV/ROVERETO/%s' % (lastFOLDER)
            path_download2 = '/home/magi/VIRTUAL/MAGIS/DOWNLOADS/NGSINFO/CONFORMITA_COV/%s' % (lastFOLDER)
            try:
                os.makedirs(path_download)
            except FileExistsError:
                pass
            try:
                os.makedirs(path_download2)
            except FileExistsError:
                pass

            # if not os.path.exists(path_download):
            #     os.makedirs(path_download)
            # if not os.path.exists(path_download2):
            #     os.makedirs(path_download2)
            
            file_download2 = (os.path.join(path_download2,'COV_CONFORMITA_'+lastFOLDER))
            self.thread_filewrite('SAMPLE'+'\t'+'MEDIA'+'\t'+'COV 10X'+'\t'+'COV 25X'+'\t'+'Target 25%'+'\t'+'Sex'+'\t'+'10X<95%'+'\n', file_download2, 'w', self.lock)
            # out_file_download2 = open(file_download2,"w")
            # out_file_download2.write()
        elif self.dest == 'b':
            dest = 'bolzano'
            path_download = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/CONFORMITA_COV/BOLZANO/%s' % (lastFOLDER)
            try:
                os.makedirs(path_download)
            except FileExistsError:
                pass
            # if not os.path.exists(path_download):
            #     os.makedirs(path_download)
        elif self.dest == 's':
            dest = 'sanfelice'
            path_download = '/home/magi/VIRTUAL/SANFELICE/DOWNLOADS/NGSINFO/CONFORMITA_COV/%s' % (lastFOLDER)
            try:
                os.makedirs(path_download)
            except FileExistsError:
                pass
            # if not os.path.exists(path_download):
            #     os.makedirs(path_download)
        elif self.dest == 'z':
            dest = 'ricerca'
            path_download = '/home/magi/VIRTUAL/RICERCA/DOWNLOADS/NGSINFO/CONFORMITA_COV/%s' % (lastFOLDER)
            try:
                os.makedirs(path_download)
            except FileExistsError:
                pass
            # if not os.path.exists(path_download):
            #     os.makedirs(path_download)

        file_download = (os.path.join(path_download,'COV_CONFORMITA_'+lastFOLDER))
        self.thread_filewrite('SAMPLE'+'\t'+'MEDIA'+'\t'+'COV 10X'+'\t'+'COV 25X'+'\t'+'Target 25%'+'\t'+'Sex'+'\t'+'10X<95%'+'\n', file_download, "w", self.lock)
        #out_file_download = open(file_download,"w")
        #out_file_download.write('SAMPLE'+'\t'+'MEDIA'+'\t'+'COV 10X'+'\t'+'COV 25X'+'\t'+'Target 25%'+'\t'+'Sex'+'\t'+'10X<95%'+'\n')
        

        if self.genome_type == 'geno37':
            raise("Geno37 is deprecated!")        

        hgmd = pd.read_csv(config.HGMD38, sep='\t',header=None,names=['#CHROM','START','END','hgmd'],encoding='latin',engine='python')
		
        sample_x = str(self.sample.name)
        #log_file2=join(folder_name,'log',sample_x+'_coverage.log')
        #writer2 = Writer(sys.stdout,log_file2)
        #sys.stdout = writer2
        bedfile = os.path.join(folder_pheno,'bed_'+sample_x)
        all_genes_exons = pd.read_csv(bedfile,sep='\t',header=0,encoding='utf-8')
        #try:
        # name_buchi = name_buchi[0]
        # name_all = name_all[0]
        # name_sex = name_sex[0]

        COV_conformita = self.cov_stat(None, name_all, name_sex, coverage_folder, all_genes_exons)
        
        try: 
            self.cov_graph(None, coverage_folder, name_buchi, all_genes_exons, hgmd)
        except Exception as e:
            print("Coverage graph for sample {} was not generated ... {}".format(self.sample.name, e)) 

        #out_file_download.write(COV_conformita)
        self.thread_filewrite(COV_conformita, file_download, 'a', self.lock)

        self.thread_print("Finished writing to COV_confromita, for sample {} in file: {}.".format(self.sample.name, file_download))
        
        if self.dest == 'r':
            #out_file_download2.write(COV_conformita)
            self.thread_filewrite(COV_conformita, file_download2, 'a', self.lock)
            print ('---------------------------------------------')
        #except IndexError:
        #	print ('\nFILE BUCHI ASSENTE',sample_x,'\n')
            
        #out_file_download.close()
        
             

