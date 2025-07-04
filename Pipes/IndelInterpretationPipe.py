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
import requests
from requests.structures import CaseInsensitiveDict
from requests.auth import HTTPBasicAuth
import json
import logging

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



class IndelInterpretationPipe(ParallelPipe):

    def __init__(self):
        super().__init__()

    def process(self, **kwargs):
        self.thread_id = os.getpid()
        self._logger = logging.getLogger(__name__)
        
        self.sample = kwargs.pop("sample")
        self.panel = kwargs.pop("panel", None)
        self.genome_type = kwargs.pop("genome", "geno38")
        self.dest = kwargs["dest"]

        self._logger.info("Starting IndelInterpretationPipe for sample: {}".format(self.sample.name))
        self.run()
        self._logger.info("IndelInterpretationPipe for sample: {} finished".format(self.sample.name))

        kwargs.update({"sample": self.sample, "panel": self.panel, "genome": self.genome_type, "dest": self.dest})
        return kwargs


    def query_varsome_cnv(self, HGVS):
        #https://api.varsome.com/lookup/cnv/chr6:64813372:64997708:DEL/1038
        url = "/".join(('https://stable-api.varsome.com/lookup/cnv',HGVS)) #,'hg38?add-source-databases=gnomad_genes'))
        headers = CaseInsensitiveDict()
        headers["Accept"] = "application/json"
        headers["user-agent"] = "VarSomeApiClientPython/2.0"
        headers["Authorization"] = "Token " + "so9N8mY0?liIn3k@c470r#6WU!VRllW3MT4CFcjV"
        resp = requests.get(url, headers=headers)
        #print (HGVS,resp,str(resp)=='<Response [200]>','AAAAAAAAAAAAAAAAA')
        if str(resp)=='<Response [200]>':
            data = json.loads(resp.text)
            data_df = pd.json_normalize(data['sv_acmg_annotation']['verdict']) #['items'])
            #print (pd.json_normalize(data['sv_acmg_annotation']['classifications']).columns)
            classification = data['sv_acmg_annotation']['classifications']
            print (classification)

            verdict=data_df['verdict']
            self._logger.info("Varsome called successfully!")
            return verdict,classification #SAMPL
        else:
            verdict='Unknown'
            self._logger.error("Something wrong with quering varsome.")
            return verdict
        

    def run(self):
        folder_name = dir_tree.principal_directory.path
        folder_coverage = dir_tree.principal_directory.coverage.path
        folder_indel = dir_tree.principal_directory.indel.path
        folder_pheno = dir_tree.principal_directory.pheno.path
        folder_bam = dir_tree.principal_directory.bam.path
        input_phenotype = join(dir_tree.principal_directory.pheno.path, "phenotype")


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

        if self.genome_type != 'geno38':
            return
        
        sample_x = str(self.sample.name)
        name_all = glob.glob(os.path.join(folder_indel, sample_x) + '/*_prefinal_indel.csv')
        self._logger.debug("Found {} files for sample {}".format(len(name_all), sample_x))
        for file in name_all:
            print("FILE={}".format(file))
            SAMPLE = pd.read_csv(file,sep='\t',header=0)
            #print(SAMPLE)
            SAMPLE['verdict']=''
            SAMPLENEW = pd.DataFrame()
            if len(SAMPLE)>0:
                #print (SAMPLE[['sample_id','hgvs','chrom','GENE']])
                group_sample = SAMPLE.groupby(['GENE'])
                for index, group in group_sample:
                    #print (group.iloc[0,:])
                    GENE=group.iloc[0,:]['GENE']
                    FIRSTSTART=group.iloc[0,:]['hgvs'].split(':')[1].split('-')[0]
                    FIRSTEND=group.iloc[0,:]['hgvs'].split(':')[1].split('-')[1]
                    FIRSTCHROM=group.iloc[0,:]['hgvs'].split(':')[0]
                    FIRSTABB=group.iloc[0,:]['ABBERATION'].replace('HOM_','')
                    LASTSTART=group.iloc[-1,:]['hgvs'].split(':')[1].split('-')[0]
                    LASTEND=group.iloc[-1,:]['hgvs'].split(':')[1].split('-')[1]
                    LASTCHROM=group.iloc[-1,:]['hgvs'].split(':')[0]
                    LASTABB=group.iloc[-1,:]['ABBERATION'].replace('HOM_','')
                    #print (ABB)
                    #CNVHGVS='chr6:64813372:64997708:DEL/1038'
                    if (int(FIRSTSTART)<int(LASTEND)):
                    #print (int(LASTEND)>int(FIRSTSTART))
                        CNVHGVS=str(FIRSTCHROM)+':'+str(FIRSTSTART)+':'+str(LASTEND)+':'+str(FIRSTABB)+'/'+'1038'
                    elif (int(FIRSTSTART)>int(LASTEND)):
                        CNVHGVS=str(FIRSTCHROM)+':'+str(LASTEND)+':'+str(FIRSTSTART)+':'+str(FIRSTABB)+'/'+'1038'
                    else:
                        CNVHGVS='Unkown!!!'
                    #CNVHGVS=str(FIRSTCHROM)+':'+str(FIRSTSTART)+':'+str(LASTEND)+':'+str(FIRSTABB)+'/'+'1038'
                    #print (CNVHGVS)
                    CNV,CLASSIFICATION = self.query_varsome_cnv(CNVHGVS)
                    print (CNVHGVS,CNV,CLASSIFICATION)
                    group.loc[group.index[group['GENE']==str(GENE)],'verdict'] = str(CNV[0])
                    group.loc[group.index[group['GENE']==str(GENE)],'CNVHGVS'] = str(CNVHGVS)
                    #group.loc[group.index[group['GENE']==str(GENE)],'CLASSIFICATION'] = str(CLASSIFICATION)
                    #group['verdict'].fillna(method='ffill',inplace=True)
                    #group['CNVHGVS'].fillna(method='ffill',inplace=True)
                    SAMPLENEW = pd.concat([SAMPLENEW,group])
                    #print (SAMPLENEW)
                #print (SAMPLENEW)
                SAMPLENEW.fillna(0,inplace=True)
                SAMPLENEW.to_csv(join(folder_indel,sample_x,sample_x+'_final_indel.csv'),sep='\t',index=None)
                self._logger.debug("Wrote to {}".format(join(folder_indel,sample_x,sample_x+'_final_indel.csv')))
                ##SAMPLENEW.to_csv(join(path_django,sample_x+'_final_indel.csv'),sep='\t',index=None)
                #system(' '.join(['scp',join(folder_indel,sample_x,sample_x+'_final_indel.csv'), "bioinfo@192.168.1.120:/home/bioinfo/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/"]))
                #system(' '.join(['cp',join(folder_indel,sample_x,sample_x+'_final_indel.csv'), "/home/bioinfo/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/"]))
                ##system(' '.join(['scp',join(folder_indel,sample_x,sample_x+'_final_indel.csv'), "bioinfo@192.168.1.133:/home/bioinfo/VIRTUAL38/apimagi_prod/NGS_RESULT/INDEL/"]))

