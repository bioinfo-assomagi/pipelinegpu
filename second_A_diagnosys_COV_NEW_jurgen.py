#!/home/magi/miniconda3/envs/PY310/bin/python
#17/06/2023
# -*- coding: utf-8 -*-
"""
Created on Wed MAY 18 2023

@author: Giuseppe Marceddu
"""

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
from os import listdir, system
from os.path import isfile, join
import time
pd.options.mode.chained_assignment = None
import first_diagnosys_core_V3_jurgen as step
from cutCDS_jurgen import cutCDS,cutCDS37,create_vertical

path = os.getcwd()
geno37 = join('/home/magi/','dataset/GENOME/37/all_chr37.fa')
geno38 = join('/home/magi/','dataset/GENOME/38/all_chr38.fa')
##########################################################################################
dbsnp144_37 = join('/home/magi/','dataset/dbsnp144/37/common_all_20150605_2.vcf')
dbsnp144_38 = join('/home/magi/','dataset/dbsnp144/38/common_all_20150603_2.vcf')
indel_37 = join('/home/magi/','dataset/dbsnp144/37/common_all_20150605_indel.vcf')
indel_38 = join('/home/magi/','dataset/dbsnp144/38/common_all_20150603_indel.vcf')
clinvar = join('/home/magi/','dataset/dbsnp144/38/clinvar_20140929_2.vcf')
clinvar_indel = join('/home/magi/','dataset/dbsnp144/38/clinvar_20140929_indel.vcf')
#############################################################################################
BUCHIARTIFICIALI = join('/home/magi/', 'PROJECT/diagnosys/bin/BUCHIARTIFICIALI.txt')
##############################################################################################
def InputPar():
    ####Introducing arguments
	parser = argparse.ArgumentParser(prog= 'MAGI EUREGIO DIAGNOSYS',description='Pipe from BAM to VCF ',
			epilog='Need to SAMTOOLS v1.9 and in global path or put your path [--samtools] [--bcftools]')

	parser.add_argument('-p','--path', metavar='PATH', default=path,
	                    help='[Default = Pipe run in same folder when launch command]')

	parser.add_argument('-proj','--project',metavar='PROJECT NAME',
	                    help='Insert Family name to create Folder')

	parser.add_argument('-g','--genome', metavar='choose Assembling', choices=['geno37','geno38'],
	                    default='geno38',
	                    help='Choices: geno37, geno38 Run Analysis with genome GRCh38 or GRCh37')

	parser.add_argument('-q','--quality',  metavar='QUALITYFILTER',default=18,
	                    help=' Choose a quality threshold [Default = 18]')

	parser.add_argument('-N','--threads', metavar='THREADS', default='32',
	                    help='Number of threads to run [32]')

	parser.add_argument('-s','--samtools', metavar='SAMTOOLS', default='samtools',
	                    help='Insert SAMTOOLS path')

	parser.add_argument('-v','--bcftools', metavar='BCFTOOLS', default='bcftools',
	                    help='Insert BCFTOOLS path')

	parser.add_argument('-t','--samstat', metavar='SAMSTAT', default='samstat',
	                    help='Insert SAMSTAT path')

	parser.add_argument('-gk','--gatk', metavar='GATK', default='gatk', help='Insert gatk path')

	parser.add_argument('-ovr','--over', metavar='New/OLD project', choices=['True','False'],
	                    default='False',
	                    help='Choices: Attention option False overwrite old data [Default = False] ')

	parser.add_argument('-d','--dest', metavar='destination', choices=['b','r','s','z','p'],
	                    required=True,
	                    help='Choices destination: b = bolzano; r = rovereto; s = sanfelice; z = ricerca; p = privato,  - required ')

	parser.add_argument('-pan','--panel',metavar='PANEL',required=True,
				help='Pannelli Consigliati: 1:CLM2; 2:OSSNEW; 3:OCULARE; 4:OCULARE2; 5:SIFSR; 6:ALLGENE;'
					'7:trusightone; 8:GENETEST1; 9:GENETEST2; 10:CANCER; 11:MALFORMATIONI VASCOLARI;'
					'12:ANOMALIE VASCOLARI; 13:INFERTILITA; 14:OBESITA; 15:GENODERMATOSI\nRequired')

#	parser.add_argument('-o','--name', metavar='Output ', required=True,help='Choose an output filename (Required)')
	return parser.parse_args()

class Writer:

    def __init__(self, stdout, filename):
        self.stdout = stdout
        self.logfile = open(filename, 'a')

    def write(self, text):
        self.stdout.write(text)
        self.logfile.write(text)

    def close(self):
        self.stdout.close()
        self.logfile.close()

def PrintLog(command,folder):
    ###Print all commands in a log file
	path = join(folder,'Commands2.log')
	cl = open(path,'a')
	cl.write(command)
	cl.write('\n')
	cl.close()
##########################################################################
def count_unbalance(param,genome,vertical,verticalX,folder,bam,phenotype, vertical_macro):
	print ('Trovato:',bam)
	a = bam.split("/")
	if '-' in a[-1]:
		x = a[-1].split("-")
		x = (x[0]+'.'+x[1])
	else:
		x = a[-1].split(".")
		x = (x[0]+'.'+x[1])

	x1 = x.split("_")
	sample = x1[0]
	print ('-------->'+sample+'<------------')

	file_to_count = join(folder,'temp/',sample+'_to_count')
	file_vcf = join(folder,'temp/',sample+'_samt.vcf')

	system(' '.join(['samtools','mpileup',
				'-Q 0 -q 0 -d10000000 -L 100000 -A', bam, '>', file_to_count]))

	print ('load count!!!')
	a = pd.read_csv(file_to_count,sep='\t',header=None,quoting=csv.QUOTE_NONE,encoding='utf-8',low_memory=False,
                              on_bad_lines='skip',names=['CHROM','POS','info','DEPTH','CALL','quality'],chunksize=40*100024)
	folder_to_count = join(folder,'temp/','to_count/')
	if not os.path.exists(folder_to_count):
		os.makedirs(folder_to_count)
	folder_to_macroarea = join(folder,'temp/','to_macroarea/')
	if not os.path.exists(folder_to_macroarea):
		os.makedirs(folder_to_macroarea)

	i=0
	for chunk in a:
		i+=1
		chunk['sample'] = sample
		count = count_coverage(folder,sample,chunk)
		name = join(folder_to_count,sample+'_'+str(i))
		count_ = ''
		count_sex = ''
		count_ = filter_for_panel(vertical,folder,name,sample)
		count_sex = filter_for_SEX(verticalX,folder,name,sample)
		try: filter_for_macroarea(vertical_macro,folder,sample,i)
		except: filter_for_macroarea(verticalX,folder,sample,i)
		print ('filter_disease!!!')
		filter_disease(folder,name,phenotype, count_,sample)
		print (i)

	print ('load vcf!!!')
########################################################################
	print ('------------------------')
	return str(sample)
########################################################################
def count_coverage(folder,sample,TO_COUNT):
	TO_COUNT=TO_COUNT[['CHROM','POS','DEPTH','CALL']]

	TO_COUNT['.']=TO_COUNT['CALL'].str.count(r'[.]')
	TO_COUNT[',']=TO_COUNT['CALL'].str.count(r'[,]')
	TO_COUNT['*']=TO_COUNT['CALL'].str.count(r'[*]')
	TO_COUNT['$']=TO_COUNT['CALL'].str.count(r'[$]')
	TO_COUNT['C']=TO_COUNT['CALL'].str.count(r'[cC]')
	TO_COUNT['G']=TO_COUNT['CALL'].str.count(r'[gG]')
	TO_COUNT['T']=TO_COUNT['CALL'].str.count(r'[tT]')
	TO_COUNT['A']=TO_COUNT['CALL'].str.count(r'[aA]')
	TO_COUNT['ins']=TO_COUNT['CALL'].str.count(r'\+[0-9]+[ACGTNacgtn]+')
	TO_COUNT['del']=TO_COUNT['CALL'].str.count(r'\-[0-9]+[ACGTNacgtn]+')
	TO_COUNT['other']=TO_COUNT['CALL'].str.count(r'[*><$^]')

	TO_COUNT['sum'] =TO_COUNT['.']+TO_COUNT[',']+TO_COUNT['.']+TO_COUNT['C']+TO_COUNT['T']+TO_COUNT['A']+TO_COUNT['G']+TO_COUNT['ins']+TO_COUNT['del']

	TO_COUNT['#CHROM'] = TO_COUNT['CHROM']
	TO_COUNT = TO_COUNT[['#CHROM','POS','DEPTH','C','G','T','A','ins','del','sum']]

	TO_COUNT['DEPTH'].fillna(0,inplace=True)
	TO_COUNT['sum'].fillna(0,inplace=True)

	TO_COUNT['sum']=TO_COUNT['sum'].astype(int)
	TO_COUNT['POS']=TO_COUNT['POS'].astype(int)
	TO_COUNT['DEPTH']=TO_COUNT['DEPTH'].astype(int)

	TO_COUNT['selection'] = 1
	TO_COUNT_final=TO_COUNT[['#CHROM','POS','DEPTH','C','G','T','A','ins','del','sum','selection']]
	print ('COUNT OK!!!')
	print ('all count:',len(TO_COUNT_final))

	name_filter_intermedio = (sample+'_unbalance_inter.txt')
	file_filtered_intermedio = join(folder,'temp/',name_filter_intermedio)

	TO_COUNT_final.to_csv(file_filtered_intermedio,sep='\t',index=False)
	return

##########################################################
######################################################################
def filter_for_macroarea(vertical,folder,sample,i):
	print('filter_for_macroarea')
	name_filter_intermedio = (sample+'_unbalance_inter.txt')
	file_filtered_intermedio = join(folder,'temp/',name_filter_intermedio)
	COUNT = pd.read_csv(file_filtered_intermedio,sep='\t',header=0)

	vertical['filt'] = 1
	print( '-------FILTERING-------')
	FILTER = pd.merge(COUNT,vertical,on=['#CHROM','POS'],how='left')

	FILTER['DEPTH'].fillna(0,inplace=True)
	FILTER['sum'].fillna(0,inplace=True)
	FILTER['selection'].fillna(0,inplace=True)
	FILTER['filt'].fillna(0,inplace=True)

	FILTER = FILTER[FILTER['filt']==1]

	FILTER['C%'] = (FILTER['C'].astype(float)/FILTER['sum'].astype(float))
	FILTER['G%'] = (FILTER['G'].astype(float)/FILTER['sum'].astype(float))
	FILTER['T%'] = (FILTER['T'].astype(float)/FILTER['sum'].astype(float))
	FILTER['A%'] = (FILTER['A'].astype(float)/FILTER['sum'].astype(float))
	FILTER['ins%'] = (FILTER['ins'].astype(float)/FILTER['sum'].astype(float))
	FILTER['del%'] = (FILTER['del'].astype(float)/FILTER['sum'].astype(float))

	FILTER['C%'] = FILTER['C%'].map('{:,.3f}'.format)
	FILTER['G%'] = FILTER['G%'].map('{:,.3f}'.format)
	FILTER['T%'] = FILTER['T%'].map('{:,.3f}'.format)
	FILTER['A%'] = FILTER['A%'].map('{:,.3f}'.format)
	FILTER['ins%'] = FILTER['ins%'].map('{:,.3f}'.format)
	FILTER['del%'] = FILTER['del%'].map('{:,.3f}'.format)
	FILTER.drop_duplicates(subset=['#CHROM','POS','DEPTH'], keep='last')
	FILTER =  FILTER[['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum','DEPTH','GENE',
                                'exone','length','strand','refseq','hgmd',
                                'selection','filt']].sort_values(by=['#CHROM','POS'],ascending=[True,True])

	FILTER['sum'] = FILTER['sum'].astype(int)
	FILTER['POS'] = FILTER['POS'].astype(int)
	FILTER['DEPTH'] = FILTER['DEPTH'].astype(int)
	FILTER['selection'] = FILTER['selection'].astype(int)
	FILTER['filt'] = FILTER['filt'].astype(int)
	FILTER['exone'] = FILTER['exone'].astype(int)
	FILTER['length'] = FILTER['length'].astype(int)
	FILTER['strand'] = FILTER['strand'].astype(int)
	FILTER2 = FILTER[FILTER['filt']==1]

	print ('len FILTER:',len(FILTER2))
	print ('FILTER OK!!!')

	folder_to_macroarea = join(folder,'temp/','to_macroarea')
	FILTER2.to_csv(join(folder_to_macroarea,sample+'_'+str(i)+'_macroarea'), sep='\t', index=False)
	return FILTER2

def filter_for_panel(vertical,folder,name,sample):
	name_filter_intermedio = (sample+'_unbalance_inter.txt')
	file_filtered_intermedio = join(folder,'temp/',name_filter_intermedio)
	COUNT = pd.read_csv(file_filtered_intermedio,sep='\t',header=0)
	vertical['filt'] = 1
	name_filter = (sample+'_unbalance_filter.txt')
	file_filtered = join(folder,'temp/',name_filter)
	print( '-------FILTERING-------')
	FILTER = pd.merge(COUNT,vertical,on=['#CHROM','POS'],how='left')
	FILTER['DEPTH'].fillna(0,inplace=True)
	FILTER['sum'].fillna(0,inplace=True)
	FILTER['selection'].fillna(0,inplace=True)
	FILTER['filt'].fillna(0,inplace=True)
	FILTER = FILTER[FILTER['filt']==1]
	FILTER['C%'] = (FILTER['C'].astype(float)/FILTER['sum'].astype(float))
	FILTER['G%'] = (FILTER['G'].astype(float)/FILTER['sum'].astype(float))
	FILTER['T%'] = (FILTER['T'].astype(float)/FILTER['sum'].astype(float))
	FILTER['A%'] = (FILTER['A'].astype(float)/FILTER['sum'].astype(float))
	FILTER['ins%'] = (FILTER['ins'].astype(float)/FILTER['sum'].astype(float))
	FILTER['del%'] = (FILTER['del'].astype(float)/FILTER['sum'].astype(float))
	FILTER['C%'] = FILTER['C%'].map('{:,.3f}'.format)
	FILTER['G%'] = FILTER['G%'].map('{:,.3f}'.format)
	FILTER['T%'] = FILTER['T%'].map('{:,.3f}'.format)
	FILTER['A%'] = FILTER['A%'].map('{:,.3f}'.format)
	FILTER['ins%'] = FILTER['ins%'].map('{:,.3f}'.format)
	FILTER['del%'] = FILTER['del%'].map('{:,.3f}'.format)
	FILTER.drop_duplicates(subset=['#CHROM','POS','DEPTH'], keep='last')
	FILTER =  FILTER[['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum','DEPTH','GENE',
                                'exone','length','strand','refseq','hgmd',
                                'selection','filt']].sort_values(by=['#CHROM','POS'],ascending=[True,True])
	FILTER['sum'] = FILTER['sum'].astype(int)
	FILTER['POS'] = FILTER['POS'].astype(int)
	FILTER['DEPTH'] = FILTER['DEPTH'].astype(int)
	FILTER['selection'] = FILTER['selection'].astype(int)
	FILTER['filt'] = FILTER['filt'].astype(int)
	FILTER['exone'] = FILTER['exone'].astype(int)
	FILTER['length'] = FILTER['length'].astype(int)
	FILTER['strand'] = FILTER['strand'].astype(int)
	FILTER2 = FILTER[FILTER['filt']==1]
	print ('len FILTER:',len(FILTER2))
	print ('FILTER OK!!!')
	return FILTER2
######################################################
def filter_for_SEX(verticalx,folder,name,sample):
	name_filter_intermedio = (sample+'_unbalance_inter.txt')
	file_filtered_intermedio = join(folder,'temp/',name_filter_intermedio)
	COUNT = pd.read_csv(file_filtered_intermedio,sep='\t',header=0)
	vertical = verticalx[verticalx['GENE'].isin(['AMELX','AMELY','SRY'])]
	vertical['filt'] = 1
	name_filter = (sample+'_unbalance_filter.txt')
	file_filtered = join(folder,'temp/',name_filter)
	print ('-------FILTERING-------')
	FILTER = pd.merge(COUNT,vertical,on=['#CHROM','POS'],how='left')
	FILTER['DEPTH'].fillna(0,inplace=True)
	FILTER['sum'].fillna(0,inplace=True)
	FILTER['selection'].fillna(0,inplace=True)
	FILTER['filt'].fillna(0,inplace=True)
	FILTER = FILTER[FILTER['DEPTH']>=5]
	FILTER = FILTER[FILTER['filt']==1]
	FILTER['C%'] = (FILTER['C'].astype(float)/FILTER['sum'].astype(float))
	FILTER['G%'] = (FILTER['G'].astype(float)/FILTER['sum'].astype(float))
	FILTER['T%'] = (FILTER['T'].astype(float)/FILTER['sum'].astype(float))
	FILTER['A%'] = (FILTER['A'].astype(float)/FILTER['sum'].astype(float))
	FILTER['ins%'] = (FILTER['ins'].astype(float)/FILTER['sum'].astype(float))
	FILTER['del%'] = (FILTER['del'].astype(float)/FILTER['sum'].astype(float))
	FILTER['C%'] = FILTER['C%'].map('{:,.1f}'.format)
	FILTER['G%'] = FILTER['G%'].map('{:,.1f}'.format)
	FILTER['T%'] = FILTER['T%'].map('{:,.1f}'.format)
	FILTER['A%'] = FILTER['A%'].map('{:,.1f}'.format)
	FILTER['ins%'] = FILTER['ins%'].map('{:,.1f}'.format)
	FILTER['del%'] = FILTER['del%'].map('{:,.1f}'.format)
	FILTER.drop_duplicates(subset=['#CHROM','POS','DEPTH'], keep='last')
	FILTER.drop('selection', axis=1,inplace=True)
	FILTER.drop('filt', axis=1,inplace=True)
	FILTER['sum'] = FILTER['sum'].astype(int)
	FILTER['POS'] = FILTER['POS'].astype(int)
	FILTER['DEPTH'] = FILTER['DEPTH'].astype(int)
	FILTER['exone'] = FILTER['exone'].astype(int)
	FILTER['length'] = FILTER['length'].astype(int)
	FILTER['strand'] = FILTER['strand'].astype(int)
	FILTER_SEX =  FILTER[['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum','DEPTH','GENE',
                                'exone','length','strand','refseq','hgmd']].sort_values(by=['#CHROM','POS'],ascending=[True,True])
	
	FILTER_SEX['sample'] = str(sample)
	print ('len FILTER SEX:',len(FILTER_SEX))
	print ('FILTER SEX OK!!!')
	FILTER_SEX.to_csv(name+'_SEX', sep='\t',header=False,index=False)
	return FILTER_SEX
########################################################################################
########################################################################################
def filter_disease(folder,name,phenotype,COUNT,sample):
	folder_coverage = join(folder,'coverage/',sample)
	result = join(folder_coverage,sample+'_all')
	phenotype_ = phenotype[phenotype['sample'].astype(str) == str(sample)]
	a = phenotype_[['malattia','gene']].drop_duplicates()
	x1 = pd.DataFrame(a['gene'])
	b = COUNT #[COUNT['GENE'].isin(x1['gene'])]
	b.drop('selection', axis=1,inplace=True)
	b.drop('filt', axis=1,inplace=True)
	b.to_csv(name+'_disease', sep='\t',header=False,index=False)
	return
######################################################################################
def filter_macroarea_files(vertical_macro,folder,sample,buchiartificiali):
	print  ('start filter_macroarea_files!!!')
	folder_to_macroarea = join(folder,'temp/','to_macroarea')
	path_macro = join(folder_to_macroarea,sample+'_*_macroarea')
	print(path_)
	files_macro = glob.glob(path_macro)
	macro_count = pd.DataFrame(columns=['#CHROM', 'POS', 'C%', 'G%', 'T%', 'A%', 'ins%', 'del%', 'sum', 'DEPTH','GENE',	'exone', 'length', 'strand', 'refseq', 'hgmd', 'selection', 'filt'])
	for file_macro in files_macro:
		#concat files
		print(file_macro)
		macro = pd.read_csv(file_macro, sep='\t')
		macro_count = macro_count._append(macro, ignore_index=True)
	print(macro_count)
	if True:
		macro_count =  macro_count.sort_values(by=['#CHROM','POS'],ascending=[True,True])
		print(macro_count)
		b=macro_count
		for index,row in buchiartificiali.iterrows():
			x = row['#CHROM']
			y = row['START']
			z = row['END']
			mask1 =  ((b['#CHROM'] == x) & (b['POS'] >= y) & (b['POS'] <= z))
			b.loc[mask1,'DEPTH'] = 0

		b['DEPTH'].fillna(0,inplace=True)
		b['sum'].fillna(0,inplace=True)
		b['filt'].fillna(0,inplace=True)
		b['sum'] = b['sum'].astype(int)
		b['POS'] = b['POS'].astype(int)
		b['DEPTH'] = b['DEPTH'].astype(int)
		b['filt'] = b['filt'].astype(int)
		b['exone'] = b['exone'].astype(int)
		b['length'] = b['length'].astype(int)
		b['strand'] = b['strand'].astype(int)

		COUNT = b
		COUNT['sample'] = sample
		COUNT.fillna(0,inplace=True)
		print (len(COUNT))
		print ('Sequence: ',len(COUNT))

		folder_coverage = join(folder,'coverage/',sample)
		if not os.path.exists(folder_coverage):
			os.makedirs(folder_coverage)

		result = join(folder_coverage,sample+'_all_macro')
		result_buchi = join(folder_coverage,sample+'_buchi_macro')
		result_non_buchi = join(folder_coverage,sample+'_marked_macro')
		result_not_marked = join(folder_coverage,sample+'_only_0_macro')

		cov = len(COUNT[COUNT['DEPTH'] >= 10])
		# cov = len(COUNT[COUNT['DEPTH'] >= 20])
		buchi = len(COUNT[COUNT['DEPTH'] < 10])
		# buchi = len(COUNT[COUNT['DEPTH'] < 20])
		tot = cov+buchi
		COUNT['sample'] = str(sample)

		buchi = COUNT[COUNT['DEPTH'] < 10]
		marked = COUNT[COUNT['DEPTH'] >= 10]
		# buchi = COUNT[COUNT['DEPTH'] < 10]
		# marked = COUNT[COUNT['DEPTH'] >= 10]
		not_marked = COUNT[COUNT['DEPTH'] == 0]

		buchi.to_csv(result_buchi,sep='\t',index=False)
		marked.to_csv(result_non_buchi,sep='\t',index=False)
		not_marked.to_csv(result_not_marked,sep='\t',index=False)
		COUNT.to_csv(result,sep='\t',index=False)

		# print ('Len COV > 10: ',cov)
		print ('Len BUCHI: ',len(buchi))
		try: print ('% Sequence COVERED: ', '{:,.1f}'.format(float(cov)/float(tot)*100),'%')
		except ZeroDivisionError: print (0)

		print ('Len disease: ',len(COUNT))
	return

def filter_disease_files(BED,vertical,folder,phenotype, COUNT,sample,buchiartificiali):
	print  ('start final filter!!!')
	COUNT =  COUNT[['#CHROM','POS','C%','G%','T%','A%','ins%',
				'del%','sum','DEPTH']].sort_values(by=['#CHROM','POS'],ascending=[True,True])
	phenotype_ = phenotype[phenotype['sample'].astype(str) == str(sample)]
	a = phenotype_[['malattia','gene']].drop_duplicates()
	x1=pd.DataFrame(BED['GENE'].drop_duplicates())
	#x1 = pd.DataFrame(a['gene'])
	x2= pd.DataFrame({'GENE':pd.Series(['CTD-3074O7.11','TM4SF2'])})
	x = x1._append(x2)
	FILTER = vertical[vertical['GENE'].isin(x['GENE'])] # so merge two times with the vertical?

	b = pd.merge(COUNT,FILTER,on=['#CHROM','POS'],how='right')
	for index,row in buchiartificiali.iterrows():
		x = row['#CHROM']
		y = row['START']
		z = row['END']
		mask1 =  ((b['#CHROM'] == x) & (b['POS'] >= y) & (b['POS'] <= z))
		b.loc[mask1,'DEPTH'] = 0

	b['DEPTH'].fillna(0,inplace=True)
	b['sum'].fillna(0,inplace=True)
	try: b['filt'].fillna(0,inplace=True)
	except: b['filt']=0
	b['sum'] = b['sum'].astype(int)
	b['POS'] = b['POS'].astype(int)
	b['DEPTH'] = b['DEPTH'].astype(int)
	b['filt'] = b['filt'].astype(int)
	b['exone'] = b['exone'].astype(int)
	b['length'] = b['length'].astype(int)
	b['strand'] = b['strand'].astype(int)

	COUNT = b
	COUNT['sample'] = sample
	COUNT.fillna(0,inplace=True)
	COUNT.drop(['filt'], axis=1,inplace=True)
	print (len(COUNT))
	print ('Sequence: ',len(COUNT))

	folder_coverage = join(folder,'coverage/',sample)
	if not os.path.exists(folder_coverage):
		os.makedirs(folder_coverage)

	result = join(folder_coverage,sample+'_all')
	result_buchi = join(folder_coverage,sample+'_buchi')
	result_non_buchi = join(folder_coverage,sample+'_marked')
	result_not_marked = join(folder_coverage,sample+'_only_0')

	cov = len(COUNT[COUNT['DEPTH'] >= 10])
	buchi = len(COUNT[COUNT['DEPTH'] < 10])
	# cov = len(COUNT[COUNT['DEPTH'] >= 20])
	# buchi = len(COUNT[COUNT['DEPTH'] < 20])
	tot = cov+buchi
	COUNT['sample'] = str(sample)

	# buchi = COUNT[COUNT['DEPTH'] < 20]
	# marked = COUNT[COUNT['DEPTH'] >= 20]
	buchi = COUNT[COUNT['DEPTH'] < 10]
	marked = COUNT[COUNT['DEPTH'] >= 10]
	not_marked = COUNT[COUNT['DEPTH'] == 0]

	buchi.to_csv(result_buchi,sep='\t',index=False)
	marked.to_csv(result_non_buchi,sep='\t',index=False)
	not_marked.to_csv(result_not_marked,sep='\t',index=False)
	COUNT.to_csv(result,sep='\t',index=False)

	print( 'Len COV > 10: ',cov)
	print ('Len BUCHI: ',len(buchi))
	try: print ('% Sequence COVERED: ', '{:,.1f}'.format(float(cov)/float(tot)*100),'%')
	except ZeroDivisionError: print (0)

	print ('Len disease: ',len(COUNT))
	return
##############################################################################
##############################################################################
if __name__=="__main__":
	lock = Lock()
	args=InputPar()

	if args.over == 'False':
		folder_name = step.principal_folder(args,step.create_folder(),over='False')
	if args.over == 'True':
		folder_name = step.principal_folder(args,step.create_folder(),over='True')
		step.path_creation(args,folder_name)
		print ('NOW YOU HAVE TO COPY BAM FILES IN FOLDER (BAM)\nAND RUN COMMAND AGAIN\nWHITOUT (--OVER)')
		sys.exit

	folder_pheno = join(folder_name,'pheno')
	input_phenotype = join(folder_name,'pheno/phenotype')
	print('here')

	print_args=vars(args)#take args in a dict and print in a log
	for k,v in print_args.items():
		strategy_out='='.join([str(k),str(v)])
		path = join(folder_name,'log')
		PrintLog(strategy_out,path)
#####################################################################################
####################################################################################
	if args.genome == 'geno38':
		genome = geno38
		dbsnp = dbsnp144_38
	elif args.genome == 'geno37':
		genome = geno37
		dbsnp = dbsnp144_37
	path_ = join(folder_name,'bam/*.bam')
	path_bai = join(folder_name,'bam/*.bai')
	print(path_)
	files_bai = glob.glob(path_bai)
	for bai2 in files_bai:
		a = bai2.split("/")
		if '-' in a[-1]:
			x = a[-1].split("-")
			x = (x[0]+'.'+x[1])
		else:
			x = a[-1].split(".")
			x = (x[0]+'.'+x[1])
		x1 = x.split("_")
		name = x1[0]
		bai =  join(folder_name,'bam/',name+'_final.bai')
		os.rename(bai2,bai)

	files_bam = glob.glob(path_)
	for bam2 in files_bam:
		#_samples_ = ['E116.2018','E190.2018','E191.2018','E194.2018','E235.2018']
		#_samples_ = ['E243.2018']
		#_samples_ = ['74.2018','76.2018','96.2018','306.2018','309.2018','310.2018',
		#		'327.2018','346.2018','347.2018','348.2018','350.2018','351.2018']
		# _samples_ =['E68.2021']
		a = bam2.split("/")
		if '-' in a[-1]:
			x = a[-1].split("-")
			x = (x[0]+'.'+x[1])
		else:
			x = a[-1].split(".")
			x = (x[0]+'.'+x[1])

		x1 = x.split("_")
		name = x1[0]

		# if name in _samples_:
		bam =  join(folder_name,'bam/',name+'_final.bam')
		os.rename(bam2,bam)
		print(name)
		log_file2=join(folder_name,'log',name+'_vcf.log')
		writer2 = Writer(sys.stdout,log_file2)
		sys.stdout = writer2
		phenotype = pd.read_csv(input_phenotype, sep='\t', header=0) #,encoding='utf-8')
		phenotype['sample'] = phenotype['sample'].astype('str')
		phenotype['sample'].replace(r'\.202$','.2020',inplace=True,regex=True)
		phenotype['sample'] = phenotype['sample'].astype('str')
		GENELIST = []
		GENELIST = list(phenotype['gene'][phenotype['sample'].astype(str)==str(name)])
		print ('NAME:',name)
		vertical = ''
		verticalX = ''
	####################################################################################
		#try:
		print( GENELIST,'###########')
		if args.genome == 'geno38':
				vertical,verticalX,BED = cutCDS(GENELIST,name,folder_pheno,args.dest)
		elif args.genome == 'geno37':
	                        vertical = cutCDS37(GENELIST,name,folder_pheno,args.dest)

		sospetto = phenotype['malattia'][phenotype['sample'].astype(str)==str(name)].unique()[0]
		print(sospetto)

		if sospetto == 'CONNETTIVOPATIE':
			print('connettivopatie')
			bed_ = join(folder_pheno,'bed_'+str(name))
			bed = pd.read_csv(bed_, sep='\t')
			region = pd.DataFrame({'#CHROM':pd.Series(['chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9']),
			'START':pd.Series([86584322,86585077,86585812,86585652,86586188,86586587,86586797,86587759,86588201,86588817,86589432,86590377,86591910,86592604,86593110]),
				'END':pd.Series([86584355,86585246,86585827,86585734,86586271,86586641,86587104,86587887,86588314,86588888,86589504,86590420,86591966,86592701,86593167]),
				'GENE':pd.Series(['int','int','int','int','int','int','int','int','int','int','int','int','int','int','int']),
				'length':pd.Series([34,170,16,83,84,55,308,129,114,72,73,44,57,98,58]),
				'exone':pd.Series([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),'strand':pd.Series([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),
				'refseq':pd.Series(['N','N','N','N','N','N','N','N','N','N','N','N','N','N','N']),
				'hgmd':pd.Series(['int','int','int','int','int','int','int','int','int','int','int','int','int','int','int'])})
			print(bed)
			print(region)
			bed = pd.concat([bed,region])
			bed.to_csv(bed_,sep='\t', index=False)
			bedx_ = join(folder_pheno,'bedX_'+str(name))
			bedx=pd.read_csv(bedx_,sep='\t')
			vertical,verticalX = create_vertical(bed,bedx,str(name),folder_pheno)
			print(vertical)
			vertical_ = join(folder_pheno,'vertical_'+str(name))
			vertical = pd.read_csv(vertical_,sep='\t')
			verticalx_ = join(folder_pheno,'verticalX_'+str(name))
			verticalX = pd.read_csv(verticalx_,sep='\t')
		elif sospetto == 'CROMATINOPATIE':
			print('CROMATINOPATIE, in ')
			bed_ = join(folder_pheno,'bed_'+str(name))
			bed = pd.read_csv(bed_, sep='\t')
			region = pd.DataFrame({'#CHROM':pd.Series(['chr11','chr12','chr15','chr15','chr15','chr15','chr16','chr16','chr22','chr3','chr3','chr3','chr5','chr5','chr7','chr7','chr7','chr7','chr7','chr7','chr3','chr22','chr22','chr22','chr22','chr22','chr22','chrX']),
			'START':pd.Series([119077253,112910723,66995991,66995821,67001018,96875310,30134513,55515763,19770411,20202572,20202683,8591535,78280735,78280904,140624418,150693467,150700214,150700331,150700720,150700214,20216072,19748403,19753887,19754000,19754265,19754305,19444417,136649004]),
				'END':pd.Series([119077347,112910816,66996065,66995870,67001101,96875443,30134555,55515815,19770570,20202743,20202743,8591605,78280783,78280957,140624528,150693561,150700516,150700413,150700808,150700494,20216167,19748627,19754144,19754071,19754311,19754311,19444466,136649122]),
				'GENE':pd.Series(['CBL','PTPN11','int','SMAD6','int','NR2F2','MAPK3','int','int','int','int','int','int','ARSB','BRAF','NOS3','int','int','int','int','SGOL1','TBX1','TBX1','int','int','int','UFD1L','ZIC3']),
				'length':pd.Series([95,94,75,50,84,134,43,53,160,172,61,71,49,54,111,95,303,83,89,281,96,225,258,72,47,7,50,119]),
				'exone':pd.Series([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
				'strand':pd.Series([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),
				'refseq':pd.Series(['N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N']),
				'hgmd':pd.Series(['int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int'])})
			#print(bed)
			#print(region)
			bed = pd.concat([bed,region])
			bed.to_csv(bed_,sep='\t', index=False)
			bedx_ = join(folder_pheno,'bedX_'+str(name))
			bedx = pd.read_csv(bedx_,sep='\t')
			vertical,verticalX = create_vertical(bed,bedx,str(name),folder_pheno)
			print(vertical)
			vertical_ = join(folder_pheno,'vertical_'+str(name))
			vertical = pd.read_csv(vertical_,sep='\t')
			verticalx_ = join(folder_pheno,'verticalX_'+str(name))
			verticalX = pd.read_csv(verticalx_,sep='\t')


		if not vertical.empty:
			print('ifnot')
			buchiartificiali = pd.read_csv(BUCHIARTIFICIALI,sep='\t',header=0)
			if 'OCULAR' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_ONA20509.2020')
			elif 'VASCULAR' in args.panel:_vertical_macro=join(args.path,'bin','VERTICAL','vertical_VASCNA20509.2020')
			elif 'NEUROLOGY' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_NEURNA20509.2020')
			elif 'INTEGRACARDIOSTANCHEZZA' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_INTEGRACARDIOSTANCHEZZA.2023')
			elif 'MIXED' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_MIX120509.2020')
			elif 'LYMPHOBESITY' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_LIMPHNA20509.2020')
			elif 'INFERTILIT' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_INA20509.2020')
			elif 'CHERATOCONO' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_CHERATOCONO.2020')
			elif 'PCDH19' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_PCDH19')
			elif 'GENEOBNEW' in args.panel: _vertical_macro=join(args.path,'bin','VERTICAL','vertical_GENEOBNEW')
			#elif 'CANCER' in args.panel: _vertical_macro==join(args.path,'bin','VERTICAL','vertical_CANC20509.2022')
			elif 'CUSTOM' in args.panel: _vertical_macro=join(folder_pheno,'vertical_'+name)
			elif 'CANCER' in args.panel:
				print ('CANCERRRRRRRRRRRRRRRRRRRRRRRR',join(args.path,'bin','VERTICAL','vertical_CANC20509.2022'))
				_vertical_macro=join(args.path,'bin','VERTICAL','vertical_CANC20509.2022')
			else: _vertical_macro = ''  #'' #print('VERTICAL not found')

			try: vertical_macro =  pd.read_csv(_vertical_macro, sep='\t')
			except: print ('ERRORE!!!') #verti
			sample = count_unbalance(args,genome,vertical,verticalX,folder_name,bam,phenotype,vertical_macro)

			cols = ['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum','DEPTH',
		                       'GENE','exone','length','strand','refseq','hgmd','sample']
			folder_to_count = join(folder_name,'temp/','to_count/')

			folder_to_count_file = join(folder_name,'temp/','to_count/',sample+'*_disease')
			glob_ = glob.glob(folder_to_count_file)
			glob_ = str(glob_)[1:-1].replace(',',' ')
			glob_ = str(glob_)[1:-1].replace('\'','')
			system(' '.join(['cat',glob_,'>',join(folder_name,'temp/',sample+'_final_disease')]))
			name = join(folder_name,'temp/',sample+'_final_disease')
			count = pd.read_csv(name,sep='\t',header=None,names=cols)
			count.fillna(0,inplace=True)

			folder_to_count_sex = join(folder_name,'temp/','to_count/',sample+'*_SEX')
			glob_sex = glob.glob(folder_to_count_sex)
			glob_sex = str(glob_sex)[1:-1].replace(',',' ')
			glob_sex = str(glob_sex)[1:-1].replace('\'','')
			system(' '.join(['cat',glob_sex,'>',join(folder_name,'temp/',sample+'_final_sex')]))
			namesex = join(folder_name,'temp/',sample+'_final_sex')
			count_sex = pd.read_csv(namesex,sep='\t',header=None,names=cols)
			count_sex.fillna(0,inplace=True)

			filter_disease_files(BED,vertical,folder_name,phenotype,count,sample,buchiartificiali)

			file_to_count = join(folder_name,'temp/',sample+'_to_count')
			a = pd.read_csv(file_to_count,sep='\t',header=None,quoting=csv.QUOTE_NONE,encoding='utf-8',low_memory=False,
									on_bad_lines='skip',names=['CHROM','POS','info','DEPTH','CALL','quality'],chunksize=40*100024)
			try: filter_macroarea_files(vertical_macro,folder_name,sample,buchiartificiali)
			except: filter_macroarea_files(vertical,folder_name,sample,buchiartificiali)
			print('outplot')
			print(join(folder_name,'coverage/',sample,sample+'_final_sex'))
			count_sex.to_csv(join(folder_name,'coverage/',sample,sample+'_final_sex'),sep='\t',index=False)
			print('to_Csv')

			system(' '.join(['rm -r',join(folder_to_count)]))
			print('rm')
			#else: pass
		else: pass
