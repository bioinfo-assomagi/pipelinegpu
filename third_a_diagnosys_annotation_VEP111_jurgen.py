#!/home/magi/miniconda3/envs/PY310/bin/python
#27/06/2015

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 2022
@author: Giuseppe Marceddu
"""

import sqlite3
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
from os import listdir, system
from os.path import isfile, join
import time
import numpy as np
from io import StringIO
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from os import getenv

pd.options.mode.chained_assignment = None
import first_diagnosys_core_V3_jurgen as step

import warnings
from pandas.errors import SettingWithCopyWarning,PerformanceWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=PerformanceWarning)

path = os.getcwd()
geno37 = join('/home/magi/','dataset/GENOME/37/all_chr37.fa')
geno38 = join('/home/magi/','dataset/GENOME/38/all_chr38.fa')

eredita_37 = join('/home/magi/','dataset/MAGIS/37/eredita.ridotto.txt')
eredita_38 = join('/home/magi/','dataset/MAGIS/38/eredita.ridotto.txt')

HGMD_37 = join('/home/magi/','dataset/HGMD/37/HGMD_pro_2014_2.bed')
HGMD_38 = join('/home/magi/','dataset/HGMD/38/LiftOverHGMD_pro_2014_2.bed')
#########################################################################
appris = join('/home/magi/', 'PROJECT/diagnosys/bin/APPRIS_PRINCIPAL')
eccezioni = join('/home/magi/', 'PROJECT/diagnosys/bin/APPRIS_ECCEZIONI')
#########################################################################
#GATK = join('/home/magi/','tools/GenomeAnalysisTK36.jar')

#VEPS = join('/home/magi/','tools/ensembl-vep-release-99/vep')
#VEPS = join('/home/magi/','tools/ensembl-vep-release-104/vep')
VEPS = join('/home/magi/','tools/ensembl-vep-release-111/vep')
# ########################################################################
# dbNSFP38_gz = join('/home/magi/','.vep/Plugins/38_version4/dbNSFP4.0a.gz')
# dbscSNV11_gz = join('/home/magi/','.vep/Plugins/38_version4/dbscSNV1.1_GRCh38.txt.gz')

dbNSFP38_gz = join('/home/magi/','.vep/Plugins/38_version4/dbNSFP4.1a.gz')
dbscSNV11_gz = join('/home/magi/','.vep/Plugins/38_version4/SNV/dbscSNV1.1_GRCh38.txt.gz')
########################################################################
def InputPar():
    ####Introducing arguments
	parser = argparse.ArgumentParser(prog= 'MAGI EUREGIO DIAGNOSYS',description='Pipe from VCF to ANNOTATION ',
			epilog='Need to VEP of ENSEMBLE version 95')

	parser.add_argument('-p','--path', metavar='PATH', default=path,
	                    help='[Default = Pipe run in same folder when launch command]')

	parser.add_argument('-proj','--project',metavar='PROJECT NAME',
	                    help='Insert Family name to create Folder')

	parser.add_argument('-g','--genome', metavar='choose Assembling', choices=['geno37','geno38'],
	                    default='geno38',
	                    help='Choices: geno37, geno38 -Run Analysis with genome GRCh38 or GRCh37')

	parser.add_argument('-d','--dest', metavar='destination', choices=['b','r','s','z','p'],
	                    required=True,
	                    help='Choices destination: b = bolzano; r = rovereto; s = sanfelice; z = ricerca; p = privato,  - required ')

	parser.add_argument('-ovr','--over', metavar='New/OLD project', choices=['True','False'],
	                    default='False',
	                    help='Choices: Attention option False overwrite old data [Default = False] ')

	parser.add_argument('-pan','--panel',metavar='PANEL',required=True,
				help='Pannelli Consigliati: 1:CLM2; 2:OSSNEW; 3:OCULARE; 4:OCULARE2; 5:SIFSR; 6:ALLGENE;'
					'7:trusightone; 8:GENETEST1; 9:GENETEST2; 10:CANCER; 11:MALFORMATIONI VASCOLARI;'
					'12:ANOMALIE VASCOLARI; 13:INFERTILITA; 14:OBESITA; 15:GENODERMATOSI\nRequired')

	parser.add_argument('-vep','--VEP', metavar='variant effect predictor', default=VEPS,
                        help='Insert VEP PATH')

	return parser.parse_args()

class Writer():
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
	path = join(folder,'Commands3.log')
	cl = open(path,'a')
	cl.write(command)
	cl.write('\n')
	cl.close()


def sample(name):
	x = name.split("/")
	x = (x[-1])
	x = x.split("_")
	sample = str(x[0])
	return str(sample)


#def annotation1(param,folder_in,folder_out):
#	if param.genome == 'geno38':
#		system(' '.join([param.VEP,'--cache --refseq --offline --assembly GRCh38 --fasta', geno38,
#				'--fork 8 --force --vcf --format vcf --pick -plugin dbNSFP,%s,genename,rs_dbSNP142',
#				'-i',folder_in,'-o', folder_out]) % (dbNSFP38_gz))

#	elif param.genome == 'geno37':
#		system(' '.join([param.VEP,'--cache --refseq --offline --assembly GRCh37 --fasta', geno37,
#				'--fork 8 --force --pick --vcf --format vcf -plugin dbNSFP,%s,genename,rs_dbSNP141',
#				'-i',folder_in,'-o', folder_out]) % (dbNSFP37_gz))

#--fork 10 TODO: find the best trade off between vep fork and input files parallelisation - how many files can fork process - 1 sample and 150 fork (try different combinations)
def annotation2(param,folder_in,folder_out):
		system(' '.join([param.VEP,'--cache --refseq --offline --use_given_ref --assembly GRCh38 --fasta', geno38,
				'--fork 8 --force_overwrite --vcf --format vcf --everything --af_1kg',
				'--buffer_size 500', '--force','--xref_refseq','-i',folder_in,'-o', folder_out]))

#--fork 10
def annotation3(param,folder_in,folder_out):
	#print('annotation3')
	if param.genome == 'geno38':
		system(' '.join([param.VEP,'--cache --refseq --offline --use_given_ref --assembly GRCh38 --fasta', geno38,
				'--fork 8 --force_overwrite --vcf --format vcf --everything --af_1kg',
				'--buffer_size 500 --force --xref_refseq --exclude_predicted --use_transcript_ref',
				'--plugin dbNSFP,%s,'
				'Eigen-PC-raw_coding_rankscore,CADD_raw_rankscore,DANN_rankscore,MetaSVM_rankscore,'
				'MetaLR_rankscore,FATHMM_converted_rankscore,MutationTaster_converted_rankscore,'
				'MutationAssessor_rankscore,Polyphen2_HDIV_rankscore,Polyphen2_HVAR_rankscore,'
				'SIFT_converted_rankscore,LRT_converted_rankscore,MutPred_rankscore,'
				'PROVEAN_converted_rankscore,VEST4_rankscore,REVEL_rankscore,'
				'SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,'
				'PROVEAN_pred,MetaSVM_pred,MetaLR_pred,M-CAP_pred,fathmm-MKL_coding_pred,'
				'GERP++_RS_rankscore,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian_rankscore,'
				'phyloP30way_mammalian_rankscore,phastCons100way_vertebrate_rankscore,'
				'phastCons30way_mammalian_rankscore,SiPhy_29way_logOdds_rankscore,'
				'clinvar_clnsig,clinvar_review,Interpro_domain,'
				'gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_nhomalt,gnomAD_exomes_controls_AF,'
				'gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,'
				'gnomAD_exomes_controls_AC,'
				'clinvar_review,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,gnomAD_exomes_controls_nhomalt,',
				'--plugin dbscSNV,%s,GRCh38,ada_score,rf_score',
				'-i',folder_in,'-o', folder_out]) % (dbNSFP38_gz,dbscSNV11_gz))

def vcf(lock,param,folder_temp,folder_in,folder_out):
	print('vcf')
	file_list = glob.glob(folder_in+'*')
	name_list = []
	for vcf in file_list:
		x = sample(vcf)
		name_list.append(x)
	for name in set(name_list):
		cds_samtools = pd.read_csv(join(folder_in,(sample(name)+'_samt_only_CDS.vcf')),sep='\t',header=0, encoding='utf8')
		cds_gatk = pd.read_csv(join(folder_in,(sample(name)+'_unfied_only_CDS.vcf')),sep='\t',header=0, encoding='utf8')
		all_samtools = pd.read_csv(join(folder_in,(sample(name)+'_samt_all.vcf')),sep='\t',header=0, encoding='utf8')
		all_gatk = pd.read_csv(join(folder_in,(sample(name)+'_unfied_all.vcf')),sep='\t',header=0, encoding='utf8')
		ALL_ = join(folder_in,name+'_temp_all_annot.csv')
		CDS_ = join(folder_temp,name+'_temp_CDS_annot.csv')
		all_samtools['refseq'] = all_samtools['refseq'].astype(str)
		all_samtools['hgmd'] = all_samtools['hgmd'].astype(str)
		all_samtools['GENE'] = all_samtools['GENE'].astype(str)
		all_gatk['refseq'] = all_gatk['refseq'].astype(str)
		all_gatk['hgmd'] = all_gatk['hgmd'].astype(str)
		all_gatk['GENE'] = all_gatk['GENE'].astype(str)

		CDS = pd.merge(cds_samtools,cds_gatk,on=['#CHROM','POS','ID','GENE',
								'GENE','exone','length','strand',
								'refseq','hgmd'],how='outer')
		ALL = pd.merge(all_samtools,all_gatk,on=['#CHROM','POS','ID','GENE',
								'GENE','exone','length','strand',
								'refseq','hgmd'],how='outer')
		samtools = name+'_x'
		gatk = name+'_y'

		CDS[name+'_samt'] = CDS[samtools]
		CDS[name+'_gatk'] = CDS[gatk]
		ALL[name+'_samt'] = ALL[samtools]
		ALL[name+'_gatk'] = ALL[gatk]
		CDS['QUAL_x'].fillna(-999,inplace=True)
		CDS['QUAL'] = CDS['QUAL_x'].astype(int)
		CDS['INFO'] = CDS['INFO_x']
		CDS['FILTER'] = CDS['FILTER_y']
		CDS.loc[CDS['FILTER'].isnull(),'FILTER'] = '.'
		CDS['REF'] = CDS['REF_y']
		CDS.loc[CDS['REF'].isnull(),'REF'] = CDS['REF_x']
		CDS['ALT'] = CDS['ALT_y']
		CDS.loc[CDS['ALT'].isnull(),'ALT'] = CDS['ALT_x']
		CDS['ALT'] = np.where(CDS['ALT'].str.split(',').str.get(0) == '*',
					CDS['ALT'].str.split(',').str.get(1),
					CDS['ALT'].str.split(',').str.get(0))
		CDS = CDS[CDS['ALT']!='*']
		maskalt = CDS['ALT'] != '*'
		CDS.loc[maskalt,'ALT'] = CDS['ALT']
		mask1 = CDS['INFO'].isnull()
		CDS.loc[mask1,'INFO'] = CDS['INFO_y']
		CDS['FORMAT'] = 'GT:PL:DP'
		CDS.drop(['QUAL_x','INFO_x','FORMAT_x',samtools,gatk], axis=1,inplace=True)
		ALL['QUAL_x'].fillna(-999,inplace=True)
		ALL['QUAL'] = ALL['QUAL_x'].astype(int)
		ALL['FILTER'] = ALL['FILTER_y']
		ALL.loc[ALL['FILTER'].isnull(),'FILTER'] = '.'
		ALL['INFO'] = ALL['INFO_x']
		ALL['REF'] = ALL['REF_y']
		ALL.loc[ALL['REF'].isnull(),'REF'] = ALL['REF_x']
		ALL['ALT'] = ALL['ALT_y']
		ALL.loc[ALL['ALT'].isnull(),'ALT'] = ALL['ALT_x']
		ALL['ALT'] = np.where(ALL['ALT'].str.split(',').str.get(0) == '*',
							ALL['ALT'].str.split(',').str.get(1),
							ALL['ALT'].str.split(',').str.get(0))
		ALL = ALL[ALL['ALT']!='*']

		maskalt = ALL['ALT'] != '*'
		ALL.loc[maskalt,'ALT'] = ALL['ALT']
		mask1 = ALL['INFO'].isnull()
		ALL.loc[mask1,'INFO'] = ALL['INFO_y']
		ALL['FORMAT'] = 'GT:PL:DP'
		ALL.drop(['QUAL_x','INFO_x','FORMAT_x',samtools,gatk], axis=1,inplace=True)
		CDS.drop(['QUAL_y','INFO_y','FORMAT_y'], axis=1,inplace=True)
		ALL.drop(['QUAL_y','INFO_y','FORMAT_y'], axis=1,inplace=True)
		ALL['ALT'] = ALL['ALT'].str.split(',').str.get(0)
		CDS['ALT'] = CDS['ALT'].str.split(',').str.get(0)
		samtools = name+'_samt'
		gatk = name+'_gatk'
		CDS['POS'] = CDS['POS'].astype(int)
		ALL['POS'] = ALL['POS'].astype(int)
		ALL['DEPTH'] = ALL['INFO'].str.split('DP=').str.get(1).str.split(';').str.get(0)
		# ALL['DEPTH'].fillna(-999,inplace=True)

		# ALL['DEPTH'] = np.where(ALL['INFO'].str.split(';').str.get(0) == 'INDEL',
		# 			ALL['INFO'].str.split(';').str.get(3).str.split('=').str.get(1),
		# 			ALL['INFO'].str.split(';').str.get(0).str.split('=').str.get(1))

		# print(ALL['DEPTH'])
		# ALL['DEPTH'] = ALL['DEPTH'].str.replace(",",'.')
		# print(ALL['DEPTH'])
		if ((param.panel == 'ALLGENE')|(param.panel == 'trusightone')|(param.panel == 'CANCER')):
			ALL = ALL[( pd.to_numeric(ALL['DEPTH'])>=10)]
			# ALL = ALL[( pd.to_numeric(ALL['DEPTH'])>=20)]
			ALL = ALL[ pd.to_numeric(ALL['QUAL'])>=18]
		else:
			ALL = ALL[( pd.to_numeric(ALL['DEPTH'])>=1)|( pd.to_numeric(ALL['DEPTH'])==1)]
			#ALL = ALL[ALL['QUAL']>=18]

		CDS = CDS[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',samtools,gatk,
					'GENE','exone','length','strand','refseq','hgmd']]
		ALL = ALL[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',samtools,gatk,
					'GENE','exone','length','strand','refseq','hgmd']]

		CDS.iloc[:,:].fillna('.',inplace=True)
		ALL.iloc[:,:].fillna('.',inplace=True)
		CDS[samtools].fillna('.',inplace=True)
		ALL[samtools].fillna('.',inplace=True)
		CDS[gatk].fillna('.',inplace=True)
		ALL[gatk].fillna('.',inplace=True)
		CDS.to_csv( CDS_, sep='\t', header=True, index=False,encoding='utf-8')
		ALL.to_csv(ALL_, sep='\t', header=True, index=False,encoding='utf-8')
		vcf_annotate_all = join(folder_out,name+'_all_annot.vcf')
		vcf_annotate_CDS = join(folder_out,name+'_CDS_annot.vcf')
		vcf_annotate_CDS2 = join(folder_out,name+'_CDS_annot_dbNSFP.vcf')

#####################ANNOTATION FUNCTION##########################
####		annotation1(param,ALL_,vcf_annotate_all)	   #
		print ('---->'+name+' ('+str(len(CDS))+')'+'<----')
		######annotation2(param, CDS_,vcf_annotate_CDS)
		annotation3(param, CDS_,vcf_annotate_CDS) #todo decommenta
#################################################################
	return
###################################################################################################################
###################################################################################################################
def final_annotation(param,folder_annot,folder_final,HGMD,EREDITA,folder_COV,APPRIS,ECCEZIONI):
	file_list = glob.glob(folder_annot+'*_CDS_annot.vcf')
	for name in set(file_list):

		sample_x = sample(name)
		sample_x = str(sample_x)
		# if sample_x=='PRE356.2021': #todo:commenta
		if True:
			COVERAGE = pd.read_csv(join(folder_COV,sample_x,sample_x+'_all'),sep='\t',header=0)
			COVERAGE = COVERAGE[['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum']]
			CDS = pd.read_csv(name,sep='\t',header=None,comment='#',
					names=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',
					sample_x+'_samt',sample_x+'_gatk','GENE','exone',
								'length','strand','refseq','hgmd'])
			more_annot2 = pd.merge(CDS,COVERAGE,on=['#CHROM','POS'],how='left')
			more_annot3 = pd.merge(more_annot2,EREDITA,on=['GENE'],how='left')
			result = more_annot3
			# print('hereeeeeeeeeeeeeee')
			# print(result.loc[result['POS']==38273375,'INFO'][351], 'INFO')
			# print(result.loc[result['POS']==38285768,'INFO'][352], 'INFO')
			# print(result.loc[result['POS']==38273375,'INFO'][351].str.split('|,'), 'INFO')
			# print(result.loc[result['POS']==38285768,'INFO'][352].str.split('|,'), 'INFO')
			result['INFO3'] = result['INFO'].str.split('|,').str.get(0)
			result.loc[result['POS']==38273375,'INFO3'] = result.loc[result['POS']==38273375,'INFO'].str.split('|').str.get(0) #todo

			INFOFIRST = pd.DataFrame(result['INFO'].str.split('CSQ=').str.get(0))
			INFOFIRST['INFOFIRST'] = INFOFIRST['INFO']
			INFOFIRST.drop('INFO',axis=1,inplace=True)

			INFO = result['INFO'].str.split('CSQ=').str.get(1)
			#INFO2 = pd.DataFrame(INFO.apply(lambda x: str(x).split('|,')))
			INFO2 = pd.DataFrame(INFO.apply(lambda x: str(x).split(',')))
			INFO2['INFO2'] = INFO2['INFO']
			INFO2.drop('INFO',axis=1,inplace=True)


			INFOFIRST.reset_index(inplace=True)
			INFO2.reset_index(inplace=True)
			result.reset_index(inplace=True)

			result2 = pd.merge(result,INFOFIRST,on='index',how='left')
			result3 = pd.merge(result2,INFO2,on='index',how='left')
			result3['PRINCIPAL'] = None

			REFSEQPRINCIPAL = APPRIS['refseq'].str.split('.').str.get(0)+'.'
			REFSEQPRINCIPALECCEZIONI = ECCEZIONI['refseq'].str.split('.').str.get(0)+'.'
			ECCEZIONI.dropna(subset=['refseq'],inplace=True)
			REFSEQPRINCIPALECCEZIONI.dropna(inplace=True)
			for index,row in result3.iterrows():
				for x in row.INFO2:
					for refseq in REFSEQPRINCIPALECCEZIONI:
						if refseq in x:
							#print (str(refseq))
							if (ECCEZIONI['GENE'][ECCEZIONI['refseq'].str.split('.').str.get(0)+'.' == refseq]==row['GENE']).item():
								result3.loc[index,'PRINCIPAL'] = str(row['INFOFIRST'])+'CSQ='+x
							else: pass
					for refseq in REFSEQPRINCIPAL:
						if refseq in x:
							if (APPRIS['GENE'][APPRIS['refseq'].str.split('.').str.get(0)+'.' == refseq]==row['GENE']).item():
								result3.loc[index,'PRINCIPAL'] = str(row['INFOFIRST'])+'CSQ='+x
							else: pass

			result3['PRINCIPAL2'] = np.where(result3['PRINCIPAL'].isnull(),
							result3['INFO3'],np.nan)

			result3['INFO'] = result3['PRINCIPAL']
			result3['INFO'].fillna(result3['PRINCIPAL2'],inplace=True)
			result3.drop('PRINCIPAL',axis=1,inplace=True)
			result3.drop('INFOFIRST',axis=1,inplace=True)
			result3.drop('INFO2',axis=1,inplace=True)
			result3.drop('INFO3',axis=1,inplace=True)
			result3.drop('PRINCIPAL2',axis=1,inplace=True)
			#result3[result3['INFO'].str.split('|').str.get(10).str.contains("XM_")]
			#result = result3
	####################################################################################
			if param.genome == 'geno37':
				final = adapt_annotation2(param,folder_final,sample_x,result3,HGMD)  ##
			elif param.genome == 'geno38':
				final = adapt_annotation2(param,folder_final,sample_x,result3,HGMD)  ##
	###################################################################################
	return

def adapt_annotation2(param,folder,sample,result,HGMD):
	name = str(sample)
	samtools = name+'_samt'
	gatk = name+'_gatk'
	cols = ['sample','HGVS','CHROM','POS','GENE','ALTGENE','ID','types','QUAL','DEPTH','mapquality','INHERITANCE',
		'samtools_geno','gatk_geno','unbalance','consequence','impact','clin_sign','strand',
		'refseq','hgmd','num_exon','num_intron','HGVS_c','HGVS_p','CDS_position',
		'amino_acids','codons','variation_cosmic','sift','polyphen','GMAF','ExAC_MAF','MAX_MAF','Adj_MAF','AFR_MAF','AMR_MAF','EAS_MAF',
		'EUR_MAF','SAS_MAF','AA_MAF','EA_MAF','pubmed',
		'hgmd_mutation','hgmd_function','hgmd_phenotype','hgmd_pubmed',

		'CADD_rankscore','DANN_rankscore','EigenPC_rankscore','FATHMM_rankscore','FATHMM_pred','GERP_rankscore','Interpro_domain',
		'LRT_rankscore','LRT_pred','MCAP_pred','MetaLR_pred','MetaLR_rankscore','MetaSVM_pred','MetaSVM_rankscore',
		'MutPred_rankscore','MutationAssessor_pred','MutationAssessor_rankscore','MutationTaster_rankscore','MutationTaster_pred',
		'PROVEAN_rankscore','PROVEAN_pred','Polyphen2HDIV_pred','Polyphen2HDIV_rankscore','Polyphen2HVAR_pred',
		'Polyphen2HVAR_rankscore','REVEL_rankscore','SIFT_rankscore','SIFT_pred','SiPhy29way_rankscore',
		'VEST3_rankscore','clinvar_clnsig','fathmmMKL_pred',
		'phastCons100way_rankscore','phastCons20way_rankscore','phyloP100way_rankscore','phyloP20way_rankscore',

		'ada_score','rf_score',

		'clinvar_MedGen_id','clinvar_OMIM_id','clinvar_Orphanet_id','clinvar_clnsig','clinvar_review',
		'gnomAD_exomes_POPMAX_AF','gnomAD_exomes_POPMAX_nhomalt','gnomAD_exomes_controls_AF',
		'gnomAD_exomes_controls_nhomalt','gnomAD_genomes_controls_nhomalt','gnomAD_genomes_POPMAX_AF','gnomAD_genomes_POPMAX_nhomalt',
		'gnomAD_genomes_controls_POPMAX_AF','gnomAD_genomes_controls_POPMAX_nhomalt','gnomAD_exomes_controls_AC',
		'gnomAD_genomes_controls_AC'
		]

######################################################################################################
	if len(result) != 0:
		for index,row in result.iterrows():
			try:
				# if ((row['G%'] != 0.0) & (row['ALT'] == 'G')):
				# 	result.loc[index,'unbalance'] = 'G='+unicode(row['G%'])
				# elif ((row['C%'] != 0.0) & (row['ALT'] == 'C')):
				# 	result.loc[index,'unbalance'] = 'C='+unicode(row['C%'])
				# elif ((row['T%'] != 0.0) & (row['ALT'] == 'T')):
				# 	result.loc[index,'unbalance'] = 'T='+unicode(row['T%'])
				# elif ((row['A%'] != 0.0) & (row['ALT'] == 'A')):
				# 	result.loc[index,'unbalance'] = 'A='+unicode(row['A%'])
				# if ((row['ins%'] != 0.0)):
				# 	result.loc[index,'unbalance'] = 'ins='+unicode(row['ins%'])
				# if ((row['del%'] != 0.0)):
				# 	result.loc[index,'unbalance'] = 'del='+unicode(row['del%'])
				if ((row['G%'] != 0.0) & (row['ALT'] == 'G')):
					result.loc[index,'unbalance'] = 'G='+str(row['G%'])
				elif ((row['C%'] != 0.0) & (row['ALT'] == 'C')):
					result.loc[index,'unbalance'] = 'C='+str(row['C%'])
				elif ((row['T%'] != 0.0) & (row['ALT'] == 'T')):
					result.loc[index,'unbalance'] = 'T='+str(row['T%'])
				elif ((row['A%'] != 0.0) & (row['ALT'] == 'A')):
					result.loc[index,'unbalance'] = 'A='+str(row['A%'])
				if ((row['ins%'] != 0.0)):
					result.loc[index,'unbalance'] = 'ins='+str(row['ins%'])
				if ((row['del%'] != 0.0)):
					result.loc[index,'unbalance'] = 'del='+str(row['del%'])


				result.loc[index,'DEPTH2'] = row['sum']
			except KeyError:
				print ('TROVA VARIANTE TERMINATO!!!')
######################################################################################################
		result['samtools'] = result[samtools].str.split(':').str.get(0)
		result['gatk'] = result[gatk].str.split(':').str.get(0)

		# result['DEPTH'] = np.where(result['INFO'].str.split(';').str.get(0) == 'INDEL',
		# 			result['INFO'].str.split(';').str.get(3).str.split('=').str.get(1),
		# 			result['INFO'].str.split(';').str.get(0).str.split('=').str.get(1))
		result['DEPTH'] = result['INFO'].str.split('DP=').str.get(1).str.split(';').str.get(0)

		result['mapquality'] = result['INFO'].str.split('CSQ=').str.get(0).str.split('MQ=').str.get(1).str.split(';').str.get(0)

		result['consequence'] = result['INFO'].str.split('|').str.get(1)
		result['impact'] = result['INFO'].str.split('|').str.get(2)
		result['symbol'] = result['INFO'].str.split('|').str.get(3)
		#result['transcript_id'] = result['INFO'].str.split('|').str.get(6)
		result['num_exon'] = result['INFO'].str.split('|').str.get(8)
		result['num_intron'] = result['INFO'].str.split('|').str.get(9)
		result['HGVS_c'] = result['INFO'].str.split('|').str.get(10)
		result['HGVS_p'] = result['INFO'].str.split('|').str.get(11)
		result['cDNA_position'] = result['INFO'].str.split('|').str.get(12)
		result['CDS_position'] = result['INFO'].str.split('|').str.get(13)
		result['protein_position'] = result['INFO'].str.split('|').str.get(14)
		result['amino_acids'] = result['INFO'].str.split('|').str.get(15)
		result['codons'] = result['INFO'].str.split('|').str.get(16)

		result['variation'] = result['INFO'].str.split('|').str.get(17).str.split('&').str.get(0)
		try: result['variation2'] = result['INFO'].str.split('|').str.get(17).str.split('&').str.get(1)
		except: result['variation2'] = ''

		result['variation'].fillna('',inplace=True)
		result['variation2'].fillna('',inplace=True)
		try:
			result['ID'] = np.where(result['variation'].str.contains('rs'),result['variation'],
                                        np.where(result['variation2'].str.contains('rs'),result['variation2'],np.nan))
		except:
			print ('ERROR IN NP.WHERE!!!!')
			result['ID'] = np.nan
		result['cosmic'] = result['INFO'].str.split('|').str.get(17).str.split('&').str.get(1)
		result['cosmic'].fillna('',inplace=True)

		try:
			result['variation_cosmic'] = np.where(result['cosmic'].str.contains('COSM'),result['cosmic'],np.nan)
		except:
			print ('ERROR COSMIC IN NP.WHERE!!!!')
			result['variation_cosmic'] = np.nan
		result['distance'] = result['INFO'].str.split('|').str.get(18)
		result['strand'] = result['INFO'].str.split('|').str.get(19)
		result['hgnc_id'] = result['INFO'].str.split('|').str.get(23)
		result['canonical'] = result['INFO'].str.split('|').str.get(24)

		result['gene_pheno'] = result['INFO'].str.split('|').str.get(38)
		result['sift'] = result['INFO'].str.split('|').str.get(39)
		result['polyphen'] = result['INFO'].str.split('|').str.get(40)
		result['domain'] = result['INFO'].str.split('|').str.get(41)
		result['GMAF'] = result['INFO'].str.split('|').str.get(44)
		result['AFR_MAF'] = result['INFO'].str.split('|').str.get(51)
		result['AMR_MAF'] = result['INFO'].str.split('|').str.get(52)
		result['EAS_MAF'] = result['INFO'].str.split('|').str.get(54)
		result['EUR_MAF'] = result['INFO'].str.split('|').str.get(56)
		result['SAS_MAF'] = result['INFO'].str.split('|').str.get(58)
		result['AA_MAF'] = result['INFO'].str.split('|').str.get(48)
		result['EA_MAF'] = result['INFO'].str.split('|').str.get(59)
		#result['ExAC_MAF'] = result['INFO'].str.split('|').str.get(47)
		#result['Adj_MAF'] = result['INFO'].str.split('|').str.get(48)
		result['ExAC_MAF'] = result['INFO'].str.split('|').str.get(50) #GnomadAF
		result['MAX_MAF'] = result['INFO'].str.split('|').str.get(70)  #MAX_MAF
		result['Adj_MAF'] = result['INFO'].str.split('|').str.get(71)  #MAX_AF_POPS
		result['clin_sign'] = result['INFO'].str.split('|').str.get(72)
		result['somatic'] = result['INFO'].str.split('|').str.get(73)
		result['PHENO'] = result['INFO'].str.split('|').str.get(74)
		result['pubmed'] = result['INFO'].str.split('|').str.get(75)
		result['pubmed'] = result['pubmed'].str.split('&').str.get(0)+' '+result['pubmed'].str.split('&').str.get(1)+' '+result['pubmed'].str.split('&').str.get(2)
		try:result['CADD_rankscore'] = result['INFO'].str.split('|').str.get(81).str.split(',').str.get(0)
		except: result['CADD_rankscore'] = result['INFO'].str.split('|').str.get(81)
		try: result['DANN_rankscore'] = result['INFO'].str.split('|').str.get(82).str.split(',').str.get(0)
		except: result['DANN_rankscore'] = result['INFO'].str.split('|').str.get(82)
		try: result['EigenPC_rankscore'] = result['INFO'].str.split('|').str.get(83).str.split(',').str.get(0)
		except: result['EigenPC_rankscore'] = result['INFO'].str.split('|').str.get(83)
		try: result['FATHMM_rankscore'] = result['INFO'].str.split('|').str.get(84).str.split(',').str.get(0)
		except: result['FATHMM_rankscore'] = result['INFO'].str.split('|').str.get(84)
		result['FATHMM_pred'] = result['INFO'].str.split('|').str.get(85)
		try: result['GERP_rankscore'] = result['INFO'].str.split('|').str.get(86).str.split(',').str.get(0)
		except: result['GERP_rankscore'] = result['INFO'].str.split('|').str.get(86)
		result['Interpro_domain'] = result['INFO'].str.split('|').str.get(87)
		try: result['LRT_rankscore'] = result['INFO'].str.split('|').str.get(88).str.split(',').str.get(0)
		except: result['LRT_rankscore'] = result['INFO'].str.split('|').str.get(88)
		result['LRT_pred'] = result['INFO'].str.split('|').str.get(89)
		result['MCAP_pred'] = result['INFO'].str.split('|').str.get(90)
		result['MetaLR_pred'] = result['INFO'].str.split('|').str.get(91)
		try: result['MetaLR_rankscore'] = result['INFO'].str.split('|').str.get(92).str.split(',').str.get(0)
		except:result['MetaLR_rankscore'] = result['INFO'].str.split('|').str.get(92)
		result['MetaSVM_pred'] = result['INFO'].str.split('|').str.get(93)
		try: result['MetaSVM_rankscore'] = result['INFO'].str.split('|').str.get(94).str.split(',').str.get(0)
		except:result['MetaSVM_rankscore'] = result['INFO'].str.split('|').str.get(94)
		try: result['MutPred_rankscore'] = result['INFO'].str.split('|').str.get(95).str.split(',').str.get(0)
		except:result['MutPred_rankscore'] = result['INFO'].str.split('|').str.get(95)
		result['MutationAssessor_pred'] = result['INFO'].str.split('|').str.get(96)
		try: result['MutationAssessor_rankscore'] = result['INFO'].str.split('|').str.get(97).str.split(',').str.get(0)
		except: result['MutationAssessor_rankscore'] = result['INFO'].str.split('|').str.get(97)
		try: result['MutationTaster_rankscore'] = result['INFO'].str.split('|').str.get(98).str.split(',').str.get(0)
		except:result['MutationTaster_rankscore'] = result['INFO'].str.split('|').str.get(98)
		result['MutationTaster_pred'] = result['INFO'].str.split('|').str.get(99)
		try: result['PROVEAN_rankscore'] = result['INFO'].str.split('|').str.get(100).str.split(',').str.get(0)
		except:result['PROVEAN_rankscore'] = result['INFO'].str.split('|').str.get(100)
		result['PROVEAN_pred'] = result['INFO'].str.split('|').str.get(101)
		result['Polyphen2HDIV_pred'] = result['INFO'].str.split('|').str.get(102)
		try: result['Polyphen2HDIV_rankscore'] = result['INFO'].str.split('|').str.get(103).str.split(',').str.get(0)
		except:result['Polyphen2HDIV_rankscore'] = result['INFO'].str.split('|').str.get(103)
		result['Polyphen2HVAR_pred'] = result['INFO'].str.split('|').str.get(104)
		try: result['Polyphen2HVAR_rankscore'] = result['INFO'].str.split('|').str.get(105).str.split(',').str.get(0)
		except:result['Polyphen2HVAR_rankscore'] = result['INFO'].str.split('|').str.get(105)
		try: result['REVEL_rankscore'] = result['INFO'].str.split('|').str.get(106).str.split(',').str.get(0)
		except:result['REVEL_rankscore'] = result['INFO'].str.split('|').str.get(106)
		try: result['SIFT_rankscore'] = result['INFO'].str.split('|').str.get(107).str.split(',').str.get(0)
		except:result['SIFT_rankscore'] = result['INFO'].str.split('|').str.get(107)
		result['SIFT_pred'] = result['INFO'].str.split('|').str.get(108)
		try: result['SiPhy29way_rankscore'] = result['INFO'].str.split('|').str.get(109).str.split(',').str.get(0)
		except:result['SiPhy29way_rankscore'] = result['INFO'].str.split('|').str.get(109)
		try: result['VEST3_rankscore'] = result['INFO'].str.split('|').str.get(110).str.split(',').str.get(0)
		except:result['VEST3_rankscore'] = result['INFO'].str.split('|').str.get(110)
		try: result['clinvar_MedGen_id'] = result['INFO'].str.split('|').str.get(111).str.split(',').str.get(0)
		except:result['clinvar_MedGen_id'] = result['INFO'].str.split('|').str.get(111)
		try: result['clinvar_OMIM_id'] = result['INFO'].str.split('|').str.get(112).str.split(',').str.get(0)
		except:result['clinvar_OMIM_id'] = result['INFO'].str.split('|').str.get(112)
		try: result['clinvar_Orphanet_id'] = result['INFO'].str.split('|').str.get(113).str.split(',').str.get(0)
		except:result['clinvar_Orphanet_id'] = result['INFO'].str.split('|').str.get(113)
		try: result['clinvar_clnsig'] = result['INFO'].str.split('|').str.get(114).str.split(',').str.get(0)
		except:result['clinvar_clnsig'] = result['INFO'].str.split('|').str.get(114)
		try: result['clinvar_review'] = result['INFO'].str.split('|').str.get(115).str.split(',').str.get(0)
		except:result['clinvar_review'] = result['INFO'].str.split('|').str.get(115)
		try: result['fathmmMKL_pred'] = result['INFO'].str.split('|').str.get(116).str.split(',').str.get(0)
		except: result['fathmmMKL_pred'] = result['INFO'].str.split('|').str.get(116)
		try: result['gnomAD_exomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(117).str.split(',').str.get(0)
		except:result['gnomAD_exomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(117)
		try: result['gnomAD_exomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(118).str.split(',').str.get(0)
		except:result['gnomAD_exomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(118)
		try: result['gnomAD_exomes_controls_AC'] = result['INFO'].str.split('|').str.get(119).str.split(',').str.get(0)
		except:result['gnomAD_exomes_controls_AC'] = result['INFO'].str.split('|').str.get(119)
		try: result['gnomAD_exomes_controls_AF'] = result['INFO'].str.split('|').str.get(120).str.split(',').str.get(0)
		except:result['gnomAD_exomes_controls_AF'] = result['INFO'].str.split('|').str.get(120)
		try: result['gnomAD_exomes_controls_nhomalt'] = result['INFO'].str.split('|').str.get(121).str.split(',').str.get(0)
		except:result['gnomAD_exomes_controls_nhomalt'] = result['INFO'].str.split('|').str.get(121)
		try: result['gnomAD_genomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(122).str.split(',').str.get(0)
		except:result['gnomAD_genomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(122)
		try: result['gnomAD_genomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(123).str.split(',').str.get(0)
		except:result['gnomAD_genomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(123)
		try: result['gnomAD_genomes_controls_AC'] = -999 #result['INFO'].str.split('\|').str.get(115).str.split(',').str.get(0)
		except:result['gnomAD_genomes_controls_AC'] =  -999 #result['INFO'].str.split('\|').str.get(115)
		try: result['gnomAD_genomes_controls_POPMAX_AF'] =  -999 #result['INFO'].str.split('\|').str.get(116).str.split(',').str.get(0)
		except:result['gnomAD_genomes_controls_POPMAX_AF'] =  -999 #result['INFO'].str.split('\|').str.get(116)
		try: result['gnomAD_genomes_controls_POPMAX_nhomalt'] =  -999 #result['INFO'].str.split('\|').str.get(117).str.split(',').str.get(0)
		except:result['gnomAD_genomes_controls_POPMAX_nhomalt'] =  -999 #result['INFO'].str.split('\|').str.get(117)
		try: result['gnomAD_genomes_controls_nhomalt'] =  -999 #result['INFO'].str.split('\|').str.get(118).str.split(',').str.get(0)
		except:result['gnomAD_genomes_controls_nhomalt'] =  -999#result['INFO'].str.split('\|').str.get(118)
		try: result['phastCons100way_rankscore'] = result['INFO'].str.split('\|').str.get(124).str.split(',').str.get(0)
		except:result['phastCons100way_rankscore'] = result['INFO'].str.split('\|').str.get(124)
		try: result['phastCons20way_rankscore'] = result['INFO'].str.split('\|').str.get(125).str.split(',').str.get(0)
		except:result['phastCons20way_rankscore'] = result['INFO'].str.split('\|').str.get(125)
		try: result['phyloP100way_rankscore'] = result['INFO'].str.split('\|').str.get(126).str.split(',').str.get(0)
		except:result['phyloP100way_rankscore'] = result['INFO'].str.split('\|').str.get(126)
		try: result['phyloP20way_rankscore'] = result['INFO'].str.split('\|').str.get(127).str.split(',').str.get(0)
		except:result['phyloP20way_rankscore'] = result['INFO'].str.split('\|').str.get(127)
		try: result['ada_score'] = result['INFO'].str.split('\|').str.get(128).str.split(',').str.get(0)
		except:result['ada_score'] = result['INFO'].str.split('\|').str.get(128)
		try: result['rf_score'] = result['INFO'].str.split('\|').str.get(129).str.split(',').str.get(0)
		except:result['rf_score'] = result['INFO'].str.split('\|').str.get(129)

		mask1a = result['samtools'] == '0/1'
		mask1a2 = result['samtools'] =='1/2'
		mask1b = result['gatk'] == '0/1'
		mask2a = result['samtools'] == '0/0'
		mask2b = result['gatk'] == '0/0'
		mask3a = result['samtools'] == '1/1'
		mask3b = result['gatk'] == '1/1'
		result.loc[mask1a, 'samtools_geno'] = 'het'
		result.loc[mask1a2, 'samtools_geno'] = 'het'
		result.loc[mask1b, 'gatk_geno'] = 'het'
		result.loc[mask2a, 'samtools_geno'] = 'homo_wild'
		result.loc[mask2b, 'gatk_geno'] = 'homo_wild'
		result.loc[mask3a, 'samtools_geno'] = 'homo'
		result.loc[mask3b, 'gatk_geno'] = 'homo'

		HGMD['POS'] = HGMD['END']
		result2 = pd.merge(result,HGMD,on=['#CHROM','POS'],how='left')
		result2.drop_duplicates(subset=['#CHROM','POS'],inplace=True)
		result2['hgmd'] = result2['hgmd_x']
		try:
			result2['hgmd_mutation'] = result2['hgmd_y'].str.split('|').str.get(2)
			result2['hgmd_function'] = result2['hgmd_y'].str.split('|').str.get(3)
			result2['hgmd_phenotype'] = result2['hgmd_y'].str.split('|').str.get(4)
			result2['hgmd_pubmed'] = result2['hgmd_y'].str.split('|').str.get(5)
		except AttributeError:
			result2['hgmd_mutation'] = ''
			result2['hgmd_function'] = ''
			result2['hgmd_phenotype'] = ''
			result2['hgmd_pubmed'] = ''
			print ('AttributeError!!!')

		result2.drop('hgmd_y',axis=1,inplace=True)
		result2.drop('hgmd_x',axis=1,inplace=True)
		result2['length'] = result2['length'].map('{:,.0f}'.format)

		for index, row in result2.iterrows():
			count_ref = len(row['REF'])
			count_alt = len(row['ALT'])
			result2.loc[index,'count_ref'] = int(count_ref)
			result2.loc[index,'count_alt'] = int(count_alt)
			
		result2['START'] = np.where((result2['count_ref']>1),
					result2['POS']+1,result2['POS'])

		result2['END'] = np.where((result2['count_ref'] > result2['count_alt']),
					result2['POS']+(result2['count_ref']-1),
					np.where((result2['count_ref'] < result2['count_alt']),
					result2['START']+((result2['count_alt']-1)-(result2['count_ref']-1)),
					result2['POS']))

		result2['count_ref'] = result2['count_ref'].astype(int)
		result2['count_alt'] = result2['count_alt'].astype(int)
		for index, row in result2.iterrows():
			if (row['count_ref'] > row['count_alt'])&(row['count_ref']>1):
				row['REF_2'] = row['REF'][row['count_alt']:]
				row['ALT_2'] = '-'
				result2.loc[index,'REF_2'] = row['REF_2']
				result2.loc[index,'ALT_2'] = row['ALT_2']
				result2.loc[index,'types'] = 'DELETION'
			elif (row['count_ref'] < row['count_alt'])&(row['count_alt']>1):
				row['REF_2'] = '-'
				row['ALT_2'] = row['ALT'][row['count_ref']:]
				result2.loc[index,'REF_2'] = row['REF_2']
				result2.loc[index,'ALT_2'] = row['ALT_2']
				result2.loc[index,'types'] = 'INSERTION'
			elif (row['count_ref'] == row['count_alt'])&(row['count_ref'] == 1):
				row['REF_2'] = row['REF']
				row['ALT_2'] = row['ALT']
				result2.loc[index,'REF_2'] = row['REF_2']
				result2.loc[index,'ALT_2'] = row['ALT_2']
				result2.loc[index,'types'] = 'SVN'

		result2['START'] = result2['START'].astype(int)
		result2['END'] = result2['END'].astype(int)
		result2['HGVS'] = result2['#CHROM']+':'+result2['START'].astype(str)+'-'+result2['END'].astype(str)+':'+result2['REF_2']+'/'+result2['ALT_2']

		result2['HGVS_c'].fillna('',inplace=True)
		result2['HGVS_p'].fillna('',inplace=True)


		for index, row in result2.iterrows():
			if row['strand'] == str(-1):
				try: seq_x = Seq(row['REF_2']) #, IUPAC.unambiguous_dna)
				except: seq_x = Seq(row['REF']) #, IUPAC.unambiguous_dna)
				try: seq_y = Seq(row['ALT_2']) #, IUPAC.unambiguous_dna)
				except: seq_x = Seq(row['ALT']) #, IUPAC.unambiguous_dna)
				seq_ref = seq_x.reverse_complement()
				seq_alt = seq_y.reverse_complement()
				row['REF_2'] = seq_ref
				row['ALT_2'] = seq_alt
			#try:
			row['HGVS_c'] = row['HGVS_c'].replace('N>',str(row['REF_2']+'>'))
			DEL = 'del'+'N'*(len(row['REF_2']))
			row['HGVS_c'] = row['HGVS_c'].replace(DEL,'del'+str(row['REF_2']))
			#except:  #row['HGVS_c'] = ''
			result2.loc[index,'HGVS_c'] = row['HGVS_c']

		mask1 = result2['HGVS_p'].str.contains('p.%3D')
		result2.loc[mask1,'HGVS_p'] = '(p.%3D)'

		result2['HGVS_c'] = result2['HGVS_c'].replace(r'NN*','',regex=True)
		result2['HGVS_c'] = result2['HGVS_c'].replace(r'delN+','',regex=True)
		result2['HGVS_c'] = result2['HGVS_c'].replace(r'del-','',regex=True)
		result2['HGVS_c'] = result2['HGVS_c'].replace(r'^M_','NM_',regex=True)
		result2['HGVS_c'] = result2['HGVS_c'].replace(r'^R_','NR_',regex=True)

		result2['refseq'].fillna(result2['HGVS_c'].str.split('.').str.get(0),inplace=True)
		result2['variation_cosmic'].fillna('RISN',inplace=True)

		mask1 = result2['variation_cosmic'].str.contains('RISN')
		result2.loc[mask1,'variation_cosmic'] = ''
		mask2 = result2['variation_cosmic'].str.contains('rs')
		result2.loc[mask2,'variation_cosmic'] = ''
		mask3 = result2['variation_cosmic'].str.contains('RPGR')
		result2.loc[mask3,'variation_cosmic'] = ''

		result2['hgmd_mutation'].replace('null;null','',inplace=True)
		result2['CHROM'] = result2['#CHROM']
		result2['ALTGENE'] = result2['GENE']
		result2['GENE'] = result2['symbol']

		#result2['GENE'] = result2['GENE'].fillna(result2['symbol'])
		#result2['GENE'] = result2['symbol']

		result2['sample'] = str(name)
		#result2b = result2[result2['QUAL'] >= 18]
		result2b = result2
		result2b['hgmd_mutation'].fillna('',inplace=True)
		result2b.sort_values(by=['#CHROM','POS'],ascending=[True,True],inplace=True)

		result2b['DEPTH'] = result2b['DEPTH'].str.split(',').str.get(0)
		result2b['DEPTH2'].fillna(result2b['DEPTH'],inplace=True)
		result2b['DEPTH2'].fillna(-999,inplace=True)
		result2b['DEPTH'] = result2b['DEPTH2'].astype(int)

		result2b['mapquality'].fillna(-999,inplace=True)
		result2b['mapquality2'] = result2b['mapquality'].astype(float)
		result2b['mapquality'] = result2b['mapquality2'].astype(int)
		result2b['GENE'].replace('DFNB31','WHRN',inplace=True)
		result3 = result2b[cols]
		final = join(folder,name+'_CDS_annot.csv')
		result3.to_csv(final,sep='\t',index=False,encoding='utf-8')
		print ('Annotation CDS len:', len(result3),'->',str(name))

	else:
		result3 = pd.DataFrame(columns=cols)
		final = join(folder,name+'_CDS_annot.csv')
		result3.to_csv(final,sep='\t',index=False,encoding='utf-8')
		print ('Annotation CDS len:', len(result3),'->',str(name))

	return result3


def filter_disease(param,phenotype,folder_name,folder_final,folder_annotation):
	file_list = glob.glob(folder_final+'*_CDS_annot.csv')
	for name in set(file_list):
		CDS = pd.read_csv(name,sep='\t',header=0)
#########################################################
		if ((param.panel == 'ALLGENE')|(param.panel == 'trusightone')|(param.panel == 'CANCER')):
			CDS = CDS[(CDS['DEPTH']>=10)]
			# CDS = CDS[(CDS['DEPTH']>=20)]
			CDS = CDS[CDS['QUAL']>=18]
		else:
			CDS = CDS[(CDS['DEPTH']>=10)] #|(CDS['DEPTH']==1)]
#########################################################
		sample_x = sample(name)
		z = str(sample_x).split("_")
		sample_ = str(z[0])
		MISSING = join(folder_name,'pheno','MISSING_'+sample_)
		try:
			_missing_ = pd.read_csv(MISSING,sep='\t',header=0)
		except:
			_missing_ = None
			print ('FILE MISSING NON ESISTE!!!!',sample_)

		phenotype_ = phenotype[phenotype['sample'].astype(str) == sample_]
		a = phenotype_[['malattia','gene']].drop_duplicates()
		x1 = pd.DataFrame(a['gene'])
		if _missing_ is not None:
			for gene in a['gene']:
				for index, gene_miss in _missing_.iterrows():
					if str(gene) == str(gene_miss.MISSING):
						x2= pd.DataFrame({'gene':pd.Series([str(str(gene_miss.HGNC_SYMBOL))])})
						x1 = x1.append(x2)
					else: x1=x1

		b = CDS
		b['strand'].fillna(1,inplace=True)
		b['strand'] = b['strand'].astype(int)
		b.sort_values(by=['GENE'],inplace=True)

		b['sample'] = sample_

		result = join(folder_final,sample_+'_pheno_annot.csv')
		# if param.dest == 'r':
		# 	dest = 'rovereto'
		# 	path_django = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/annot'
		# elif param.dest == 'b':
		# 	dest = 'bolzano'
		# 	path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/annot'
		# elif param.dest == 's':
		# 	dest = 'sanfelice'
		# 	path_django = '/home/magi/VIRTUAL/SANFELICE/NGS_RESULT/annot'
		# elif param.dest == 'z':
		# 	dest = 'ricerca'
		# 	path_django = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/annot'
		# elif param.dest == 'p':
		# 	dest = 'privato'
		# 	path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/annot'
		#
		# result_django = join(path_django,sample_+'_pheno_annot.csv')

		#b.to_csv(result_django,sep='\t',index=False,encoding='utf-8')
		b.to_csv(result,sep='\t',index=False,encoding='utf-8')
		print ('Len disease annot: ',len(b),'->',str(sample_))

		# if param.dest == 'r':
		# 	dest = 'rovereto'
		# 	django_folder = '/home/magi/VIRTUAL/MAGIS/MAGIS/templates/ngs/stat/%s/' % (str(sample_))
		# elif param.dest == 'b':
		# 	dest = 'bolzano'
		# 	django_folder = '/home/magi/VIRTUAL/EUREGIO/EUREGIO/templates/ngs/stat/%s/' % (str(sample_))
		# elif param.dest == 's':
		# 	dest = 'sanfelice'
		# 	django_folder = '/home/magi/VIRTUAL/SANFELICE/SANFELICE/templates/ngs/stat/%s/' % (str(sample_))
		# elif param.dest == 'z':
		# 	dest = 'ricerca'
		# 	django_folder = '/home/magi/VIRTUAL/RICERCA/RICERCA/templates/ngs/stat/%s/' % (str(sample_))
		#
		# if not os.path.exists(django_folder):
		# 	os.makedirs(django_folder)
		# system(' '.join(['cp',join(folder_annotation,sample_+'*.html'), django_folder]))
	return


if __name__=="__main__":
	lock = Lock()
	args=InputPar()

	folder_name = step.principal_folder(args,step.create_folder(),over='False')
	input_phenotype = join(folder_name,'pheno/phenotype')

	input_vcf =  join(folder_name,'vcf/')
	output_annotation = join(folder_name,'annotation/')
	output_final = join(folder_name,'final/')
	folder_temp = join(folder_name,'temp/')
	folder_COV = join(folder_name,'coverage/')

	print_args=vars(args) #take args in a dict and print in a log
	for k,v in print_args.items():
		strategy_out='='.join([str(k),str(v)])
		path = join(folder_name,'log')
		PrintLog(strategy_out,path)
#############################################################
	vcf(lock,args,folder_temp,input_vcf,output_annotation)
#############################################################
	if args.genome == 'geno37':
		eredita = pd.read_csv(eredita_37,sep='\t',header=0)
		hgmd = pd.read_csv(HGMD_37,sep='\t',header=None,names=['#CHROM','START','END','hgmd'])
#########################################################################################################
		final = final_annotation(args,output_annotation,output_final,hgmd,eredita,folder_COV)
		phenotype = pd.read_csv(input_phenotype, sep='\t', header=0)
#######################################################################################################
		filter_disease(args,phenotype,folder_name,output_final,output_annotation)

	elif args.genome == 'geno38':
		print('annot')
		eredita = pd.read_csv(eredita_38,sep='\t',header=0)
		hgmd = pd.read_csv(HGMD_38,sep='\t',header=None,names=['#CHROM','START','END','hgmd'],encoding='latin')
		APPRIS = pd.read_csv(appris,sep='\t',header=0)
		ECCEZIONI = pd.read_csv(eccezioni,sep='\t',header=0)
#####################################################################################################################
		final = final_annotation(args,output_annotation,output_final,hgmd,eredita,folder_COV,APPRIS,ECCEZIONI)	##
		phenotype = pd.read_csv(input_phenotype, sep='\t', header=0,encoding='utf-8')			##
#####################################################################################################################
		filter_disease(args,phenotype,folder_name,output_final,output_annotation)