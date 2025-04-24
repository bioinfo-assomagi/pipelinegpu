#!/home/magi/miniconda3/envs/PY310/bin/python
#10/04/2023

# -*- coding: utf-8 -*-
"""
Created on Wed AUG 18 2021
@author: Elisa Sorrentino
"""
import pandas as pd
import os
from os import listdir, system
from os.path import isfile, join
import time
import numpy as np
from bioservices import HGNC
import re
import config

import mygene
#from pyGeno.Genome import *
import myvariant

gtf_path = r'/home/magi/PROJECT/bedfiles/mane_phase14_rs.ucsc_seqids.gtf'
#######################################################################
#hgmd = join('/home/magi/PROJECT/diagnosys/bin/geneToNM_hgmd.txt')
#hgmd_df = pd.read_csv(hgmd,sep='\t',header=0)
#hgmd_df['hgmd'] = 'hgmd'
########################################################################
#hgmd = join('/home/magi/', 'PROJECT/somatic/bin/APPRIS_PRINCIPAL')
#hgmd_df = pd.read_csv(hgmd,sep='\t',header=0)
#hgmd_df['refseq2'] = hgmd_df['refseq'].str.split('.').str.get(0)
#hgmd_df['refseq'] = hgmd_df['refseq2']
#hgmd_df['hgmd'] = 'hgmd'
#hgmd_df = hgmd_df[['GENE','refseq','hgmd']]
###############################################################
hgmd = join('/home/magi/', 'PROJECT/diagnosys/bin/APPRIS_PRINCIPAL')
hgmd_df = pd.read_csv(hgmd,sep='\t',header=0)
hgmd_df['refseq2'] = hgmd_df['refseq'].str.split('.').str.get(0)
hgmd_df['refseq'] = hgmd_df['refseq2']
hgmd_df['APPRIS2'] = hgmd_df['APPRIS'].str.split('.').str.get(0)
hgmd_df['APPRIS'] = hgmd_df['APPRIS2']
hgmd_df['hgmd'] = 'hgmd'
hgmd_df = hgmd_df[['GENE','refseq','hgmd']]
######################################################################
eccezioni = join('/home/magi/', 'PROJECT/diagnosys/bin/APPRIS_ECCEZIONI')
eccezioni_df = pd.read_csv(eccezioni,sep='\t',header=0)
eccezioni_df['refseq2'] = eccezioni_df['refseq'].str.split('.').str.get(0)
eccezioni_df['refseq'] = eccezioni_df['refseq2']
#eccezioni_df['APPRIS2'] = eccezioni_df['APPRIS'].str.split('.').str.get(0)
#eccezioni_df['APPRIS'] = eccezioni_df['APPRIS2']
eccezioni_df['hgmd'] = 'other'
eccezioni_df = eccezioni_df[['GENE','refseq','hgmd']]
########################################################################
chainfrom37to38 = join('/home/magi/dataset/CHAIN_LIFTOVER','hg19ToHg38.over.chain.gz')
chainfrom38to37 = join('/home/magi/dataset/CHAIN_LIFTOVER','hg38ToHg19.over.chain.gz')
mg = mygene.MyGeneInfo()
from pyliftover import LiftOver
lo = LiftOver(chainfrom38to37)
lo37 = LiftOver(chainfrom37to38)
#######################################################################################################################
CHROM = []
hgvs38 = []
hgvs37 = []
START = []
END = []
START37 = []
END37 = []
geno = []
exone = []
tipo = []
gene = []
STRAND = []
variazione = []
annotazione = []
CDSSTART=[]
CDSEND=[]
refseq = []
INFO =[]
#######################################################################################################################
def create_null_coverage(x,y):
	RANGE= pd.DataFrame({'POS':range(x,y+1,1)})
	return RANGE

def create_vertical2(bed_df):
	bed_df['POS'] = [[x for x in range(bed_df.at[i, 'START'], bed_df.at[i, 'END'] + 1)] for i in bed_df.index]
	vertical = bed_df.drop(['START', 'END'], axis=1).explode('POS').reset_index()
	vertical.drop_duplicates(subset=['#CHROM','POS', 'GENE'], inplace=True)
	return vertical[['#CHROM','POS','GENE','exone','length','strand','refseq','hgmd']]


def create_vertical(DATA,DATAX,FILENAME,FOLDER):
	df1 = pd.DataFrame()
	print ('START VERTICAL!!!'+str(FILENAME))
	for row_index,row in DATA.iterrows():
		START = int(row['START'])
		END = int(row['END'])
		gene = row['GENE']
		chrom = row['#CHROM']
		exone = row['exone']
		length = row['length']
		strand = row['strand']
		refseq = row['refseq']
		hgmd = 	row['hgmd']

		result = create_null_coverage(START,END)
		result['GENE'] = gene
		result['#CHROM'] = chrom
		result['exone'] = exone
		result['length'] = length
		result['strand'] = strand
		result['refseq'] = refseq
		result['hgmd'] = hgmd

		df1 = pd.concat([df1,result])
	df1.POS = df1['POS'].astype(int)
	df1.length = df1['length'].astype(int)
	df1.strand = df1['strand'].astype(int)
	df1.exone = df1['exone'].astype(int)
	vertical = df1.sort_values(by=['#CHROM','POS','hgmd'],ascending=[True,True,True])
	vertical.drop_duplicates(subset=['#CHROM','POS','GENE'],inplace=True)
	vertical = vertical[['#CHROM','POS','GENE','exone','length','strand','refseq','hgmd']]

	df2 = pd.DataFrame()
	print ('START VERTICALX!!!'+str(FILENAME))
	for row_index,rowx in DATAX.iterrows():
		START = int(rowx['START'])
		END = int(rowx['END'])
		gene = rowx['GENE']
		chrom = rowx['#CHROM']
		exone = rowx['exone']
		length = rowx['length']
		strand = rowx['strand']
		refseq = rowx['refseq']
		hgmd = 	rowx['hgmd']

		result = create_null_coverage(START,END)
		result['GENE'] = gene
		result['#CHROM'] = chrom
		result['exone'] = exone
		result['length'] = length
		result['strand'] = strand
		result['refseq'] = refseq
		result['hgmd'] = hgmd

		df2 = pd.concat([df2,result])
  
	df2.POS = df2['POS'].astype(int)
	df2.length = df2['length'].astype(int)
	df2.strand = df2['strand'].astype(int)
	df2.exone = df2['exone'].astype(int)
	verticalX = df2.sort_values(by=['#CHROM','POS','hgmd'],ascending=[True,True,True])
	verticalX.drop_duplicates(subset=['#CHROM','POS','GENE'],inplace=True)
	verticalX = verticalX[['#CHROM','POS','GENE','exone','length','strand','refseq','hgmd']]

	print ('END VERTICAL!!!'+str(FILENAME))
	return vertical,verticalX

def liftover_from37_to38(CHROM,POS):
	lift = lo37.convert_coordinate(CHROM,POS)
	return lift[0][1]

def liftover_from38_to37(CHROM,POS):
	lift = lo.convert_coordinate(CHROM,POS)
	return lift[0][1]
############################################################################################################################
#############################################################################################################################
def cutCDS(GENE,sample_obj,FOLDER,DEST):
			SAMPLE = sample_obj.name # TODO: fix naming

			# lastFOLDER = FOLDER.split('/')[-2]
			#print lastFOLDER
			# lastFolder is equal to principal_directory name
			# if DEST == 'r':
			# 	dest = 'ROVERETO'
			# 	download_folder = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % (dest,lastFOLDER)
			# 	if not os.path.exists(download_folder):
			# 		os.makedirs(download_folder)
			# elif DEST == 'b':
			# 	dest = 'BOLZANO'
			# 	download_folder = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % (dest,lastFOLDER)
			# 	if not os.path.exists(download_folder):
			# 		os.makedirs(download_folder)
			# elif DEST == 's':
			# 	dest = 'SANFELICE'
			# 	download_folder = '/home/magi/VIRTUAL/SANFELICE/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % (dest,lastFOLDER)
			# 	if not os.path.exists(download_folder):
			# 		os.makedirs(download_folder)
			# elif DEST == 'z':
			# 	dest = 'RICERCA'
			# 	download_folder = '/home/magi/VIRTUAL/RICERCA/DOWNLOADS/NGSINFO/BED_INFO/%s/' % (lastFOLDER)
			# 	if not os.path.exists(download_folder):
			# 		os.makedirs(download_folder)

			mg = ''
			mv = ''
			mv = myvariant.MyVariantInfo()
			mg = mygene.MyGeneInfo()
			a = []

			CHROM = []
			hgvs38 = []
			hgvs37 = []
			START = []
			END = []
			START37 = []
			END37 = []
			geno = []
			exone = []
			tipo = []
			gene = []
			STRAND = []
			variazione = []
			annotazione = []
			CDSSTART=[]
			CDSEND=[]
			refseq = []
			INFO =[]
			ufficial_gene=[]
			snplist = []
			a = []
			TABLE = ''
			gene_DF = pd.DataFrame()
			PREsnp_DF3 = pd.DataFrame()
			TABLE_SNP = pd.DataFrame()

			a = mg.querymany(list(GENE), scopes="symbol",fields=["symbol","exons"], species="human", as_dataframe=True,returnall=True)
			try:
				gene_DF1 = a['out']
				gene_DF1['symbol'].fillna('ciao',inplace=True)
				gene_DF1 = gene_DF1[gene_DF1['symbol']!='ciao']
				gene_DF1['info'] =  'OfficialName'
			except: 
				gene_DF1 = pd.DataFrame()

			# Write missing
			file_ = (join(FOLDER,'MISSING_'+str(SAMPLE)))
			file_download = (join(config.DOWNLOAD_FOLDER, 'MISSING_'+str(SAMPLE)))
			missing = a['missing']
			out_file = open(file_,"w")
			out_file_download = open(file_download,"w")
			out_file.write('MISSING'+'\t'+'HGNC_SYMBOL'+'\n')
			out_file_download.write('MISSING'+'\t'+'HGNC_SYMBOL'+'\n')

			regex=re.compile('^rs\d*')
			regex2=re.compile('^chr\d*:g\.\d*:\w*>\w*')
			regex2B=re.compile('^chr\d*:g\.\d+')
			regex3=re.compile('^NM_\d*\.?(\d+)?:c.\d*\w*>\w*')
			regex4=re.compile('^chr\d*:g\.\d*_\d*:\w*')
			regex4B=re.compile('^chr\d*:g\.\d+_\d+')
			regex5=re.compile('^chrX:g\.\d*:\w*>\w*')
			regex5B=re.compile('^chrX:g\.\d+')

			# Deal with missing 
			if len(missing)>0:
			#try:
				for index,row in missing.iterrows():
					genemiss = row['query']
					#print (genemiss)
					if re.match(regex, str(genemiss)):
						print ('TROVATOOOOO11111111!!!!')
						snplist.append(str(genemiss))
					elif re.match(regex2, str(genemiss)):
						print ('TROVATOOOOO22222222AAAAAAAAAAA!!!!')
						_chr_=genemiss.split(':')[0]
						_pos38_=genemiss.split(':')[1][2:]
						_geno_=genemiss.split(':')[2]
						PREsnp_DF3.loc[index,'#CHROM'] = _chr_
						PREsnp_DF3.loc[index,'START'] = int(int(_pos38_)-1)
						PREsnp_DF3.loc[index,'END'] = int(int(_pos38_)+1)
						POS37=liftover_from38_to37(_chr_,int(_pos38_))
						genemiss37 = str(_chr_)+':g.'+str(POS37)+str(_geno_)
						snplist.append(str(genemiss37))
					elif re.match(regex2B, str(genemiss)):
						print (genemiss,'TROVATOOOOO22222BBBBBBBBBBBBBB!!!!')
						_chr_=genemiss.split(':')[0]
						_pos38_=genemiss.split(':')[1][2:].split('_')[0]
						_geno_='' #genemiss.split(':')[2]
						PREsnp_DF3.loc[index,'#CHROM'] = _chr_
						PREsnp_DF3.loc[index,'START'] = int(int(_pos38_)-1)
						PREsnp_DF3.loc[index,'END'] = int(int(_pos38_)+1)
						POS37=liftover_from38_to37(_chr_,int(_pos38_))
						genemiss37 = str(_chr_)+':g.'+str(POS37)+str(_geno_)
						snplist.append(str(genemiss37))
					elif re.match(regex4, str(genemiss)):
						print ('TROVATOOOO44444444444444!!!!')
						_chr_=genemiss.split(':')[0]
						_pos38START_=genemiss.split(':')[1][2:].split('_')[0]
						_pos38END_=genemiss.split(':')[1][2:].split('_')[1]
						_geno_=genemiss.split(':')[2]
						PREsnp_DF3.loc[index,'#CHROM'] = _chr_
						PREsnp_DF3.loc[index,'START'] = int(int(_pos38START_)-1)
						PREsnp_DF3.loc[index,'END'] = int(int(_pos38END_)+1)
						POS37=liftover_from38_to37(_chr_,int(_pos38START_))
						genemiss37 = str(_chr_)+':g.'+str(POS37)+str(_geno_)
						snplist.append(str(genemiss37))
					elif re.match(regex4, str(genemiss)):
						print ('TROVATOOOO4444444444BBBBBBBBBBBBBB!!!!')
						_chr_=genemiss.split(':')[0]
						_pos38START_=genemiss.split(':')[1][2:].split('_')[0]
						_pos38END_=genemiss.split(':')[1][2:].split('_')[1]
						_geno_='' #genemiss.split(':')[2]
						PREsnp_DF3.loc[index,'#CHROM'] = _chr_
						PREsnp_DF3.loc[index,'START'] = int(int(_pos38START_)-1)
						PREsnp_DF3.loc[index,'END'] = int(int(_pos38END_)+1)
						POS37=liftover_from38_to37(_chr_,int(_pos38START_))
						genemiss37 = str(_chr_)+':g.'+str(POS37)+str(_geno_)
						#print (genemiss37)
						snplist.append(str(genemiss37))
					elif re.match(regex5, str(genemiss)):
						print ('TROVATOOOOO55555555555555555!!!!')
						_chr_=genemiss.split(':')[0]
						_pos38_=genemiss.split(':')[1][2:]
						_geno_=genemiss.split(':')[2]
						PREsnp_DF3.loc[index,'#CHROM'] = _chr_
						PREsnp_DF3.loc[index,'START'] = int(int(_pos38_)-1)
						PREsnp_DF3.loc[index,'END'] = int(int(_pos38_)+1)
						POS37=liftover_from38_to37(_chr_,int(_pos38_))
						genemiss37 = str(_chr_)+':g.'+str(POS37)+str(_geno_)
						#print (genemiss37)
						snplist.append(str(genemiss37))
					elif re.match(regex5B, str(genemiss)):
						print (genemiss,'TROVATOOOOO55555555BBBBBBBBBBB!!!!')
						_chr_=genemiss.split(':')[0]
						_pos38_=genemiss.split(':')[1][2:].split('_')[0]
						_geno_='' #genemiss.split(':')[2]
						PREsnp_DF3.loc[index,'#CHROM'] = _chr_
						PREsnp_DF3.loc[index,'START'] = int(int(_pos38_)-1)
						PREsnp_DF3.loc[index,'END'] = int(int(_pos38_)+1)
						POS37=liftover_from38_to37(_chr_,int(_pos38_))
						genemiss37 = str(_chr_)+':g.'+str(POS37)+str(_geno_)
						snplist.append(str(genemiss37))
					else:
						#print (genemiss,'AAAAAAAAAAAAAAAAAAAAAAA')
						if not genemiss == 'MITOCONDRIO':
							h = HGNC()
							search = h.search(genemiss)
							print (search,genemiss,'bbbbbbbbbbbb')
							b = search['response']['docs'][0]['symbol']
							out_file.write(str(genemiss)+'\t')
							out_file.write(str(b)+'\n')
							#out_file.write(''+'\n')
							out_file_download.write(str(genemiss)+'\t')
							out_file_download.write(str(b)+'\n')
							ufficial_gene.append(str(b))
			else: ufficial_gene=[]
			#except: ufficial_gene=[]

			try:
				a2 = mg.querymany(list(ufficial_gene), scopes="symbol",fields=["symbol","exons"], species="human", as_dataframe=True,returnall=True)
				gene_DF2 = a2['out']
				gene_DF2['symbol'].fillna('ciao',inplace=True)
				gene_DF2 = gene_DF2[gene_DF2['symbol']!='ciao']
				gene_DF2['info'] =  'NotOfficialName'
			except: gene_DF2 = pd.DataFrame()

			try:
				snp = mv.querymany(list(snplist), scopes='_id,dbsnp.rsid', fields='dbsnp,dbnsfp,gnomad_exome,gnomad_genome',species="human",as_dataframe=True,returnall=True)
				snp_DF3a =  snp['out']
				snp_DF3 = snp_DF3a[['dbsnp.rsid','dbsnp.gene.symbol','dbsnp.chrom',
										'dbsnp.hg19.start','dbsnp.hg19.end',
										'dbsnp.vartype','dbsnp.ref','dbsnp.alt']]
				snp_DF3['dbsnp.class'] = ''
				snp_DF3['dbsnp.gene.strand']=-1
				snp_DF3.reset_index(inplace=True)
			except:
				snp_DF3 = pd.DataFrame()

			for index,snp in snp_DF3.iterrows():
				#print (snp,'aaaaaaaaaaaaaaaaaaa')
				_liftCHROM_ = snp['dbsnp.chrom']
				liftCHROM = 'chr'+str(_liftCHROM_)
				liftPOS = snp['dbsnp.hg19.start']
				try: POS38 = liftover_from37_to38(liftCHROM,liftPOS)
				except: POS38 = -999
				snp_DF3.loc[index,'#CHROM'] = liftCHROM
				snp_DF3.loc[index,'START'] = int(int(POS38)-1)
				snp_DF3.loc[index,'END'] = int(int(POS38)+1)
				snp_DF3.loc[index,'GENE'] = snp['dbsnp.gene.symbol']
				snp_DF3.loc[index,'length'] = 3
				if snp['dbsnp.gene.strand'] == '+': _strand_=1
				elif snp['dbsnp.gene.strand'] == '-': _strand_=-1
				else: _strand_=1
				snp_DF3.loc[index,'strand'] = _strand_
				snp_DF3.loc[index,'exone'] = 0
				snp_DF3.loc[index,'refseq'] = 'unknown'	#snp['dbsnp.gene.rnas']
				snp_DF3.loc[index,'hgmd'] = snp['dbsnp.rsid'] #'SNP'
				snp_DF3.loc[index,'REF']=snp['dbsnp.ref'] #'SNP'
				snp_DF3.loc[index,'ALT']=snp['dbsnp.alt'] #'SNP'
				snp_DF3.loc[index,'START']=snp_DF3.loc[index,'START'].astype(int)
				snp_DF3.loc[index,'END']=snp_DF3.loc[index,'END'].astype(int)
				_CHROM_=snp_DF3.loc[index,'#CHROM']
				_END_=snp_DF3.loc[index,'END'].astype(int)
				_REF_=snp_DF3.loc[index,'REF']
				_ALT_=snp_DF3.loc[index,'ALT']

				try:
					_G_ = mg.query(snp['dbsnp.gene.symbol'], scopes="symbol",fields=["symbol","refseq"], species="human", as_dataframe=True)
					#print (snp,_G_,'bbbbbbbbbbbbbbbbb')
					_ref_ =(_G_['refseq.rna'][0][0].split('.')[0])
					if _ref_=='N': _ref_ =(_G_['refseq.rna'][0].split('.')[0])
				except: _ref_ = 'unknown'

				try:
					PREsnp_DF3[['GENE','exone','length','strand','refseq','hgmd']]=np.nan
					#print (PREsnp_DF3)
					snp_DF3.loc[index,'refseq'] = (_ref_)
					TABLE_SNP = pd.concat([PREsnp_DF3,snp_DF3[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']]]).sort_values(by=['#CHROM','START','hgmd'])
					TABLE_SNP = TABLE_SNP[TABLE_SNP['#CHROM']!='chrnan']
				except:
					snp_DF3.loc[index,'refseq'] = snp['dbsnp.gene.symbol']
					TABLE_SNP = snp_DF3[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['#CHROM','START','hgmd'])
					TABLE_SNP = TABLE_SNP[TABLE_SNP['#CHROM']!='chrnan']

			if len(TABLE_SNP)>0:
				TABLE_SNP['strand'].fillna(-1,inplace=True)
				TABLE_SNP['exone'].fillna(0,inplace=True)
				TABLE_SNP['hgmd'].fillna('notfound',inplace=True)
				TABLE_SNP['length'].fillna((TABLE_SNP['END'].astype(int)-TABLE_SNP['START'].astype(int)+1),inplace=True)

				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>94001180)&(TABLE_SNP['START']<94115089)),'ABCA4',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>94001180)&(TABLE_SNP['START']<94115089)),'NM_000350',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>94001180)&(TABLE_SNP['START']<94115089)),-1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>215794439)&(TABLE_SNP['START']<216074139)),'USH2A',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>215794439)&(TABLE_SNP['START']<216074139)),'NM_206933',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>215794439)&(TABLE_SNP['START']<216074139)),-1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>88101180)&(TABLE_SNP['START']<88101185)),'CEP290',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>88101180)&(TABLE_SNP['START']<88101185)),'NM_025114',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>88101180)&(TABLE_SNP['START']<88101185)),-1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>15988232)&(TABLE_SNP['START']<15988239)),'PROM1',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>15988232)&(TABLE_SNP['START']<15988239)),'NM_006017',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>15988232)&(TABLE_SNP['START']<15988239)),-1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>63488520)&(TABLE_SNP['START']<63488539)),'ACE',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>63488520)&(TABLE_SNP['START']<63488539)),'NM_000789',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>63488520)&(TABLE_SNP['START']<63488539)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>233760223)&(TABLE_SNP['START']<233760243)),'UGT1A1',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>233760223)&(TABLE_SNP['START']<233760243)),'NM_000463',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>233760223)&(TABLE_SNP['START']<233760243)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>101126415)&(TABLE_SNP['START']<101126435)),'SERPINE1',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>101126415)&(TABLE_SNP['START']<101126435)),'NM_000602',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>101126415)&(TABLE_SNP['START']<101126435)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>7870849)&(TABLE_SNP['START']<7870871)),'MTRR',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>7870849)&(TABLE_SNP['START']<7870871)),'NM_002454',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>7870849)&(TABLE_SNP['START']<7870871)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>8745633)&(TABLE_SNP['START']<8745655)),'CAV3',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>8745633)&(TABLE_SNP['START']<8745655)),'NM_033337',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>8745633)&(TABLE_SNP['START']<8745655)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>34370492)&(TABLE_SNP['START']<34370514)),'KCNE2',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>34370492)&(TABLE_SNP['START']<34370514)),'NM_172201',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>34370492)&(TABLE_SNP['START']<34370514)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>97305353)&(TABLE_SNP['START']<97305375)),'DPYD',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>97305353)&(TABLE_SNP['START']<97305375)),'NM_000110',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>97305353)&(TABLE_SNP['START']<97305375)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>99672905)&(TABLE_SNP['START']<99672927)),'CYP3A5',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>99672905)&(TABLE_SNP['START']<99672927)),'NM_000777',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>99672905)&(TABLE_SNP['START']<99672927)),1,TABLE_SNP['strand'])
				TABLE_SNP['GENE'] = np.where(((TABLE_SNP['START']>169549800)&(TABLE_SNP['START']<169549822)),'F5',TABLE_SNP['GENE'])
				TABLE_SNP['refseq'] = np.where(((TABLE_SNP['START']>169549800)&(TABLE_SNP['START']<169549822)),'NM_000130',TABLE_SNP['refseq'])
				TABLE_SNP['strand'] = np.where(((TABLE_SNP['START']>169549800)&(TABLE_SNP['START']<169549822)),1,TABLE_SNP['strand'])


				TABLE_SNP.sort_values(by=['#CHROM','START','hgmd'],ascending=False, inplace=False)
				TABLE_SNP.reset_index(inplace=True)
				TABLE_SNP.drop(['index'],axis=1,inplace=True)
				try:
					TABLE_SNP.loc[:,'START'] = TABLE_SNP['START'].astype(int)
					TABLE_SNP.loc[:,'END'] = TABLE_SNP['END'].astype(int)
					TABLE_SNP.loc[:,'strand'] = TABLE_SNP['strand'].astype(int)
					TABLE_SNP.loc[:,'exone'] = TABLE_SNP['exone'].astype(int)
					TABLE_SNP.loc[:,'length'] = TABLE_SNP['length'].astype(int)
				except:
					TABLE_SNP = pd.DataFrame()

			out_file.close()
			out_file_download.close()

			if (len(gene_DF1)==0 & len(gene_DF2)):
				if len(TABLE_SNP)>0:
					print ('-----------')
				else: print ('NESSUN DATO TROVATO!!!')
			else:
				gene_DF = pd.concat([gene_DF1,gene_DF2])
				try:
					gene_DF.dropna(subset=['exons'],axis=0,how='all',inplace=True)
					gene_DF.reset_index(inplace=True)
				except: gene_DF=pd.DataFrame()

			if (len(gene_DF1)==0 & len(gene_DF2)):
				if len(TABLE_SNP)>0:
					print ('-----------')
				else: print ('NESSUN DATO TROVATO!!!')
			else:
				gene_DF = pd.concat([gene_DF1,gene_DF2])
				gene_DF.dropna(subset=['exons'],axis=0,how='all',inplace=True)
				gene_DF.reset_index(inplace=True)

			if len(gene_DF)>=1:
					i=0
					z = dict(gene_DF['exons'])
					for DICT in z.keys():
						refseq_rna = z[DICT]
						for RNA in refseq_rna:
							try:
								chromosome = RNA['chr']
								if len(chromosome.split('_')) > 1:
									chromosome = RNA[1]['chr']
									strand = RNA[1]['strand']
									cdsstart = RNA[1]['cdsstart']
									cdsend = RNA[1]['cdsend']
									exons_number = len(RNA[1]['position'])
									if int(strand) == -1:
										_exon_ = max(range(exons_number))+1
									elif int(strand) == 1:
										_exon_ = min(range(exons_number))+1
									for exons in RNA[1]['position']:
										_CHROM_ = 'chr'+str(chromosome)
										start = exons[0]
										end = exons[1]
										CHROM.append('chr'+str(chromosome))
										START.append(str(start))
										END.append(str(end))
										CDSSTART.append(str(cdsstart))
										CDSEND.append(str(cdsend))
										refseq.append(str(RNA['transcript']))
										geno.append('ND')
										exone.append(str(_exon_))
										tipo.append('FRAGMENT')
										gene.append(str(gene_DF['symbol'][i]))
										INFO.append(str(gene_DF['info'][i]))
										STRAND.append(str(strand))
										variazione.append(str(gene_DF['symbol'][i])+':'+str(RNA['transcript'])+':ex'+str(_exon_))
										annotazione.append('')
										if int(strand) == -1:
											_exon_ = _exon_-1
										elif int(strand) == 1:
											_exon_ = _exon_+1
								else:
									chromosome = RNA[0]['chr']
									strand = RNA[0]['strand']
									cdsstart = RNA[0]['cdsstart']
									cdsend = RNA[0]['cdsend']
									exons_number = len(RNA[0]['position'])

									if int(strand) == -1:
										_exon_ = max(range(exons_number))+1
									elif int(strand) == 1:
										_exon_ = min(range(exons_number))+1
									for exons in RNA[0]['position']:
										_CHROM_ = 'chr'+str(chromosome)
										start = exons[0]
										end = exons[1]
										CHROM.append('chr'+str(chromosome))
										START.append(str(start))
										END.append(str(end))
										CDSSTART.append(str(cdsstart))
										CDSEND.append(str(cdsend))
										refseq.append(str(RNA['transcript']))
										geno.append('ND')
										exone.append(str(_exon_))
										tipo.append('FRAGMENT')
										gene.append(str(gene_DF['symbol'][i]))
										INFO.append(str(gene_DF['info'][i]))
										STRAND.append(str(strand))
										variazione.append(str(gene_DF['symbol'][i])+':'+str(RNA['transcript'])+':ex'+str(_exon_))
										annotazione.append('')
										if int(strand) == -1:
											_exon_ = _exon_-1
										elif int(strand) == 1:
											_exon_ = _exon_+1
							except KeyError:
								chromosome = RNA['chr']
								strand = RNA['strand']
								cdsstart = RNA['cdsstart']
								cdsend = RNA['cdsend']
								exons_number = len(RNA['position'])
								if int(strand) == -1:
									_exon_ = max(range(exons_number))+1
								elif int(strand) == 1:
									_exon_ = min(range(exons_number))+1
									
								for exons in RNA['position']:
									_CHROM_ = 'chr'+str(chromosome)
									start = exons[0]
									end = exons[1]
									CHROM.append('chr'+str(chromosome))
									START.append(str(start))
									END.append(str(end))
									CDSSTART.append(str(cdsstart))
									CDSEND.append(str(cdsend))
									refseq.append(str(RNA['transcript']))
									geno.append('ND')
									exone.append(str(_exon_))
									tipo.append('FRAGMENT')
									gene.append(str(gene_DF['symbol'][i]))
									INFO.append(str(gene_DF['info'][i]))
									STRAND.append(str(strand))
									variazione.append(str(gene_DF['symbol'][i])+':'+str(RNA['transcript'])+':ex'+str(_exon_))
									annotazione.append('')
									if int(strand) == -1:
										_exon_ = _exon_-1
									elif int(strand) == 1:
										_exon_ = _exon_+1
							else:
								print ('GENE NON TROVATO!!!',str(gene_DF['symbol'][i]))
								pass
						i+=1
			else:
						CHROM.append('chr'+str(0))
						START.append(str(0))
						END.append(str(0))
						CDSSTART.append(str(0))
						CDSEND.append(str(0))
						refseq.append(str('Unknwon'))
						geno.append('ND')
						exone.append(str(0))
						tipo.append('FRAGMENT')
						gene.append(str('Unknown'))
						INFO.append(str('Unknown'))
						STRAND.append(str(0))
						variazione.append(str('Unknown'))
						annotazione.append('')

			TABLE = pd.DataFrame({'#CHROM':pd.Series(CHROM),'cdsstart':pd.Series(CDSSTART),
					'cdsend':pd.Series(CDSEND),'GENE':pd.Series(gene),'refseq':pd.Series(refseq),
					'START':pd.Series(START),'END':pd.Series(END),'INFO':pd.Series(INFO),
					'exone':pd.Series(exone),'geno':pd.Series(geno),'tipo':pd.Series(tipo),
					'strand':pd.Series(STRAND),'variazione':pd.Series(variazione),'annotazione':pd.Series(annotazione)})
			TABLE['NUMCHROM'] = TABLE['#CHROM'].apply(lambda x: len(x.split('_')))
			TABLE = TABLE[TABLE['NUMCHROM']==1]
			TABLE.drop(['NUMCHROM'], axis=1,inplace=True)
			TABLE2_A = pd.merge(TABLE,hgmd_df,on=['GENE','refseq'],how='left')
			TABLE2_A.dropna(subset=['hgmd'],inplace=True)
			TABLE2_B = pd.merge(TABLE,eccezioni_df,on=['GENE','refseq'],how='left')
			TABLE2_B.dropna(subset=['hgmd'],inplace=True)
			TABLE2 = pd.concat([TABLE2_A,TABLE2_B]).sort_values(by=['hgmd','GENE','exone'])
			if len(TABLE2) == 0:
				TABLE2 = TABLE
				TABLE2['hgmd'] = 'hgmd'
			TABLE2.reset_index(inplace=True)
			TABLE2.drop(['index'],axis=1,inplace=True)
			TABLE2['START'] = TABLE2['START'].astype(int)
			TABLE2['END'] = TABLE2['END'].astype(int)
			TABLE2['cdsstart'] = TABLE2['cdsstart'].astype(int)
			TABLE2['cdsend'] = TABLE2['cdsend'].astype(int)
			mask1 = TABLE2['START'] < TABLE2['cdsstart']
			mask2 = TABLE2['END'] < TABLE2['cdsstart']
			mask3 = TABLE2['START'] > TABLE2['cdsend']
			mask4 = TABLE2['END'] > TABLE2['cdsend']
			print (TABLE2)
###############################################################################
			TABLE2.loc[mask1,'CODING_START'] = 'NoNcoding'
			TABLE2.loc[mask2,'CODING_END'] = 'NoNcoding'
			TABLE2.loc[mask3,'CODING_START'] = 'NoNcoding'
			TABLE2.loc[mask4,'CODING_END'] = 'NoNcoding'
			TABLE2['CODING_START'].fillna('coding',inplace=True)
			TABLE2['CODING_END'].fillna('coding',inplace=True)

			pattern = 'NR_'
			maskNR = TABLE2['refseq'].str.contains(pattern)
			TABLE2.loc[maskNR,'CODING_START'] = 'coding'
			TABLE2.loc[maskNR,'CODING_END'] = 'coding'

			mask1 = TABLE2['CODING_START']=='NoNcoding'
			mask2 = TABLE2['CODING_END']=='NoNcoding'
			TABLE2.loc[mask1,'START'] = TABLE2['cdsstart']
			TABLE2.loc[mask2,'END'] = TABLE2['cdsend']
###############################################################################
			TABLE2.sort_values(by=['hgmd','GENE','exone'],inplace=True)
			TABLE_GROUP = TABLE2.groupby(['GENE'])

			df = pd.DataFrame()
			dfx = pd.DataFrame()
			for index, group in TABLE_GROUP:
				if len(group) <= 2:
					group.drop_duplicates(subset=['#CHROM','START','END'],keep='last',inplace=True)
				elif len(group) > 2:
					if (((group['CODING_START']=='NoNcoding').all())&((group['CODING_END']=='NoNcoding').all())):
						if group['hgmd'].isin(['hgmd']).any():
							group = group[group['hgmd'].isin(['hgmd'])]
							group.drop_duplicates(subset=['#CHROM','START','END'],keep='last',inplace=True)
						else:
							group.drop_duplicates(subset=['#CHROM','START','END'],keep='first',inplace=True)
					else:
						group.drop_duplicates(subset=['#CHROM','START','END'],keep='first',inplace=True) # keeps exons from the hgmd transcript only
				else: pass
				df = pd.concat([df,group])
				dfx = pd.concat([df,group])
			#print df[df['GENE']=='GRIN2A']
###############################################################################
			TABLE2 = df
			TABLE2['START'] = TABLE2['START'].astype(int)
			TABLE2['END'] = TABLE2['END'].astype(int)
			TABLE2['strand'] = TABLE2['strand'].astype(int)
			TABLE2['exone'] = TABLE2['exone'].astype(int)
			TABLE2['length'] = (TABLE2['END'])-(TABLE2['START'])
			TABLE2['length'] = TABLE2['length'].astype(int)

			#TABLE2['START'] = TABLE2['START']-15
			#TABLE2['END'] = TABLE2['END']+15
			TABLE2['START'] = TABLE2['START']-5
			TABLE2['END'] = TABLE2['END']+5
###############################################################################################################
			group = TABLE2.groupby(['GENE'])
			size = pd.DataFrame(group.size()).reset_index()
			size['size'] = size[0]
			gene_unique = size[size['size']==1]
			for _gene_ in list(gene_unique['GENE']): # the exons that are the only exons for a gene will be marked as coding in both start and end
				mask = TABLE2['GENE'] == _gene_
				TABLE2.loc[mask,'CODING_START'] = 'coding'
				TABLE2.loc[mask,'CODING_END'] = 'coding'

			#pattern = 'NM_'
			#TABLE2 = TABLE2[TABLE2['refseq'].str.contains(pattern)]
			TABLE2.drop_duplicates(subset=['#CHROM','START','END'],inplace=True)
			TABLE2B = TABLE2[~((TABLE2['CODING_START']=='NoNcoding')&(TABLE2['CODING_END']=='NoNcoding'))]
			TABLE3B = TABLE2B #[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']]

			_col_ = ['#CHROM','START','END','GENE','exone','refseq',
				'hgmd','strand','variazione','length',
				'INFO','annotazione','cdsend','cdsstart',
				'geno','tipo','CODING_START','CODING_END']

			if TABLE3B['GENE'].isin(['CBS']).any():
				print ("TROVATO!!!"),
				MORETABLE = TABLE3B[TABLE3B['refseq']=='NM_001178008'].sort_values(by=['GENE','START'],ascending=[True,True])
				MORETABLE.drop_duplicates(subset=['GENE','exone'],keep='last',inplace=True)
				_TABLE_ = TABLE3B[TABLE3B['GENE']!='CBS']
				TABLE3 = pd.concat([_TABLE_,MORETABLE])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
				TABLE3B = pd.concat([_TABLE_,MORETABLE])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
			else:
				TABLE3 = TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

			if TABLE3B['GENE'].isin(['CFTR']).any():
				mask1 = ((TABLE3B['GENE']=='CFTR') & (TABLE3B['exone']==10))
				TABLE3B.loc[mask1,'START'] = 117548600
				TABLE3 =  TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
				TABLE3B = TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
			else:
				TABLE3 = TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

			TABLE4_38MINUS = TABLE3[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['#CHROM','START','END'])
			if 'HRURF' in list(GENE): # add this region (exon) for HRURF if it is present in the original gene list
				TABLE_PLUS = pd.DataFrame([['chr8','22130597','22130712','HRURF',1,106,-1,'NM_001394132','hgmd']],
					   columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])
				TABLE4_38A = pd.concat([TABLE4_38MINUS,TABLE_PLUS])
			else:TABLE4_38A=TABLE4_38MINUS
			#print (list(GENE))
			if 'MITOCONDRIO' in list(GENE): # add these regions if MITOCONDRION is present in the original gene list
				#print ('WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW')
				TABLE_MITOCONDRIO = pd.DataFrame([['chrMT','1','16569','MITOCONDRIO',0,16569,1,'NC_012920','hgmd']],
				   	columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])
				TABLE4_38A = pd.concat([TABLE4_38MINUS,TABLE_MITOCONDRIO]).sort_values(by=['#CHROM','START','END'])
			else:TABLE4_38A=TABLE4_38MINUS

			#TABLE4_38A = TABLE3.sort_values(by=['#CHROM','START','END'])
			TABLE4_38A.drop_duplicates(subset=['#CHROM','START'],keep='last',inplace=True)
			TABLE4_38A.drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)

			if len(TABLE_SNP) > 0:
				TABLE_SNP.drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)
				TABLE4_38B = pd.concat([TABLE4_38A,TABLE_SNP]).reset_index()
			else:
				TABLE4_38B = TABLE4_38A

			TABLE3 = TABLE4_38B.sort_values(by=['GENE','START'])
			TABLE4 = TABLE4_38B[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['GENE','START'])
			file_ = (join(FOLDER,'TOTAL_BED_'+str(SAMPLE)))
			TABLE3.to_csv(file_,sep='\t',index=False)
			file_ = (join(FOLDER,'bed_'+str(SAMPLE)))
			sample_obj.BED = file_ # TODO: fix this
			file_download = (join(config.DOWNLOAD_FOLDER,'bed_'+str(SAMPLE)))
			TABLE4.to_csv(file_,sep='\t',index=False)
			TABLE4.to_csv(file_download,sep='\t',index=False)
#######################################################################################################################################################
			TABLE2X = dfx
			TABLE2X['START'] = TABLE2X['START'].astype(int)
			TABLE2X['END'] = TABLE2X['END'].astype(int)
			TABLE2X['strand'] = TABLE2X['strand'].astype(int)
			TABLE2X['exone'] = TABLE2X['exone'].astype(int)
			TABLE2X['length'] = (TABLE2X['END'])-(TABLE2X['START'])
			TABLE2X['length'] = TABLE2X['length'].astype(int)
			TABLE2X['START'] = TABLE2X['START']-15
			TABLE2X['END'] = TABLE2X['END']+15
#######################################################################################################################################################
			group = TABLE2X.groupby(['GENE'])
			size = pd.DataFrame(group.size()).reset_index()
			size['size'] = size[0]
			gene_unique = size[size['size']==1]
			for _gene_ in list(gene_unique['GENE']):
				mask = TABLE2X['GENE'] == _gene_
				TABLE2X.loc[mask,'CODING_START'] = 'coding'
				TABLE2X.loc[mask,'CODING_END'] = 'coding'

			TABLE2X.drop_duplicates(subset=['#CHROM','START','END'],inplace=True)
			TABLE2BX = TABLE2X[~((TABLE2X['CODING_START']=='NoNcoding')&(TABLE2X['CODING_END']=='NoNcoding'))]
			TABLE3BX = TABLE2BX #[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']]

			_col_ = ['#CHROM','START','END','GENE','exone','refseq',
				'hgmd','strand','variazione','length',
				'INFO','annotazione','cdsend','cdsstart',
				'geno','tipo','CODING_START','CODING_END']

			if TABLE3BX['GENE'].isin(['CBS']).any():
				print ("TROVATO!!!"),
				MORETABLEX = TABLE3BX[TABLE3BX['refseq']=='NM_001178008'].sort_values(by=['GENE','START'],ascending=[True,True])
				MORETABLEX.drop_duplicates(subset=['GENE','exone'],keep='last',inplace=True)
				_TABLE_X = TABLE3BX[TABLE3BX['GENE']!='CBS']
				TABLE3X = pd.concat([_TABLE_X,MORETABLEX])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
				TABLE3BX = pd.concat([_TABLE_X,MORETABLEX])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
			else:
				TABLE3X = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

			if TABLE3BX['GENE'].isin(['CFTR']).any():
				mask1 = ((TABLE3BX['GENE']=='CFTR') & (TABLE3BX['exone']==10))
				TABLE3BX.loc[mask1,'START'] = 117548600
				TABLE3X =  TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
				TABLE3BX = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
			else:
				TABLE3X = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

			if 'HRURF' in list(GENE):
				TABLE_PLUS = pd.DataFrame([['chr8','22130597','22130712','HRURF',1,106,-1,'NM_001394132','hgmd']],
					   columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])
				TABLE3X = pd.concat([TABLE3BX[TABLE3BX['GENE']!='HRURF'],TABLE_PLUS])
				TABLE3BX = pd.concat([TABLE3BX[TABLE3BX['GENE']!='HRURF'],TABLE_PLUS])
			else:
				TABLE3X = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

			TABLE4_38AX = TABLE3X[_col_].sort_values(by=['#CHROM','START','END'])
			TABLE4_38AX[_col_].drop_duplicates(subset=['#CHROM','START'],keep='last',inplace=True)
			TABLE4_38AX[_col_].drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)

			if len(TABLE_SNP) > 0:
				TABLE_SNP.drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)
				TABLE4_38BX = pd.concat([TABLE4_38AX,TABLE_SNP]).reset_index()
			else:
				TABLE4_38BX = TABLE4_38AX

			TABLE3X = TABLE4_38BX.sort_values(by=['GENE','START'])
			TABLE4X = TABLE4_38BX[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['GENE','START'])
			TABLE4X_XY = pd.DataFrame([['chrX','11296910','11296920','AMELX',5,3,1,'NM_001142','hgmd'],
					   ['chrY','6869900','6869920','AMELY',5,9,-1,'NM_001143','hgmd'],
					   ['chrY','2787170','2787270','SRY',1,94,-1,'NM_003140','hgmd']],
					   columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])

			TABLE4X = pd.concat([TABLE4X,TABLE4X_XY])
			file_ = (join(FOLDER,'TOTAL_BEDX_'+str(SAMPLE)))
			TABLE3X.to_csv(file_,sep='\t',index=False)
			file_ = (join(FOLDER,'bedX_'+str(SAMPLE)))
			file_download = (join(config.DOWNLOAD_FOLDER,'bed_'+str(SAMPLE)))

			TABLE4X.to_csv(file_,sep='\t',index=False)
			TABLE4X.to_csv(file_download,sep='\t',index=False)
############################################################################################################
			# TABLE4 is the BED file
			vertical,verticalX = create_vertical(TABLE4,TABLE4X,SAMPLE,FOLDER)
			file_2 = (join(FOLDER,'vertical_'+str(SAMPLE)))
			vertical.to_csv(file_2,sep='\t',index=False)
			file_2X = (join(FOLDER,'verticalX_'+str(SAMPLE)))
			verticalX.to_csv(file_2X,sep='\t',index=False)
			print ('FINITO!!! TUTTO OK!!!')

			sample_obj.vertical = file_2
			sample_obj.verticalX = file_2X
			sample_obj.bedX = file_
			sample_obj.bed = (join(FOLDER,'bed_'+str(SAMPLE)))
			sample_obj.saveJSON()

			return vertical,verticalX,TABLE4
		#except: return None
############################################################################################################
############################################################################################################

def cutCDS_mane(gene_list, sample_obj, folder, dest):
    """ Function that builds the bed file from the MANE gtf file. The bedX, which is identical to bed except the three AMELX, AMELY, and SRY genes, is also built. Both beds are saved as text files.
    It additionally builds the respective vertical files for each bed. Verticals are also saved as text files.

    Args:
        gene_list (list): a gene list that the bed file will contain (redundant, since the genes  can be extracted from the phenotype file by having the sample_id).
        sample_id (str): the sample id which can be used to extract the genelist.
        folder (str): the path of the pheno folder, where the bed files will be placed.
        dest (str): r, d, or z, specifying where the server where this sample is stored (bolzano or rovereto) - not used.

    Returns:
        pandas.DataFrame: The bed is returned as a dataframe
        pandas.DataFrame: vertical of the bed is returned as a dataframe
        pandas.DataFrame: verticalX of the bedX is returned as a dataframe


    Todo:
        * It is important to split this function into smaller pieces. The creation of verticals and the file saving should be done by seperate functions.
    """
    
    sample_id = sample_obj.name

    # Define output paths
    bed_out_path = os.path.join(folder, "bed_{}".format(sample_id))
    bedX_out_path = os.path.join(folder, "bedX_{}".format(sample_id))
    vertical_out_path = os.path.join(folder, "vertical_{}".format(sample_id))
    verticalX_out_path = os.path.join(folder, "verticalX_{}".format(sample_id))

    # gene_list = pd.read_csv(gene_list_path, header=None, names=['GENE'])

    #sample_id = "VerRGC20509.2024"

    # # Read the gene list
    # phenotype = pd.read_csv("phenotype", sep="\t")
    # gene_list = list(phenotype[phenotype["sample"] == sample_id]["gene"])
    # gene_list.append("CBS")

    #gene_list = list(["ORAI1", "SHANK3"])

    colnames_bed_mane = ['CHROM', 'SOURCE', 'TYPE', 'START', 'END', 'SCORE', 'STRAND', 'FRAME', 'INFO']
    bed_mane = pd.read_csv(gtf_path, sep='\t', comment='#', names=colnames_bed_mane)

    bed_mane.drop(columns=['SOURCE', 'SCORE', 'FRAME'], inplace=True)
    bed_mane = bed_mane[bed_mane['TYPE'].str.contains('CDS')]

    bed_mane['GENE'] = bed_mane['INFO'].str.extract('gene "(.*?)"')
    bed_mane['TRANSCRIPT_ID'] = bed_mane['INFO'].str.extract('transcript_id "(.*?)"')
    bed_mane['EXON_NUMBER'] = bed_mane['INFO'].str.extract('exon_number "(.*?)"')
    bed_mane['TAG'] = bed_mane['INFO'].str.extract('tag "(.*?)"')

    new_bed = bed_mane[bed_mane['GENE'].isin(gene_list)] # TODO: check if all genes were found in the gtf

    new_bed['STRAND'] = new_bed['STRAND'].replace({'+': 1, '-': -1})

    new_bed['LENGTH'] = new_bed.apply(lambda row: row['END'] - row['START'], axis=1)
    new_bed['LENGTH'] = new_bed['LENGTH'] + 1                                                 #+1 per aggiustare lunghezza esone...
    new_bed.drop(columns=['INFO', 'TYPE'], inplace=True)


    # Filter rows with TAG "MANE Select"
    mane_select = new_bed[new_bed['TAG'] == 'MANE Select']

    # Merge and filter to find unique rows
    merged = new_bed.merge(mane_select[['CHROM', 'START', 'END', 'GENE']],
                        on=['CHROM', 'START', 'END', 'GENE'],
                        how='left', indicator=True)

    filtered_df = merged[(merged['_merge'] == 'left_only') | (merged['TAG'] == 'MANE Select')]
    filtered_df.drop(columns=['_merge'], inplace=True)

    filtered_df['EXON_NUMBER'] = filtered_df['EXON_NUMBER'].astype(int)
    filtered_df.sort_values(by=['GENE', 'EXON_NUMBER'], inplace=True)


    ###

    filtered_df["TAG"] = filtered_df["TAG"].str.replace("\s", "_", regex=True)
    filtered_df["CHROM"] = filtered_df["CHROM"].str.replace("_.*", "", regex=True)

    filtered_df = filtered_df.rename(columns={"CHROM": "#CHROM", "TRANSCRIPT_ID": "refseq", "LENGTH": "length", "STRAND": "strand", "EXON_NUMBER": "exone", "TAG":"hgmd"})
    bed = filtered_df.loc[:,['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']]
    bed = bed[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['GENE','START'])

    

    # Build bedx
    bedX = bed[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']]
    sex_exons = pd.DataFrame([
            ['chrX','11298548','11298973','AMELX',5, 425, 1,'NM_001142.2','MANE_Select'],
            ['chrY','6868037','6868462','AMELY', 6, 425, -1,'NM_001143.2','MANE_Select'],
            ['chrY','2786992','2787603','SRY', 1, 611, -1,'NM_003140.3','MANE_Select']],
            columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])
    bedX = pd.concat([bedX, sex_exons])
	
	# Convert to numeric
    bed[['START', 'END', 'strand', 'exone', 'length']] = bed[['START', 'END', 'strand', 'exone', 'length']].apply(pd.to_numeric)
    bedX[['START', 'END', 'strand', 'exone', 'length']] = bed[['START', 'END', 'strand', 'exone', 'length']].apply(pd.to_numeric)  

	# Include flanking regions
    bed['START'] = bed['START'] - 5
    bed['END'] = bed['END'] + 5
    bedX['START'] = bedX['START'] - 15
    bedX['END'] = bedX['END'] + 15
	
    vertical, verticalX = create_vertical(bed, bedX, sample_id, folder)


    # Save each file
    bed.to_csv(bed_out_path, sep='\t', index=False)
    bedX.to_csv(bedX_out_path, sep='\t', index=False)
    vertical.to_csv(vertical_out_path, sep='\t', index=False)
    verticalX.to_csv(verticalX_out_path, sep='\t', index=False)
	
    sample_obj.bed = bed_out_path
    sample_obj.bedX = bedX_out_path
    sample_obj.vertical = vertical_out_path
    sample_obj.verticalX = verticalX_out_path
    sample_obj.saveJSON()
	
    print ('FINITO!!! TUTTO OK!!!')
    
    return vertical, verticalX, bed



