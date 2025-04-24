import config
import pandas as pd
import numpy as np
import os
import mygene
#import myvariant
import re
from pyliftover import LiftOver
from bioservices import HGNC

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

hgmd = os.path.join(config.HOME_DIR, 'PROJECT/diagnosys/bin/APPRIS_PRINCIPAL')
hgmd_df = pd.read_csv(hgmd,sep='\t',header=0)
df_orig = pd.read_csv(hgmd,sep='\t',header=0)

mg = mygene.MyGeneInfo()
chainfrom37to38 = os.path.join('/home/magi/dataset/CHAIN_LIFTOVER','hg19ToHg38.over.chain.gz')
chainfrom38to37 = os.path.join('/home/magi/dataset/CHAIN_LIFTOVER','hg38ToHg19.over.chain.gz')
lo = LiftOver(chainfrom38to37)
lo37 = LiftOver(chainfrom37to38)

hgmd_df['refseq2'] = hgmd_df['refseq'].str.split('.').str.get(0)
hgmd_df['refseq'] = hgmd_df['refseq2']
hgmd_df['APPRIS2'] = hgmd_df['APPRIS'].str.split('.').str.get(0)
hgmd_df['APPRIS'] = hgmd_df['APPRIS2']
hgmd_df['hgmd'] = 'hgmd'
hgmd_df = hgmd_df[['GENE','refseq','hgmd']]

eccezioni = os.path.join(config.HOME_DIR, 'PROJECT/diagnosys/bin/APPRIS_ECCEZIONI')
eccezioni_df = pd.read_csv(eccezioni,sep='\t',header=0)
eccezioni_df['refseq2'] = eccezioni_df['refseq'].str.split('.').str.get(0)
eccezioni_df['refseq'] = eccezioni_df['refseq2']
#eccezioni_df['APPRIS2'] = eccezioni_df['APPRIS'].str.split('.').str.get(0)
#eccezioni_df['APPRIS'] = eccezioni_df['APPRIS2']
eccezioni_df['hgmd'] = 'other'
eccezioni_df = eccezioni_df[['GENE','refseq','hgmd']]


def create_null_coverage(start, end, step=1):
	return pd.DataFrame({'POS': range(start, end, step)})

""" 
Split each chromosome position into unique rows.
Input
    DATA: ,
	DATAX: ,
	FILENAME: ,
	FOLDER: 
"""
def create_vertical(DATA, DATAX, FILENAME, FOLDER):
	df1 = pd.DataFrame()
	print ('START VERTICAL!!!' + str(FILENAME))
	for row_index, row in DATA.iterrows():
		START = int(row['START'])
		END = int(row['END'])
		gene = row['GENE']
		chrom = row['#CHROM']
		exone = row['exone']
		length = row['length']
		strand = row['strand']
		refseq = row['refseq']
		hgmd = 	row['hgmd']

		result = create_null_coverage(START, END)
		result['GENE'] = gene
		result['#CHROM'] = chrom
		result['exone'] = exone
		result['length'] = length
		result['strand'] = strand
		result['refseq'] = refseq
		result['hgmd'] = hgmd

		df1 = pd.concat([df1, result])
		
	df1.POS = df1['POS'].astype(int)
	df1.length = df1['length'].astype(int)
	df1.strand = df1['strand'].astype(int)
	df1.exone = df1['exone'].astype(int)
	vertical = df1.sort_values(by=['#CHROM','POS','hgmd'],ascending=[True,True,True])
	vertical.drop_duplicates(subset=['#CHROM','POS','GENE'],inplace=True)
	vertical = vertical[['#CHROM','POS','GENE','exone','length','strand','refseq','hgmd']]

	df2 = pd.DataFrame()
	print ('START VERTICALX!!!' + str(FILENAME))
	for row_index, rowx in DATAX.iterrows():
		START = int(rowx['START'])
		END = int(rowx['END'])
		gene = rowx['GENE']
		chrom = rowx['#CHROM']
		exone = rowx['exone']
		length = rowx['length']
		strand = rowx['strand']
		refseq = rowx['refseq']
		hgmd = 	rowx['hgmd']

		result = create_null_coverage(START, END)
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
	
	return vertical, verticalX

def liftover_from37_to38(CHROM,POS):
	lift = lo37.convert_coordinate(CHROM,POS)
	return lift[0][1]

def liftover_from38_to37(CHROM,POS):
	lift = lo.convert_coordinate(CHROM,POS)
	return lift[0][1]

# input_phenotype = os.path.join('', '/home/bioinfo/PROJECT/diagnosys/RESULT/17_Dec_2023_CLM2/pheno/phenotype')
# phenotype_df = pd.read_csv(input_phenotype, sep='\t', header=0, dtype={'sample': str})

# gene_list = list(phenotype_df['gene'][phenotype_df['sample']==str(name)])

# print(phenotype_df)
# print(gene_list)