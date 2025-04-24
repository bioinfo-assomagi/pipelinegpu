import pandas as pd
import numpy as np

from Pipe import Pipe


def make_mafdecisional(samplefinalnew,sample_x):
	staralg = samplefinalnew	#[['gnomAD_exomes_POPMAX_AF','MAX_MAF','Adj_MAF','ExAC_MAF']]
	staralg['Adj_MAF2'] = staralg['Adj_MAF']
	staralg['Adj_MAF2'].fillna('unknown',inplace=True)
	staralg['Adj_MAF2'] = staralg['Adj_MAF2'].astype(str)
	staralg['Adj_MAF2'] = staralg['Adj_MAF2'].str.split('&').str.get(0)
	staralg['Adj_MAF2'].replace('gnomAD_ASJ','exclude',inplace=True)
	staralg['Adj_MAF2'].replace('gnomAD_FIN','exclude',inplace=True)
	staralg['Adj_MAF2'].replace('gnomAD_OTHER','exclude',inplace=True)
	staralg['Adj_MAF2'].replace('gnomAD_OTH','exclude',inplace=True)
	staralg['MAX_MAF'].fillna(-999,inplace=True)
	mask1 = staralg['Adj_MAF2'].str.contains('gnomAD')
	staralg.loc[mask1,'MAX_MAF2'] = staralg['MAX_MAF']
	mask2 = staralg['MAX_MAF2'] == 1
	staralg.loc[mask2,'MAX_MAF2'] = np.nan
	staralg['decisionINFO'] = None
	staralg['decisionmaf'] = staralg['gnomAD_exomes_POPMAX_AF']
	staralg['decisionINFO'] = np.where(staralg['decisionmaf'].notnull(),'POPMAX',staralg['decisionmaf'])
	staralg['decisionmaf'].fillna(staralg['MAX_MAF2'],inplace=True)
	staralg['decisionINFO'] = np.where(staralg['decisionmaf']==staralg['MAX_MAF2'],'POP',staralg['decisionINFO'])
	staralg['decisionmaf'].fillna(staralg['ExAC_MAF'],inplace=True)
	staralg['decisionINFO'].replace('nan','ALL',inplace=True)
	print (staralg[['decisionmaf','MAX_MAF2','gnomAD_exomes_POPMAX_AF']])
	return staralg

""" Pipe which, having a VEP annotated VCF, is responsible for selecting which MAF to use. """
class MAFDecisionPipe(Pipe):

    def __init__(self):
        pass

    def process(self, **kwargs):
        self.sample = kwargs["sample"]
        

        return kwargs