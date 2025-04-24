import numpy as np
import pandas as pd


import os

from Pipe import Pipe


""" Takes as an input the SangerPredictionPipe output. The sanger_probs field of a Sample object. """
def evaluate_contamination(sampledata):
	hmean = 0.471707
	hstd = 0.067168
	x2 = 0.340057441652
	y2= 0.603356263732

	#print sampledata[['types','QUAL','DEPTH','samtools_geno']].dtypes
	try:
		sampledata_correct = sampledata[(sampledata['samtools_geno']=='het') & (sampledata['gatk_geno']=='het') &\
			(sampledata['types']=='SVN') & (sampledata['QUAL']>=18) & (sampledata['DEPTH']>20) & (sampledata['unbalance'] != 'del=nan')]
		sampledata_correct.loc[:,'UNBAL'] = sampledata_correct['unbalance'].str.split('=').str.get(1).astype(float)
		xmean = sampledata_correct['UNBAL'].mean()
		xstd = sampledata_correct['UNBAL'].std()
		sampledata_correct['zscore'] = (sampledata_correct['UNBAL']-hmean)/hstd
		count = len(sampledata_correct)
		menocount = len(sampledata_correct[sampledata_correct['UNBAL']<x2])
		piucount = len(sampledata_correct[sampledata_correct['UNBAL']>y2])
		totcount = menocount+piucount

		#perccount = (float(totcount)/float(count))*100
		try: perccount = (float(totcount)/float(count))*100
		except: perccount = 0

		menozscore = len(sampledata_correct[sampledata_correct['zscore']<-2])
		piuzscore =  len(sampledata_correct[sampledata_correct['zscore']>+2])
		totzscore = menozscore+piuzscore
		#perczscore = (float(totzscore)/float(count))*100
		try: perczscore = (float(totzscore)/float(count))*100
		except: perczscore = 0

		#print totcount,count,perccount,perczscore
	  	#sampledata['xmean'] = '{:,.2f}'.format(xmean)
	  	#sampledata['xstd'] = '{:,.2f}'.format(xstd)
	  	#sampledata['count'] = '{:,.2f}'.format(count)
	  	#sampledata['menocount'] = '{:,.2f}'.format(menocount)
	  	#sampledata['piucount'] = '{:,.2f}'.format(piucount)
	  	#sampledata['totcount'] = '{:,.2f}'.format(totcount)
	  	#sampledata['perccount'] = '{:,.2f}'.format(perccount)
	  	#sampledata['menozscore'] = '{:,.2f}'.format(menozscore)
	  	#sampledata['piuzscore'] = '{:,.2f}'.format(piuzscore)
	  	#sampledata['totzscore'] = '{:,.2f}'.format(totzscore)
	  	#sampledata['perczscore'] = '{:,.2f}'.format(perczscore)
		sampledata['probcontaminazione'] = '{:,.2f}'.format(perczscore)

	except:
		sampledata['probcontaminazione'] = 0

	#print sampledata['probcontaminazione']
	return sampledata

class ContaminationPipe(Pipe):
    """ 
	This pipe is responsible for doing the contamination analysis. The idea behind relies on the fact that
	DNA sample can be contaminated from external factors during sample handling, logistics, etc. 
	It compares the distribution of allele frequencies in the present sample, coming from the results
	of the coverage analysis (i.e. ``samtools mpileup``, ``count_unbalance()``, etc.)
	"""

    def __init__(self):
        pass

    def process(self, **kwargs):
        
        self.sample = kwargs["sample"]
        sanger_probs = self.sample.sanger_probs
        sample_probs_df = pd.read_csv(sanger_probs, sep="\t", header=0)
        contamination_data = evaluate_contamination(sample_probs_df)
        
        # TODO: think of a way to make this pipes modular. Removal of one should not break the flow.
        # Overwrite the sanger_probs would do for now
        contamination_data.to_csv(self.sample.sanger_probs, sep="\t")

        return kwargs