#!/usr/bin/env python3



# import pandas as pd
# import numpy as np
# import os
# import dir_tree

# from Pipe import Pipe
# import utils


# from sklearn.linear_model import LogisticRegression
# from sklearn.svm import SVC
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn import tree


# def funcSVN(DF):
#     hmean = 0.471707
#     hstd = 0.067168
#     x2 = 0.340057441652
#     y2= 0.603356263732

#     sample = DF['sample_id'].drop_duplicates()
#     DF = DF[['hgvs','conferma','samtools','gatk','UNBAL','qual',
#             'depth','internal_maf','tipo','sample_id','gene']]

#     DF['UNBAL'] = DF['UNBAL'].astype(float)
#     mask1 = DF['UNBAL'] < x2
#     mask2 = (DF['UNBAL'] > y2) & (DF['UNBAL'] < 0.9)
#     DF.loc[mask1,'UNBALCLASS'] = 'NO'
#     DF.loc[mask2,'UNBALCLASS'] = 'NO'
#     DF['UNBALCLASS'].fillna('SI',inplace=True)
#     DFQUAL = DF[['conferma','gatk']]
#     DFforstats = DF[['conferma','gatk','qual']]
# ########LABEL ENCODING##############################################
#     to_be_encoded_cols = DFQUAL.columns.values
#     utils.label_encode(DFforstats, to_be_encoded_cols)
# ########TEST TRAINING################################################
#     y_col = 'conferma'
#     train_test_ratio = 0.7
#     X_ALL,Y_ALL,X_train,Y_train,X_test,Y_test = utils.get_train_test(DFforstats, y_col, train_test_ratio)
# ########DECISIONAL TREE PROCESS AND GRAPH#########################################
#     _treeall_ = tree.DecisionTreeClassifier(criterion='entropy',max_depth=2,random_state=0)
#     _treeall_.fit(X_ALL,Y_ALL)
#     #print _treeall_.predict([[0,20]])
#     #print _treeall_.predict_proba([[0,20]])
#     #print _treeall_.score([[0,20]],[[1]])
#     return _treeall_.fit(X_ALL,Y_ALL)


# class SangerPredictionPipe(Pipe):

#     def __init__(self):
#         pass
    
#     def addprediction(self, _treeall_):
#         pheno_annot_path = self.sample.pheno_annot
#         sampledata = pd.read_csv(pheno_annot_path,sep='\t',header=0)
#         for index,data in sampledata.iterrows():
#             if data['types'] == 'SVN':
#                 if data['QUAL']: qual=int(data['QUAL'])
#                 else: qual=0
#                 if data['gatk_geno'] == 'homo':gatk = 1
#                 elif data['gatk_geno'] == 'het':gatk = 1
#                 else: gatk=0
#                 if data['DEPTH'] >= 11:
#                     gatk=1
#                 else:
#                     gatk=0
#                     qual=0

#                 predprob = _treeall_.predict_proba([[gatk,qual]])
#                 sampledata.loc[index,'probsanger'] = '{:,.2f}'.format(predprob[0][1])
#             else:
#                 sampledata.loc[index,'probsanger'] = 0
#         return sampledata
	

#     def process(self, **kwargs):

#         self.sample = kwargs["sample"]

#         folder_name = dir_tree.principal_directory.path

#         path_sklearn = os.path.join('/', folder_name.split('/')[1], folder_name.split('/')[2], folder_name.split('/')[3], folder_name.split('/')[4])
#         SVN = pd.read_csv(os.path.join(path_sklearn,'bin/SKLEARN','basedatiSKLEARNSVN_SENZA_varDubbie.csv'),sep='\t',header=0)
#         _treeall_ = funcSVN(SVN)
#         samplepred = self.addprediction(_treeall_)
#         self.sample.sanger_probs = os.path.join(dir_tree.principal_directory.temp_dir.path, "sanger_probs.csv")
#         samplepred.to_csv(self.sample.sanger_probs, sep="\t") 
        
#         self.sample.saveJSON()
        

#         return kwargs










import glob
import os
from pathlib import Path
import pandas as pd
import numpy as np
from joblib import load


def addprediction(sampledata, samplename, model, p_threshold=0.9):
    """
    Per ogni variante SVN con QUAL>0, calcola la probabilità di conferma e la 
    trasforma in un call binario usando la soglia p_threshold.
    probsanger: se 1 nosanger, se 0 sanger
    """
    sampledata['QUAL'] = pd.to_numeric(sampledata['QUAL'], errors='coerce')
    sampledata['DEPTH'] = pd.to_numeric(sampledata['DEPTH'], errors='coerce')

    
    mask = (sampledata['types'] == 'SVN') & (sampledata['QUAL'] > 0)  # Maschera sulle righe da processare

    feats = sampledata.loc[mask, ['QUAL', 'DEPTH']].copy()
    feats['qual'] = np.log(feats['QUAL'])
    feats['depth'] = feats['DEPTH']
    feats = feats[['qual', 'depth']]
    X = feats[['qual', 'depth']].values                               # qui prendi l'array numpy, non il DataFrame

    sampledata['probsanger'] = 0                                      # Se ci sono righe valide, calcola probabilità e poi threshold
    if X.shape[0] > 0:
        probs = model.predict_proba(X)[:, 1] 
        calls = (probs >= p_threshold).astype(int)
        sampledata.loc[mask, 'probsanger'] = calls


    return sampledata


def evaluate_contamination(sampledata, hmean: float = 0.471707, hstd:  float = 0.067168, low_th: float = 0.340057441652,high_th: float = 0.603356263732):
    """
    Lo script filtra il DataFrame alle varianti SNP eterozigoti confermate da Samtools e GATK con QUAL elevato, profondità adeguata e unbalance valido
    Dalla colonna “unbalance” estrae il valore numerico e lo converte in un campo UNBAL in formato float
    Per ogni variante calcola uno zscore normalizzando UNBAL rispetto a una media e deviazione standard predefinite
    Vengono quindi contati i siti con zscore estremi o con UNBAL al di fuori delle soglie attese
    Infine ne ricava la percentuale di varianti anomale, la formatta con due decimali e la assegna alla colonna probcontaminazione
    La percentuale di contaminazione è: (# righe con |zscore| > 2) / (# righe totali) * 100
    """
    mask = (
        (sampledata['gatk_geno']    == 'het') &
        (sampledata['types']        == 'SVN') &
        (pd.to_numeric(sampledata['QUAL'], errors='coerce') >= 60) & # 60 of phred quality score for gatk correpsond to (error)1/1000000 = 0.0001%	(accuracy 1- error) 99.9999%
        (pd.to_numeric(sampledata['DEPTH'], errors='coerce') > 20) &
        sampledata['unbalance'].notna() &
        (sampledata['unbalance'] != 'del=nan')
    )

    unbal = (sampledata.loc[mask, 'unbalance'].str.split('=', n=1).str[1].astype(float))
    z = (unbal - hmean) / hstd
    tot = len(z)
    num_extreme = (z.abs() > 2).sum()
    percz = (num_extreme / tot * 100) if tot else 0.0
    sampledata['probcontaminazione'] = f"{percz:,.2f}"

    return sampledata



def make_mafdecisional(sampledata):
    """
    Returns a DataFrame with two new columns:
      - 'decisionmaf': the chosen MAF value
      - 'decisionINFO': source of the chosen value ('POPMAX', 'POP', 'ExAC', or 'ALL')
    """
    sampledata = sampledata.copy()
    # Prepare Adj_MAF2 by stripping suffixes and excluding certain gnomAD populations
    exclude_pops = ['gnomAD_ASJ', 'gnomAD_FIN', 'gnomAD_OTHER', 'gnomAD_OTH']
    sampledata['Adj_MAF2'] = (sampledata['Adj_MAF'].fillna('unknown').astype(str).str.split('&', n=1).str[0].replace({pop: 'exclude' for pop in exclude_pops}))

    # Compute MAX_MAF2: only for rows where Adj_MAF2 contains 'gnomAD', and drop exact-1 values
    sampledata['MAX_MAF2'] = np.where(sampledata['Adj_MAF2'].str.contains('gnomAD', na=False),sampledata['MAX_MAF'].replace(1, np.nan),np.nan)

    # Build decisionmaf by priority: POPMAX > MAX_MAF2 > ExAC_MAF
    sampledata['decisionmaf'] = (sampledata['gnomAD_exomes_POPMAX_AF'].combine_first(sampledata['MAX_MAF2']).combine_first(sampledata['ExAC_MAF']))

    # Annotate source of decision with vectorized select
    conditions = [sampledata['gnomAD_exomes_POPMAX_AF'].notna(),sampledata['decisionmaf'] == sampledata['MAX_MAF2']]
    choices = ['POPMAX', 'POP']
    sampledata['decisionINFO'] = np.select(conditions, choices, default='ALL')

    return sampledata
##############################################################################

if __name__ == "__main__":
    directory = "/home/alessandro/PROJECT/RESULT/job22_VARIANTSRED_R_q25_90_OCULARE_OCULARE/final"
    pattern = "*_pheno_annot.csv"
    file_list = glob.glob(os.path.join(directory, pattern))

    model = load('/home/alessandro/PROJECT/SKLEARN/logistic_model.pkl')

    for file in file_list:
        filepath = Path(file)
        sample_id = filepath.name.split('_')[0]
        pheno_annot = pd.read_csv(file, sep='\t')

        samplepred = addprediction(pheno_annot, sample_id, model)
        PREsamplefinalnew = evaluate_contamination(samplepred)
        samplefinalnew = make_mafdecisional(PREsamplefinalnew)

        # out_cols = ['sample', 'HGVS', 'types', 'DEPTH', 'QUAL', 'probsanger', 'probcontaminazione']
        samplefinalnew.to_csv('/home/alessandro/PROJECT/RESULT/test.csv', index=False)





    
	# if args.genome == 'geno37':
	# 	pass
	# elif args.genome == 'geno38':
	# 	print('in')
	# 	folder_list = glob.glob(folder_final+'*')
		
	# 	for folder_sample in set(folder_list):
	# 		#print folder_sample
	# 		if 'pheno_annot' in folder_sample:
	# 			#print folder_sample
	# 			sample_x = folder_sample.split('/')[-1]
	# 			sample_x = str(sample_x)
	# 			print( '-----'+sample_x+'------')
	# 			sampledata = pd.read_csv(folder_sample,sep='\t',header=0)
	# 			print ('Start Data',len(sampledata))
	# 			_treeall_ = funcSVN(SVN)
	# 			#########################################################
	# 			samplepred = addprediction(sampledata,sample_x,_treeall_)
				##################################
				# PREsamplefinalnew = evaluate_contamination(samplepred,sample_x)
				##################################
				# samplefinalnew = make_mafdecisional(PREsamplefinalnew,sample_x)
	# 			##################################
	# 			# if args.dest == 'r':
	# 			# 	dest = 'rovereto'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/MAGIS/NGS_RESULT/annot'
	# 			# elif args.dest == 'b':
	# 			# 	dest = 'bolzano'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/EUREGIO/NGS_RESULT/annot'
	# 			# elif args.dest == 's':
	# 			# 	dest = 'sanfelice'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/SANFELICE/NGS_RESULT/annot'
	# 			# elif args.dest == 'z':
	# 			# 	dest = 'ricerca'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/RICERCA/NGS_RESULT/annot'

	# 			sample_ = sample_x.split('_')[0]
	# 			result1 = join(folder_final,sample_+'_pheno_predict.csv')
	# 			# result_django1 = join(path_django,sample_+'_pheno_predict.csv')
	# 			result2 = join(folder_final,sample_+'_other_annot.csv')
	# 			# result_django2 = join(path_django,sample_+'_other_annot.csv')

	# 			print ('End Data',len(samplefinalnew))
	# 			#print samplefinalnew.columns
	# 			if len(samplefinalnew) > 0: samplefinalnew1 = samplefinalnew[cols1]
	# 			else: samplefinalnew1 = pd.DataFrame(columns=cols1)
	# 			#samplefinalnew2 = samplefinalnew[cols2]
	# 			samplefinalnew2 = samplefinalnew[((samplefinalnew['consequence']!='synonymous_variant') & (samplefinalnew['consequence']!='5_prime_UTR_variant') & (samplefinalnew['consequence']!='upstream_gene_variant') & (samplefinalnew['consequence']!='downstream_gene_variant') & (samplefinalnew['consequence']!='3_prime_UTR_variant') & (samplefinalnew['consequence']!='intergenic_variant') & (samplefinalnew['consequence']!='non_coding_transcript_exon_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='intron_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='regulatory_region_variant') & (samplefinalnew['consequence']!='splice_region_variant&intron_variant') & (samplefinalnew['consequence']!='splice_region_variant&intron_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='splice_region_variant&synonymous_variant'))][cols2]
	# 			samplefinalnew2bis = samplefinalnew2.fillna(-999)
	# 			samplefinalnew1['sample'] = sample_
	# 			samplefinalnew2bis['sample'] = sample_
	# 			samplefinalnew1.to_csv(result1,sep='\t',index=False,encoding='utf-8')
	# 			#samplefinalnew1.to_csv(result_django1,sep='\t',index=False,encoding='utf-8')
	# 			samplefinalnew2bis.to_csv(result2,sep='\t',index=False,encoding='utf-8')
	# 			#samplefinalnew2bis.to_csv(result_django2,sep='\t',index=False,encoding='utf-8')

	# 			#system(' '.join(['scp',result1,"remo@192.168.4.230:/home/remo/venv/apimagi_dev/NGS_RESULT/annot"]))
	# 			#system(' '.join(['scp',result2,"remo@192.168.4.230:/home/remo/venv/apimagi_dev/NGS_RESULT/annot"]))
	# 			#system(' '.join(['scp',result1,"remo@192.168.4.230:/home/remo/venv/apimagi_prod/NGS_RESULT/annot"]))
	# 			#system(' '.join(['scp',result2,"remo@192.168.4.230:/home/remo/venv/apimagi_prod/NGS_RESULT/annot"]))

	# 		if 'family_annot_pre1' in folder_sample:
	# 			#print folder_sample
	# 			sample_x = folder_sample.split('/')[-1]
	# 			sample_x = str(sample_x)
	# 			print ('-----'+sample_x+'------')
	# 			sampledata = pd.read_csv(folder_sample,sep='\t',header=0)
	# 			print ('Start Data',len(sampledata))
	# 			_treeall_ = funcSVN(SVN)
	# 			###############################################################
	# 			samplepred = addprediction(sampledata,sample_x,_treeall_)
	# 			###############################################################
	# 			PREsamplefinalnew = evaluate_contamination(samplepred,sample_x)
	# 			##################################
	# 			samplefinalnew = make_mafdecisional(PREsamplefinalnew,sample_x)
	# 			###############################################################
	# 			# if args.dest == 'r':
	# 			# 	dest = 'rovereto'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/MAGIS/NGS_RESULT/annot'
	# 			# elif args.dest == 'b':
	# 			# 	dest = 'bolzano'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/EUREGIO/NGS_RESULT/annot'
	# 			# elif args.dest == 's':
	# 			# 	dest = 'sanfelice'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/SANFELICE/NGS_RESULT/annot'
	# 			# elif args.dest == 'z':
	# 			# 	dest = 'ricerca'
	# 			# 	path_django = '/home/bioinfo/VIRTUAL/RICERCA/NGS_RESULT/annot'
	# 			sample_ = sample_x.split('_')[0]
	# 			result1 = join(folder_final,sample_+'_family_annot_pre.csv')
	# 			# result_django1 = join(path_django,sample_+'_family_annot.csv')
	# 			print ('End Data',len(samplefinalnew))
	# 			#print samplefinalnew.columns
	# 			try:
	# 				samplefinalnew1 = samplefinalnew[cols3]	#[((samplefinalnew['consequence']!='intron_variant') & (samplefinalnew['consequence']!='synonymous_variant') & (samplefinalnew['consequence']!='5_prime_UTR_variant') & (samplefinalnew['consequence']!='upstream_gene_variant') & (samplefinalnew['consequence']!='downstream_gene_variant') & (samplefinalnew['consequence']!='3_prime_UTR_variant') & (samplefinalnew['consequence']!='intergenic_variant') & (samplefinalnew['consequence']!='non_coding_transcript_exon_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='intron_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='regulatory_region_variant') & (samplefinalnew['consequence']!='splice_region_variant&intron_variant') & (samplefinalnew['consequence']!='splice_region_variant&intron_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant') & (samplefinalnew['consequence']!='splice_region_variant&synonymous_variant'))][cols3]
	# 			except: samplefinalnew1 = pd.DataFrame(columns=cols3)
	# 			samplefinalnew1['strand'] = samplefinalnew1['strand'].astype(int)
	# 			samplefinalnew1bis = samplefinalnew1.fillna(-999)
	# 			samplefinalnew1bis['sample'] = sample_
	# 			samplefinalnew1bis.to_csv(result1,sep='\t',index=False,encoding='utf-8')
	# 			#samplefinalnew1bis.to_csv(result_django1,sep='\t',index=False,encoding='utf-8')
	# 			#system(' '.join(['scp',result1, "remo@192.168.4.230:/home/remo/venv/apimagi_dev/NGS_RESULT/annot"]))
	# 			#system(' '.join(['scp',result1, "remo@192.168.4.230:/home/remo/venv/apimagi_prod/NGS_RESULT/annot"]))
