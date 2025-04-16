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




# def PrintLog(command,folder):
#     ###Print all commands in a log file
# 	path = join(folder,'Commands4.log')
# 	cl = open(path,'a')
# 	cl.write(command)
# 	cl.write('\n')
# 	cl.close()

# def sample(name):
# 	x = name.split("/")
# 	x = (x[-1])
# 	x = x.split("_")
# 	sample = str(x[0])
# 	return str(sample)
# #############################################################################
# def label_encode(df, columns):
#     for col in columns:
#         le = LabelEncoder()
#         col_values_unique = list(df[col].unique())
#         le_fitted = le.fit(col_values_unique)
#         col_values = list(df[col].values)
#         le.classes_
#         col_values_transformed = le.transform(col_values)
#         df[col] = col_values_transformed
# #############################################################################
# def get_train_test(df, y_col, ratio):
# 	mask = np.random.rand(len(df)) < ratio
# 	df_train = df[mask]
# 	df_test = df[~mask]
# 	Y_ALL = df[y_col].values
# 	Y_train = df_train[y_col].values
# 	Y_test = df_test[y_col].values
# 	del df_train[y_col]
# 	del df_test[y_col]
# 	del df[y_col]
# 	X_ALL = df.values
# 	X_train = df_train.values
# 	X_test = df_test.values
# 	return X_ALL,Y_ALL,X_train, Y_train, X_test, Y_test
# #############################################################################
# def plot_decision_region(X,y, classifier,test_idx=None,resolution=0.02):
# 	markers = ('o','x','v','^','s')
# 	colors = ('orange','cyan','lightgreen','green','blue')
# 	colors2 = ('red','green','lightgreen','blue','blue')
# 	cmap = ListedColormap(colors[:len(np.unique(y))])
# 	cmap2 = ListedColormap(colors2[:len(np.unique(y))])
# 	x1_min,x1_max= X[:,0].min() - 0.5, X[:,0].max() +0.5
# 	x2_min,x2_max = X[:,1].min() - 0.5, X[:,1].max() +0.5
# 	#x3_min,x3_max = X[:,2].min() - 0.5, X[:,2].max() +0.5
# 	xx1, xx2, = np.meshgrid(np.arange(x1_min,x1_max,resolution), np.arange(x2_min,x2_max,resolution))
# 	Z = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
# 	Z = Z.reshape(xx1.shape)
# 	###################################################
# 	plt.figure(1,figsize=(10,8))
# 	#plt.pcolormesh(xx1,xx2,Z,alpha=0.4,cmap=plt.cm.Paired)
# 	plt.contourf(xx1,xx2,Z,alpha=0.6,cmap=cmap)
# 	plt.xlim(xx1.min(),xx1.max())
# 	plt.ylim(xx2.min(),xx2.max())
# 	#plt.xticks(())
# 	#plt.yticks(())
# 	X_test, y_test = X[test_idx,:],y[test_idx]
# 	for idx, cl in enumerate(np.unique(y)):
# 		plt.scatter(x=X[y==cl,0],y=X[y==cl,1],alpha=0.4,c=cmap2(idx), marker=markers[idx],label=cl)
# 	if test_idx:
# 		x_test, y_test = X[test_idx,:], y[test_idx]
# 		plt.scatter(X_test[:,0],X_test[:,1],c='',alpha=1.0,linewidth=1,
# 				marker='o',s=90,label='test set')

# 	plt.xlabel('PCA')
# 	plt.ylabel('Conferma Sanger')
# 	plt.legend(loc='best')
# 	plt.show()
# #############################################################################
# def funcSVN(DF):
# 	hmean = 0.471707
# 	hstd = 0.067168
# 	x2 = 0.340057441652
# 	y2= 0.603356263732

# 	sample = DF['sample_id'].drop_duplicates()
# 	DF = DF[['hgvs','conferma','samtools','gatk','UNBAL','qual',
# 			'depth','internal_maf','tipo','sample_id','gene']]

# 	DF['UNBAL'] = DF['UNBAL'].astype(float)
# 	mask1 = DF['UNBAL'] < x2
# 	mask2 = (DF['UNBAL'] > y2) & (DF['UNBAL'] < 0.9)
# 	DF.loc[mask1,'UNBALCLASS'] = 'NO'
# 	DF.loc[mask2,'UNBALCLASS'] = 'NO'
# 	DF['UNBALCLASS'].fillna('SI',inplace=True)
# 	DFQUAL = DF[['conferma','gatk']]
# 	DFforstats = DF[['conferma','gatk','qual']]
# ########LABEL ENCODING##############################################
# 	to_be_encoded_cols = DFQUAL.columns.values
# 	label_encode(DFforstats, to_be_encoded_cols)
# ########TEST TRAINING################################################
# 	y_col = 'conferma'
# 	train_test_ratio = 0.7
# 	X_ALL,Y_ALL,X_train,Y_train,X_test,Y_test = get_train_test(DFforstats, y_col, train_test_ratio)
# ########DECISIONAL TREE PROCESS AND GRAPH#########################################
# 	_treeall_ = tree.DecisionTreeClassifier(criterion='entropy',max_depth=2,random_state=0)
# 	_treeall_.fit(X_ALL,Y_ALL)
# 	#print _treeall_.predict([[0,20]])
# 	#print _treeall_.predict_proba([[0,20]])
# 	#print _treeall_.score([[0,20]],[[1]])
# 	return _treeall_.fit(X_ALL,Y_ALL)
#############################################################################
# def evaluate_contamination(sampledata,samplename):
# 	hmean = 0.471707
# 	hstd = 0.067168
# 	x2 = 0.340057441652
# 	y2= 0.603356263732

# 	#print sampledata[['types','QUAL','DEPTH','samtools_geno']].dtypes
# 	try:
# 		sampledata_correct = sampledata[(sampledata['samtools_geno']=='het') & (sampledata['gatk_geno']=='het') &\
# 			(sampledata['types']=='SVN') & (sampledata['QUAL']>=18) & (sampledata['DEPTH']>20) & (sampledata['unbalance'] != 'del=nan')]
# 		sampledata_correct.loc[:,'UNBAL'] = sampledata_correct['unbalance'].str.split('=').str.get(1).astype(float)
# 		xmean = sampledata_correct['UNBAL'].mean()
# 		xstd = sampledata_correct['UNBAL'].std()
# 		sampledata_correct['zscore'] = (sampledata_correct['UNBAL']-hmean)/hstd
# 		count = len(sampledata_correct)
# 		menocount = len(sampledata_correct[sampledata_correct['UNBAL']<x2])
# 		piucount = len(sampledata_correct[sampledata_correct['UNBAL']>y2])
# 		totcount = menocount+piucount

# 		#perccount = (float(totcount)/float(count))*100
# 		try: perccount = (float(totcount)/float(count))*100
# 		except: perccount = 0

# 		menozscore = len(sampledata_correct[sampledata_correct['zscore']<-2])
# 		piuzscore =  len(sampledata_correct[sampledata_correct['zscore']>+2])
# 		totzscore = menozscore+piuzscore
# 		#perczscore = (float(totzscore)/float(count))*100
# 		try: perczscore = (float(totzscore)/float(count))*100
# 		except: perczscore = 0

# 		#print totcount,count,perccount,perczscore
# 	  	#sampledata['xmean'] = '{:,.2f}'.format(xmean)
# 	  	#sampledata['xstd'] = '{:,.2f}'.format(xstd)
# 	  	#sampledata['count'] = '{:,.2f}'.format(count)
# 	  	#sampledata['menocount'] = '{:,.2f}'.format(menocount)
# 	  	#sampledata['piucount'] = '{:,.2f}'.format(piucount)
# 	  	#sampledata['totcount'] = '{:,.2f}'.format(totcount)
# 	  	#sampledata['perccount'] = '{:,.2f}'.format(perccount)
# 	  	#sampledata['menozscore'] = '{:,.2f}'.format(menozscore)
# 	  	#sampledata['piuzscore'] = '{:,.2f}'.format(piuzscore)
# 	  	#sampledata['totzscore'] = '{:,.2f}'.format(totzscore)
# 	  	#sampledata['perczscore'] = '{:,.2f}'.format(perczscore)
# 		sampledata['probcontaminazione'] = '{:,.2f}'.format(perczscore)

# 	except:
# 		sampledata['probcontaminazione'] = 0

# 	#print sampledata['probcontaminazione']
# 	return sampledata


# def make_mafdecisional(samplefinalnew,sample_x):
# 	staralg = samplefinalnew	#[['gnomAD_exomes_POPMAX_AF','MAX_MAF','Adj_MAF','ExAC_MAF']]
# 	staralg['Adj_MAF2'] = staralg['Adj_MAF']
# 	staralg['Adj_MAF2'].fillna('unknown',inplace=True)
# 	staralg['Adj_MAF2'] = staralg['Adj_MAF2'].astype(str)
# 	staralg['Adj_MAF2'] = staralg['Adj_MAF2'].str.split('&').str.get(0)
# 	staralg['Adj_MAF2'].replace('gnomAD_ASJ','exclude',inplace=True)
# 	staralg['Adj_MAF2'].replace('gnomAD_FIN','exclude',inplace=True)
# 	staralg['Adj_MAF2'].replace('gnomAD_OTHER','exclude',inplace=True)
# 	staralg['Adj_MAF2'].replace('gnomAD_OTH','exclude',inplace=True)
# 	staralg['MAX_MAF'].fillna(-999,inplace=True)
# 	mask1 = staralg['Adj_MAF2'].str.contains('gnomAD')
# 	staralg.loc[mask1,'MAX_MAF2'] = staralg['MAX_MAF']
# 	mask2 = staralg['MAX_MAF2'] == 1
# 	staralg.loc[mask2,'MAX_MAF2'] = np.nan
# 	staralg['decisionINFO'] = None
# 	staralg['decisionmaf'] = staralg['gnomAD_exomes_POPMAX_AF']
# 	staralg['decisionINFO'] = np.where(staralg['decisionmaf'].notnull(),'POPMAX',staralg['decisionmaf'])
# 	staralg['decisionmaf'].fillna(staralg['MAX_MAF2'],inplace=True)
# 	staralg['decisionINFO'] = np.where(staralg['decisionmaf']==staralg['MAX_MAF2'],'POP',staralg['decisionINFO'])
# 	staralg['decisionmaf'].fillna(staralg['ExAC_MAF'],inplace=True)
# 	staralg['decisionINFO'].replace('nan','ALL',inplace=True)
# 	print (staralg[['decisionmaf','MAX_MAF2','gnomAD_exomes_POPMAX_AF']])
# 	return staralg

# def evaluate_sex(samplefinalnew,samplename):
# 	print ('samplename')
# 	pass

from joblib import load
import os
import glob
import pandas as pd
from pathlib import Path
import numpy as np







import glob
import os
from pathlib import Path
import pandas as pd
import numpy as np
from joblib import load


def addprediction(sampledata, samplename, model):
    sampledata['QUAL'] = pd.to_numeric(sampledata['QUAL'], errors='coerce')
    sampledata['DEPTH'] = pd.to_numeric(sampledata['DEPTH'], errors='coerce')

    # Crea una maschera per selezionare le righe di interesse (tipicamente, variant type SVN con QUAL > 0)
    mask = (sampledata['types'] == 'SVN') & (sampledata['QUAL'] > 0)
    
    # Prepara il DataFrame delle feature per la predizione
    # Nota: il modello è stato addestrato con colonne in minuscolo: 'qual' (logaritmo) e 'depth'
    features = sampledata.loc[mask, ['QUAL', 'DEPTH']].copy()
    features['qual'] = np.log(features['QUAL'])
    features['depth'] = features['DEPTH']
    features = features[['qual', 'depth']]  # Riordina le colonne in modo che corrispondano
       
    # Inizializza la colonna 'probsanger' a 0 per tutte le righe
    sampledata['probsanger'] = 0
    
    # Se esistono righe selezionate, calcola le probabilità tramite il modello
    if not features.empty:
        # Il modello è una pipeline: le feature in input verranno trasformate e poi passate al logistic
        probs = model.predict_proba(features)[:, 1]  
        sampledata.loc[mask, 'probsanger'] = np.round(probs, 2)
    
    # Rinomina le colonne in minuscolo (se necessario) per mantenere la consistenza
    sampledata.rename(columns={'QUAL': 'qual', 'DEPTH': 'depth'}, inplace=True)
    return sampledata

##############################################################################
if __name__ == "__main__":
    directory = "/home/alessandro/PROJECT/RESULT/job22_VARIANTSRED_R_q25_90_OCULARE_OCULARE/final"
    pattern = "*_pheno_annot.csv"
    file_list = glob.glob(os.path.join(directory, pattern))

    file_pkl = '/home/alessandro/PROJECT/SKLEARN/logistic_model.pkl'
    model = load(file_pkl)

    if file_list:
        for file in file_list:
            filepath = Path(file)
            sample_id = filepath.name.split('_')[0]  # esempio: '415.2024_pheno_annot.csv' → '415.2024'
            pheno_annot = pd.read_csv(file, sep='\t')
            samplepred = addprediction(pheno_annot, sample_id, model)
            # Stampo le colonne desiderate:
            # Uso le colonne originali "QUAL" e "DEPTH" (invece di "qual" e "depth") per evitare KeyError.
            print(samplepred[['sample', 'HGVS', 'GENE', 'qual', 'depth', 'probsanger']])
            samplepred[['sample','HGVS','types','depth' ,'qual','probsanger']].to_csv('/home/alessandro/PROJECT/RESULT/test.csv', index = False)





    
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
	# 			samplefinalnew = make_mafdecisional(PREsamplefinalnew,sample_x)
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
