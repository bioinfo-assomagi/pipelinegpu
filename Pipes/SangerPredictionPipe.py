import pandas as pd
import numpy as np
import os
import dir_tree

from Pipe import Pipe
import utils


from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn import tree

# test git ale
def funcSVN(DF):
    hmean = 0.471707
    hstd = 0.067168
    x2 = 0.340057441652
    y2= 0.603356263732

    sample = DF['sample_id'].drop_duplicates()
    DF = DF[['hgvs','conferma','samtools','gatk','UNBAL','qual',
            'depth','internal_maf','tipo','sample_id','gene']]

    DF['UNBAL'] = DF['UNBAL'].astype(float)
    mask1 = DF['UNBAL'] < x2
    mask2 = (DF['UNBAL'] > y2) & (DF['UNBAL'] < 0.9)
    DF.loc[mask1,'UNBALCLASS'] = 'NO'
    DF.loc[mask2,'UNBALCLASS'] = 'NO'
    DF['UNBALCLASS'].fillna('SI',inplace=True)
    DFQUAL = DF[['conferma','gatk']]
    DFforstats = DF[['conferma','gatk','qual']]
########LABEL ENCODING##############################################
    to_be_encoded_cols = DFQUAL.columns.values
    utils.label_encode(DFforstats, to_be_encoded_cols)
########TEST TRAINING################################################
    y_col = 'conferma'
    train_test_ratio = 0.7
    X_ALL,Y_ALL,X_train,Y_train,X_test,Y_test = utils.get_train_test(DFforstats, y_col, train_test_ratio)
########DECISIONAL TREE PROCESS AND GRAPH#########################################
    _treeall_ = tree.DecisionTreeClassifier(criterion='entropy',max_depth=2,random_state=0)
    _treeall_.fit(X_ALL,Y_ALL)
    #print _treeall_.predict([[0,20]])
    #print _treeall_.predict_proba([[0,20]])
    #print _treeall_.score([[0,20]],[[1]])
    return _treeall_.fit(X_ALL,Y_ALL)


class SangerPredictionPipe(Pipe):

    def __init__(self):
        pass
    
    def addprediction(self, _treeall_):
        pheno_annot_path = self.sample.pheno_annot
        sampledata = pd.read_csv(pheno_annot_path,sep='\t',header=0)
        for index,data in sampledata.iterrows():
            if data['types'] == 'SVN':
                if data['QUAL']: qual=int(data['QUAL'])
                else: qual=0
                if data['gatk_geno'] == 'homo':gatk = 1
                elif data['gatk_geno'] == 'het':gatk = 1
                else: gatk=0
                if data['DEPTH'] >= 11:
                    gatk=1
                else:
                    gatk=0
                    qual=0

                predprob = _treeall_.predict_proba([[gatk,qual]])
                sampledata.loc[index,'probsanger'] = '{:,.2f}'.format(predprob[0][1])
            else:
                sampledata.loc[index,'probsanger'] = 0
        return sampledata
	

    def process(self, **kwargs):

        self.sample = kwargs["sample"]

        folder_name = dir_tree.principal_directory.path

        path_sklearn = os.path.join('/', folder_name.split('/')[1], folder_name.split('/')[2], folder_name.split('/')[3], folder_name.split('/')[4])
        SVN = pd.read_csv(os.path.join(path_sklearn,'bin/SKLEARN','basedatiSKLEARNSVN_SENZA_varDubbie.csv'),sep='\t',header=0)
        _treeall_ = funcSVN(SVN)
        samplepred = self.addprediction(_treeall_)
        self.sample.sanger_probs = os.path.join(dir_tree.principal_directory.temp_dir.path, "sanger_probs.csv")
        samplepred.to_csv(self.sample.sanger_probs, sep="\t") 
        
        self.sample.saveJSON()
        

        return kwargs

