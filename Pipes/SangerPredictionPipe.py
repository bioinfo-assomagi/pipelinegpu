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



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sanger Prediction Pipe con modello SVM pre-addestrato.
Carica il modello SVM salvato (svm_model.pkl) e lo usa per elaborare
i dati del campione e predire le probabilità di conferma Sanger.
"""


import os
import pandas as pd
import numpy as np
import pickle
import dir_tree  # Modulo della tua pipeline per gestire la struttura delle cartelle
from Pipe import Pipe  # Classe base delle pipe nella tua pipeline

def addprediction(df, svm_model):
    """
    Calcola in maniera vettorializzata le probabilità di conferma Sanger utilizzando il modello SVM pre-addestrato.
    
    Parametri:
      df : DataFrame con le annotazioni fenotipiche (si assume contenga almeno le colonne 'tipo', 'gatk', 'qual', 'depth')
      svm_model : modello SVM addestrato con probability=True
      
    Restituisce:
      df : DataFrame aggiornato con la colonna 'probsanger'
    """
    df = df.copy()
    df.columns = df.columns.str.lower().str.strip()
    
    # Converte 'qual' e 'depth' in formato numerico
    df['qual'] = pd.to_numeric(df['qual'], errors='coerce').fillna(0).astype(int)
    df['depth'] = pd.to_numeric(df['depth'], errors='coerce').fillna(0).astype(int)
    
    # Mappa il campo 'gatk': "homo" o "het" → 1, altrimenti 0.
    df['processed_gatk'] = df['gatk'].astype(str).str.lower().map(lambda x: 1 if x in ['homo', 'het'] else 0)
    
    # Se la profondità è inferiore a 11, imposta 'processed_gatk' e 'qual' a 0.
    shallow_mask = df['depth'] < 11
    df.loc[shallow_mask, 'processed_gatk'] = 0
    df.loc[shallow_mask, 'qual'] = 0
    
    # Inizializza la colonna dei risultati
    df['probsanger'] = 0.0
    
    # Applica la predizione solo sulle righe in cui 'tipo' è "SVN"
    svn_mask = df['tipo'].str.upper() == 'SVN'
    if svn_mask.any():
        features = df.loc[svn_mask, ['processed_gatk', 'qual']].values
        pred_probs = svm_model.predict_proba(features)[:, 1]  # probabilità della classe 1
        df.loc[svn_mask, 'probsanger'] = [f"{p:.2f}" for p in pred_probs]
        
    return df

class SangerPredictionPipe(Pipe):
    def __init__(self):
        super().__init__()
    
    def process(self, **kwargs):
        """
        Esegue il processamento del campione:
          - Carica il modello SVM pre-addestrato dal file "svm_model.pkl".
          - Legge il file di annotazione del campione (pheno_annot).
          - Calcola le probabilità di conferma Sanger.
          - Salva il risultato in "sanger_probs.csv" nella cartella temporanea definita in dir_tree.
        """
        self.sample = kwargs["sample"]
        
        # Costruisce il percorso del modello pre-addestrato (assunto nella stessa cartella dello script)
        model_path = os.path.join(os.path.dirname(__file__), "svm_model.pkl")
        with open(model_path, 'rb') as f:
            svm_model = pickle.load(f)
        
        # Legge il file di annotazione del campione (pheno_annot)
        sampledata = pd.read_csv(self.sample.pheno_annot, sep='\t', header=0)
        sample_pred = addprediction(sampledata, svm_model)
        
        # Definisce il percorso di output per il file dei risultati e lo salva
        output_path = os.path.join(dir_tree.principal_directory.temp_dir.path, "sanger_probs.csv")
        sample_pred.to_csv(output_path, sep="\t", index=False)
        self.sample.sanger_probs = output_path
        
        self.sample.saveJSON()  # Aggiorna lo stato del sample
        return kwargs

if __name__ == "__main__":
    # Esempio di utilizzo stand-alone della pipe.
    class Sample:
        def __init__(self, pheno_annot):
            self.pheno_annot = pheno_annot
            self.sanger_probs = None
        def saveJSON(self):
            print("Sample salvato.")
            print("File pheno_annot:", self.pheno_annot)
            print("File sanger_probs:", self.sanger_probs)
    
    # Percorso del file CSV di input contenente le annotazioni fenotipiche
    sample_input_path = "/home/alessandro/PROJECT/SKLEARN/SKLEARN.csv"
    sample_obj = Sample(pheno_annot=sample_input_path)
    
    pipe = SangerPredictionPipe()
    pipe.process(sample=sample_obj)