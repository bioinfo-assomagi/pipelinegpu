#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script per addestrare il modello SVM e salvarlo su disco.
Il modello viene addestrato sul dataset fornito e salvato in "svm_model.pkl".
"""

import os
import pandas as pd
import numpy as np
import pickle
from sklearn.svm import SVC

def train_svm_model(training_file, output_model_path):
    # Legge il dataset di training
    df = pd.read_csv(training_file, sep='\t', header=0)
    df.columns = df.columns.str.lower().str.strip()

    # Seleziona le colonne necessarie per il modello: 'conferma', 'gatk', 'qual'
    required_cols = ['conferma', 'gatk', 'qual']
    df = df[required_cols].copy()
    
    # Converte 'qual' in numerico, sostituendo errori con 0
    df['qual'] = pd.to_numeric(df['qual'], errors='coerce').fillna(0).astype(int)
    
    # Mappa il campo 'gatk': "homo" o "het" (case-insensitive) → 1, altrimenti 0.
    df['gatk'] = df['gatk'].astype(str).str.lower().map(lambda x: 1 if x in ['homo', 'het'] else 0)
    
    # Codifica la colonna 'conferma' in valori numerici (ad esempio: "NO" → 0, "SI" → 1) utilizzando factorize
    df['conferma'], uniques = pd.factorize(df['conferma'])
    
    # Costruisce il dataset delle feature (X) e l'etichetta (y)
    X = df[['gatk', 'qual']].values
    y = df['conferma'].values
    
    # Inizializza il modello SVM con probability=True
    svm_model = SVC(probability=True, random_state=0)
    svm_model.fit(X, y)
    
    # Salva il modello tramite pickle
    with open(output_model_path, 'wb') as f:
        pickle.dump(svm_model, f)
    print(f"Modello SVM addestrato e salvato in: {output_model_path}")

if __name__ == "__main__":
    # Specifica il percorso del file di training.
    # Modifica questo percorso in base al tuo ambiente; ad es., potrebbe essere:
    # training_file = "/home/alessandro/PROJECT/SKLEARN/basedatiSKLEARNSVN_SENZA_varDubbie.csv"
    training_file = "/home/alessandro/PROJECT/SKLEARN/SKLEARN_training.csv"
    output_model_path = "svm_model.pkl"
    train_svm_model(training_file, output_model_path)