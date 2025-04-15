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
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns


def train_logistic_model(training_file, output_model_path, test_size=0.2, random_state=0):
    # Legge il dataset di training (modifica il separatore se necessario)
    df = pd.read_csv(training_file, header=0)  # se il file è un CSV standard con separatore virgola
    print("Columns found in training file:", df.columns.tolist())
    
    # Normalizza i nomi delle colonne (minuscolo e senza spazi)
    df.columns = df.columns.str.lower().str.strip()
    
    # Lista delle colonne richieste
    required_cols = ['conferma', 'qual', 'depth']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise KeyError(f"Le seguenti colonne mancano nel file di training: {missing}")
    df = df[required_cols].copy()
    
    # Converte 'qual' in numerico (errori convertiti in 0)
    df['qual'] = pd.to_numeric(df['qual'], errors='coerce').fillna(0)
    # Sostituisce valori inferiori a 0.001 con 0.001 e applica la trasformazione logaritmica ln(1 + x)
    df['qual'] = np.log1p(df['qual'].clip(lower=0.001))
    
    # Converte 'depth' in numerico (errori convertiti in 0)
    df['depth'] = pd.to_numeric(df['depth'], errors='coerce').fillna(0).astype(int)
    
    # Codifica la colonna 'conferma' in valori numerici (vengono assegnati codici in base all'ordine)
    df['conferma'], uniques = pd.factorize(df['conferma'])
    print("Mapping di 'conferma':", dict(enumerate(uniques)))
    
    # Costruisce il dataset delle feature utilizzando le colonne 'qual' (logaritmizzata) e 'depth'
    X = df[['qual', 'depth']].values
    y = df['conferma'].values
    
    # Suddivide in training e test set per valutare il modello
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y)
    
    # Inizializza il modello LogisticRegression
    logistic_model = LogisticRegression(random_state=random_state, max_iter=1000)
    logistic_model.fit(X_train, y_train)
    
    # Valutazione sul test set
    y_pred = logistic_model.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    print(f"Accuratezza sul test set: {acc:.2f}")
    print("Classification Report:")
    print(classification_report(y_test, y_pred))
    
    # Calcola e visualizza la matrice di confusione
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(6, 4))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
    plt.title("Matrice di Confusione")
    plt.xlabel("Predetto")
    plt.ylabel("Reale")
    plt.tight_layout()
    plt.show()
    
    # Salva il modello tramite pickle
    with open(output_model_path, 'wb') as f:
        pickle.dump(logistic_model, f)
    print(f"Modello di regressione logistica addestrato e salvato in: {output_model_path}")
    
    return logistic_model

if __name__ == "__main__":
    # Specifica il percorso del file di training (modifica se necessario)
    training_file = "/home/alessandro/PROJECT/SKLEARN/SKLEARN.csv"
    output_model_path = "/home/alessandro/PROJECT/SKLEARN/logistic_model.pkl"
    
    # Avvia il training e la valutazione
    train_logistic_model(training_file, output_model_path)