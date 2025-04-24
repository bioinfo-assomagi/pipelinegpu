#!/usr/bin/env python3
"""
Script completo per:
 - Caricamento e preprocessing del dataset
 - Creazione e validazione incrociata di una pipeline con trasformazione polinomiale e regressione logistica
 - Visualizzazione della matrice di confusione e della mappa di contorni delle probabilità predette
 - Salvare la pipeline in un file pickle per uso futuro
"""

import os
import numpy as np
import pandas as pd
import pickle
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import train_test_split, cross_val_score, KFold
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.pipeline import Pipeline

def load_and_preprocess_data(csv_path):
    df = pd.read_csv(csv_path, header=0)
    df.columns = df.columns.str.lower().str.strip()
    
    data = df[['conferma', 'qual', 'depth']].copy()
    data['qual'] = pd.to_numeric(data['qual'], errors='coerce')
    data['depth'] = pd.to_numeric(data['depth'], errors='coerce')
    
    data.loc[data['qual'] < 0, 'qual'] = 0.001
    data['qual'] = np.log(data['qual'])
    
    data = data[data['qual'] >= 0]
    data['conferma'] = data['conferma'].map({'NO': 0, 'SI': 1}).astype('category')
    
    return data

def create_and_train_pipeline(data):
    X = data[['qual', 'depth']].values
    Y = data['conferma'].values
    
    # Creazione di una pipeline che include:
    # - Trasformazione polinomiale (grado 2) senza bias
    # - Regressione logistica
    pipeline = Pipeline([
        ('poly', PolynomialFeatures(degree=2, include_bias=False)),
        ('logistic', LogisticRegression(random_state=0, max_iter=1000))
    ])
    
    cv = KFold(n_splits=5, shuffle=True, random_state=42)                         # Configurazione della validazione incrociata a 5 fold
    scores = cross_val_score(pipeline, X, Y, cv=cv, scoring='accuracy')
    print("Punteggi di accuratezza per ogni fold:", scores)
    print("Accuratezza media:", scores.mean())
    
    pipeline.fit(X, Y)                                                            # Addestra il modello sui dati interi
    return pipeline

def plot_confusion_matrix(model, data):
    X = data[['qual', 'depth']].values
    Y = data['conferma'].values
    Y_pred = model.predict(X)
    
    cm = confusion_matrix(Y, Y_pred)
    
    plt.figure(figsize=(6,4))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title("Matrice di Confusione (Dataset Intero)")
    plt.xlabel("Predetto")
    plt.ylabel("Reale")
    plt.tight_layout()
    plt.show()

def plot_probability_contours(model, data, save_path="/home/alessandro/PROJECT/SKLEARN/probability_contours.png"):
    qual_seq = np.linspace(data['qual'].min(), data['qual'].max(), 100)
    depth_seq = np.linspace(data['depth'].min(), data['depth'].max(), 100)
    
    qual_grid, depth_grid = np.meshgrid(qual_seq, depth_seq)
    grid_points = np.c_[qual_grid.ravel(), depth_grid.ravel()]
    
    
    pred_probs = model.predict_proba(grid_points)[:, 1]                            # Calcola le probabilità della classe 1 per ogni punto della griglia
    Z = pred_probs.reshape(qual_grid.shape)
    
    plt.figure(figsize=(12,10))
    contour_levels = np.arange(0, 1.1, 0.1)
    contour = plt.contour(qual_seq, depth_seq, Z, levels=contour_levels,
                          colors='grey', linestyles='dashed')
    plt.clabel(contour, fmt='%1.1f', inline=True, fontsize=8)
    
    colors = data['conferma'].map({0: 'black', 1: 'cyan'})
    plt.scatter(data['qual'], data['depth'], c=colors, edgecolor='k', s=50,
                label='Punti del dataset')
    
    plt.xlabel("qual")
    plt.ylabel("depth")
    plt.title("Mappa a Contorni delle Probabilità Predette")
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    if save_path is not None:
        plt.savefig(save_path)
        print("Immagine salvata in:", save_path)
    
    plt.show()

def threshold_predictions(model, data, p_threshold=0.9):
    
    X_full = model.named_steps['poly'].transform(data[['qual', 'depth']].values)    # Trasforma le feature del dataset usando lo step 'poly' della pipeline
    pred_probs_full = model.named_steps['logistic'].predict_proba(X_full)[:, 1]     # Calcola le probabilità usando la regressione logistica
    
    Y_hat = np.where(pred_probs_full < p_threshold, 0, 1)                           # Converti le probabilità in classi utilizzando la soglia definita
    
    cm_threshold = confusion_matrix(data['conferma'].values, Y_hat)
    print("Matrice di Confusione (threshold = {:.1f}):".format(p_threshold))
    print(cm_threshold)

def main():
    csv_path = "/home/alessandro/PROJECT/SKLEARN/SKLEARN.csv"
    data = load_and_preprocess_data(csv_path)                                       # Carica e preprocessa i dati
    pipeline_model = create_and_train_pipeline(data)                                # Crea ed addestra la pipeline
    plot_confusion_matrix(pipeline_model, data)                                     # Visualizza la matrice di confusione sui dati interi
    plot_probability_contours(pipeline_model, data)                                 # Visualizza la mappa a contorni delle probabilità predette
    threshold_predictions(pipeline_model, data, p_threshold=0.9)                    # Calcola le predizioni utilizzando una soglia custom e visualizza la matrice di confusione
    model_output_path = "/home/alessandro/PROJECT/SKLEARN/logistic_model.pkl"
    joblib.dump(pipeline_model, model_output_path)

    print("Modello salvato in:", model_output_path)


if __name__ == "__main__":
    main()
