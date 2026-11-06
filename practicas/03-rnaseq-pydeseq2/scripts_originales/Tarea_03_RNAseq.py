#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 10:53:43 2025

@author: vjsanchez
"""
# %% [0] Import packages

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# %% [1] Set up paths and load dataframe

wd = Path.cwd()
counts = Path("featurecounts/gene_counts.txt")


# Leer el archivo y especificar nombres de columnas manualmente
counts_df = pd.read_csv('featurecounts/gene_counts.txt', 
                 sep='\t', 
                 skiprows=2,  # Salta metadatos
                 names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length',
                        'SRR3480371', 'SRR3480372', 'SRR3480373', 'SRR3480374'])

# %% [2] Calculo de estadisticas por gen y muestra


# Seleccionar solo las columnas de conteos (excluyendo metadatos)
count_columns = ['SRR3480371', 'SRR3480372', 'SRR3480373', 'SRR3480374']
counts_only = counts_df[count_columns]

# Estadísticas por columna (muestra)
print("Estadísticas por muestra:")
print(counts_only.describe())

# Estadísticas por gen (promedio across muestras)
counts_df['mean_counts'] = counts_only.mean(axis=1)
counts_df['total_counts'] = counts_only.sum(axis=1)

print("\nEstadísticas de conteos por gen:")
print(counts_df[['mean_counts', 'total_counts']].describe())

# %% [3] Representaciones graficas

# Configurar el estilo
plt.style.use('default')
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# 1. Histograma de conteos totales por gen
axes[0,0].hist(counts_df['total_counts'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
axes[0,0].set_xlabel('Conteos totales por gen')
axes[0,0].set_ylabel('Frecuencia')
axes[0,0].set_title('Distribución de conteos totales por gen')
axes[0,0].set_yscale('log')  # Escala log para mejor visualización

# 2. Boxplot de conteos por muestra
counts_only.boxplot(ax=axes[0,1])
axes[0,1].set_title('Distribución de conteos por muestra')
axes[0,1].set_ylabel('Conteos')
axes[0,1].tick_params(axis='x', rotation=45)

# 3. Distribución en escala log
log_counts = np.log10(counts_df['total_counts'] + 1)  # +1 para evitar log(0)
axes[1,0].hist(log_counts, bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
axes[1,0].set_xlabel('Log10(conteos totales + 1)')
axes[1,0].set_ylabel('Frecuencia')
axes[1,0].set_title('Distribución en escala logarítmica')

# 4. Violin plot
counts_melted = counts_only.melt(var_name='Sample', value_name='Counts')
sns.violinplot(data=counts_melted, x='Sample', y='Counts', ax=axes[1,1])
axes[1,1].set_title('Distribución de conteos por muestra (Violin plot)')
axes[1,1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.show()