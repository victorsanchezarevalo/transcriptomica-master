#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 11:21:14 2025

@author: slimbook
"""
# %%


from pathlib import Path

import pandas as pd
import numpy as np


from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

import seaborn as sns
# %%

# Ruta hacia los datos
data_path = Path('/home/slimbook/Documentos/Universidad_UFV/Asignaturas/Master_Bioinfo_UFV/Ejercicios/Tarea_03_RNAseq/featurecounts/gene_counts.txt')

# Verificar que el archivo existe
if not data_path.exists():
    raise FileNotFoundError(f"No se encontró el archivo: {data_path}")

# Columnas a cargar
columns = [
    'Geneid',
    'star_alignment/SRR3480371_Aligned.sortedByCoord.out.bam',
    'star_alignment/SRR3480372_Aligned.sortedByCoord.out.bam',
    'star_alignment/SRR3480373_Aligned.sortedByCoord.out.bam',
    'star_alignment/SRR3480374_Aligned.sortedByCoord.out.bam'
]

# Mapeo para renombrar columnas (más conciso)
column_mapping = {col: col.split('/')[-1].split('_')[0] if 'SRR' in col else col 
                  for col in columns}

# Cargar y procesar el DataFrame
counts_df = pd.read_csv(
    data_path, 
    sep='\t', 
    skiprows=1, 
    usecols=columns
).rename(columns=column_mapping)

counts_df.set_index('Geneid', inplace=True)

counts_df = counts_df.astype(int)
# Establecer Geneid como índice (opcional, útil para análisis posteriores)

# Verificar la carga
print(f"Dimensiones del DataFrame: {counts_df.shape}")
print("\nPrimeras filas:")
print(counts_df.head())
print(f"\nColumnas: {list(counts_df.columns)}")

# Verificar datos faltantes
print(f"\nDatos faltantes por columna:\n{counts_df.isnull().sum()}")

# %%

sample_names = ['control1', 'control2', 'treatment1', 'treatment2']
sample_data = pd.DataFrame(index=sample_names, columns=['condition'])

# Y luego asignar las condiciones
sample_data['condition'] = ['control', 'control', 'treatment', 'treatment']


# %%

inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=sample_data,
    design="~condition",  # compare samples based on the "condition"
    # column ("B" vs "A")
    refit_cooks=True,
    inference=inference,
)