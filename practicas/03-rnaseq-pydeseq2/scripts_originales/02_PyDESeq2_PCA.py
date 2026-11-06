#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 16:08:23 2024

@author: vjsanchez
"""
# %% [1] Import packages

# Rutas y busqueda de archivos
from pathlib import Path

# Manejo de datos
import pandas as pd
import anndata

# DESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

# PCA and graphs
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# Logging processes
from loguru import logger
# Basic loguru configuration
logger.add("debug.log", level="DEBUG", rotation="500 MB")

# %% [2]  Functions

def plot_pca(adata):
    """
    Realiza y muestra un gráfico PCA para los datos de expresión génica
    filtrados por un punto temporal específico en un objeto AnnData,
    utilizando la capa 'vst_counts' para los datos de expresión.
    
    Parameters:
    - adata (AnnData): Objeto AnnData que contiene los datos de expresión y metadatos.
    
    Descripción de procesos:
    - Verifica la existencia de la capa 'vst_counts' en el objeto AnnData.
    - Realiza un análisis PCA utilizando los datos de la capa 'vst_counts'.
    - Genera y muestra un gráfico de PCA utilizando Seaborn y Matplotlib.
    
    Ejemplo de uso:
    >>> plot_pca_by_time(adata)
    """
    # Verificar que la capa 'vst_counts' existe en el objeto AnnData
    if 'vst_counts' not in adata.layers:
        raise KeyError("La capa 'vst_counts' no se encuentra en el objeto AnnData. "
                       "Verifica que la capa exista o utiliza una capa válida.")
   
    counts_time = adata.layers['vst_counts']
    
    # Realizar el PCA
    pca = PCA(n_components=2)
    pca_results = pca.fit_transform(counts_time)
    
    # Convertir los resultados del PCA a un DataFrame
    pca_df = pd.DataFrame(pca_results, columns=['PC1', 'PC2'], index=adata.obs.index)
    pca_df['condition'] = adata.obs['condition'].values
    pca_df['sample'] = adata.obs.index  # Añadir etiquetas de muestra
    
    # Configurar el estilo de seaborn
    sns.set_style("whitegrid")
    
    # Crear el gráfico de PCA con Seaborn
    plt.figure(figsize=(12, 8))
    
    # Crear scatter plot con seaborn
    scatter = sns.scatterplot(
        data=pca_df,
        x='PC1',
        y='PC2',
        hue='condition',
        s=100,  # Tamaño de los puntos
        alpha=0.8,
        edgecolor='black',
        linewidth=0.5
    )
    
    # Añadir etiquetas de muestra a cada punto
    for i, row in pca_df.iterrows():
        plt.annotate(
            row['sample'], 
            (row['PC1'], row['PC2']),
            xytext=(5, 5),  # Offset de las etiquetas
            textcoords='offset points',
            fontsize=8,
            alpha=0.7
        )
    
    # Personalizar el gráfico
    plt.title('PCA RNA-seq data', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=12)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', fontsize=12)
    
    # Personalizar la leyenda
    plt.legend(title='Condition', title_fontsize=12, fontsize=10, loc='best')
    
    # Ajustar el layout
    plt.tight_layout()
    
    # Mostrar el gráfico
    plt.show()
    
    # Imprimir información adicional sobre la varianza explicada
    print(f"Varianza explicada por PC1: {pca.explained_variance_ratio_[0]:.2%}")
    print(f"Varianza explicada por PC2: {pca.explained_variance_ratio_[1]:.2%}")
    print(f"Varianza total explicada: {sum(pca.explained_variance_ratio_):.2%}")
    
    return pca_df, pca


# %% [3] Set working directory

# Set working directory
working_dir = Path.home()

# Access to data folder
data_dir = Path(
    working_dir / 'Documentos/Universidad_UFV/Asignaturas/Master_Bioinfo_UFV/Ejercicios/Tarea 3_Expresion_diferencial/')


# %% [4] Reading DataFrame

data_df = pd.read_csv('counts_final.tsv', sep='\t')
data_df.head()


# %%


#Lectura del objeto AnnData
adata = anndata.read(data_dir / "gene_counts.h5ad")
adata

# %% [5] Getting Variant Stabilizing Transformation

logger.info('Generating DeseqDataSet object')

inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    adata=adata,
    design_factors='condition',
    refit_cooks=True,
    inference=inference,
)

logger.info('Getting Variant Stabilizing Transformation')

# Variance stabilizing transformation for PCA and visualization heatmpas
dds.vst(use_design=False, fit_type='mean')
vst_counts = dds.layers["vst_counts"]
vst_counts_df = pd.DataFrame(vst_counts,
                             index=dds.obs_names,
                             columns=dds.var_names)
# Save to CSV
#vst_counts_df.to_csv(DATA / 'vst_counts.csv')
vst_counts_df.head()


# %% [6] Generate PCAs

logger.info('Generating PCA plots')
# Realizar el PCA y crear gráficos para cada tiempo
plot_pca(adata)
