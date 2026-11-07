# Práctica 9 — Análisis de scRNA-seq: PBMC 3k

**Máster en Bioinformática — Módulo de Transcriptómica**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/09-scrna-seq-pbmc/ScrRNAseq_PBMC3k_Colab.ipynb)

Primera de las dos sesiones de transcriptómica de célula única. Se analiza el dataset de referencia **PBMC 3k de 10x Genomics** (~2.700 PBMC de un donante sano, química v1) con **scanpy**, siguiendo el flujo clásico de análisis de scRNA-seq.

## Contenido

1. Descarga y carga de la salida de Cell Ranger (formato Matrix Market disperso).
2. Anatomía del objeto `AnnData` (X, obs, var, layers, obsm).
3. Control de calidad por célula (genes detectados, conteos, % mitocondrial) y por gen.
4. Normalización por tamaño de librería (`normalize_total`) y `log1p`.
5. Selección de genes altamente variables (conexión con la práctica 6).
6. PCA → grafo de vecinos → UMAP → clustering de **Leiden**.
7. Marcadores diferenciales por clúster (Wilcoxon) y **anotación de tipos celulares** con marcadores canónicos (CD3/IL7R, MS4A1, CD14, FCGR3A, NKG7, FCER1A, PPBP...).
8. Ejercicios: efecto de la resolución de Leiden; interpretación de los QC.

## Datos

- Descarga automática desde el CDN de 10x Genomics (~7,6 MB), con espejo de respaldo en `datos/pbmc3k_filtered_gene_bc_matrices.tar.gz` (el notebook prueba ambos y admite subida manual).

## Requisitos

- Google Colab (CPU estándar, ~2 GB de RAM). Instalación incluida en el notebook: `scanpy`, `leidenalg`.
- Duración estimada de ejecución completa: ~2–3 minutos tras la instalación.

## Referencias

- Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 19, 15 (2018).
- Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 9, 5233 (2019).
- Zheng, G. X. Y. *et al.* Massively parallel digital transcriptional profiling of single cells. *Nature Communications* 8, 14049 (2017).
