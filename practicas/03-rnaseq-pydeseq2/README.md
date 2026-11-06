# Práctica 3 — Expresión Diferencial con PyDESeq2 y PCA

## Análisis del knockdown de UPF1 en células madre embrionarias humanas (GSE81214)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/03-rnaseq-pydeseq2/RNAseq_PyDESeq2_Colab.ipynb)

## Objetivo

Continuación de la [Práctica 2](../02-rnaseq-upf1/): partiendo de la matriz de conteos de
**featureCounts**, realizar el análisis de expresión diferencial entre hESCs control y
hESCs con *knockdown* de UPF1 usando **PyDESeq2** (implementación Python de DESeq2),
evaluar la estructura de los datos con **VST + PCA**, y **validar biológicamente** el
experimento comprobando la reducción de expresión de UPF1.

## Estructura

- **PARTE I** — Exploración de la matriz de conteos: estadísticas y visualización.
- **PARTE II** — PyDESeq2: modelo binomial negativa, test de Wald, volcano plot y
  validación del knockdown de UPF1 con conteos normalizados.
- **PARTE III** — Variance Stabilizing Transformation (VST) y PCA de las muestras.

> **Nota didáctica:** se trabaja con la anotación reducida de 6 genes del chr19
> (UPF1, EML2, ITPKC, JUND, LPAR2, SIRT2) y 2 réplicas por condición. Las preguntas del
> notebook discuten explícitamente las limitaciones estadísticas de este diseño.

## Contenido

| Archivo | Descripción |
|---|---|
| `RNAseq_PyDESeq2_Colab.ipynb` | Notebook completo para Google Colab |
| `datos/gene_counts.txt` | Matriz de conteos de featureCounts (salida de la Práctica 2) |
| `scripts_originales/` | Scripts de referencia en Python (versión Spyder del análisis) |

## Prerequisitos

Haber completado la Práctica 2 (o usar la matriz de conteos precomputada incluida en
`datos/`). El notebook también permite subir tu propio `gene_counts.txt` (`MODO = "SUBIR"`).

## Resultados esperados

Con los datos de la práctica deberías obtener:

- **UPF1 (ENSG00000005007): log2FC ≈ −3.0, padj ≈ 1e-118** — reducción de ~87% en conteos
  normalizados. El knockdown se valida claramente.
- **LPAR2** (padj ≈ 0.004) e **ITPKC** (padj ≈ 0.017) también significativos, aunque con
  cambios pequeños (|log2FC| < 0.5).
- En el **PCA** sobre datos VST, PC1 explica ~99% de la varianza y separa perfectamente
  control vs knockdown.

Pequeñas variaciones son normales según la versión de pydeseq2.

## Referencias

- Love, Huber & Anders (2014). *Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2*. Genome Biology 15:550.
- Muzellec et al. (2023). *PyDESeq2: a python package for bulk RNA-seq differential
  expression analysis*. Bioinformatics 39(9):btad547.
