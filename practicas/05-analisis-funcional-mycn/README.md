# Práctica 5 — Análisis funcional: knockout de MYCN en neuroblastoma

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/05-analisis-funcional-mycn/Analisis_Funcional_MYCN_Colab.ipynb)

Análisis funcional completo de un RNA-seq de células de neuroblastoma con **knockout de MYCN** (CRISPR, 5 sgMYCN) frente a controles (5 sgAAVS1).

## Contenido

- `Analisis_Funcional_MYCN_Colab.ipynb` — notebook de Google Colab.
- `datos/degs.csv` — resultados del análisis de expresión diferencial (KO vs Control), 16.235 genes.
- `datos/normalized_counts.csv` — matriz de cuentas normalizadas (genes × 10 muestras).

## Qué se hace en el notebook

1. **Carga y limpieza**: símbolos a mayúsculas, eliminación de filas sin `gene_id`, agregación de símbolos duplicados.
2. **Comprobación de coherencia** entre `degs.csv` y `normalized_counts.csv` + PCA exploratorio.
3. **ORA** (Enrichr) de genes UP y DOWN (FDR < 0.05, |log2FC| > 1) con GO BP 2023, KEGG 2021 y MSigDB Hallmark 2020.
4. **GSEA pre-ranked** con métrica `Rank = -log10(padj) · sign(log2FC) + 0.01 · log2FC` (evita empates).
5. **GSEA clásico** desde la matriz de expresión (`permutation_type="phenotype"` con fallback a `"gene_set"`).
6. **Heatmap de leading-edge genes** de la ruta MYC/E2F más significativa.
7. **ssGSEA**: matriz rutas × muestras, heatmap y t-test ruta a ruta.
8. Exportación de todos los resultados en un ZIP.

## Notas sobre los datos

- El gen **MYCN no aparece** en las matrices (filtrado en pasos previos del pipeline); el efecto del knockout se evalúa a través de sus genes diana (Hallmark MYC targets, E2F targets).
- `normalized_counts.csv` contiene muchas filas sin `gene_id` y símbolos duplicados (RNAs repetitivos como `Y_RNA`); el notebook los gestiona automáticamente.
