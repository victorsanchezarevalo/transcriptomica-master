# Práctica 4 — Expresión diferencial: TNF en células T3M4

**Máster en Bioinformática — Módulo de Transcriptómica**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/04-tnf-t3m4-pydeseq2/04-tnf-t3m4-pydeseq2/TNF_T3M4_PyDESeq2_Colab.ipynb)

## Objetivo

Analizar los cambios transcriptómicos inducidos por **TNF** en células **T3M4** de cáncer de páncreas, identificando genes diferencialmente expresados y explorando las vías biológicas afectadas.

## Diseño experimental

- **Control (`Untreated`)**: `shCtrl_Untreated_rep1`, `shCtrl_Untreated_rep2`, `shCtrl_Untreated_rep3`.
- **Tratamiento (`TNF`)**: `shCtrl_TNF_rep1`, `shCtrl_TNF_rep2`, `shCtrl_TNF_rep3`.
- **Réplicas biológicas**: 3 por condición.

## Datos

- `datos/counts_final_example.csv`: matriz featureCounts con **60.155 genes × 6 muestras**, ya restringida al diseño del ejercicio.
- La matriz `counts_final.tsv` original contiene 24 muestras, incluyendo `shBPTF` y una cuarta réplica; para esta práctica se seleccionan únicamente las 6 muestras `shCtrl` indicadas.

## Contenido del notebook

1. Carga y validación de la matriz de conteos.
2. Creación de un objeto `AnnData` con metadatos de muestras y genes.
3. Normalización con PyDESeq2/DESeq2.
4. Transformación VST.
5. PCA y heatmap de correlación entre muestras.
6. Contraste **TNF vs Untreated**.
7. Filtrado por **FDR < 0,05** y **|log2FC| > 1**.
8. MA plot, volcano plot y heatmap de los 50 genes con menor FDR.
9. Comprobación de genes canónicos de respuesta a TNF.
10. Sección opcional de enriquecimiento funcional con Enrichr.
11. Exportación de resultados en ZIP.

## Nota de validación

Con los datos aportados, el subconjunto de 6 muestras coincide exactamente con las columnas correspondientes de la matriz completa. Tras filtrar genes con suma de conteos ≥ 10 quedan **22.515 genes**. En los resultados de referencia aportados se obtienen **367 genes diferencialmente expresados** (266 sobre-expresados y 101 sub-expresados con FDR < 0,05 y |log2FC| > 1). Los genes canónicos de respuesta a TNF —por ejemplo `TNFAIP3`, `IL1B`, `BIRC3`, `NFKBIA`, `CXCL8`, `CCL2`, `ICAM1` y `RELB`— aparecen sobre-expresados, lo que respalda la coherencia biológica del experimento.

## Uso en Colab

1. Abrir el cuaderno con el botón **Open in Colab**.
2. Ejecutar **Entorno de ejecución → Ejecutar todas**.
3. Si el archivo de datos todavía no está subido al repositorio, cambiar `MODO = "URL"` a `MODO = "SUBIR"` en la celda 1 y subir `counts_final_example.csv`.
4. Al final se descargará un ZIP con tablas y figuras.
