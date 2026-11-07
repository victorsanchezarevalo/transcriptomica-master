# Práctica 6 — WGCNA: redes de coexpresión en el modelo 5xFAD de Alzheimer

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/06-wgcna-5xfad/WGCNA_5xFAD_Colab.ipynb)

Análisis de redes de coexpresión génica con **PyWGCNA** sobre un RNA-seq del modelo murino de Alzheimer **5xFAD**: 192 muestras (corteza e hipocampo; genotipos 5xFADWT y 5xFADHEMI; edades 4, 8, 12 y 18 meses; ambos sexos).

## Contenido

- `WGCNA_5xFAD_Colab.ipynb` — notebook de Google Colab.
- `datos/expressionList.csv.gz` — matriz de expresión (192 muestras × ~22.000 genes, TPM). Descomprimir antes de usar en local.
- `datos/sampleInfo.csv` — metadatos de las muestras (`Age`, `Tissue`, `Sex`, `Genotype`).
- `datos/geneList` — equivalencia Ensembl → símbolo de gen (evita la consulta a BioMart).

## Qué se hace en el notebook

1. Creación del objeto `WGCNA` y exploración del AnnData.
2. `preprocess()`: filtrado de genes poco expresados y muestras con muchos faltantes.
3. `findModules()`: elección del soft power, TOM, clustering y fusión de módulos.
4. Metadatos: `updateSampleInfo()`, `setMetadataColor()`, `updateGeneInfo()`.
5. `analyseWGCNA()` y heatmaps módulo–rasgo personalizados (edad como variable continua, combinación edad × sexo).
6. Genes hub con `top_n_hub_genes()`.
7. Redes de coexpresión con `CoexpressionModulePlot()`.
8. Enriquecimiento funcional (GO, KEGG ratón, Reactome) vía Enrichr.
9. Interacciones proteína–proteína contra STRING (`request_PPI`, `PPI_network`).
10. Guardado del objeto (`saveWGCNA`) y exportación de resultados en ZIP.

## Notas

- Los pasos de enriquecimiento y PPI necesitan internet (Enrichr y STRING).
- `findModules()` es la parte pesada: varios minutos en Colab.
- Los datos proceden del tutorial oficial de PyWGCNA (mortazavilab/PyWGCNA), artículo del modelo 5xFAD.
