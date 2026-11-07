# Práctica 8 — Deconvolución inmune: citometría digital con TumorDecon

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/08-deconvolution-tumordecon/Deconvolucion_TumorDecon_Colab.ipynb)

Estimación de la **composición celular del microambiente tumoral** a partir de RNA-seq bulk con **TumorDecon** (Le et al., 2020), sobre 592 tumores de cáncer colorrectal de **TCGA COADREAD** y la matriz de signaturas **LM22** (22 tipos de células inmunes).

## Contenido

- `Deconvolucion_TumorDecon_Colab.ipynb` — notebook de Google Colab.
- `datos/coadread_tcga_reducida.txt.gz` — matriz de expresión TCGA COADREAD (formato cBioPortal) **reducida** a 5.411 genes (los de LM22 + los 5.000 más variables) × 592 pacientes. La reducción no afecta a CIBERSORT/DeconRNASeq (resultados idénticos a la matriz completa, verificado) y hace viable el alojamiento en GitHub.
- `datos/LM22.txt` — matriz de signaturas de 22 tipos de células inmunes (547 genes).
- `datos/kmeans_signature_matrix_qval.txt` — matriz de signaturas alternativa (16 subtipos, IDs Ensembl) para la sección opcional.

## Qué se hace en el notebook

1. Instalación de TumorDecon + parche de compatibilidad (paquete de 2020: `collections.Mapping` y `DataFrame.append`).
2. **CIBERSORT** (ν-SVR) sobre los 592 pacientes.
3. **DeconRNASeq** (mínimos cuadrados no negativos) + comparación de concordancia entre métodos por tipo celular.
4. **ssGSEA** y **SingScore** (métodos de rangos) + su comparación.
5. Visualización: boxplot de la cohorte, barchart por paciente, pair plot, clustermap y agrupación de tipos celulares (`combine_celltypes`).
6. Sección opcional con una matriz de signaturas alternativa (conversión Ensembl→Hugo automática).
7. Exportación de resultados en ZIP.

## Notas

- Tiempos aproximados en Colab: CIBERSORT ~2–3 min, ssGSEA ~1–2 min, el resto es inmediato.
- Los métodos de rangos (ssGSEA/SingScore) usan la matriz reducida como fondo de ranking; con el transcriptoma completo los scores variarían ligeramente (se explica en el notebook).
