# Transcriptómica — Máster en Bioinformática
## Universidad Francisco de Vitoria

Material docente del **módulo de Transcriptómica**: prácticas, datos y scripts de análisis.

## Prácticas

| # | Práctica | Descripción | Abrir |
|---|----------|-------------|-------|
| 1 | [Microarrays LLA](practicas/01-microarrays-lla/) | Expresión diferencial **BCR/ABL vs NEG** en Leucemia Linfoide Aguda con limma (InMoose) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/01-microarrays-lla/Analisis_microarrays_LLA_Colab.ipynb) |
| 2 | [RNA-seq UPF1](practicas/02-rnaseq-upf1/) | QC, alineamiento con **STAR**, cuantificación con **featureCounts** y pseudoalineamiento con **Salmon** + pymportx (GSE81214, UPF1 en hESCs) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/02-rnaseq-upf1/RNAseq_UPF1_Colab.ipynb) |
| 3 | [PyDESeq2 y PCA](practicas/03-rnaseq-pydeseq2/) | Expresión diferencial con **PyDESeq2** y **PCA** sobre VST: validación del knockdown de UPF1 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/03-rnaseq-pydeseq2/RNAseq_PyDESeq2_Colab.ipynb) |
| 4 | [TNF en T3M4](practicas/04-tnf-t3m4-pydeseq2/) | Expresión diferencial **TNF vs Untreated** en células T3M4 de cáncer de páncreas con **PyDESeq2**, PCA y visualizaciones | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/04-tnf-t3m4-pydeseq2/04-tnf-t3m4-pydeseq2/TNF_T3M4_PyDESeq2_Colab.ipynb) |
| 5 | [Análisis funcional MYCN](practicas/05-analisis-funcional-mycn/) | **ORA, GSEA (pre-ranked y clásico) y ssGSEA** con gseapy sobre un knockout de **MYCN** en neuroblastoma | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/05-analisis-funcional-mycn/Analisis_Funcional_MYCN_Colab.ipynb) |
| 6 | [WGCNA 5xFAD](practicas/06-wgcna-5xfad/) | **Redes de coexpresión con PyWGCNA** en el modelo murino de Alzheimer 5xFAD: módulos, relación módulo–rasgo, genes hub, redes y enriquecimiento | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/06-wgcna-5xfad/WGCNA_5xFAD_Colab.ipynb) |
| 7 | [decoupler bulk](practicas/07-decoupler-bulk/) | **Actividad de TFs y rutas con decoupler** (CollecTRI + ULM, PROGENy + MLM, Hallmark + ORA) en HSCs tratadas con **TGFβ** (GSE151251) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/07-decoupler-bulk/Decoupler_Bulk_Colab.ipynb) |
| 8 | [Deconvolución inmune](practicas/08-deconvolution-tumordecon/) | **Citometría digital con TumorDecon** (CIBERSORT, DeconRNASeq, ssGSEA, SingScore) con la matriz **LM22** sobre 592 tumores colorrectales de **TCGA COADREAD** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/08-deconvolution-tumordecon/Deconvolucion_TumorDecon_Colab.ipynb) |
| 9 | [scRNA-seq PBMC 3k](practicas/09-scrna-seq-pbmc/) | **Análisis de célula única con scanpy** (QC, normalización, HVG, PCA, UMAP, Leiden, marcadores y anotación) sobre 2.700 PBMC de 10x Genomics | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/09-scrna-seq-pbmc/ScrRNAseq_PBMC3k_Colab.ipynb) |
| 10 | [Integración y pseudobulk](practicas/10-integracion-pseudobulk/) | **Integración con Harmony** y **DE pseudobulk con PyDESeq2** en PBMC de 8 donantes ± **IFN-β** (Kang 2018, GSE96583); la falacia de las pseudoréplicas | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/10-integracion-pseudobulk/Integracion_Pseudobulk_Colab.ipynb) |

## Cómo usar este material (estudiantes)

1. Pulsa el botón **Open in Colab** de la práctica.
2. En Colab: **Entorno de ejecución → Ejecutar todas**.
3. La primera celda de cada notebook instala los paquetes necesarios y **puede tardar
   varios minutos** (en la práctica 1, `inmoose` se compila desde el código fuente y
   puede tardar 5–10 minutos). Es normal: espera a que aparezca el mensaje de
   confirmación antes de continuar.
4. Los datos se descargan automáticamente (desde este repositorio o desde GEO/TCGA/10x
   según la práctica); no necesitas subir ningún archivo.
5. Al final del notebook se genera un ZIP con todos los resultados (CSV + figuras) que
   se descarga a tu ordenador. **Guárdalo**: Colab borra los archivos al cerrar la sesión.

## Ejecución en local (opcional)

Si prefieres trabajar en tu equipo (Spyder, Jupyter…), cada carpeta de práctica incluye
un `README.md` con las instrucciones y los datos necesarios. La práctica 1 además
incluye un entorno conda listo para usar:

```bash
conda env create -f practicas/01-microarrays-lla/inmoose_env.yaml
conda activate inmoose_env
```

## Estructura del repositorio

```
transcriptomica-master/
├── README.md                  ← este archivo
└── practicas/
    ├── 01-microarrays-lla/
    ├── 02-rnaseq-upf1/
    ├── 03-rnaseq-pydeseq2/
    ├── 04-tnf-t3m4-pydeseq2/
    ├── 05-analisis-funcional-mycn/
    ├── 06-wgcna-5xfad/
    ├── 07-decoupler-bulk/
    ├── 08-deconvolution-tumordecon/
    ├── 09-scrna-seq-pbmc/
    └── 10-integracion-pseudobulk/
```

Cada carpeta de práctica contiene su `README.md` con la guía, el notebook de Colab
y, cuando aplica, los datos y scripts auxiliares. Ejemplo (práctica 1):

```
01-microarrays-lla/
├── README.md          ← guía de la práctica
├── Analisis_microarrays_LLA_Colab.ipynb
├── Limma_inmoose.py   ← versión script (Spyder)
├── inmoose_env.yaml   ← entorno conda
└── datos/
    ├── exprs.csv.gz   ← matriz de expresión comprimida (12 625 sondas × 128 muestras)
    └── pData.csv      ← metadatos clínicos y subtipos moleculares
```

## Licencia

Material docente con fines educativos. Se recomienda distribuir bajo licencia
[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.es).

Los datos proceden de fuentes públicas: dataset `ALL` de Bioconductor (Chiaretti et
al., 2004, *Blood*; plataforma Affymetrix HG-U95Av2), series de **GEO** (GSE81214,
GSE151251, GSE96583 — Kang et al., 2018, *Nature Biotechnology*), **TCGA**
(COADREAD, vía GDC), la matriz de referencia **LM22** de CIBERSORT (Newman et al.,
2015, *Nature Methods*) y el dataset **PBMC 3k** de 10x Genomics.
