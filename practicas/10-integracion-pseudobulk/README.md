# Práctica 10 — Integración de muestras y análisis pseudobulk: PBMC + IFN-β

**Máster en Bioinformática — Módulo de Transcriptómica**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/10-integracion-pseudobulk/Integracion_Pseudobulk_Colab.ipynb)

Segunda sesión de célula única. Se analiza el dataset de **Kang et al. 2018** (GSE96583): PBMC de **8 donantes** multiplexados genéticamente, cultivados 6 h sin (`ctrl`) o con interferón-β (`stim`). La asignación célula → donante (demuxlet) y la anotación de tipos celulares se toman de los metadatos de los autores.

## Contenido

1. Construcción de un `AnnData` multimuestra desde los ficheros crudos de GEO (dos corridas 10x).
2. Control de calidad y visualización del efecto de condición y donante.
3. Integración con **Harmony** (batch = donante) y comparación antes/después; discusión de qué es batch y qué es biología.
4. **Ejercicio de anotación**: reproducir los tipos celulares de los autores con los marcadores canónicos de la práctica 9 (dotplot sobre los clústeres de Leiden) y cuantificar el acuerdo.
5. **Análisis composicional**: proporciones de tipos celulares por donante y condición (la respuesta al IFN-β es transcripcional, no composicional); mención a scCODA para el análisis riguroso.
6. **Pseudobulk**: agregación de conteos por donante × condición × tipo celular (la réplica real es el donante, no la célula).
7. Expresión diferencial con **PyDESeq2** en monocitos CD14+ (stim vs ctrl) — recupera las prácticas 3–4. Volcán con los ISG clásicos (IFIT, IFITM, ISG15/20, MX, OAS, DDX58...).
8. Comparación con el test "ingenuo" célula a célula: la falacia de las pseudoréplicas (Squair et al. 2021).
9. (Opcional) Análisis funcional con **decoupler** + Hallmark sobre el estadístico de DESeq2 — conexión con la práctica 7. Top vías: respuesta a interferón α/γ.

## Datos

- Descarga automática desde GEO (NCBI): dos matrices `.mtx.gz` (~29 MB cada una), barcodes, genes y la tabla de metadatos de los autores. No requiere ficheros del repositorio.

## Requisitos

- Google Colab (CPU estándar, ~4 GB de RAM). Instalación incluida en el notebook: `scanpy`, `leidenalg`, `harmonypy==0.0.10` (la versión compatible con el wrapper de scanpy), `pydeseq2`, `decoupler`.
- Duración estimada de ejecución completa: ~3–5 minutos tras la instalación y la descarga.

## Referencias

- Kang, H. M. *et al.* Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nature Biotechnology* 36, 89–94 (2018).
- Korsunsky, I. *et al.* Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods* 16, 1289–1296 (2019).
- Squair, J. W. *et al.* Confronting false discoveries in single-cell differential expression. *Nature Communications* 12, 5692 (2021).
- Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15, 550 (2014).
