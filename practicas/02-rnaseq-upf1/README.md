# Práctica 2 — RNA-seq: Control de Calidad, Trimming y Mapeo o Pseudoalineamiento

## Análisis del Rol de UPF1 en la Diferenciación de Células Madre Embrionarias Humanas

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/02-rnaseq-upf1/RNAseq_UPF1_Colab.ipynb)

## Contexto biológico

El *Nonsense-Mediated RNA Decay* (NMD) es un mecanismo de control de calidad celular que
degrada ARNm con codones de terminación prematuros. **UPF1** es la proteína central de este
proceso, actuando como helicasa y regulando la estabilidad del ARNm. En células madre
embrionarias humanas (hESCs), UPF1 no solo elimina transcritos defectuosos, sino que también
regula la expresión de genes importantes para el mantenimiento de la pluripotencia y la
diferenciación celular.

## Objetivo del ejercicio

Analizar el impacto transcriptómico de la pérdida de función de UPF1 en células madre
embrionarias humanas mediante técnicas de RNA-seq, identificando genes diferencialmente
expresados y evaluando su relevancia biológica.

## Diseño experimental

| Característica | Detalle |
|---|---|
| Condiciones | Control: hESCs + siRNA *scrambled* / Tratamiento: hESCs + siRNA contra UPF1 |
| Plataforma | Illumina HiSeq 2500 |
| Tipo de lecturas | 50 bp single-end, *stranded* |
| Replicados | 2 réplicas biológicas por condición |
| Dataset | [GSE81214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81214) (GEO/NCBI) |

| ID muestra | Condición | Descripción | Nº lecturas (originales) |
|---|---|---|---|
| SRR3480371 | Control 1 | siRNA control | 42 059 631 |
| SRR3480372 | Control 2 | siRNA control | 54 401 313 |
| SRR3480373 | Tratamiento 1 | siRNA UPF1 | 52 991 800 |
| SRR3480374 | Tratamiento 2 | siRNA UPF1 | 64 635 061 |

> **Para la práctica en clase** se utilizan **submuestras** de estas lecturas
> (`datos/*.subset.fastq.gz`, ~7 000-12 000 lecturas por muestra) junto con una anotación
> reducida (`datos/subset_gencode.v38.annotation.gtf`, GENCODE v38) que cubre 6 genes del
> cromosoma 19: **UPF1**, EML2, ITPKC, JUND, LPAR2 y SIRT2. Esto permite ejecutar el
> pipeline completo en minutos en Google Colab. El análisis de los datos completos
> seguiría exactamente los mismos pasos.

## Estructura de la práctica

- **PARTE I** — Preparación y control de calidad: formato FASTQ, FastQC.
- **PARTE II** — Alineamiento con STAR (genoma de referencia: cromosoma 19 humano).
- **PARTE III** — Cuantificación génica con featureCounts.
- **PARTE IV** — Cuantificación alternativa con Salmon (quasi-mapping) y agregación a nivel
  de gen con [pymportx](https://github.com/victorsanchezarevalo/Pymportx).
- **PARTE V** — Comparación de métodos de cuantificación (STAR+featureCounts vs Salmon+pymportx).

Todas las preguntas de las actividades están integradas en el notebook, en celdas de
markdown identificadas con ❓.

## Contenido

| Archivo | Descripción |
|---|---|
| `RNAseq_UPF1_Colab.ipynb` | Notebook completo para Google Colab (recomendado) |
| `datos/*.subset.fastq.gz` | Submuestras de las 4 muestras SRA (50 bp single-end) |
| `datos/subset_gencode.v38.annotation.gtf` | Anotación GENCODE v38 reducida (6 genes del chr19, incluido UPF1) |
| `enunciado.md` | Enunciado original de la tarea con las instrucciones de entrega |

## Instrucciones de entrega (estudiantes)

Debes entregar **dos archivos**:

1. **Documento PDF:** respuestas teóricas y análisis de resultados — claro, bien
   estructurado, con gráficos y tablas cuando sean relevantes y conclusiones fundamentadas.
2. **Scripts:** el notebook ejecutado de principio a fin, con código funcional, comentarios
   exhaustivos (propósito, lógica, parámetros, interpretación) y salidas visibles.

Se evaluará especialmente: claridad de las respuestas, presentación del código,
funcionalidad (debe ejecutarse sin errores) y completitud (todas las preguntas respondidas).
