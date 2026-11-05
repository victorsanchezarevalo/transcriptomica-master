# Práctica 1 — Análisis de expresión diferencial con microarrays: BCR/ABL vs NEG en LLA

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TU_USUARIO/transcriptomica-master/blob/main/practicas/01-microarrays-lla/Analisis_microarrays_LLA_Colab.ipynb)

## Objetivo principal

Identificar el perfil de expresión génica diferencial entre los subtipos moleculares
**BCR/ABL** y **NEG** en Leucemia Linfoide Aguda mediante análisis de microarrays
usando **limma** en **InMoose**.

## Objetivos específicos

- **Análisis diferencial:** detectar genes significativos (FDR < 0.05, |log2FC| > 0.5),
  clasificarlos en sobre- e infra-expresados y priorizar candidatos.
- **Control de calidad estadística:** corrección de Benjamini-Hochberg, suavizado
  bayesiano empírico (eBayes) y detección robusta de columnas.
- **Visualización:** volcano plot, MA plot y heatmaps con clustering.
- **Reproducibilidad:** pipeline automatizado de código abierto, exportación en CSV.
- **Preparación posterior:** conjuntos de genes UP/DOWN para análisis de enriquecimiento
  funcional y estudios de firma génica.

## Datos

Dataset público `ALL` (Chiaretti et al., 2004, *Blood* — paquete `ALL` de Bioconductor),
plataforma Affymetrix HG-U95Av2:

- `datos/exprs.csv` — matriz de expresión: 12 625 sondas × 128 muestras.
- `datos/pData.csv` — metadatos: subtipo molecular (`mol.biol`), datos clínicos y citogenéticos.
- Comparación de la práctica: **74 muestras NEG vs 37 BCR/ABL** (111 muestras).

## Contenido

| Archivo | Descripción |
|---|---|
| `Analisis_microarrays_LLA_Colab.ipynb` | Notebook completo para Google Colab (recomendado) |
| `Limma_inmoose.py` | Mismo análisis como script con celdas `# %%` para Spyder |
| `inmoose_env.yaml` | Entorno conda para ejecución local |
| `datos/` | Datos de expresión y metadatos |

## Resultados esperados

Con los umbrales por defecto (FDR < 0.05 y |log2FC| > 0.5) deberías obtener
aproximadamente **760 genes con FDR significativo**, de los cuales unos **270
sobre-expresados** y **65 infra-expresados** superan también el umbral de fold change.
Pequeñas variaciones son normales según la versión de inmoose.

## Instrucciones

### Opción A — Google Colab (recomendada)

1. Abre el notebook con el botón **Open in Colab** de arriba.
2. **Entorno de ejecución → Ejecutar todas**.
3. La celda de instalación de `inmoose` tarda 5-10 minutos la primera vez.
4. Descarga el ZIP de resultados al final.

### Opción B — Local con conda (Spyder)

```bash
conda env create -f inmoose_env.yaml
conda activate inmoose_env
spyder
```

Abre `Limma_inmoose.py`, ajusta las rutas de `expr_path` y `pdata_path` a la carpeta
`datos/` y ejecuta por celdas (`# %%`).
