# Transcriptómica — Máster en Bioinformática

Material docente del **módulo de Transcriptómica**: prácticas, datos y scripts de análisis.

## Prácticas

| # | Práctica | Descripción | Abrir |
|---|----------|-------------|-------|
| 1 | [Microarrays LLA](practicas/01-microarrays-lla/) | Expresión diferencial **BCR/ABL vs NEG** en Leucemia Linfoide Aguda con limma (InMoose) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/01-microarrays-lla/Analisis_microarrays_LLA_Colab.ipynb) |

## Cómo usar este material (estudiantes)

1. Pulsa el botón **Open in Colab** de la práctica.
2. En Colab: **Entorno de ejecución → Ejecutar todas**.
3. La primera celda instala el paquete `inmoose` y **puede tardar entre 5 y 10 minutos**
   (se compila desde el código fuente). Es normal: espera a que aparezca el mensaje de
   confirmación.
4. Los datos se descargan automáticamente desde este repositorio; no necesitas subir
   ningún archivo.
5. Al final del notebook se genera un ZIP con todos los resultados (CSV + figuras) que
   se descarga a tu ordenador. **Guárdalo**: Colab borra los archivos al cerrar la sesión.

## Ejecución en local (opcional)

Si prefieres trabajar en tu equipo (Spyder, Jupyter…), cada práctica incluye un entorno
conda listo para usar:

```bash
conda env create -f practicas/01-microarrays-lla/inmoose_env.yaml
conda activate inmoose_env
```

## Estructura del repositorio

```
transcriptomica-master/
├── README.md                  ← este archivo
└── practicas/
    └── 01-microarrays-lla/
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
Los datos de expresión proceden del dataset público `ALL` de Bioconductor
(Chiaretti et al., 2004, *Blood*), plataforma Affymetrix HG-U95Av2.
