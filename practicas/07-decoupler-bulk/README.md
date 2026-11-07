# Práctica 7 — Análisis funcional con decoupler: actividad de TFs y rutas en bulk RNA-seq

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/victorsanchezarevalo/transcriptomica-master/blob/main/practicas/07-decoupler-bulk/Decoupler_Bulk_Colab.ipynb)

Inferencia de **actividad de factores de transcripción y rutas de señalización** con `decoupler` (2.x) sobre el dataset **GSE151251**: células madre hematopoyéticas (HSCs) humanas, 3 controles frente a 3 muestras tratadas con **TGFβ**.

## Contenido

- `Decoupler_Bulk_Colab.ipynb` — notebook de Google Colab.
- `datos/counts.txt` — matriz de conteos crudos (salida de featureCounts: ~64.000 genes × 6 muestras).

## Qué se hace en el notebook

1. Filtrado de genes por conteos (criterio `filterByExpr` de edgeR, con `dc.pp.filter_by_expr`).
2. Expresión diferencial TGFβ vs control con **PyDESeq2** (volcano plot).
3. **Actividad de TFs**: red **CollecTRI** + método **ULM** (barplot, red de TFs, volcano por TF).
4. **Actividad de rutas**: modelo **PROGENy** + método **MLM** (barplot y scatter de genes diana de TGFb).
5. **Enriquecimiento** con **Hallmark (MSigDB)** + ORA (barplot y running score tipo GSEA).
6. *(Opcional)* Interacciones **ligando–receptor** con LIANA + ULM.

## Notas

- El notebook usa la **API de decoupler 2.x** (`dc.mt.ulm`, `dc.op.collectri`, `dc.pl.barplot`…), distinta de la 1.x de los tutoriales antiguos (`dc.run_ulm`, `dc.get_collectri`…).
- La descarga de los recursos (CollecTRI, PROGENy, Hallmark) y LIANA necesita internet.
- Resultado esperado de control: la ruta **TGFb de PROGENy** y el Hallmark **TGF_BETA_SIGNALING** deben aparecer como los más activados.
