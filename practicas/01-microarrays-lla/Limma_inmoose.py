"""
Created on Mon Nov 24 16:06:11 2025
@author: vjsanchez

Análisis de Datos de Microarrays con InMoose + Limma
Comparación: BCR/ABL vs NEG (Leucemia Linfoide Aguda)
VERSIÓN FINAL Y COMPLETA con Heatmap CORREGIDO y PALETA ACTUALIZADA
"""

# %%[0] Importación de librerías 
# ============================================================================
# CELDA 1: Importación de librerías y Configuración de parámetros
# ============================================================================

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
import matplotlib.pyplot as plt
import seaborn as sns 
from matplotlib.patches import Patch

# Suprimir advertencias para output más limpio
warnings.filterwarnings('ignore')

# Importar limma desde inmoose (implementación Python de limma de R)
from inmoose import limma
import patsy  # Para crear matrices de diseño con fórmulas R-style

print("✓ Librerías importadas correctamente")

# %% [1] Configuración de parámetros y rutas
# ============================================================================
# CONFIGURACIÓN: Define rutas de datos y parámetros estadísticos
# ============================================================================

# --- Rutas a los archivos de datos ---
expr_path = Path("/home/vjsanchez/Documentos/Universidad_UFV/Asignaturas/Master_Bioinfo_UFV/Microarray/Datos/exprs.csv")
pdata_path = Path("/home/vjsanchez/Documentos/Universidad_UFV/Asignaturas/Master_Bioinfo_UFV/Microarray/Datos/pData.csv")

# --- Directorio de salida para resultados ---
OUT_DIR = Path("/home/vjsanchez/Documentos/Universidad_UFV/Asignaturas/Master_Bioinfo_UFV/Microarray/Datos/limma_results_LLA")
OUT_DIR.mkdir(exist_ok=True)  # Crear directorio si no existe

# --- Umbrales estadísticos para significancia ---
P_VALUE_THRESHOLD = 0.05    # FDR < 5% (corrección múltiple de Benjamini-Hochberg)
LOG_FC_THRESHOLD = 0.5      # |log2 Fold Change| > 0.5 (≈ 1.4x cambio en expresión)

print(f"✓ Configuración: FDR < {P_VALUE_THRESHOLD}, log2FC > |{LOG_FC_THRESHOLD}|")

# %% [3] Carga de datos
# ============================================================================
# CARGA Y PREPARACIÓN: Leer datos de expresión y metadatos
# ============================================================================

print("\n📊 CARGANDO Y ALINEANDO DATOS LLA...")

# Cargar matriz de expresión (genes × muestras)
# index_col=0 usa la primera columna como nombres de genes
datos_expresion_full = pd.read_csv(expr_path, index_col=0) 

# Cargar metadatos (información de muestras: grupos, características)
# index_col=0 usa la primera columna como IDs de muestras
design_matrix_full = pd.read_csv(pdata_path, index_col=0) 

# --- ALINEAMIENTO CRÍTICO: Asegurar que muestras coincidan en ambos archivos ---
# Intersección de nombres de muestras entre expresión y metadatos
common_samples = datos_expresion_full.columns.intersection(design_matrix_full.index)

# Filtrar ambos datasets para mantener solo muestras comunes
datos_expresion_full = datos_expresion_full[common_samples]
design_matrix_full = design_matrix_full.loc[common_samples]

print(f"✓ {len(common_samples)} muestras alineadas correctamente")

# %% [4] Filtrado de muestras BCR/ABL vs NEG
# ============================================================================
# FILTRADO: Seleccionar solo los grupos de interés para el análisis
# ============================================================================

print("\n🔬 FILTRANDO MUESTRAS PARA BCR/ABL vs NEG...")

# Definir grupos de interés para la comparación
grupos_interes = ['BCR/ABL', 'NEG']

# Filtrar metadatos para mantener solo muestras de estos grupos
design_lla = design_matrix_full[design_matrix_full['mol.biol'].isin(grupos_interes)].copy()

# Renombrar columna para evitar problemas con puntos en nombres
design_lla = design_lla.rename(columns={'mol.biol': 'mol_biol'}) 

# Obtener IDs de muestras filtradas
muestras_lla = design_lla.index

# Filtrar matriz de expresión para mantener solo muestras de interés
datos_expresion = datos_expresion_full[muestras_lla]

# Convertir a matriz numpy para limma (requerimiento de inmoose)
expr_matrix = datos_expresion.values.astype(float) 

# --- PREPARAR VARIABLE CATEGÓRICA ---
# Convertir a categoría y definir orden (NEG como referencia)
design_lla['mol_biol'] = design_lla['mol_biol'].astype('category').cat.reorder_categories(
    ['NEG', 'BCR/ABL'], ordered=True  # NEG será el grupo de referencia
)

print(f"✓ {len(muestras_lla)} muestras filtradas: {design_lla['mol_biol'].value_counts().to_dict()}")

# %% [5] Preparando la matriz de diseño y análisis con Limma
# ============================================================================
# MODELADO: Crear matriz de diseño y ejecutar análisis diferencial con limma
# ============================================================================

print("\n⚙️ PREPARANDO MATRIZ DE DISEÑO CON PATSY...")

# Fórmula estilo R: Treatment coding con NEG como referencia
# C() indica variable categórica, Treatment() define grupo de referencia
formula = "~ C(mol_biol, Treatment(reference='NEG'))"

# Crear matriz de diseño usando patsy (similar a model.matrix en R)
design_df = patsy.dmatrix(formula, design_lla, return_type="dataframe")

# Convertir a DesignMatrix de patsy (compatible con limma)
design = patsy.DesignMatrix(design_df.values, design_df.design_info)
coef_names = design.design_info.column_names  # Nombres de coeficientes

print(f"✓ Matriz de diseño creada: {coef_names}")

# --- ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON LIMMA ---
print("\n🧬 INICIANDO ANÁLISIS LIMMA...")

# 1. AJUSTAR MODELO LINEAL PARA CADA GEN
# lmFit ajusta: Expresión_gen = β₀ + β₁*grupo + ε para cada gen
fit = limma.lmFit(expr_matrix, design)

# 2. DEFINIR CONTRASTE DE INTERÉS (BCR/ABL vs NEG)
# Buscar automáticamente el nombre del coeficiente para BCR/ABL
bcr_name = [c for c in coef_names if 'BCR' in c or 'BCR/ABL' in c or 'BCR_ABL' in c][0]

# Crear matriz de contrastes: 1 para BCR/ABL, 0 para otros coeficientes
contrast_df = pd.DataFrame(0, index=coef_names, columns=["BCRABL_vs_NEG"])
contrast_df.loc[bcr_name, "BCRABL_vs_NEG"] = 1

# 3. APLICAR CONTRASTES AL MODELO
fit2 = limma.contrasts_fit(fit, contrast_df)

# 4. MODERACIÓN BAYESIANA (eBayes) - PASO CLAVE
# Estabiliza varianzas usando información de todos los genes
fit2 = limma.eBayes(fit2)

# 5. EXTRAER RESULTADOS
# topTable: tabla ordenada por significancia
# adjust_method="fdr_bh": corrección FDR de Benjamini-Hochberg
results_df = pd.DataFrame(limma.topTable(fit2, coef="BCRABL_vs_NEG", number=np.inf, adjust_method="fdr_bh"))
results_df.index = datos_expresion.index  # Restaurar nombres de genes

print("✓ Análisis limma completado")

# --- PROCESAR Y CLASIFICAR RESULTADOS ---
print("\n📈 ANALIZANDO RESULTADOS...")

def detect_cols_robust(res_df):
    """
    Función robusta para detectar automáticamente nombres de columnas
    en diferentes versiones de limma/inmoose
    """
    # Buscar columna de p-valores ajustados
    p_candidates = ['adj.P.Val', 'FDR', 'padj']
    p_col = next((col for col in p_candidates if col in res_df.columns), None)
    
    if p_col is None:
        # Búsqueda más flexible
        p_candidates_generic = [col for col in res_df.columns if any(x in col.lower() for x in ['adj', 'fdr', 'padj', 'p.adj'])]
        p_col = p_candidates_generic[0] if p_candidates_generic else None
    
    if p_col is None:
        raise KeyError("No se pudo identificar la columna de p-valores ajustados.")

    # Buscar columna de log Fold Change
    log_candidates = ['logFC', 'log2FC']
    log_col = next((col for col in log_candidates if col in res_df.columns), None)
    
    if log_col is None:
        log_candidates_generic = [col for col in res_df.columns if 'FC' in col or 'log' in col.lower()]
        log_col = log_candidates_generic[0] if log_candidates_generic else None
    
    if log_col is None:
        raise KeyError(f"No se pudo identificar la columna de Log Fold Change. Columnas disponibles: {res_df.columns.tolist()}")
        
    return p_col, log_col

# Detectar nombres de columnas automáticamente
try:
    p_col, log_col = detect_cols_robust(results_df)
    print(f"✓ Columnas detectadas: p-valor='{p_col}', logFC='{log_col}'")
except KeyError as e:
    raise e

# --- CLASIFICAR GENES POR SIGNIFICANCIA ---
# Genes significativos (FDR < umbral)
genes_significativos = results_df[results_df[p_col] < P_VALUE_THRESHOLD]

# Subclasificar en sobre e infra expresados
genes_sobre = genes_significativos[genes_significativos[log_col] > LOG_FC_THRESHOLD]
genes_infra = genes_significativos[genes_significativos[log_col] < -LOG_FC_THRESHOLD]

# Añadir columna de categoría para visualizaciones
results_df['Categoria'] = 'No significativo'
results_df.loc[(results_df[p_col] < P_VALUE_THRESHOLD) & (results_df[log_col] > LOG_FC_THRESHOLD), 'Categoria'] = 'Sobre-expresado'
results_df.loc[(results_df[p_col] < P_VALUE_THRESHOLD) & (results_df[log_col] < -LOG_FC_THRESHOLD), 'Categoria'] = 'Infra-expresado'

# Paleta de colores para visualizaciones
paleta = {'No significativo': 'gray', 'Sobre-expresado': 'red', 'Infra-expresado': 'blue'}

print(f"✓ Genes clasificados: {len(genes_sobre)} sobre-expresados, {len(genes_infra)} infra-expresados")

# %% [6] MA Plot y Volcano Plot 
# ============================================================================
# VISUALIZACIÓN: Gráficos para explorar resultados de expresión diferencial
# ============================================================================

print("\n📉 CREANDO MA PLOT...")
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(10, 6))

# MA Plot: logFC vs Intensidad promedio
for categoria, color in paleta.items():
    subset = results_df[results_df['Categoria'] == categoria]
    ax.scatter(subset['AveExpr'], subset[log_col], c=color, alpha=0.6, s=15, 
               label=f'{categoria} (n={len(subset)})')

# Líneas de referencia
ax.axhline(0, color='black', linestyle='-', alpha=0.5)  # Línea de no cambio
ax.axhline(LOG_FC_THRESHOLD, color='black', linestyle='--', alpha=0.5)  # Umbral superior
ax.axhline(-LOG_FC_THRESHOLD, color='black', linestyle='--', alpha=0.5)  # Umbral inferior

ax.set_xlabel('Expresión Media (AveExpr)', fontsize=12)
ax.set_ylabel('log2 Fold Change', fontsize=12)
ax.set_title('MA Plot - Expresión Diferencial (BCR/ABL vs NEG)', fontsize=14, fontweight='bold')
ax.legend()
plt.tight_layout()
fig.savefig(OUT_DIR / "MA_Plot_BCRABL_vs_NEG.png", dpi=150)
plt.show()

print("\n🎨 CREANDO VOLCANO PLOT...")
fig, ax = plt.subplots(figsize=(10, 8))

# Volcano Plot: logFC vs -log10(p-value)
for categoria, color in paleta.items():
    subset = results_df[results_df['Categoria'] == categoria]
    ax.scatter(subset[log_col], -np.log10(subset[p_col]), c=color, alpha=0.6, s=15, 
               label=f'{categoria} (n={len(subset)})')

# Líneas de umbrales estadísticos
ax.axhline(-np.log10(P_VALUE_THRESHOLD), color='black', linestyle='--', alpha=0.7)  # Umbral FDR
ax.axvline(LOG_FC_THRESHOLD, color='black', linestyle='--', alpha=0.7)  # Umbral FC positivo
ax.axvline(-LOG_FC_THRESHOLD, color='black', linestyle='--', alpha=0.7)  # Umbral FC negativo

ax.set_xlabel('log2 Fold Change (BCR/ABL vs NEG)', fontsize=12)
ax.set_ylabel(f'-log10({p_col})', fontsize=12)
ax.set_title('Volcano Plot - Expresión Diferencial (BCR/ABL vs NEG)', fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(OUT_DIR / "volcano_BCRABL_vs_NEG.png", dpi=150)
plt.show()

# %% [7] Heatmaps
# ============================================================================
# HEATMAPS: Visualización de patrones de expresión en genes significativos
# ============================================================================

print("\n🔥 CREANDO HEATMAPS ESENCIALES...")

if len(genes_significativos) > 0:
    
    # Paleta de colores para grupos (anotaciones)
    grupo_colores = {'BCR/ABL': '#FF8C00', 'NEG': '#008000'}  # Naranja y Verde
    
    # ========================================================================
    # PREPARACIÓN DE DATOS PARA HEATMAPS
    # ========================================================================
    
    # Seleccionar top genes más significativos (máximo 40)
    N_TOP_GENES = min(40, len(genes_significativos))
    top_genes_indices = genes_significativos.nsmallest(N_TOP_GENES, p_col).index
    
    # Extraer datos de expresión para genes seleccionados
    heatmap_data = datos_expresion.loc[top_genes_indices]
    
    # Ordenar muestras por grupo para mejor visualización
    design_heatmap = design_lla.copy()
    design_heatmap['mol_biol'] = pd.Categorical(design_heatmap['mol_biol'], 
                                               categories=['NEG', 'BCR/ABL'], ordered=True)
    design_heatmap = design_heatmap.sort_values('mol_biol')
    muestras_ordenadas = design_heatmap.index.tolist()
    heatmap_data = heatmap_data[muestras_ordenadas]
    
    # Normalización Z-score por filas (genes)
    # Esto permite comparar patrones de expresión entre genes
    heatmap_zscore = heatmap_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    heatmap_zscore = heatmap_zscore.dropna()  # Eliminar genes con varianza cero
    
    print(f"✓ Datos preparados: {len(heatmap_zscore)} genes × {len(muestras_ordenadas)} muestras")
    
    # ========================================================================
    # HEATMAP 1: EXPRESIÓN GENERAL (Clustermap principal)
    # ========================================================================
    
    print("\n📊 Creando Heatmap Principal...")
    
    # Preparar anotaciones de columnas (colores por grupo)
    col_colors = design_heatmap['mol_biol'].map(grupo_colores)
    
    # Crear clustermap con clustering
    g = sns.clustermap(
        heatmap_zscore,
        cmap='RdBu_r',      # Esquema rojo-azul (rojo=sobre, azul=infra)
        center=0,           # Centrar en 0 (expresión media)
        vmin=-2, vmax=2,    # Límites para mejor contraste visual
        row_cluster=True,   # Agrupar genes por similitud de expresión
        col_cluster=False,  # No agrupar muestras (mantener orden por grupo)
        col_colors=col_colors,  # Barra de colores para grupos
        yticklabels=True,   # Mostrar nombres de genes
        xticklabels=False,  # No mostrar nombres de muestras (muchas)
        figsize=(12, 10),
        dendrogram_ratio=0.12,  # Tamaño relativo del dendrograma
        cbar_pos=(0.02, 0.87, 0.03, 0.1),  # Posición de la barra de color
        linewidths=0,       # Sin líneas entre celdas
        cbar_kws={'label': 'Z-score (σ)'}  # Etiqueta barra de color
    )
    
    # Mejorar etiquetas y formato
    g.ax_heatmap.set_ylabel('Genes', fontsize=11, fontweight='bold')
    g.ax_heatmap.set_xlabel('Muestras', fontsize=11, fontweight='bold')
    g.ax_heatmap.tick_params(axis='y', labelsize=7)  # Tamaño fuente genes
    
    # Crear leyenda para grupos
    legend_grupos = [
        Patch(facecolor=grupo_colores['NEG'], edgecolor='black', 
              label=f"NEG (n={(design_heatmap['mol_biol'] == 'NEG').sum()})"),
        Patch(facecolor=grupo_colores['BCR/ABL'], edgecolor='black', 
              label=f"BCR/ABL (n={(design_heatmap['mol_biol'] == 'BCR/ABL').sum()})")
    ]
    
    # Añadir leyenda al gráfico
    g.ax_col_dendrogram.legend(handles=legend_grupos, title='Subtipo Molecular', 
                              loc='center left', bbox_to_anchor=(0.03, 0.5), 
                              frameon=True, fontsize=10, title_fontsize=11)
    
    # Título principal
    g.fig.suptitle(f'Top {len(heatmap_zscore)} Genes Diferencialmente Expresados\nBCR/ABL vs NEG', 
                   fontsize=14, fontweight='bold', y=0.97)
    
    # Guardar heatmap
    plt.savefig(OUT_DIR / "heatmap_principal_BCRABL_vs_NEG.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✓ Heatmap principal guardado")
    
    # ========================================================================
    # HEATMAP 2: GENES UP vs DOWN SEPARADOS
    # ========================================================================
    
    print("\n↕️ Creando Heatmap de Genes UP vs DOWN...")
    
    # Seleccionar top genes sobre e infra expresados por separado
    N_MAX_UP_DOWN = 25
    genes_up_top = genes_sobre.nsmallest(min(N_MAX_UP_DOWN, len(genes_sobre)), p_col).index
    genes_down_top = genes_infra.nsmallest(min(N_MAX_UP_DOWN, len(genes_infra)), p_col).index
    
    if len(genes_up_top) > 0 and len(genes_down_top) > 0:
        # Combinar genes UP y DOWN (UP primero para orden claro)
        genes_combined = list(genes_up_top) + list(genes_down_top)
        heatmap_combined = datos_expresion.loc[genes_combined, muestras_ordenadas]
        
        # Aplicar normalización Z-score
        heatmap_combined_z = heatmap_combined.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
        heatmap_combined_z = heatmap_combined_z.dropna()
        
        # Contar genes finales después de limpieza
        num_genes_up = len(genes_up_top.intersection(heatmap_combined_z.index))
        num_genes_down = len(genes_down_top.intersection(heatmap_combined_z.index))
        
        # --- Mismo formato que el primer heatmap para consistencia ---
        col_colors = design_heatmap['mol_biol'].map(grupo_colores)
        
        # Crear clustermap sin clustering para mantener orden UP/DOWN
        g = sns.clustermap(
            heatmap_combined_z,
            cmap='RdBu_r',
            center=0,
            vmin=-2.5, vmax=2.5,  # Rango ligeramente mayor para mejor contraste
            row_cluster=False,     # Mantener orden original (UP primero, DOWN después)
            col_cluster=False,     # Mantener orden por grupo
            col_colors=col_colors, # Mismas anotaciones de grupo
            yticklabels=True,
            xticklabels=False,
            figsize=(14, 10),
            dendrogram_ratio=0.1,
            cbar_pos=(0.02, 0.83, 0.03, 0.12),
            linewidths=0,
            cbar_kws={'label': 'Z-score (σ)'}
        )
        
        # Mismo formato de etiquetas
        g.ax_heatmap.set_ylabel('Genes', fontsize=11, fontweight='bold')
        g.ax_heatmap.set_xlabel('Muestras', fontsize=11, fontweight='bold')
        g.ax_heatmap.tick_params(axis='y', labelsize=7)
        
        # LEYENDA COMBINADA (grupos + regulación)
        legend_grupos = [
            Patch(facecolor=grupo_colores['NEG'], edgecolor='black', 
                  label=f"NEG (n={(design_heatmap['mol_biol'] == 'NEG').sum()})"),
            Patch(facecolor=grupo_colores['BCR/ABL'], edgecolor='black', 
                  label=f"BCR/ABL (n={(design_heatmap['mol_biol'] == 'BCR/ABL').sum()})")
        ]
        
        legend_reg = [
            Patch(facecolor='red', edgecolor='black', label=f'↑ UP (n={num_genes_up})'),
            Patch(facecolor='blue', edgecolor='black', label=f'↓ DOWN (n={num_genes_down})')
        ]
        
        # Leyenda unificada
        g.ax_col_dendrogram.legend(
            handles=legend_grupos + legend_reg, 
            title='Leyenda', 
            loc='center left', 
            bbox_to_anchor=(0.03, 0.8), 
            frameon=True, 
            fontsize=10,
            title_fontsize=11
        )
        
        # Título informativo
        g.fig.suptitle('Top Genes Sobre e Infra-expresados\n' +
                       f'BCR/ABL vs NEG | ↑ UP: {num_genes_up} genes | ↓ DOWN: {num_genes_down} genes', 
                       fontsize=14, fontweight='bold', y=0.97)
        
        plt.savefig(OUT_DIR / "heatmap_up_down_BCRABL_vs_NEG.png", dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✓ Heatmap UP vs DOWN guardado (formato consistente)")
    
    # ========================================================================
    # RESUMEN DE HEATMAPS
    # ========================================================================
    
    print("\n📊 RESUMEN DE HEATMAPS GENERADOS:")
    print(f"  1. Clustermap Principal: {len(heatmap_zscore)} genes más significativos")
    if len(genes_up_top) > 0 and len(genes_down_top) > 0:
        print(f"  2. Heatmap UP/DOWN: {num_genes_up} genes UP + {num_genes_down} genes DOWN")
    print(f"\n✓ Todos los heatmaps guardados en: {OUT_DIR}")
    
else:
    print("⚠️ No hay genes significativos para crear heatmaps")

print("\n--- HEATMAPS COMPLETADOS ---")

# %% [8] Guardar resultados
# ============================================================================
# EXPORTACIÓN: Guardar resultados tabulares y resumen final
# ============================================================================

print("\n💾 GUARDANDO RESULTADOS...")

# Mostrar top 10 genes más significativos
print("\n🏆 TOP 10 GENES MÁS SIGNIFICATIVOS:")
top_genes = results_df.nsmallest(10, p_col)
for i, (gen, row) in enumerate(top_genes.iterrows(), 1):
    direccion = "↑ SOBRE" if row[log_col] > 0 else "↓ INFRA"
    print(f"  {i:2d}. {gen}: logFC={row[log_col]:6.2f} {direccion} {p_col}={row[p_col]:.2e}")

# --- GUARDAR ARCHIVOS DE RESULTADOS ---
# Resultados completos (todos los genes)
results_df.to_csv(OUT_DIR / "results_BCRABL_vs_NEG.csv")

# Subconjuntos por significancia
genes_significativos.to_csv(OUT_DIR / "sig_BCRABL_vs_NEG.csv")      # Todos los significativos
genes_sobre.to_csv(OUT_DIR / "up_BCRABL_vs_NEG.csv")               # Solo sobre-expresados
genes_infra.to_csv(OUT_DIR / "down_BCRABL_vs_NEG.csv")             # Solo infra-expresados

print(f"✓ Archivos guardados en: {OUT_DIR}")
print("\n--- ANÁLISIS COMPLETADO EXITOSAMENTE ---")