# Enunciado original — Tarea RNA-seq

Tarea RNA_seq: Control de Calidad, Trimming y Mapeo o Pseudoalineamiento

Análisis del Rol de UPF1 en la Diferenciación de Células Madre Embrionarias Humanas

Contexto Biológico: El Nonsense-Mediated RNA Decay (NMD) es un mecanismo de control de calidad celular que degrada ARNm con codones de terminación prematuros. UPF1 es la proteína central de este proceso, actuando como helicasa y regulando la estabilidad del ARNm. En células madre embrionarias humanas (hESCs), UPF1 no solo elimina transcritos defectuosos, sino que también regula la expresión de genes importantes para el mantenimiento de la pluripotencia y la diferenciación celular.

Objetivo del Ejercicio: Analizar el impacto transcriptómico de la pérdida de función de UPF1 en células madre embrionarias humanas mediante técnicas de RNA-seq, identificando genes diferencialmente expresados y evaluando su relevancia biológica.

Diseño Experimental

Condiciones:

    Control: hESCs tratadas con siRNA control (scrambled)

    Tratamiento: hESCs tratadas con siRNA contra UPF1 (knockdown)

Datos Disponibles:

    Plataforma: Illumina HiSeq 2500

    Tipo de lecturas: 50bp single-end, stranded

    Replicados: 2 réplicas biológicas por condición

    Dataset: GSE81214 (GEO/NCBI)

Muestras:

ID Muestra
	

Condición
	

Descripción
	

Número de Lecturas

SRR3480371
	

Control 1
	

siRNA control
	

42,059,631

SRR3480372
	

Control 2
	

siRNA control
	

54,401,313

SRR3480373
	

Tratamiento 1
	

siRNA UPF1
	

52,991,800

SRR3480374
	

Tratamiento 2
	

siRNA UPF1
	

64,635,061

 

PARTE I: Preparación y Control de Calidad

Actividad 1.1: Descarga y Preparación de Datos

    Descarga las secuencias FASTQ desde SRA usando fastq-dump

    Examina las primeras líneas de cada archivo FASTQ

    Pregunta: ¿Qué información contiene cada línea del formato FASTQ?

Actividad 1.2: Control de Calidad con FastQC

    Ejecuta FastQC en todas las muestras

    Examina los reportes HTML generados

    Preguntas:

        ¿Cuál es la calidad promedio de las lecturas?

        ¿Observas algún sesgo en la composición de bases?

        ¿Hay presencia de adaptadores o contaminantes?

        ¿Recomendarías algún filtrado adicional?

PARTE II: Alineamiento y Mapeo con STAR

Actividad 2.1: Preparación del Genoma de Referencia

    Genera el índice del genoma con STAR usando el cromosoma 19 humano

    Pregunta: ¿Por qué es importante incluir la anotación GTF durante la indexación?

Actividad 2.2: Alineamiento con STAR

    Ejecuta STAR para cada muestra con los parámetros optimizados

    Examina los archivos Log.final.out de cada muestra

    Preguntas:

        ¿Cuál es el porcentaje de lecturas mapeadas en cada muestra?

        ¿Hay diferencias significativas en las tasas de mapeo entre condiciones?

        ¿Qué porcentaje de lecturas se mapea a regiones splice junction?

Actividad 2.3: Análisis de la Calidad del Alineamiento

    Indexa los archivos BAM generados por STAR usando samtools

    Pregunta: ¿Por qué STAR genera archivos BAM ordenados directamente?

PARTE III: Cuantificación de la Expresión Génica con featureCounts

Actividad 3.1: Conteo con featureCounts

    Ejecuta featureCounts para cuantificar la expresión génica

    Examina el archivo de salida gene_counts.txt

    Preguntas:

        ¿Qué información contienen las primeras 6 columnas del archivo?

        ¿Cuántos genes se detectaron en total?

        ¿Qué porcentaje de lecturas se asignó exitosamente a genes?

Actividad 3.2: Exploración de los Datos de Conteo

    Carga el archivo de conteos en Python

    Calcula estadísticas básicas (suma, mediana, distribución)

    Preguntas:

        ¿Cuál es la distribución de conteos por gen?

        ¿Hay diferencias evidentes entre las muestras control y tratamiento?

PARTE IV: Cuantificación Alternativa con Salmon

Actividad 4.1: Preparación del Índice de Transcriptoma

    Genera el índice del transcriptoma con Salmon usando las secuencias de cDNA del cromosoma 19 humano

    Pregunta: ¿Cuáles son las ventajas de la cuantificación a nivel de transcripto vs. alineamiento genómico?

Actividad 4.2: Cuantificación con Salmon

    Ejecuta Salmon para cada muestra en modo quasi-mapping

    Examina los archivos quant.sf de cada muestra

    Preguntas:

        ¿Cuál es el porcentaje de lecturas mapeadas en cada muestra?

        ¿Hay diferencias en las tasas de mapeo comparado con STAR?

        ¿Cómo se comparan los tiempos de ejecución entre STAR y Salmon?

Actividad 4.3: Importación de Datos de Salmon a Nivel de Gen

    Utiliza pymportx para agregar los datos de transcriptos a nivel de gen

    Pregunta: ¿Por qué es recomendable usar pymportx en lugar de sumar directamente los conteos de transcriptos?

PARTE V: Comparación de Métodos de Cuantificación

Actividad 5.1: Comparación STAR vs. Salmon

    Compara los conteos obtenidos con featureCounts (STAR) vs. tximport (Salmon)

    Genera gráficos de correlación entre ambos métodos

    Preguntas:

        ¿Qué nivel de correlación observas entre ambos métodos?

        ¿Hay genes que muestran diferencias importantes entre métodos?

        ¿Cuáles pueden ser las causas de estas diferencias?

Instrucciones para la Entrega de Actividades
Contenido Requerido

Responde completamente a todas las preguntas planteadas en las actividades asignadas.
Formato de Entrega

Debes entregar dos archivos:

    Documento PDF: Contiene las respuestas teóricas y análisis de resultados
    Archivos y Scripts: Contiene el código implementado con comentarios detallados

Requisitos Específicos
Para el Documento PDF:

    Respuestas claras y bien estructuradas
    Análisis e interpretación de resultados
    Gráficos y tablas cuando sean relevantes
    Conclusiones fundamentadas

Scripts:

    Código completamente funcional
    Comentarios exhaustivos que expliquen:
        Propósito de cada función/bloque de código
        Lógica implementada
        Parámetros utilizados
        Interpretación de resultados
    Estructura ordenada con títulos y subtítulos en markdown
    Salidas de ejecución visibles

Criterios de Evaluación

Se evaluará especialmente:

    Claridad en las respuestas: Explicaciones coherentes y bien argumentadas
    Presentación del código: Organización, legibilidad y documentación
    Funcionalidad: El código debe ejecutarse sin errores
    Completitud: Todas las preguntas deben estar respondidas

Recomendaciones

    Revisa que ambos archivos estén completos antes de la entrega
    Verifica que el script se ejecute correctamente de principio a fin
    Utiliza un lenguaje técnico apropiado pero comprensible
    Incluye referencias bibliográficas si es necesario
