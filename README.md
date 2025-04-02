---
title: "Análisis de datos ómicos (M0-157) Primera prueba de evaluación continua"
subtitle: "Código R para la exploración de los datos debidamente comentado"
author: "Gabriel Peña Peña"
date: "`r Sys.Date()`"
output:
---

# Selección de Dataset

Para esta PEC seleccione el dataset [2024-Cachexia](https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv), disponible en el repositorio de GitHub del proyecto [metaboData](https://github.com/nutrimetabolomics/metaboData). Todos los conjuntos propuestos están basados en datos de metabolitos, por lo que la elección se centró en el contenido biológico del estudio. En este caso, el dataset se relaciona con la cachexia, un síndrome metabólico complejo asociado a la pérdida severa de masa muscular, común en enfermedades crónicas como el cáncer. El cual me pareció una opción interesante por su relevancia clínica y su uso para el análisis exploratorio propuesto.

# SummarizedExperiment 

## Breve explicación

La clase SummarizedExperiment es una estructura desarrollada en el entorno Bioconductor para organizar datos experimentales de alto rendimiento, como transcriptómica o metabolómica. A diferencia de ExpressionSet, que fue diseñada principalmente para microarrays, SummarizedExperiment tiene la capacidad de almacenar múltiples matrices de datos junto con sus metadatos asociados, tanto para las muestras como para las variables medidas. Al poseer mayor flexibilidad y modularidad permite un mejor manejo de datasets complejos y se integra mejor con multiples paquetes y herramientas bioinformáticas.

## Creación del SummarizedExperiment

Para la creación del objeto SummarizedExperiment primero se cargan las librerias necesarias

```{r message=FALSE, warning=FALSE}
# Se cargan las librerias necesarias
library(dplyr)
library(SummarizedExperiment)
```

Luego se carga el archivo CSV que contiene el dataset de metabolómica para luego separar los datos de metabolitos de los metadatos (ID y grupo de pérdida muscular),

```{r}
# Se cargar el dataset
human_cachexia <- read.csv("human_cachexia.csv", header = TRUE)

# Se separan metadatos y datos de metabolitos
metadata <- human_cachexia %>%
  select(Patient.ID, Muscle.loss)

metabolitos <- human_cachexia %>%
  select(-Patient.ID, -Muscle.loss)
```

Despues se preparan los metadatos convertidos en dataframe y la matriz de metabolitos. Tambien se renombran las filas y columnas para garantizar que coincidan las muestras entre ambos objetos.

```{r}
# Se crea el colData con los metadatos
colData <- data.frame(MuscleLoss = metadata$Muscle.loss)
rownames(colData) <- metadata$Patient.ID

# Se crea la matriz de datos de metabolitos
colnames(metabolitos) <- colnames(human_cachexia)[-c(1, 2)]
rownames(metabolitos) <- metadata$Patient.ID 
metabolitos <- t(as.matrix(metabolitos))
```

Finalmente se construye el objeto SummarizedExperiment, que almacenara los datos experimentales y la información asociada a cada muestra, lo que facilitara los posteriores análisis.

```{r}
# Se crea el objeto SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(metabolitos = metabolitos),
  colData = colData
)

# Se visualiza el objeto SummarizedExperiment.
se
```

```{r}
# Se exporta el objeto SummarizedExperiment en formato binario (.Rda).
save(se, file = "summarized_experiment_2024_cachexia.Rda")

# Se exportan los datos de metabolitos en formato csv
write.csv(assays(se)$metabolitos, file = "datos_metabolitos.csv", row.names = TRUE)

# Se exportan los datos  metadatos en formato csv
write.csv(as.data.frame(colData(se)), file = "metadatos_pacientes.csv", row.names = TRUE)
```

# Análisis exploratorio del dataset

## Estructura general

Se inspecciona el contenido del objeto SummarizedExperiment creado verificando su estructura, dimensiones, primeras y ultimas filas de datos y metadatos, ademas se observa la distribución de la variable MuscleLoss, que inciia si el pacinete es sano o presenta la nefemedad, tambien se realiza un resumen estadisticos de los datos de metabolitos. Esto permite verificar que los datos se ha cargado correctamente y tener una idea general de su estructura.

```{r}
# Se visualiza estructura del objeto SummarizedExperiment
str(se)
```

```{r}
# Se visualizan dimensiones de la matriz de metabolitos y metadatos
dim(assays(se)$metabolitos)
dim(colData(se))
```

```{r}
# Se visualizan las primeras filas de los metabolitos
head(assays(se)$metabolitos)
```

```{r}
# Se visualizan las primeras y ultimas filas de los metadatos
head(colData(se))
tail(colData(se))
```

```{r}
# Resumen de los datos de metabolitos
summary(t(assays(se)$metabolitos))
```

Se observa cómo están distribuidos los grupos de pérdida muscular en la muestra, entregando una idea del balance de clases.

```{r}
# Se visualiza la distribución de la variable 'MuscleLoss'
tabla_1 <- table(colData(se)$MuscleLoss)
tabla_1
```
```{r}
# Se visualiza la distribución porcentual de la variable 'MuscleLoss'
tabla_2  <- round(prop.table(tabla_1) * 100, 2)
tabla_2
```

## Concentración de metabolitos

Se representan las concentraciones del primer metabolito entre pacientes y según la condición MuscleLoss, utilizando histogramas y boxplots. Además, se visualiza la distribución de concentraciones de todos los metabolitos por paciente, lo que permite explorar el rango, la variabilidad y obtener una primera aproximación sobre aquellos pacientes que podrían presentar concentraciones metabólicas anómalas.

```{r}
# Distribución de la concentración del primer metabolito
hist(assays(se)$metabolitos[1,], 
     main="Distribución concentración del primer metabolito", 
     xlab="Concentración",
     ylab="Frecuencia")
```


```{r}
# Distribución concentración del primer metabolito según grupo 'MuscleLoss'
boxplot(assays(se)$metabolitos[1,] ~ colData(se)$MuscleLoss,
        main="Distribución del primer metabolito por condición",
        xlab="Condición", ylab="Concentración")
```

```{r}
# Boxplot general de todos los metabolitos
boxplot(assays(se)$metabolitos, 
        main="Distribución global de metabolitos",
        outline = FALSE, las = 2, cex.axis = 0.7)
```

## Análisis de correlación entre metabolitos

Se calculan y exploran las posibles relaciones, tanto directas como inversas, entre los distintos metabolitos. Esto permite identificar grupos de metabolitos con comportamientos similares entre pacientes, lo que podria sugerir relaciones funcionales.

```{r}
# Se calcula la matriz de correlación entre metabolitos
cor_matrix <- cor(t(assays(se)$metabolitos))

# Se visualización como heatmap
heatmap(cor_matrix, 
        main="Correlación entre metabolitos")
```

## Análisis de variabilidad de metabolitos

Se evalúa la variabilidad de los metabolitos entre pacientes mediante el cálculo del coeficiente de variación (CV), una medida que relaciona la desviación estándar con la media. Esto permite comparar metabolitos con diferentes escalas de concentración e identificar cuáles presentan mayor o menor dispersión entre muestras. Posteriormente, se representa gráficamente, mediante un mapa de calor, el comportamiento de los 20 metabolitos con mayor CV, lo que permite visualizar patrones de variación entre pacientes.

```{r}
# Se calcula la media de cada metabolito
media_metabolitos <- apply(assays(se)$metabolitos, 1, mean)

# Se calcula la desviación estándar de cada metabolito
sd_metabolitos <- apply(assays(se)$metabolitos, 1, sd)

# Se calcula el coeficiente de variación
cv_metabolitos <- sd_metabolitos / media_metabolitos

# Se visualizan los resultados
head(cv_metabolitos)
```

### Represetación gráfica

```{r}
# Se seleccionan los 20 metabolitos más variables
top_var_idx <- order(cv_metabolitos, decreasing = TRUE)[1:20]
top_var_data <- metabolitos[top_var_idx, ]

# Se escalan los datos para visualización
top_var_data_scaled <- t(scale(t(top_var_data)))

# Se genera un mapa de valor
heatmap(top_var_data_scaled,
        main = "Metabolitos más variables",
        Colv = NA,
        scale = "none")
```

## Analisis ANOVA para comparar grupos

Se realiza un análisis ANOVA para evaluar si existen diferencias estadísticamente significativas en la concentración de cada metabolito entre los grupos con y sin pérdida muscular. A partir de este análisis, se extraen los p-valores y se organiza una tabla con los metabolitos ordenados según su nivel de significancia. Finalmente, se visualizan los 10 metabolitos más significativos.

```{r}
# Se extrae la matriz de concentraciones de metabolitos
matriz_metabolitos <- assays(se)$metabolitos

# Se extrae la condición
grupo_muscular <- colData(se)$MuscleLoss

# Se aplica ANOVA para cada metabolito, se extrae el p-valor
pvalores <- apply(matriz_metabolitos, 1, function(x) {
  modelo <- aov(x ~ grupo_muscular) 
  summary(modelo)[[1]][["Pr(>F)"]][1] 
})

# Se crear una tabla con metabolitos y sus p-valores
tabla_anova <- data.frame(
  Metabolito = rownames(matriz_metabolitos),
  p_valor = pvalores
)
rownames(tabla_anova) <- NULL 

# Se ordena la tabla por p-valor mas significativo
tabla_anova_ordenada <- tabla_anova[order(tabla_anova$p_valor), ]
rownames(tabla_anova_ordenada) <- NULL 

# Se visualizan los 10 metabolitos más significativos
head(tabla_anova_ordenada, 10)
```

### Represetación gráfica

```{r}
# Se muestran los 10 metabolitos más significativos
top10 <- head(tabla_anova_ordenada, 10)

# Se genera un grafico con los p-valores
barplot(-log10(top10$p_valor),
        names.arg = top10$Metabolito,
        las = 2,
        main = "Top 10 metabolitos con diferencias mas signifcativas",
        ylab = "-log10(p-valor)")
```

