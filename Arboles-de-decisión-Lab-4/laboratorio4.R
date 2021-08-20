######################################################################################################
###########################    Laboratorio 2 - Analisis de Datos    ##################################
######################################################################################################
###########################    Autores:     Gary Simken             ##################################
###########################                 Sebasti?n Orellana      ##################################
###########################    Fecha:       18 - Junio - 2021       ##################################
######################################################################################################

######################################################################################################
##########################            Bibliotecas          ###########################################
library(modeest)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggcorrplot)
library(reshape2)
library(RColorBrewer)
library(missForest)
library(factoextra)
library(cluster)
library(arulesViz)
library(NbClust)
library("C50")
library("caret")
library(tidyverse)
library(waffle)

# Se limpia el ambiente de variables y los graficos antiguos
rm(list=ls())
if(length(dev.list())!=0){
  dev.off(dev.list()["RStudioGD"])
}




######################################################################################################
##########################            Data-set             ###########################################

#1. Sample code number: id number
#2. Clump Thickness: 1 - 10             CT
#3. Uniformity of Cell Size: 1 - 10     UCS
#4. Uniformity of Cell Shape: 1 - 10    UCSH
#5. Marginal Adhesion: 1 - 10           MA
#6. Single Epithelial Cell Size: 1 - 10 SECS
#7. Bare Nuclei: 1 - 10                 BN
#8. Bland Chromatin: 1 - 10             BC
#9. Normal Nucleoli: 1 - 10             NN
#10. Mitoses: 1 - 10                    M
#11. Class: (2 for benign, 4 for malignant) C

# Se definen los nombres de las columnas, estos son los mismos de provistos por la base de datos
columnas = c("id",
             "clumpThickness",
             "uniformityCellSize",
             "uniformityCellShape",
             "marginalAdhesion",
             "singleEpithCellSize",
             "bareNuclei",
             "blandChromatin",
             "normalNucleoli",
             "mitoses",
             "class"
)
# Se importa la base de dato mediante la url del repositorio UCI Machine Learning
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data"
df = read.csv(url, 
              header = F, 
              sep=",", 
              col.names = columnas)

# Se setea una semilla para las funciones que necesites utilizar nnumeros aleatorios
set.seed(1234)




######################################################################################################
##########################        Limpieza de Datos        ###########################################

# Se eliminan las variables if y class, dado que son datos no relevantes en esta parte del estudio,
# en el caso de id es un simple identificador por lo que no aporta y en el caso de clase, esta es 
# eliminada porque el Clustering es un m??todo no supervisado.
df <- subset(df, select = -c(id))
df.inicial<-df
df.inicial$class = factor(df.inicial$class, levels = c(2,4), labels = c("Benigno", "Maligno"))

# Resumen de los cuartiles de cada variable
summary(df)

# Para los datos "missing" se procede a cambiar los signos de interrogacisn por un NA, que es el valor
# default en los enteros cuando no existe un dato
df$bareNuclei[df$bareNuclei == "?"] = NA

# Por otro lado convertiremos todos los datos a tipo numerico para que no existan problemas al 
# momento de analizarlos
df[ , ] = apply(df[ , ], 2,function(x) as.numeric(as.character(x)) )

# Para reemplazar los valores nulos por datos que aporten al data-set se utiliza el metodo
# RandomForest
forestImp <- missForest(df)

# Valores imputados
df.forestImp <- forestImp$ximp
df.forestImp$bareNuclei <- round(df.forestImp$bareNuclei)

# Se guardan los nuevos datos
df<- df.forestImp
df$class<-df.inicial$class

# Se procede a realizar un grafico de cajas por cada variable para identificar los posibles
# ouliers que puedan existir. Para esto se conviernten los datos a tipo long para ser graficados
datos.long<-melt(df)

# Se crea un grafico de cajas
bp <- ggboxplot(datos.long,
                x = "variable", 
                y = "value",
                fill = "variable")

# Se muestra el grafico
print(bp)

# Se Cambian los outliers encontrados por la media de la correspondiente variable
df$marginalAdhesion[df$marginalAdhesion > 8] <- mean(df$marginalAdhesion)
df$singleEpithCellSize[df$singleEpithCellSize > 7] <- mean(df$singleEpithCellSize)
df$blandChromatin[df$blandChromatin > 9] <- mean(df$blandChromatin)
df$normalNucleoli[df$normalNucleoli > 8] <- mean(df$normalNucleoli)
df$mitoses[df$mitoses > 1] <- mean(df$mitoses)

#y regraficamos
# Se conviernten los datos a tipo long para ser graficados
datos.long<-melt(df)

# Se crea un grafico de cajas
bp <- ggboxplot(datos.long,
                x = "variable", 
                y = "value",
                fill = "variable")

print(bp)


# II. CONJUNTOS DE ENTRENAMIENTO Y PRUEBA
# 
# Para la construccion del arbol, se generan un conjunto de entrenamiento y un conjuno de prueba

# Antes verificaremos la relacion en la cual esta distribuida la clase a lo largo del conjunto completo.
set.seed(20)  #SENS=0.96,ACCURACY=0.9522
#set.seed(30)
#set.seed(2) #SENS=0.9781,ACCURACY=0.9713
#set.seed(284) #SENS=0.9922                    ,ACCURACY=.9474
table(df$class)

training.index = createDataPartition(df$class, p = 0.7)$Resample1
df.training = df[training.index, ]
df.test = df[-training.index, ]

# Ahora verificaremos si la praticion se ha realizado correctamente, tanto para el conjunto de prueba, 
# como para el conjunto de entrenamiento
prop.table(table(df.training$class))
prop.table(table(df.test$class))




#________________________________________________________________________________________________________#
# III. AJUSTE DEL MODELO

# Para este caso se utiliza el algoritmo C5.0 para la generacion del arbol, el cual es aplicado para cada
# una de las variables, excepto para la variable explicativa que en este caso corresponde a la variable
# o columna "class" (indice 10 en el data frame)

tree = C5.0(class ~ ., df.training)
tree.rules = C5.0(x = df.training[, -10], y = df.training$class, rules = T)
tree.pred.class = predict(tree, df.test, type = "class")
tree.pred.prob = predict(tree, df.test, type = "prob")


tree.pred.class

#________________________________________________________________________________________________________#
# III. EVALUACION DEL MODELO

conf.matrix.tree = confusionMatrix(table(df.test$class, tree.pred.class))
print(conf.matrix.tree)

head(tree.pred.prob)

#________________________________________________________________________________________________________#
# III. EVALUACION DEL MODELO

plot(tree)
summary(tree)

#________________________________________________________________________________________________________#
# III. BOOSTING

tree_b = C5.0(class ~ ., df.training, trials = 5)
tree_b.rules = C5.0(x = df.training[, -10], y = df.training$class, rules = T, trials = 5)
tree_b.pred.class = predict(tree_b, df.test, type = "class")
tree_b.pred.prob = predict(tree_b, df.test, type = "prob")

conf.matrix.tree.b = confusionMatrix(table(df.test$class, tree_b.pred.class))
print(conf.matrix.tree.b)
