######################################################################################################
###########################    Laboratorio 1 - Analisis de Datos    ##################################
######################################################################################################
###########################    Autores:     Gary Simken             ##################################
###########################                 Sebasti?n Orellana      ##################################
###########################    Fecha:       27 - Mayo - 2021        ##################################
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

# Se eliminan los graficos y variables antiguas
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



set.seed(1234)
######################################################################################################
##########################        Limpieza de Datos        ###########################################
#Se eliminar??n las variables de clase e id, dado que son datos no relevantes en esta parte del estudio,
#en el caso de id es un simple identificador por lo que no aporta y en el caso de 
#clase, esta es eliminada porque el Clustering es un m??todo no supervisado.
df <- subset(df, select = -c(id))
df.inicial<-df
df <- subset(df, select = -c(class))

#Para un vistaso previo aplicaremos un summary
summary(df)

#Para los datos "missing" se proceder?? a cambiar estos por un NA que es el valor default en los enteros cuando no
#existe un dato
df$bareNuclei[df$bareNuclei == "?"] = NA

#Por otro lado convertiremos todos los datos a numericos para poder trabajarlos como tal
df[ , ] = apply(df[ , ], 2,function(x) as.numeric(as.character(x)) )

#Para poder tratar estos datos se utilizara el metodo basado en Random Forest 

forestImp <- missForest(df)
summary(forestImp)

# Valores imputados
df.forestImp <- forestImp$ximp
df.forestImp$bareNuclei <- round(df.forestImp$bareNuclei)

#Ahora con todos los valores 
df<- df.forestImp

#para revisar los Valores outliers graficaremos en caja los datos
# Se conviernten los datos a tipo long para ser graficados
datos.long<-melt(df)

# Se crea un gr?fico de cajas
bp <- ggboxplot(datos.long,
                x = "variable", 
                y = "value",
                fill = "variable")

# Se muestra el gr?fico
print(bp)

# Se Cambian los outliers encontrados por la media de la correspondiente variable
df$marginalAdhesion[df$marginalAdhesion > 5] <- mean(df$marginalAdhesion)
df$singleEpithCellSize[df$singleEpithCellSize > 4] <- mean(df$singleEpithCellSize)
df$blandChromatin[df$blandChromatin > 4] <- mean(df$blandChromatin)
df$normalNucleoli[df$normalNucleoli > 5] <- mean(df$normalNucleoli)
df$mitoses[df$mitoses > 2] <- mean(df$mitoses)

# Se conviernten los datos a tipo long para ser graficados
datos.long<-melt(df)

# Se crea un gr?fico de cajas
bp <- ggboxplot(datos.long,
                x = "variable", 
                y = "value",
                fill = "variable")

# Se muestra el gr?fico
print(bp)

#Con esto se concluye la trata de datos, deshaciendonos de outliers y valores nulos


######################################################################################################
##########################       Distribuci?n de Datos      ##########################################

# Por otro lado, se utiliza la funci?n hist() para crear un histograma de los
# datos de cada grupo. Esto se usa como primera aproximaci?n para
# observar si existe alg?n tipo de asimetr?a o si los datos se comportan
# de forma normal.
par(mfrow=c(3,3))
ign<-mapply(hist,
            df,
            main=colnames(df),
            col="lightsteelblue",
            xlab="Puntaje")

# Al graficar el histograma se puede observar que los grupos de datos no parecen
# seguir distribuciones normales, y la mayor?a de las variables pareciera tener 
# un alto grado de asimetr?a negativa, es decir, los datos se 
# concentran a la derecha. Sin embargo, para estar seguros utilizaremos un
# Shapiro-Test, cuyas hipotesis a contrastar son las siguientes:

# H0: La muestra proviene de una poblaci?n cuya distribuci?n es normal
# HA: La muestra proviene de una poblaci?n cuya distribuci?n no es normal

# Se define un alfa = 0.05 para contrastar con el resultado obtenido

# Se aplica el test a todas las variables
shapirotest <- apply(df,2,shapiro.test)

# Se muestra el resultado del test
print(shapirotest)

# Se puede ver que en cada una de las columnas el p-valor es menor a 0.05, 
# por lo que se rechaza la hipotesis nula en todas las distribuciones. Esto
# indica que las variables no siguen una distribuci?n normal. Por lo que no
# se puede aplicar en test ANOVA y se deber? optar por un test no param?trico.

#############################################################################
########################Normalizacion

#para realizar clustering, es necesario escalar los datos. Esto debido a que se est?? ocupando una misma "regla" en las mediciones.
df.scale=scale(df)

#Para realizar el clustering buscaremos las distancias entre los datos
dist.eucl = dist(df.scale, method = "euclidean")
a<-fviz_dist(dist.eucl)
print(a)
b<-fviz_nbclust(df.scale, kmeans, method = "silhouette")+labs(subtitle = "Metodo de siluetas")
print(b)
b<-fviz_nbclust(df.scale, kmeans, method = "wss")+labs(subtitle = "Metodo de codo")+geom_vline(xintercept = 2, linetype = 2)
print(b)
#Podemos ver el codo marcado en x =2 por lo que colocamos una linea en el para notarlo

b<-fviz_nbclust(df.scale, kmeans, method = "gap_stat",iter.max=50)+labs(subtitle = "gap stadistict method")
print(b)

#Tenemos que 2 de los 3 metodos nos arrojan un k optimo = 2 por lo que podemos decir que este sera nuestro k
#finalmente podemos graficarlos
km.res = kmeans(df.scale, 2, nstart = 25)
cluster<-fviz_cluster(km.res, data = df.scale, palette = "jco", ggtheme = theme_minimal())
print(cluster)


#Segun los resultados obtenidos obtenidos obtenemos que el k obtimo en este caso es 2
