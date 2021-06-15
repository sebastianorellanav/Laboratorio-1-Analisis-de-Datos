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
library(NbClust)

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
df <- subset(df, select = -c(class))

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
df$marginalAdhesion[df$marginalAdhesion > 5] <- mean(df$marginalAdhesion)
df$singleEpithCellSize[df$singleEpithCellSize > 4] <- mean(df$singleEpithCellSize)
df$blandChromatin[df$blandChromatin > 4] <- mean(df$blandChromatin)
df$normalNucleoli[df$normalNucleoli > 5] <- mean(df$normalNucleoli)
df$mitoses[df$mitoses > 2] <- mean(df$mitoses)

# Se conviernten los datos a tipo long para ser graficados
datos.long<-melt(df)

# Se crea un grafico de cajas
bp <- ggboxplot(datos.long,
                x = "variable", 
                y = "value",
                fill = "variable")

# Se muestra el grafico
print(bp)

# Con esto se concluye la limpieza y preprocesamiento de los datos, deshaciendonos de outliers y
# valores nulos





########################################################################################################
####################################       Estandarizacisn         ####################################

# Antes de realizar el proceso de clustering, es necesario escalar los datos. Esto debido a que se 
# ocupa una misma "regla" para las mediciones.
df.scale=scale(df)

# Comparacion entre los datos antes y despuis de realizar el escalamiento
# Se comvierten los datos a tipo long para graficarlos
datos.long=melt(df)

# Se crea un grafico cuartil-cuartil
qq <- ggqqplot(datos.long,
               x = "value",
               color = "variable") + facet_wrap(~ variable)

# Se muestra el grafico
print(qq)

# Se comvierten los datos escalados a tipo long para graficarlos
datos.long=melt(df.scale)

# Se crea un grafico cuartil cuartil con los datos escalados
qq <- ggqqplot(datos.long,
               x = "value",
               color = "Var2") + facet_wrap(~ Var2)

# Se muestra el grafico
print(qq)




########################################################################################################
#################################       Medidas de Distancia         ##################################

#Para realizar el clustering se analizan las distancias entre los datos
#Se calculan distintas distancias para compararlas

# Distancia Euclmdea
dist.eucl <- daisy(df.scale, 
                   metric = "euclidean", 
                   stand = FALSE)

# Distancia de Gower
dist.gower <- daisy(df.scale, 
                    metric = "gower", 
                    stand = FALSE)

# Distancia de Manhattan
dist.manhattan <- daisy(df.scale, 
                        metric = "manhattan", 
                        stand = FALSE)

#Se grafican las distancias
plot.eucl<-fviz_dist(dist.eucl,
                     gradient = list(low = "#00AFBB", 
                                     mid = "white", 
                                     high = "#FC4E07"))
print(plot.eucl)

plot.gower<-fviz_dist(dist.gower,
                      gradient = list(low = "#00AFBB", 
                                      mid = "white", 
                                      high = "#FC4E07"))
print(plot.gower)

plot.manhattan<-fviz_dist(dist.manhattan,
                          gradient = list(low = "#00AFBB", 
                                          mid = "white", 
                                          high = "#FC4E07"))
print(plot.manhattan)

#dist.matrix=as.matrix(dist.eucl)




########################################################################################################
#################################       Buscando el K Optimo          ##################################

# Metodo de la silueta
silhouette <- fviz_nbclust(df.scale, 
                           kmeans, 
                           method = "silhouette") + labs(subtitle = "Metodo de siluetas")

# Mostrar grafico
print(silhouette)

# Metodo del Codo
wss <- fviz_nbclust(df.scale, 
                    kmeans, 
                    method = "wss") +labs(subtitle = "Metodo de codo")+geom_vline(xintercept = 2, linetype = 2)

# Mostrar grafico
print(wss)

# Metodo de la brecha estadmstica
gap_stat <- fviz_nbclust(df.scale, 
                         kmeans, 
                         method = "gap_stat",
                         iter.max=20) + labs(subtitle = "gap stadistict method")

# Mostrar Grafico
print(gap_stat)

# Se utiliza la funcisn nbClust que utiliza diferentes mitodos para analizar el nzmero sptimo de clusters
res.nbclust <- NbClust(df.scale, 
                       distance = "euclidean",
                       min.nc = 2, 
                       max.nc = 9, 
                       method = "complete", 
                       index ="all")

# Se crea un grafico de barras
fviz_nbclust(res.nbclust) + theme_minimal() + ggtitle("Numero Optimo de Clusters (distancia euclidea)")

# Se realiza el mismo procedimiento pero con la distancia de manhattan
res.nbclust <- NbClust(df.scale, 
                       distance = "manhattan",
                       min.nc = 2, 
                       max.nc = 9, 
                       method = "complete", 
                       index ="all")

# Se crea un grafico de barras
fviz_nbclust(res.nbclust) + theme_minimal() + ggtitle("Numero Optimo de Clusters (distancia manhattan)")

#Tenemos que 2 de los 3 metodos nos arrojan un k optimo = 2 por lo que podemos decir que este sera nuestro k

# Se realiza el procedimiento de k-means para k=2
km.res <- kmeans(df.scale, 
                2, 
                nstart = 25)

# Se grafica el cluster
cluster_k2 <- fviz_cluster(km.res, 
                         data = df.scale, 
                         palette = "jco", 
                         ggtheme = theme_minimal(), 
                         main = "K = 2")

# Se muestra el grafico
print(cluster_k2)


# Se realiza el procedimiento de k-means para k=3
km.res <- kmeans(df.scale, 
                3, 
                nstart = 25)

# Se crea el grafico de los clusters
cluster_k3 <- fviz_cluster(km.res, 
                        data = df.scale, 
                        palette = "jco", 
                        ggtheme = theme_minimal(), 
                        main = "K = 3")

# Se muestra el grafico
print(cluster_k3)


# Se realiza el procedimiento de k-means para k=4
km.res <- kmeans(df.scale, 
                4, 
                nstart = 25)

# Se crea el grafico
cluster_k4 <- fviz_cluster(km.res, 
                         data = df.scale, 
                         palette = "jco", 
                         ggtheme = theme_minimal(), 
                         main = "K = 4")

# Se muestra el grafico
print(cluster_k4)


# Se realiza el procedimiento de k-means para k=6
km.res <- kmeans(df.scale, 
                6, 
                nstart = 25)

# Se crea el grafico
cluster_k6 <- fviz_cluster(km.res, 
                         data = df.scale, 
                         palette = "jco", 
                         ggtheme = theme_minimal(), 
                         main = "K = 6")

# Se muestra el grafico
print(cluster_k6)

# Se crea una figura con todos los graficos
ggarrange(cluster_k2, 
          cluster_k3, 
          cluster_k4, 
          cluster_k6, 
          nrow = 2, 
          ncol = 2)




clase.cluster<-as.numeric(cluster_k2$data$cluster)
clase.cluster<-replace(clase.cluster,clase.cluster==1,4)
clase<-df.inicial$class

malignos<-sum(clase==4)
malignos.c<-sum(clase.cluster==4)

benignos<-sum(clase==2)
benignos.c<-sum(clase.cluster==2)


df2 <- data.frame(Prediccion=rep(c("Original", "Cluster"), each=2),
                  Clase=rep(c("Benigno", "Maligno"),2),
                  Cantidad=c(benignos, malignos, benignos.c, malignos.c))


grafico.barra<-ggplot(data=df2, aes(x=Clase, y=Cantidad, fill=Prediccion)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=Cantidad), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()





comparacion<- clase==clase.cluster
acertacion=sum(comparacion)*100/length(clase)
a<-paste("Acertacion ",as.character(round(acertacion,2)),"%",sep="")
b<-paste("Fallo ",as.character(100-round(acertacion,2)),"%",sep="")

grafico.acertacion <- pie(c(acertacion,100-acertacion), labels = c(a,b),col=rainbow(2), main="% Casos Benignos vs Malignos")

ggarrange(grafico.barra, 
          grafico.acertacion,
          nrow = 1, 
          ncol = 2)
