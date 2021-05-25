#main v2
library(modeest)
library(ggpubr)
library(cowplot)
library(corrplot)
library(ggcorrplot)
library(reshape2)
library(RColorBrewer)

rm(list=ls())
if(length(dev.list())!=0){
  dev.off(dev.list()["RStudioGD"])
}

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

#Se definen los nombres de las columnas, estos son los mismos de provistos por la base de datos
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
#mediante el url se importan los datos
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data"
df = read.csv(url, header = F, sep=",", col.names = columnas)


#Se sacan los datos nulos del datagrama
df=df[!(df$bareNuclei == "?"),] #SN= sin nulos
#Se conbvierten todos a numericos
df[ , ] = apply(df[ , ], 2,function(x) as.numeric(as.character(x)) )

#Para nalizar Outalier, primero se graficaran los datos obtenidos
#Se elimina la columna id
df$id<-NULL
#Se convierte los enteros en clases para un mayor
df$class <- as.character(df$class)
df$class[df$class == "2"] <- "Benigno"
df$class[df$class == "4"] <- "Maligno"

# #Se define una función para comparar cada variable con la clase de cancer que se tiene (benigno o maligno)
# grafico <- function(df,aes,y_string){
#   boxplt =  ggboxplot(data = df, x = "class", y = y_string, color = "class",) + border()
# 
#   ydens = axis_canvas(boxplt, axis = "y", coord_flip = TRUE) + geom_density(data = df, aes, alpha = 0.7, size = 0.2) + coord_flip()
# 
#   boxplt = insert_yaxis_grob(boxplt, ydens, grid::unit(.2, "null"), position = "right")
#   ggdraw(boxplt)
# }
# #y se llama para cada columna
# a<-grafico(df,aes(x = clumpThickness, fill = class),"clumpThickness")
# b<-grafico(df,aes(x = uniformityCellSize, fill = class),"uniformityCellSize")
# c<-grafico(df,aes(x = uniformityCellShape, fill = class),"uniformityCellShape")
# d<-grafico(df,aes(x = marginalAdhesion, fill = class),"marginalAdhesion")
# e<-grafico(df,aes(x = singleEpithCellSize, fill = class),"singleEpithCellSize")
# f<-grafico(df,aes(x = bareNuclei, fill = class),"bareNuclei")
# g<-grafico(df,aes(x = blandChromatin, fill = class),"blandChromatin")
# h<-grafico(df,aes(x = normalNucleoli, fill = class),"normalNucleoli")
# i<-grafico(df,aes(x = mitoses, fill = class),"mitoses")
# ggarrange(a,b,c,d)
# ggarrange(e,f,g,h)
# i
# 
# 
# #Al revisar minuciosamente se llega al concenso de que hay que filtrar
# #ciertos datos
# 
# df1<-df[!(((df$uniformityCellSize>5 | df$uniformityCellShape>5 | df$marginalAdhesion>7 | df$singleEpithCellSize>7 |
#               df$bareNuclei>6 | df$blandChromatin>6 | df$normalNucleoli>7 | df$mitoses>6)
#            & df$class =="Benigno")  | (df$mitoses>6 & df$class =="Maligno")),]
# 
# df1<-df1[!is.na(df1$clumpThickness),]
# 
# # y se vuelven a graficar
# a2<-grafico(df1,aes(x = clumpThickness, fill = class),"clumpThickness")
# b2<-grafico(df1,aes(x = uniformityCellSize, fill = class),"uniformityCellSize")
# c2<-grafico(df1,aes(x = uniformityCellShape, fill = class),"uniformityCellShape")
# d2<-grafico(df1,aes(x = marginalAdhesion, fill = class),"marginalAdhesion")
# e2<-grafico(df1,aes(x = singleEpithCellSize, fill = class),"singleEpithCellSize")
# f2<-grafico(df1,aes(x = bareNuclei, fill = class),"bareNuclei")
# g2<-grafico(df1,aes(x = blandChromatin, fill = class),"blandChromatin")
# h2<-grafico(df1,aes(x = normalNucleoli, fill = class),"normalNucleoli")
# i2<-grafico(df1,aes(x = mitoses, fill = class),"mitoses")
# ggarrange(a2,b2,c2,d2)
# ggarrange(e2,f2,g2,h2)
# i2

pruebas = df
pruebas$id <- NULL
pruebas$class<- NULL

datos.long<-melt(pruebas)
p2 <- ggqqplot(
  datos.long,
  x = "value",
  color = "variable"
)
p2 <- p2 + facet_wrap(~ variable)
print(p2)

p1 <- ggboxplot(
  datos.long,
  x = "variable", y = "value",
  fill = "variable"
)
print(p1)
df<-df[!(df$marginalAdhesion>8 | df$singleEpithCellSize>7 |
          df$blandChromatin>9 | df$normalNucleoli>8 | df$mitoses>8),]

df<-df[!is.na(df1$clumpThickness),]

pruebas = df
pruebas$id <- NULL
pruebas$class<- NULL

datos.long2<-melt(pruebas)
p2 <- ggqqplot(
  datos.long2,
  x = "value",
  color = "variable"
)
p2 <- p2 + facet_wrap(~ variable)
print(p2)

p1 <- ggboxplot(
  datos.long2,
  x = "variable", y = "value",
  fill = "variable"
)
print(p1)


distribucion <- apply(pruebas,2,shapiro.test)
par(mfrow=c(3,3))
#en cada uno de las columnas el p-valor es menor a 0.05, por lo que rechazamos en todas
#Para verificar se revisara el histograma de cada columna (En caso de error, agrandar la seccion de graficos)
ign<-mapply(hist,pruebas,main=colnames(pruebas),col="lightsteelblue",xlab="Puntaje")


# #Aqui estará la cantidad de datos que se trabajará
# B = sum(df1$class == "Benigno",na.rm=TRUE)
# M = sum(df1$class == "Maligno",na.rm=TRUE)
# Tot = B+M
# PB= B/Tot
# PM= M/Tot
# cat("Totales:",Tot,"\nBenignos:",B,"Malignos:",M,"\nPorcentaje Benigno:",PB,"Porcentaje Maligno:",PM)

# # GRAFICO CIRCULAR para revisar la distribucion de estos
# valores <- c(PB, PM)
# etiquetas <- c("Benignos","Malignos")
# etiquetas <- paste(etiquetas,round(valores,2))
# etiquetas <- paste(etiquetas,"%",sep="")
# pie(valores, ,labels = etiquetas,col=rainbow(2), main="% Casos Benignos vs Malignos")

#Para calcular las medidas de centralizaciony dispersion se sacara id y class

medias = apply(pruebas, 2, mean)
medianas =  apply(pruebas, 2, median)
modas = apply(pruebas, 2, mlv,  method = "mfv")
desviaciones = apply(pruebas, 2, sd)
varianzas = apply(pruebas, 2, var)
rangos = apply(pruebas, 2, range)
#Se crea una tabla para mostrar todo
tabla = matrix(c(medias,medianas,modas,varianzas,desviaciones),ncol=9,byrow=TRUE)
#Se le asigna nombre a las columnas
colnames(tabla) =  c("clumpThickness",
                     "uniformityCellSize",
                     "uniformityCellShape",
                     "marginalAdhesion",
                     "singleEpithCellSize",
                     "bareNuclei",
                     "blandChromatin",
                     "normalNucleoli",
                     "mitoses"
)
#y a las filas
rownames(tabla) = c("media",
                    "mediana",
                    "moda",
                    "varianza",
                    "SD"
)
#Finalmente se convierte a tabla y se muestra
tabla=as.table(tabla)
tabla



#Se verifica si las variables siguen una distribución normal mediante shapiro test y un alfa estandar de 0.05
#hipotesis nula
#Los datos siguen una distribucion normal
#hipostesis alternativa
#Los datos no siguen una distribucion normal

distribucion <- apply(pruebas,2,shapiro.test)
par(mfrow=c(3,3))
#en cada uno de las columnas el p-valor es menor a 0.05, por lo que rechazamos en todas
#Para verificar se revisara el histograma de cada columna (En caso de error, agrandar la seccion de graficos)
ign<-mapply(hist,pruebas,main=colnames(pruebas),col="lightsteelblue",xlab="Puntaje")


#Efectivamente, Ninguna sigue una distribuición normal, por lo que se opta por un test no paramétrico


df.cor = cor(pruebas, method = 'spearman')
#Para luego realizar un grafico de calor con la matriz de correlacion
ggcorrplot(df.cor, 
           hc.order = TRUE, 
           type = "full",
           lab = TRUE)





