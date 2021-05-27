######################################################################################################
###########################    Laboratorio 1 - Análisis de Datos    ##################################
######################################################################################################
###########################    Autores:     Gary Simken             ##################################
###########################                 Sebastián Orellana      ##################################
###########################    Fecha:       27 - Mayo - 2021        ##################################
######################################################################################################

######################################################################################################
##########################            Bibliotecas          ###########################################
library(modeest)
library(ggpubr)
library(cowplot)
library(corrplot)
library(ggcorrplot)
library(reshape2)
library(RColorBrewer)

# Se eliminan los gráficos y variables antiguas
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




######################################################################################################
##########################        Limpieza de Datos        ###########################################

# Se eliminan los datos nulos del data-set
df = df[!(df$bareNuclei == "?"),] 

# Se convierten los datos a tipo numérico
df[ , ] = apply(df[ , ], 2,function(x) as.numeric(as.character(x)) )

# Se elimina la columna id al no ser relevante para el análisis
df$id<-NULL

# Se cambian los valores de las clases por un String para una mejor visualización
df$class <- as.character(df$class)
df$class[df$class == "2"] <- "Benigno"
df$class[df$class == "4"] <- "Maligno"

# A continuación se procede a analizar posibles outliers en los datos
pruebas = df
pruebas$id <- NULL
pruebas$class<- NULL

# Se conviernten los datos a tipo long para ser graficados
datos.long<-melt(pruebas)

# Se crea un gráfico de cajas
bp1 <- ggboxplot(datos.long,
                 x = "variable", 
                 y = "value",
                 fill = "variable")

# Se muestra el gráfico
print(bp1)

# Se Cambian los outliers encontrados por la media de la correspondiente variable
df$marginalAdhesion[df$marginalAdhesion > 5] <- mean(df$marginalAdhesion)
df$singleEpithCellSize[df$singleEpithCellSize > 4] <- mean(df$singleEpithCellSize)
df$blandChromatin[df$blandChromatin > 4] <- mean(df$blandChromatin)
df$normalNucleoli[df$normalNucleoli > 5] <- mean(df$normalNucleoli)
df$mitoses[df$mitoses > 2] <- mean(df$mitoses)

# Se elimina el id y la clase para seguir analizando las variables por separado
pruebas = df
pruebas$id <- NULL
pruebas$class<- NULL

# Se convierten a tipo long los datos sin outliers
datos.long2<-melt(pruebas)

# Se crea un gráfico cuartil-cuartil para observar si se eliminaron correctamente
# los valores atípicos
qq2 <- ggqqplot(datos.long2,
                x = "value",
                color = "variable")
qq2 <- qq2 + facet_wrap(~ variable)

# Se muestra el gráfico
print(qq2)

# En adición, se crea un grafico de cajas nuevamente para observar los datos filtrados
bp2 <- ggboxplot(datos.long2,
                 x = "variable", y = "value",
                 fill = "variable")

# Se muestra el gráfico
print(bp2)




######################################################################################################
##########################     Medidas de Centralización   ###########################################

# Una vez que se realizó la limpieza de datos se proceden a calcular las siguientes métricas:

# Se calcula la media para todas las variables
medias = apply(pruebas, 2, mean)

# Se calcula la mediana para todas las variables
medianas =  apply(pruebas, 2, median)

# Se calcula la moda para todas la variables
modas = apply(pruebas, 2, mlv,  method = "mfv")

# se calcula la desviación estandar para todas las variables
desviaciones = apply(pruebas, 2, sd)

# Se calcula la varianza para todas las variables
varianzas = apply(pruebas, 2, var)

# Se calcula el rango de cada una de las variables
rangos = apply(pruebas, 2, range)

# Se crea una tabla para mostrar los resultados
tabla = matrix(c(medias,
                 medianas,
                 modas,
                 varianzas,
                 desviaciones, 
                 rangos),
               ncol=9,
               byrow=TRUE)

# Se le asignan los nombres de las variables a las columnas
colnames(tabla) =  c("clumpThickness",
                     "uniformityCellSize",
                     "uniformityCellShape",
                     "marginalAdhesion",
                     "singleEpithCellSize",
                     "bareNuclei",
                     "blandChromatin",
                     "normalNucleoli",
                     "mitoses")

# Se le asignan los nombres de las métricas a las filas
rownames(tabla) = c("media","mediana","moda","varianza","SD")

# Se convierte la variable en una tabla
tabla=as.table(tabla)

# Se muestran los resultados
print(tabla)




######################################################################################################
##########################       Distribución de Datos      ##########################################

#Se verifica si las variables siguen una distribuciÃ³n normal mediante shapiro test y un alfa estandar de 0.05
#hipotesis nula
#Los datos siguen una distribucion normal
#hipostesis alternativa
#Los datos no siguen una distribucion normal

distribucion <- apply(pruebas,2,shapiro.test)
par(mfrow=c(3,3))
#en cada uno de las columnas el p-valor es menor a 0.05, por lo que rechazamos en todas
#Para verificar se revisara el histograma de cada columna (En caso de error, agrandar la seccion de graficos)
ign<-mapply(hist,pruebas,main=colnames(pruebas),col="lightsteelblue",xlab="Puntaje")
#Efectivamente, Ninguna sigue una distribuiciÃ³n normal, por lo que se opta por un test no paramÃ©trico

#se realiza una prueba de spearman
df.cor = cor(pruebas, method = 'spearman')
#Para luego realizar un grafico de calor con la matriz de correlacion
ggcorrplot(df.cor,hc.order = TRUE,type = "full",lab = TRUE)


#Aqui estarÃ¡ la cantidad de datos que se trabajarÃ¡
B = sum(df$class == "Benigno",na.rm=TRUE)
M = sum(df$class == "Maligno",na.rm=TRUE)
Tot = B+M
PB= B/Tot
PM= M/Tot
cat("Totales:",Tot,"\nBenignos:",B,"Malignos:",M,"\nPorcentaje Benigno:",PB,"Porcentaje Maligno:",PM)
#GRAFICO CIRCULAR para revisar la distribucion de estos
valores <- c(PB, PM)
etiquetas <- c("Benignos","Malignos")
etiquetas <- paste(etiquetas,round(valores,2))
etiquetas <- paste(etiquetas,"%",sep="")
pie(valores, ,labels = etiquetas,col=rainbow(2), main="% Casos Benignos vs Malignos")


#Se define una funciÃ³n para comparar cada variable con la clase de cancer que se tiene (benigno o maligno)
grafico <- function(df,aes,y_string){
  boxplt =  ggboxplot(data = df, x = "class", y = y_string, color = "class",) + border()

  ydens = axis_canvas(boxplt, axis = "y", coord_flip = TRUE) + geom_density(data = df, aes, alpha = 0.7, size = 0.2) + coord_flip()

  boxplt = insert_yaxis_grob(boxplt, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(boxplt)
}
#y se llama para cada columna
a<-grafico(df,aes(x = clumpThickness, fill = class),"clumpThickness")
b<-grafico(df,aes(x = uniformityCellSize, fill = class),"uniformityCellSize")
c<-grafico(df,aes(x = uniformityCellShape, fill = class),"uniformityCellShape")
d<-grafico(df,aes(x = marginalAdhesion, fill = class),"marginalAdhesion")
e<-grafico(df,aes(x = singleEpithCellSize, fill = class),"singleEpithCellSize")
f<-grafico(df,aes(x = bareNuclei, fill = class),"bareNuclei")
g<-grafico(df,aes(x = blandChromatin, fill = class),"blandChromatin")
h<-grafico(df,aes(x = normalNucleoli, fill = class),"normalNucleoli")
i<-grafico(df,aes(x = mitoses, fill = class),"mitoses")
ggarrange(a,b,c,d)
ggarrange(e,f,g,h)
i


















