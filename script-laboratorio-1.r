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
                 desviaciones),
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

# Por otro lado, se utiliza la función hist() para crear un histograma de los
# datos de cada grupo. Esto se usa como primera aproximación para
# observar si existe algún tipo de asimetría o si los datos se comportan
# de forma normal.
par(mfrow=c(3,3))
ign<-mapply(hist,
            pruebas,
            main=colnames(pruebas),
            col="lightsteelblue",
            xlab="Puntaje")

# Al graficar el histograma se puede observar que los grupos de datos no parecen
# seguir distribuciones normales, y la mayoría de las variables pareciera tener 
# un alto grado de asimetría negativa, es decir, los datos se 
# concentran a la derecha. Sin embargo, para estar seguros utilizaremos un
# Shapiro-Test, cuyas hipotesis a contrastar son las siguientes:

# H0: La muestra proviene de una población cuya distribución es normal
# HA: La muestra proviene de una población cuya distribución no es normal

# Se define un alfa = 0.05 para contrastar con el resultado obtenido

# Se aplica el test a todas las variables
shapirotest <- apply(pruebas,2,shapiro.test)

# Se muestra el resultado del test
print(shapirotest)

# Se puede ver que en cada una de las columnas el p-valor es menor a 0.05, 
# por lo que se rechaza la hipotesis nula en todas las distribuciones. Esto
# indica que las variables no siguen una distribución normal. Por lo que no
# se puede aplicar en test ANOVA y se deberá optar por un test no paramétrico.




######################################################################################################
##########################       Análisis por Clases        ##########################################

# Primero se suma la cantidad de observaciones que pertenecen a benigno y las que pertenecen a maligno
B = sum(df$class == "Benigno",na.rm=TRUE)
M = sum(df$class == "Maligno",na.rm=TRUE)

# Se calculan los porcentajes de cada clase
Tot = B+M
PB= B/Tot
PM= M/Tot

# Se muestran los resultados
cat("Totales:",Tot,"\nBenignos:",B,"Malignos:",M,"\nPorcentaje Benigno:",PB,"Porcentaje Maligno:",PM)

# Se crea un gráfico para visualizar la distribución por clases.
valores <- c(PB, PM)
etiquetas <- c("Benignos","Malignos")
etiquetas <- paste(etiquetas,round(valores,2))
etiquetas <- paste(etiquetas,"%",sep="")

# Se muestra el gráfico
par(mfrow=c(1,1))
pie(valores, labels = etiquetas,col=rainbow(2), main="% Casos Benignos vs Malignos")

# Se crea una función para comparar cada variable con la clase de cancer
# a la cual pertence la observación
grafico <- function(df,aes,y_string){
                    boxplt =  ggboxplot(data = df, 
                                        x = "class", 
                                        y = y_string, 
                                        color = "class",) + border()
                    
                    ydens = axis_canvas(boxplt, 
                                        axis = "y", 
                                        coord_flip = TRUE) + 
                                        geom_density(data = df, 
                                        aes, alpha = 0.7, 
                                        size = 0.2) + 
                                        coord_flip()
  
                    boxplt = insert_yaxis_grob(boxplt, 
                                               ydens, 
                                               grid::unit(.2, "null"), 
                                               position = "right")
                    ggdraw(boxplt)

                    }# Fin de la función

# Se llama a la función para crear los gráficos
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




######################################################################################################
##########################           Correlación            ##########################################

#Por ultimo, como se mencionó que las variables no siguen una distribución normal se decide
# utiliza el test no paramétrico de Spearman para analizar la correlación de las variables

# Se aplica el test de correlación
df.cor = cor(pruebas, method = 'spearman')

# Se grafican los resultados
ggcorrplot(df.cor,hc.order = TRUE,type = "full",lab = TRUE)



#######################################################################################################
## Nota: en caso de tener error al mostrar los gráficos se debe expandir o agrandar la ventana       ##
##       en donde se muestran los gráficos. En caso de que el grafico este en blanco hacer lo mismo. ##
#######################################################################################################

# Fin del script