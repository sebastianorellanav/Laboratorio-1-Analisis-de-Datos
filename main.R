
library(modeest)
library(ggpubr)
library(cowplot)
library(corrplot)
library(ggcorrplot)
library(RColorBrewer)
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
df_SN=df[!(df$bareNuclei == "?"),] #SN= sin nulos
#Se conbvierten todos a numericos
df_SN[ , ] = apply(df_SN[ , ], 2,function(x) as.numeric(as.character(x)) )

#Se filtran las columnas id y class para calcular, medidas de centralización y dispersión
df_sin_id_class = df_SN[ , 2:10]  
medias = apply(df_sin_id_class, 2, mean)
medianas =  apply(df_sin_id_class, 2, median)
modas = apply(df_sin_id_class, 2, mlv,  method = "mfv")
desviaciones = apply(df_sin_id_class, 2, sd)
varianzas = apply(df_sin_id_class, 2, var)
rangos = apply(df_sin_id_class, 2, range)
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

#Eliminamos el campo id, ya que es irrelevante para aplicar las pruebas de normalidad
df_SN$id <- NULL
#Se verifica si las variables siguen una distribución normal mediante shapiro test y un alfa estandar de 0.05
#hipotesis nula
#Los datos siguen una distribucion normal
#hipostesis alternativa
#Los datos no siguen una distribucion normal

distribucion <- apply(df_SN,2,shapiro.test)
#en cada uno de las columnas el p-valor es menor a 0.05, por lo que rechazamos en todas
#Para verificar se revisara el histograma de cada columna
apply(df_SN,2,hist)
#Efectivamente, Ninguna sigue una distribuición normal, por lo que se opta por un test no paramétrico
df.cor = cor(df_SN, method = 'spearman')
#Para luego realizar un grafico de calor con la matriz de correlacion
ggcorrplot(df.cor, 
           hc.order = TRUE, 
           type = "full",
           lab = TRUE)


#Se convierte los enteros en clases para un mayor
df_SN$class <- as.character(df_SN$class)
df_SN$class[df_SN$class == "2"] <- "Benigno"
df_SN$class[df_SN$class == "4"] <- "Maligno"


#Se define una función para comparar cada variable con la clase de cancer que se tiene (benigno o maligno)
grafico <- function(df,aes,y_string){
  boxplt =  ggboxplot(data = df, x = "class", y = y_string, color = "class") + border()
  
  ydens = axis_canvas(boxplt, axis = "y", coord_flip = TRUE) + geom_density(data = df, aes, alpha = 0.7, size = 0.2) + coord_flip()
  
  boxplt = insert_yaxis_grob(boxplt, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(boxplt)  
}
#y se llama para cada columna
grafico(df_SN,aes(x = clumpThickness, fill = class),"clumpThickness")
grafico(df_SN,aes(x = uniformityCellSize, fill = class),"uniformityCellSize")
grafico(df_SN,aes(x = uniformityCellShape, fill = class),"uniformityCellShape")
grafico(df_SN,aes(x = marginalAdhesion, fill = class),"marginalAdhesion")
grafico(df_SN,aes(x = singleEpithCellSize, fill = class),"singleEpithCellSize")
grafico(df_SN,aes(x = bareNuclei, fill = class),"bareNuclei")
grafico(df_SN,aes(x = blandChromatin, fill = class),"blandChromatin")
grafico(df_SN,aes(x = normalNucleoli, fill = class),"normalNucleoli")
grafico(df_SN,aes(x = mitoses, fill = class),"mitoses")





