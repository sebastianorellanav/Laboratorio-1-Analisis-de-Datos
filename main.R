
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

url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data"
df = read.csv(url, header = F, sep=",", col.names = columnas)

df_SN=df[!(df$bareNuclei == "?"),] #SN= sin nulos

df_SN[ , ] = apply(df_SN[ , ], 2,function(x) as.numeric(as.character(x)) )

df_sin_id_class = df_SN[ , 2:10]  
medias = apply(df_sin_id_class, 2, mean)
medianas =  apply(df_sin_id_class, 2, median)
modas = apply(df_sin_id_class, 2, mlv,  method = "mfv")
desviaciones = apply(df_sin_id_class, 2, sd)
varianzas = apply(df_sin_id_class, 2, var)
rangos = apply(df_sin_id_class, 2, range)

tabla = matrix(c(medias,medianas,modas,varianzas,desviaciones),ncol=9,byrow=TRUE)

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

rownames(tabla) = c("media",
                    "mediana",
                    "moda",
                    "varianza",
                    "SD"
)
tabla=as.table(tabla)
tabla

df_SN$id <- NULL
#Se verifica si las variables siguen una distribución normal
distribucion <- apply(df_SN,2,shapiro.test)
apply(df_SN,2,hist)
#Ninguna sigue una distribuición normal, se opta por un test no paramétrico
df.cor = cor(df_SN, method = 'spearman')

ggcorrplot(df.cor, 
           hc.order = TRUE, 
           type = "full",
           lab = TRUE)

#desvuacuib
#rango
#varianza


df_SN$class <- as.character(df_SN$class)
df_SN$class[df_SN$class == "2"] <- "Benigno"
df_SN$class[df_SN$class == "4"] <- "Maligno"


# GRAFICOS2 # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

# Comparación Clump Thickness #
boxplot.clumpThickness =  ggboxplot(data = df_SN, x = "class", y = "clumpThickness", color = "class") + border()

ydens = axis_canvas(boxplot.clumpThickness, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = clumpThickness, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.clumpThickness = insert_yaxis_grob(boxplot.clumpThickness, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.clumpThickness)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Uniformity of Cell Size # 
boxplot.unifCellSize =  ggboxplot(data = df_SN, x = "class", y = "uniformityCellSize", color = "class") + border()

ydens = axis_canvas(boxplot.unifCellSize, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = uniformityCellSize, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.unifCellSize = insert_yaxis_grob(boxplot.unifCellSize, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.unifCellSize)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Uniformity of Cell Shape #
boxplot.unifCellShape =  ggboxplot(data = df_SN, x = "class", y = "uniformityCellShape", color = "class") + border()

ydens = axis_canvas(boxplot.unifCellShape, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = uniformityCellShape, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.unifCellShape = insert_yaxis_grob(boxplot.unifCellShape, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.unifCellShape)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Marginal Adhesion #
boxplot.marginalAdhesion =  ggboxplot(data = df_SN, x = "class", y = "marginalAdhesion", color = "class") + border()

ydens = axis_canvas(boxplot.marginalAdhesion, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = marginalAdhesion, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.marginalAdhesion = insert_yaxis_grob(boxplot.marginalAdhesion, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.marginalAdhesion)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Single Epithelial Cell Size #
boxplot.epithCellSize =  ggboxplot(data = df_SN, x = "class", y = "singleEpithCellSize", color = "class") + border()

ydens = axis_canvas(boxplot.epithCellSize, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = singleEpithCellSize, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.epithCellSize = insert_yaxis_grob(boxplot.epithCellSize, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.epithCellSize)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Bare Nuclei #
boxplot.bareNuclei =  ggboxplot(data = df_SN, x = "class", y = "bareNuclei", color = "class") + border()

ydens = axis_canvas(boxplot.bareNuclei, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = bareNuclei, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.bareNuclei = insert_yaxis_grob(boxplot.bareNuclei, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.bareNuclei)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Bland Chromatin #
boxplot.blandChromatin =  ggboxplot(data = df_SN, x = "class", y = "blandChromatin", color = "class") + border()

ydens = axis_canvas(boxplot.blandChromatin, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = blandChromatin, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.blandChromatin = insert_yaxis_grob(boxplot.blandChromatin, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.blandChromatin)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Normal Nucleoli #
boxplot.normalNucleoli =  ggboxplot(data = df_SN, x = "class", y = "normalNucleoli", color = "class") + border()

ydens = axis_canvas(boxplot.normalNucleoli, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = normalNucleoli, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.normalNucleoli = insert_yaxis_grob(boxplot.normalNucleoli, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.normalNucleoli)

# ----------------------------------------------------------------------------------------------------------------------------------- #
# Comparacion Mitosis #
boxplot.mitoses =  ggboxplot(data = df_SN, x = "class", y = "mitoses", color = "class") + border()

ydens = axis_canvas(boxplot.mitoses, axis = "y", coord_flip = TRUE) + geom_density(data = df_SN, aes(x = mitoses, fill = class), alpha = 0.7, size = 0.2) + coord_flip()

boxplot.mitoses = insert_yaxis_grob(boxplot.mitoses, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(boxplot.mitoses)
