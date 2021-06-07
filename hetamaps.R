## HEATMAP IN R

setwd("~/4º BIOTECNOLOGIA/SEGUNDO CUATRI/PRACTICAS/SSP final")
library(gplots)

# Base function: heatmap()
## scale () es una función que centra y escala las columnas de un 
## matriz numérica. Transponemos la matriz con t () para luego centrar 
## y escale los datos de cada proteína (es decir, las filas) con scale (). 
## Finalmente, transponemos los datos a la forma original.

## datos PcBMM
FPKM_PcBMM <- read.csv("FPKM_PcBMM2.txt", header=FALSE, sep = "\t")
colnames(FPKM_PcBMM) <- c("invitro",	"Col0-10h",	"Col0-16h",	"Col0-24h",	"cyp-10h",	"cyp-16h")
rownames(FPKM_PcBMM) <- c("AIM000156/SSP1.1", "AIM006935/SSP1.1", "AIM003066/SSP2", "AIM001916/SSP3", "AIM004866/SSP4", "AIM002092/SSP4", "AIM008201/SSP5", "AIM004328/SSP5", "AIM004457/SSP6", "AIM008364/SSP6.2", "AIM006883/SSP6.2", "AIM009621/SSP7", "AIM011005/SSP7", "AIM006534/SSP8.1", "AIM007454/SSP8.1", "AIM002535/SSP8.1", "AIM003826/SSP8.3", "AIM007909/SSP8.3", "AIM001556/SSP8.3", "AIM006016/SSP8.3", "AIM003702/SSP8.3", "AIM001769/SSP8.3", "AIM003089/SSP8.3", "AIM006326/SSP9", "AIM003697/SSP9", "AIM008393/SSP10", "AIM008337/SSP10", "AIM009859/SSP10", "AIM008031/SSP11", "AIM003562/SSP11", "AIM006878/SSP11", "AIM010172/SSP11", "AIM003204/SSP12", "AIM004170/SSP12", "AIM006317/SSP14.2")
FPKM_PcBMM
# scale escala por columnas, yo quiero filas
data_PcBMM <- data.matrix(FPKM_PcBMM)
data_scaled_t_PcBMM <- scale(t(data_PcBMM)) # traspuesta y escalar por columans
data_scaled_PcBMM <- t(data_scaled_t_PcBMM) # deshago la trasposición

## DATOS Pc2127
FPKM_Pc2127 <- read.csv("FPKM_Pc2127.txt", header=FALSE, sep="\t")
colnames(FPKM_Pc2127) <- c("invitro",	"Col0-10h",	"Col0-16h",	"Col0-24h",	"cyp-10h",	"cyp-16h")
rownames(FPKM_Pc2127) <- c("AIM010795/SSP1.1", "AIM001150/SSP1.1", "AIM002806/SSP2", "AIM005985/SSP3", "AIM000735/SSP4", "AIM009771/SSP4", "AIM001939/SSP4", "AIM006101/SSP5", "AIM004951/SSP5", "AIM004613/SSP6", "AIM005081/SSP6.2", "AIM007279/SSP6.2", "AIM006980/SSP7", "AIM000922/SSP7", "AIM003574/SSP8.1", "AIM004894/SSP8.1", "AIM001564/SSP8.1", "AIM007727/SSP8.3", "AIM004417/SSP8.3", "AIM003223/SSP8.3", "AIM005450/SSP8.3", "AIM006554/SSP8.3", "AIM001120/SSP8.3", "AIM002336/SSP8.3", "AIM002271/SSP8.3", "AIM006218/SSP9", "AIM006549/SSP9", "AIM005111/SSP10", "AIM007274/SSP11", "AIM006671/SSP11", "AIM007581/SSP12", "AIM009834/SSP12", "AIM003699/SSP14.1", "AIM006209/SSP14.2")
FPKM_Pc2127

data_Pc2127 <- data.matrix(FPKM_Pc2127)
data_scaled_t_Pc2127 <- scale(t(data_Pc2127))
data_scaled_Pc2127 <- t(data_scaled_t_Pc2127)

## DATOS P0831
FPKM_P0831 <- read.csv("FPKM_P0831.txt", header=FALSE, sep="\t")
colnames(FPKM_P0831) <- c("invitro",	"Col0-10h",	"Col0-16h",	"Col0-24h",	"cyp-10h",	"cyp-16h")
rownames(FPKM_P0831) <- c("AIM010024/SSP1.1", "AIM009437/SSP1.1", "AIM006788/SSP1.2", "AIM006318/SSP2", "AIM008473/SSP3", "AIM006276/SSP4", "AIM007779/SSP4", "AIM006996/SSP5", "AIM006134/SSP5", "AIM002818/SSP6.1", "AIM004146/SSP6", "AIM007721/SSP7", "AIM001611/SSP7", "AIM007235/SSP7", "AIM004258/SSP8.1", "AIM002198/SSP8.1", "AIM004053/SSP8.3", "AIM010767/SSP8.3", "AIM001570/SSP8.3", "AIM009039/SSP8.3", "AIM008659/SSP8.3", "AIM005061/SSP9", "AIM005463/SSP10", "AIM000571/SSP10", "AIM002142/SSP11", "AIM009821/SSP14", "AIM001836/SSP14.2")
FPKM_P0831

data_P0831 <- data.matrix(FPKM_P0831)
data_scaled_t_P0831 <- scale(t(data_P0831))
data_scaled_P0831 <- t(data_scaled_t_P0831)

# Color palette

my_palette <- colorRampPalette(c("blue",
                                 "white",
                                 "red"))

## DENDOGRAMAS CON TIPOS DE METODOS DE CLUSTERING JERARQUICO

library(dplyr)
matriz_distancias <- dist(x = data_scaled_PcBMM[, c("invitro", "Col0-10h",	"Col0-16h",	"Col0-24h",	"cyp-10h",	"cyp-16h")], method = "euclidean")
set.seed(567)
hc_euclidea_completo <- hclust(d = matriz_distancias, method = "complete")
hc_euclidea_single   <- hclust(d = matriz_distancias, method = "single")
hc_euclidea_average  <- hclust(d = matriz_distancias, method = "average")
hc_euclidea_ward  <- hclust(d = matriz_distancias, method = "ward.D")

par(mfrow = c(2,2))
plot(x = hc_euclidea_completo, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage complete")
plot(x = hc_euclidea_single, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage single")
plot(x = hc_euclidea_average, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage average")
plot(x = hc_euclidea_ward, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage ward")

# Todos los dendogramas me sacan 4 clusteres claramente, y el que
# me reparte los genes de manera más homogénea es ward.D, por lo que 
# es el que elijo, aunque me valdría cualquiera.


## DIAGRAMA DE CODO PARA NUMERO DE CLUSTERS

diag_codoplot <- function(data, nc, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

diag_codoplot(data_scaled_PcBMM, nc=10) 

# Según el diagrama de codo, donde se ha estudiado la suma 
# de cuadrados dentro de cada grupo (within groups sum of squares)
# se establece que el número de clusteres óptimo es 4.

# El siguiente paso es establecer el método de clustering


## Sacamos que el mejor método de clustering es WARD.D

# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(RColorBrewer)
my_group <- as.numeric(as.factor(substr(rownames(data_scaled_PcBMM), 1 , 1)))
colSide <- brewer.pal(9, "Set1")[my_group]
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
Heatmap(data_scaled_PcBMM, 
        name = "PcBMM", #title of legend
        column_title = "PcBMM Gene expression", row_title = "Gene name",
        row_names_gp = gpar(fontsize = 7),
        clustering_method_row = "ward.D",
        cluster_rows = TRUE,
        column_order = c("invitro", "Col0-10h", "Col0-16h", "Col0-24h", "cyp-10h", "cyp-16h"),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:6),
                labels = c("1", "2", "3", "4"), 
                labels_gp = gpar(col = "white", fontsize = 10))),
        row_km = 4
        )

# ?Heatmap

Heatmap(data_scaled_Pc2127,  #title of legend
        column_title = "Pc2127",
        row_names_gp = gpar(fontsize = 7),# Text size for row names
        clustering_method_row = "ward.D",
        cluster_rows = TRUE,
        column_order = c("invitro", "Col0-10h", "Col0-16h", "Col0-24h", "cyp-10h", "cyp-16h"),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:6),
                  labels = c("4", "3", "2", "1"),
                  labels_gp = gpar(col = "white", fontsize = 10)),
                  width = unit(0.4, "cm")),
        row_km = 4,
        width = unit(6, "cm"), height = unit(8, "cm"),
        column_names_side = "top"
)

Heatmap(data_scaled_P0831, 
        name = "P0831", #title of legend
        column_title = "P0831 Gene expression", row_title = "Gene name",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        clustering_method_row = "ward.D",
        cluster_rows = TRUE,
        column_order = c("invitro", "Col0-10h", "Col0-16h", "Col0-24h", "cyp-10h", "cyp-16h"),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:6),
                labels = c("1", "2", "3", "4"), 
                labels_gp = gpar(col = "white", fontsize = 10))),
        row_km = 4
)

## HEATMAP COMUN

setwd("~/4º BIOTECNOLOGIA/SEGUNDO CUATRI/TFG")
comun <- read.csv("rel.txt", header = FALSE, sep = "\t")
comun
colnames(comun) <- c("PcBMM-Col0-10h",	"PcBMM-Col0-16h", "PcBMM-Col0-24h", "PcBMM-cyp-10h", "PcBMM-cyp-16h", "Pc2127-Col0-10h", "Pc2127-Col0-16h", "Pc2127-Col0-24h", "Pc2127-cyp-10h", "Pc2127-cyp-16h", "P0831-Col0-10h", "P0831-Col0-16h", "P0831-Col0-24h", "P0831-cyp-10h", "P0831-cyp-16h")
rownames(comun) <- c("PcBMM_AIM000156", "PcBMM_AIM003066", "PcBMM_AIM001916", "PcBMM_AIM004866", "PcBMM_AIM002092", "PcBMM_AIM004328", "PcBMM_AIM009621", "PcBMM_AIM011005", "PcBMM_AIM006534", "PcBMM_AIM003697", "PcBMM_AIM008393", "PcBMM_AIM006317")
comun
lista_genes <- c("PcBMM_AIM000156, Pc2127_AIM010795, P0831_AIM010024", "PcBMM_AIM003066, Pc2127_AIM002806, P0831_AIM006318", "PcBMM_AIM001916, Pc2127_AIM005985, P0831_AIM008473", "PcBMM_AIM004866, Pc2127_AIM000735, P0831_AIM006276", "PcBMM_AIM002092, Pc2127_AIM009771, P0831_AIM007779", "PcBMM_AIM004328, Pc2127_AIM004951, P0831_AIM006134", "PcBMM_AIM009621, Pc2127_AIM006980, P0831_AIM007721", "PcBMM_AIM011005, Pc2127_AIM000922, P0831_AIM001611", "PcBMM_AIM006534, Pc2127_AIM003574, P0831_AIM004258", "PcBMM_AIM003697, Pc2127_AIM006549, P0831_AIM005061", "PcBMM_AIM008393, Pc2127_AIM005111, P0831_AIM005463", "PcBMM_AIM006317, Pc2127_AIM006209, P0831_AIM001836")

comun_data <- data.matrix(comun)
# escalado por filas de manera independiente
data_scaled_comun_t <- scale(t(comun_data)) # traspuesta y escalar por columans
data_scaled_comun <- t(data_scaled_comun_t)
data_scaled_comun
# escalado todas las filas de la misma manera
data_scaled_total_t <- scale(t(comun_data), scale = c(8, 8, 8, 8,8, 8,8, 8,8, 8,8, 8)) # traspuesta y escalar por columans
data_scaled_total <- t(data_scaled_total_t)
data_scaled_total

matriz_distancias2 <- dist(x = data_scaled_comun[, c(c("PcBMM-Col0-10h",	"PcBMM-Col0-16h", "PcBMM-Col0-24h", "PcBMM-cyp-10h", "PcBMM-cyp-16h", "Pc2127-Col0-10h", "Pc2127-Col0-16h", "Pc2127-Col0-24h", "Pc2127-cyp-10h", "Pc2127-cyp-16h", "P0831-Col0-10h", "P0831-Col0-16h", "P0831-Col0-24h", "P0831-cyp-10h", "P0831-cyp-16h"))], method = "euclidean")
set.seed(567)
hc_euclidea_completo2 <- hclust(d = matriz_distancias2, method = "complete")
hc_euclidea_single2 <- hclust(d = matriz_distancias2, method = "single")
hc_euclidea_average2 <- hclust(d = matriz_distancias2, method = "average")
hc_euclidea_ward2 <- hclust(d = matriz_distancias2, method = "ward.D")

par(mfrow = c(2,2))
plot(x = hc_euclidea_completo2, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage complete")
plot(x = hc_euclidea_single2, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage single")
plot(x = hc_euclidea_average2, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage average")
plot(x = hc_euclidea_ward2, cex = 0.6, xlab = "", ylab = "", sub = "",
     main = "Distancia euclídea, Linkage ward")
# complete o ward??

## VALIDO NUMERO DE CLUSTERS
diag_codoplot(data_scaled_comun, nc=10) 


Heatmap(data_scaled_comun, #title of legend
        row_names_gp = gpar(fontsize = 7),
        clustering_method_row = "ward.D",
        cluster_rows = TRUE,
        column_names_gp = gpar(fontsize = 9),
        column_order = c("PcBMM-Col0-10h",	"PcBMM-Col0-16h", "PcBMM-Col0-24h", "PcBMM-cyp-10h", "PcBMM-cyp-16h", "Pc2127-Col0-10h", "Pc2127-Col0-16h", "Pc2127-Col0-24h", "Pc2127-cyp-10h", "Pc2127-cyp-16h", "P0831-Col0-10h", "P0831-Col0-16h", "P0831-Col0-24h", "P0831-cyp-10h", "P0831-cyp-16h"),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:7),
                                                          labels = c("5", "4", "3", "2", "1"),
                                                          labels_gp = gpar(col = "white", fontsize = 10)),
                                         width = unit(0.4, "cm")),
        row_km = 5,
        width = unit(6, "cm"), height = unit(8, "cm"),
        column_names_side = "top"
)


