library(data.table)
library(corrplot)
library(preputils)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(emmeans)
library(sva)
library(haven)
library(tidyr)
library(limma)
library(tidyverse)
library(pcaMethods)
library(rgl)





### PCA & separation by plate color ###

#perform pca with prcomp()
#results will be saved as a dataframe where $x seems to hold the values of the PCA 

#consider reversing the signs of the $rotation and $x sections
#plot the results

sframe <- fread("data/sframeSmoke.csv")
sframe <- fread("data/sframeSmokeCombat.csv")
sframe <- fread("data/sframeSmokeLimmaSmoke.csv")

sframe <- fread("data/sframe10to20Smoke.csv")
sframe <- fread("data/sframe10to20SmokeCombat.csv")
sframe <- fread("data/sframe10to20SmokeLimmaSmoke.csv")

sframe <- fread("data/sframe100Smoke.csv")
sframe <- fread("data/sframe100SmokeCombat.csv")
sframe <- fread("data/sframe100SmokeLimmaSmoke.csv")




#move smoking to earlier position

sframe <- sframe %>% relocate(smoking_status, .before = colnames(sframe[, 8]))



pca_test <- sframe[, 9:117]
pca_test <- pca_test[complete.cases(pca_test), ]
results <- prcomp(pca_test) #no scaling --> input data is already in percentages; scaling only sets variance between 1 & 0

results$rotation <- -1*results$rotation
results$rotation

results$x <- -1*results$x
head(results$x)

sframe <- sframe[complete.cases(sframe), ]
pca_scores <- as.data.frame(results$x)
pca_scores$batch <- sframe$Sentrix_Barcode
pca_scores$sex <- sframe$Sex
pca_scores$smoking_status <- sframe$smoking_status

pca_scores$batch <- as.factor(pca_scores$batch)
pca_scores$sex <- as.factor(pca_scores$sex)
pca_scores$smoking_status <- as.factor(pca_scores$smoking_status)

pairs(pca_scores[, 1:5],
      main = "Scatterplot Matrix",
      pch = 19,          
      col = pca_scores$sex,
)


#boxplot pc1

boxplot(PC2 ~ sex + batch, pca_scores)

boxplot(PC1 ~ batch + sex, pca_scores)

plot(results) #check which principal components are the most important by plotting





?prcomp
?pcaMethods
?plot3d



#in case of driver issues

options(rgl.printRglwidget = TRUE)

with(pca_scores, plot3d(PC1, PC2, PC3, type="s", size=3, col=as.numeric(factor(sframe$smoking_status))))
with(pca_scores, plot3d(PC1, PC2, PC3, type='p', size=10))

with(pca_scores, plot3d(PC1, PC2, PC3, type="s", size=scale(sframe$Gestational_age), col=as.numeric(factor(sframe$Sex))))
with(pca_scores, plot3d(PC1, PC2, PC3, type="s", col=as.numeric(factor(sframe$Sentrix_Position)), size=as.numeric(factor(sframe$Sex))))

with(pca_scores, plot3d(PC1, PC2, PC3, type="s", size=sframe$smoking_status, col=as.numeric(factor(sframe$Sex))))





#changing types (as test)

with(pca_scores, plot3d(PC1[sframe$smoking_status == 1], PC2[sframe$smoking_status == 1], PC3[sframe$smoking_status == 1], type="s", size=as.numeric(factor(pca_scores$sex[sframe$smoking_status == 1])), col=as.numeric(factor(pca_scores$sex[sframe$smoking_status == 1]))))
with(pca_scores, plot3d(PC1[sframe$smoking_status == 2], PC2[sframe$smoking_status == 2], PC3[sframe$smoking_status == 2], type="l", col=as.numeric(factor(pca_scores$sex[sframe$smoking_status == 2])), add=TRUE))

#again with gestational age

with(pca_scores, plot3d(PC1[sframe$smoking_status == 1], PC2[sframe$smoking_status == 1], PC3[sframe$smoking_status == 1], type="s", size=scale(sframe$Gestational_age[sframe$smoking_status == 1]), col=as.numeric(factor(pca_scores$sex[sframe$smoking_status == 1]))))
with(pca_scores, plot3d(PC1[sframe$smoking_status == 2], PC2[sframe$smoking_status == 2], PC3[sframe$smoking_status == 2], type="l", col=as.numeric(factor(pca_scores$sex[sframe$smoking_status == 2])), add=TRUE))


#alternative

with(pca_scores1, plot3d(PC1, PC2, PC3, type="s", size=as.numeric(pca_scores1$sex), col=as.numeric(pca_scores1$sex)))
with(pca_scores2, plot3d(PC1, PC2, PC3, type="l", col=as.numeric(pca_scores2$sex), add=TRUE))














