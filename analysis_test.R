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

setwd("P/Scripts")

load("../data/output/beta_QN_filtered_cgName.RData")

sframe <- read.csv("../sample_with_pheno_sexOnly_relevantColumns.csv")



#make naming scheme of sframe match beta2

sframe$Sentrix_full <- paste(sframe$Sentrix_Barcode, sframe$Sentrix_Position, sep = "_")



#check if there are unmatched elements

matches <- sframe$Sentrix_full %in% colnames(beta2)
unmatching_values <- colnames(beta2)[!matches]
unmatching_values <- which(!matches)
unmatching_values



#add the columns that match

sframe <- sframe %>%
  mutate(cpg1 = beta2["cg01086462", match(sframe$Sentrix_full, colnames(beta2))])

fwrite(sframe, "./datasets/sframe.csv")

sframe <- fread("./datasets/sframe.csv")




### 10 & more cpg sites ###

#randomize beta2

set.seed(123)

beta2 <- beta2[sample(nrow(beta2)),]
beta2


for (i in seq(from = 1, to = 100, by = 1)) {
  new_name <- rownames(beta2)[i]
  
  sframe <- sframe %>%
    mutate(!!new_name := beta2[new_name, match(sframe$Sentrix_full, colnames(beta2))])
}

sframe$cpg1 <- NULL

fwrite(sframe, "./datasets/sframe10to20.csv")
fwrite(sframe, "./datasets/sframe100.csv")


colnames(sframe)[colnames(sframe) == "cg01086462"] <- "cpg1"
colnames(sframe)[colnames(sframe) == "cpg1"] <- "cg01086462"

colnames(sframe)[colnames(sframe) == "adjusted"] <- "cpg1"
colnames(sframe)[colnames(sframe) == "cpg1"] <- "adjusted"







#test LM

colors <- c("lightcoral", "lightgreen")

group_index <- as.integer(factor(sframe$Sex))


plot(group_index, sframe$cpg1, 
       main="Title", 
       xlab="X", 
       ylab="Y",
       pch=19, 
       col=colors[group_index])

groups <- sort(unique(sframe$Sex))




#general summary

sframe$Sentrix_Barcode <- as.factor(sframe$Sentrix_Barcode)
reg <- lm(cpg1 ~ Sex, data = sframe)
summary(reg)




#emmeans alternative

emmeans(reg, "Sex")

regBatch <- lm(cpg1 ~ Sex + Sentrix_Barcode, data = sframe)
regBatch0 <- lm(cpg1 ~ Sentrix_Barcode, data = sframe)
plot(emmeans(reg, "Sex"))
plot(emmeans(regBatch, "Sex"))
plot(emmeans(regBatch, "Sentrix_Barcode"))
plot(emmeans(regBatch, "Sentrix_Barcode", "Sex"))
plot(emmeans(regBatch0, "Sentrix_Barcode"))

boxplot(cpg1 ~ Sex + Sentrix_Barcode, sframe)
boxplot(cpg1 ~ Sentrix_Barcode + Sex, sframe) 

boxplot(cpg1 ~ Sex, sframe)
boxplot(cpg1 ~ Sentrix_Barcode, sframe)



#with ggplot2

#make summary for errorbar plot

sf_summary <- sframe %>%
  group_by(Sex) %>%
  summarise(mean = mean(cpg1, na.rm = TRUE), sd = sd(cpg1, na.rm = TRUE), .groups = 'drop')


#regular plot

ggplot(sframe, aes(x = Sex, y = cpg1, color = Sex, group = Sex)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("female" = "lightcoral", "male" = "lightgreen", "NA" = "black")
  )


#errorbar plot

ggplot(sf_summary, aes(x = Sex, y = mean, color = Sex)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  scale_color_manual(values = c("female" = "lightcoral", "male" = "lightgreen", "NA" = "black")
  )


#boxplot & beeswarm plot

ggplot(sframe, aes(x = Sex, y = cpg1, color = Sex)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_boxplot() +
  geom_beeswarm() +
  scale_color_manual(values = c("female" = "lightcoral", "male" = "lightgreen", "NA" = "black")
  )
  



#other summary
  
x <- sframe$cpg1
summary_stats <- data.frame(
  Mean = mean(x, na.rm = TRUE),
  SD   = sd(x, na.rm = TRUE),
  Min  = min(x, na.rm = TRUE),
  Max  = max(x, na.rm = TRUE)
)

summary_stats














### combat ###


#get data into correct format for combat

col_test <- sframe$Sentrix_full
row_test <- colnames(sframe)[8:17]

testdf <- data.frame()
#testdf <- cbind(testdf, setNames(as.data.frame(matrix(NA, nrow = nrow(testdf), ncol = length(col_test))), col_test))
testdf[row_test, col_test] <- NA


#turn correct sentrix codes into index for sframe and assign values to testdf where they match

rownames(sframe) <- sframe$Sentrix_full
sframe <- as.data.frame(sframe)
testdf[row_test, col_test] <- t(sframe[col_test, row_test])



#test combat

?ComBat

combattest <- ComBat(testdf, sframe$Sentrix_Barcode, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
combattest <- cbind(Cpg_Site = rownames(combattest), combattest)
fwrite(combattest, "./datasets/combat_results100.csv")


#prepare to plot combat

cb <- fread("./datasets/batch_correction_results/combat_results.csv")
sframe <- fread("./datasets/sframe100.csv")

cb <- as.data.frame(cb)
sframe <- as.data.frame(sframe)


#copy cb into sframe

rownames(cb) <- cb[, 1]
sframe[, row_test] <- t(cb[row_test, col_test])



ggplot(sframe, aes(x = Sex, y = cg01086462, color = Sex)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_boxplot() +
  geom_beeswarm() +
  scale_color_manual(values = c("female" = "lightcoral", "male" = "#00BFC4", "NA" = "black")
  )




#new variables

sm <- read_sav("P_Methylated.sav")

smoke <- sm[, "Nikotin_tgl_FrühSS"]
smoke[, "Nikotinkonsum"] <- sm[, "Nikotinkonsum"]
smoke[, "Vpn"] <- sm[, "Probandennummer"]

smoke <- smoke %>% 
  mutate(Nikotin_tgl_FrühSS = if_else(Nikotin_tgl_FrühSS <= 2, 2, 1))

smoke <- smoke %>% 
  mutate(smoking_status = if_else(Nikotin_tgl_FrühSS != 1 | Nikotinkonsum != 1, 2, 1))

smoke <- smoke %>% unnest(Vpn)


sframe <- sframe %>% 
  left_join(smoke %>% select(Vpn, smoking_status), 
            by = c("vpn" = "Vpn"))


fwrite(sframe, "./datasets/sframeBatchSmoke.csv")
sframe <- fread("./datasets/sframeSmokeCombat.csv")
  

#emmeans plot

sframe$Sentrix_Barcode <- as.factor(sframe$Sentrix_Barcode)

regSmoke <- lm(cpg1 ~ Sentrix_Barcode + Sex + smoking_status, data = sframe)
regSmokeI <- lm(cpg1 ~ Sentrix_Barcode + Sex * smoking_status, data = sframe)
regSmokeI2 <- lm(cpg1 ~ Sentrix_Barcode + Sex * Gestational_age, data = sframe)

plot(emmeans(regSmoke, "smoking_status"))
plot(emmeans(regSmoke, "smoking_status", "Sex"))
plot(emmeans(regSmoke, "Sentrix_Barcode", "smoking_status"))
plot(emmeans(regSmoke, ~Sentrix_Barcode + smoking_status))
plot(emmeans(regSmoke, ~smoking_status + Sentrix_Barcode))
plot(emmeans(regSmoke, ~smoking_status + Sex))
plot(emmeans(regSmoke, ~smoking_status * Sex))

plot(emmeans(regSmokeI, "smoking_status"))
plot(emmeans(regSmokeI, "smoking_status", "Sex"))
plot(emmeans(regSmokeI, "Sentrix_Barcode", "smoking_status"))
plot(emmeans(regSmokeI, ~smoking_status + Sentrix_Barcode))
plot(emmeans(regSmokeI, ~smoking_status * Sentrix_Barcode))

plot(emmeans(regSmokeI2, ~Sex + Gestational_age))
plot(emmeans(regSmokeI2, ~Sex * Gestational_age))



#beeswarm plot

sframe <- sframe %>%
  mutate(smoking_status = recode(smoking_status,
                        `1` = "non-smoker",
                        `2` = "smoker"))

ggplot(sframe, aes(x = Sex, y = cpg1, color = smoking_status)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_boxplot() +
  geom_beeswarm() +
  scale_color_manual(values = c("non-smoker" = "darksalmon", "smoker" = "darkslategrey")
  )









### limma ###

sframe <- fread("./datasets/sframeSmoke.csv")

expression_matrix <- t(sframe[, 8:17])
colnames(expression_matrix) <- unlist(sframe[,"Sentrix_full"])


batch_corr <- removeBatchEffect(x = expression_matrix, batch = sframe$Sentrix_Barcode, design = model.matrix(~sframe$smoking_status))
batch_corr <- data.frame(Cpg_site = rownames(batch_corr), batch_corr, check.names = FALSE)
fwrite(batch_corr, "./datasets/limma_resultsSmoke.csv")

sframe <- as.data.frame(sframe)
batch_corr$Cpg_site <- NULL
sframe[, 8:17] <- t(batch_corr)


fwrite(sframe, "./datasets/sframe10to20SmokeLimmaSmoke.csv")
sframe <- fread("./datasets/sframeSmokeLimmaSmoke.csv")























