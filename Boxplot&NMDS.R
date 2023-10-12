setwd("C:/Users/hanna/Documents/Universitetet/PhD/GitHub/Cirsium/")

library(tidyverse)
library(readxl)
library(reshape2)
library(data.table)
library(janitor)

##Kalle run this chunk:
Total <- read.csv("Control_adjusted_samples.csv")
Total$Total <- as.numeric(Total$Total)
prop <- read.csv("cirsium_proportions.csv")

###Total scent###

shapiro.test(Total$Total) #W = 0.52813, p-value < 2.2e-16

kruskal.test(only_floral$Total~only_floral$Pop) #Kruskal-Wallis chi-squared = 139.47, df = 34, p-value = 9.989e-15

boxplot(Total$Total~Total$Sample.type, xlab="Tissue", ylab= "Relative abundance", main = "Relative emission per tissue", col = c(1:length(as.factor(Total$Sample.type))))

#Logtransform total scent

totallog <- Total
totallog$Total <- log(totallog$Total)

barplot(totallog$Total)

boxplot(totallog$Total~totallog$Sample.type)

#Divide into species

het <- subset(Total, Species=="HET")
ole <- subset(Total, Species=="OLE")

boxplot(het$Total~het$Sample.type, xlab="Tissue", ylab= "Relative abundance", main = "Relative emission per tissue C. heterophyllum", col = c(1:length(as.factor(het$Sample.type))))
boxplot(ole$Total~ole$Sample.type, xlab="Tissue", ylab= "Relative abundance", main = "Relative emission per tissue C. heterophyllum", col = c(1:length(as.factor(ole$Sample.type))))

#Test emissions of tissues within species

boxplot(het$Total~het$Sample.type, xlab="Tissue", ylab= "Relative abundance", main = "Relative emission per tissue C. heterophyllum", col = c(1:length(as.factor(het$Sample.type))))

boxplot(ole$Total~ole$Sample.type, xlab="Tissue", ylab= "Relative abundance", main = "Relative emission per tissue C. oleraceum", col = c(1:length(as.factor(ole$Sample.type))))


###Compounds### Kalle run this
library(vegan)
library(RColorBrewer)
library(reshape2)
library(pals)


# Find smallest non-zero variable ; Kalle run this
scentValues <- as.vector(as.matrix(prop[,12:58]))
smallestNonZero <- min(scentValues[scentValues > 0])
smallestNonZero

# Adding this value as a dummy variable to every sample; Kalle run this
prop2 <- prop
prop2$DummyVariable <- smallestNonZero

# PERMANOVA; Kalle run this
adonis2(prop2[,c(12:58,59)] ~ Sample.type, data = prop, permutations = 9999, method = "bray") 
# p = 1e-04***



# NMDS; Kalle run this
NMDS3 <- metaMDS(prop2[,c(12:58,59)], distance = "bray", k =2, try = 10, trymax = 200, autotransform = FALSE)
NMDS3 # Stress:0.1468733 


stressplot(NMDS3)

# Plot for both species and tissue type; Kalle run this
dev.off()

layout(rbind(1,2), heights=c(7,1))

ordiplot(NMDS3, display = "sites", type = "n",  main = "Tissue types")
points(NMDS3, display = "sites", pch = c(21:23)[as.factor(prop2$Sample.type)], cex=2,  bg =c("#8400D8", "#11811D")[as.factor(prop$Species)])

par(mar=c(0, 0, 0, 0))

legend("bottomleft", legend = levels(as.factor(prop2$Sample.type)), pch =  c(21:23), col ="black",
       ncol=1,bty ="n")
#legend( legend = levels(as.factor(prop2$Species)), pch =  c(5,19), cex=2, col = 1)

#col =c("#8400D8", "#11811D")[as.factor(prop$Species)],
 ##heterophyllum##

hetprop <- subset(prop, Species=="HET")

# Find smallest non-zero variable
scentValues <- as.vector(as.matrix(hetprop[,12:58]))
smallestNonZero <- min(scentValues[scentValues > 0])
smallestNonZero

# Adding this value as a dummy variable to every sample
hetprop2 <- hetprop
hetprop2$DummyVariable <- smallestNonZero

# PERMANOVA
adonis2(hetprop2[,c(12:58,59)] ~ Sample.type, data = hetprop, permutations = 9999, method = "bray") 
# p = 0.015 *

# NMDS
NMDS3 <- metaMDS(hetprop2[,c(12:58,59)], distance = "bray", k =3, try = 10, trymax = 200, autotransform = FALSE)
NMDS3 # Stress:0.03560533


stressplot(NMDS3)

# Plot for both species and tissue type
ordiplot(NMDS3, display = "sites", type = "n",  main = "Tissue types C. heterophyllum")
points(NMDS3, display = "sites", pch = 19, cex=1.5, col = as.factor(hetprop$Sample.type))
legend("bottomleft", legend = levels(as.factor(hetprop$Sample.type)), pch = 19, col =c(1:3))

##oleraceum##

oleprop <- subset(prop, Species=="OLE")

# Find smallest non-zero variable
scentValues <- as.vector(as.matrix(oleprop[,12:58]))
smallestNonZero <- min(scentValues[scentValues > 0])
smallestNonZero

# Adding this value as a dummy variable to every sample
oleprop2 <- oleprop
oleprop2$DummyVariable <- smallestNonZero

# PERMANOVA
adonis2(oleprop2[,c(12:58,59)] ~ Sample.type, data = oleprop, permutations = 9999, method = "bray") 
# p = 0.0035 **

# NMDS
NMDS3 <- metaMDS(oleprop2[,c(12:58,59)], distance = "bray", k =3, try = 10, trymax = 200, autotransform = FALSE)
NMDS3 # Stress:9.092996e-05


stressplot(NMDS3)

# Plot for both species and tissue type
ordiplot(NMDS3, display = "sites", type = "n",  main = "Tissue types C. oleraceum")
points(NMDS3, display = "sites", pch = 19, cex=1.5, col = as.factor(oleprop$Sample.type))
legend("bottomleft", legend = levels(as.factor(oleprop$Sample.type)), pch = 19, col =c(1:3))

