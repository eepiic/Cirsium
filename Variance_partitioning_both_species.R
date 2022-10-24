setwd("C:/Users/hanna/Documents/Universitetet/PhD/GitHub/Cirsium")


#load libraries
library(readxl)
library(Hmsc)
library(corrplot)

#data wrangling
library(data.table)
library(dplyr)
library(reshape)

#plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggsci)

#set a seed to make data reproducible
set.seed(1)


#####data preparation#####

#response matrix: floral scent
#needs to have scent compounds in the columns and samples in the rows

comps <- read.csv("Control_adjusted_samples.csv")

Y <- as.matrix(comps[,12:58])

#data frame with explanatory variables
#one column per factor
#the rows should be the same as for Y
#XData cannot contain NA?s
XData <- data.frame(Sample.type      = as.factor    (comps$Sample.type),
                    Species   = as.factor   (comps$Species)
) 

#data frame with random factors
studyDesign = data.frame(ID         = as.factor (comps$ID),
                         Pot.ID     = as.factor (comps$Pot.ID),
                         Syringe.no = as.factor(comps$Syringe.no.)
)

#ranLevels: random factors/ latent variables 
rL1 = HmscRandomLevel(units = levels(studyDesign$ID    ))    
rL2 = HmscRandomLevel(units = levels(studyDesign$Pot.ID       ))
rL3 = HmscRandomLevel(units = levels(studyDesign$Syringe.no       ))



#####model construction#####

#using abundance data (as opposed to occurence data)
#use square root transformed scent data to deal with extremely high count values in your scent data
#YScale: scales responses for which normal distribution is assumed 
#YScale: this is another way of dealing with extremely high count values in your scent data
#I am assuming normal distribution of errors
m = Hmsc(Y=sqrt(Y), 
         XData       = XData, 
         XFormula    = ~Sample.type          +
           Species,             
         
         YScale      = T,                               
         distr       = ("normal"),
         studyDesign = studyDesign,
         ranLevels   = list(ID = rL1, 
                            Pot.ID = rL2,
                            Syringe.no = rL3
                            
         ))



#check the model input
m$Y[1:5, 1:5]
head(m$X)
m$rL
m 
#This R package was designed for community data. 
#That is why the model calls the data in the columns "species"
#even though it will be the scent compounds in your case. 
#Don?t get confused by that.


#####fit model#####

#you can use these settings for a quick test run, 
#but you obviously have to increase sampling for your final analysis
m = sampleMcmc(m,                 
               thin = 1, 
               samples = 100,     
               transient = 10,    
               nChains = 1,       
               verbose = 20)       


#test model fit
mpost = convertToCodaObject(m)

pdf("betapost.pdf")
plot(mpost$Beta)
dev.off()
#the left hand side plots show the sampling and should cover the whole range,
#meaning there should be large fluctuations
#the right hand side plots should show normal distribution, if not, you can try 
#to increase sampling

#this tests whether the sampling covered the whole possible parameter space
#the mean should be somewhat close to the "samples" size you set above in the
#"sampleMCMC" command
summary(effectiveSize(mpost$Beta))


# variance partitioning
# recommended to group variables related to the same theme, as then the variance 
# component accounts also for co-variation within the group
m$covNames

#to compute variances for each factor level
VP = computeVariancePartitioning(m)

#to group all levels of one factor together
#example:
VP = computeVariancePartitioning(m, group = c(1,1,1,2,3), 
                                 groupnames = c("Pop","Mating","Type.x"))

plotVariancePartitioning(m, VP = VP)


#evaluate model explanatory power in terms of R2
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
MF$R2
MF$R2[is.na(MF$R2)] <- 0
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))




##### Everything below here is optional: you can use this script if you wish 
#to include the proportion of unexplained variance into your plot, and/ or
#if you wish to plot in ggplot. It involves a lot of data wrangling. #####

#extract values of the variance partitioning
VPvals1 <- as.data.frame(VP$vals)
VPvals  <- transpose(as.data.frame(VPvals1))
colnames(VPvals) <- rownames(VPvals1)
VPvals$compound  <- colnames(VPvals1)
VPvals[is.na(VPvals)] <- 0

#extract R2 for each compound (i.e. variance explained by the model)
VPvals$R2<-MF$R2

#create a new variable for the variance not explained by the model
VPvals$Unexplained   <- 1-VPvals$R2

#scale the variance explained by my explanatory variables by the R2
#that way, all the variance explained by each variance component 
#plus the unexplained variance will sum to 1


count <- 0
for (i in 1:5) {      #instead of "9" you should add the number of factors you have 
  res<-VPvals[,i]*VPvals$R2
  count <- count + 1
  VPvals[count] <- res
}


#calculate the mean variance explained by each of my variance components
means <- rep(NA, 5)
count <- 0
for (i in 1:5) {
  mean<-mean(VPvals[,i])
  count <- count + 1
  means[count] <- mean
}

#the means start with the component that is on the left in the VPvals df
means
mean(VPvals$Unexplained)

#check that the scaled variances plus the unexplained variance sum to 1
sum(means)+mean(VPvals$Unexplained)



#change format of the table so I can read it into ggplot
VPvals<-VPvals[,-7] #remove R2 column
VPvals<-melt(VPvals, id="compound")



#arrange variance plot by variable of choice
df <- VPvals %>% 
  filter(variable == "Sample.type") %>% 
  arrange(desc(value))
VPvals$compound <- factor(VPvals$compound, levels=df$compound)

ggplot(VPvals, aes(x=compound, y=value, fill=variable)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x="Compound",y="Variance proportion") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) + 
  scale_fill_npg(name = "Source of variation",
                 labels=c("Sample type (34%)",
                          "Species (10%)",
                          "Random: ID (13%)",
                          "Random: Pot ID (15%)",
                          "Random: Syringe (5%)",
                          "Unexplained variance (22%)")) + 
  scale_y_continuous(labels = scales::percent) 


dev.copy(tiff,"plot.tiff",  width=4500, height=2200, res=600)
dev.off()

