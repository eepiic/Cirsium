
setwd("C:/Users/hanna/Documents/Universitetet/PhD/GitHub/Cirsium/")


install.packages("devtools")
install.packages("bouquet")
install.packages("tidyverse")
install.packages("readxl")
install.packages("reshape2")
install.packages("data.table")
install.packages("janitor")

library(devtools)
library(bouquet)
library(tidyverse)
library(readxl)
library(reshape2)
library(data.table)
library(janitor)

View(metadata)
View(chemtable)
View(sampletable)


file1 <- read.csv("scents_master.csv")
file1 <- subset(file1, select=-c(1))

####Create the metadata table####

meta <- subset(file1, select=c(1:11))

write.csv(meta, "Filter_script_metadata.csv")

####Create chemtable####

#Put ctrls before samples
dat <- subset(file1, select=-c(1:7))

dat_list <- t(dat)
dat_list <- as.data.frame(dat_list)

chemtable <- as.data.frame(rownames(dat_list)) #list of compounds
colnames(chemtable) <- c("names")


####create sample table###

samples <- subset(file1, select=c(11:50))
rownames(samples) <- samples[,1]
samples <- samples[,-1]


####Run filter function####

filter <- filter_ambient_date(chemtable, samples, meta)
prune <- prune_sampletable_dates(samples,chemtable,meta,filter) #wont remove sample 8913, have to do that manually in csv file

write.csv(prune, "Ctrl_adjusted_ng_per_flowerORdrymass_per_hour.csv")

###Create big master sheet###

data <- read.csv("Ctrl_adjusted_ng_per_flowerORdrymass_per_hour.csv")
meta <- read.csv("Filter_script_metadata.csv")
names(data)[1] <- "ID"

total <- left_join(meta, data, by="ID")
total <- subset(total, select=-c(X, Corr.ctrl, date, type))

write.csv(total, "Master_scent_ng_per_flower_OR_drymass_per_hour_none_removed.csv")
