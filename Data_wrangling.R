setwd("C:/Users/hanna/Documents/Universitetet/PhD/GitHub/Cirsium/")

```{r}
library(tidyverse)
library(readxl)
library(dplyr)
```

```{r}



##getting data loaded in long format
path <- file.choose("Second_integration")
a <- path %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path = path, skip=4, sheet = .x,), .id = "sheet") %>% 
  rename(ID=sheet, Peak=`...1`) %>% select(Apex_RT, ID, Peak, Area) %>% 
  mutate(file="201103")
#View(a)

write_csv(a, "2nd_cirsium_data.csv")

#gets occurrence frequencies
counts<-a %>% group_by(Peak) %>% summarize(N=length(ID))
view(counts)

write_csv(counts, "2nd_cirsium_counts.csv")

#cleaning up names

new_names <- read_csv("Compound_dictionary.csv")

#need to check all of these names for consistency

compound_list<-c("B.Phellandrene",
                 "A.Phellandrene",
                 "Eucalyptol",
                 "Cis.Beta.Ocimene",
                 "Trans.Beta.Ocimene",
                 "3.Hexenyl.acetate",
                 "3.Hexen.1.ol",
                 "Unknown.terpenoid.1",
                 "Cis.3.Hexenyl.isovalerate",
                 "Copaene",
                 "Benzaldehyde",
                 "Linalool",
                 "B.Cubebene",
                 "B.Elemene",
                 "Caryophyllene",
                 "2.Amino.1.phenylethanol",
                 "B.Farnesene",
                 "A.Caryophyllene",
                 "4.Oxoisophorone",
                 "G.Muurolene",
                 "Unknown.terpenoid.2",
                 "Lilac.alcohol.1",
                 "A.Farnesene.Z.E",
                 "B.Selinene",
                 "Lilac.alcohol.2",
                 "A.Farnesene",
                 "Methyl.salicylate",
                 "Lilac.alcohol.3",
                 "Phenylethyl.acetate",
                 "Lilac.alcohol.4",
                 "Benzyl.alcohol",
                 "Cinnamaldehyde",
                 "Phenylethyl.alcohol",
                 "Isoamyl.benzoate",
                 "3.Phenyl.1.propanol.acetate",
                 "Cis.jasmone",
                 "Caryophyllene.oxide",
                 "4.Methoxybenzaldehyde",
                 "Benzenepropanol",
                 "Benzyl.tiglate",
                 "3.Hexen.1.ol.benzoate",
                 "Cinnamyl.acetate",
                 "Eugenol",
                 "Jasmine.lactone",
                 "Cinnamyl.alcohol",
                 "Benzophenone",
                 "Benzyl.benzoate",
                 "Phenethyl.benzoate"
)


                 
scents_clean <- left_join(a, new_names, by="Peak") %>% 
  relocate(`New_name`, .after=ID) %>% select(ID, Apex_RT, `New_name`, Area, file) %>% rename(Peak=`New_name`)


write_csv(scents_clean, "clean_long_scent_data_pre_filter.csv")

scents <- scents_clean %>% 
  filter(Peak %in% compound_list) 

scents <- distinct(scents)

write_csv(scents, "clean_long_cirsium_data.csv")


##Create big master sheet

scents <- read_csv("clean_long_cirsium_data.csv")
master <- read.csv("SPME_sampling.csv")



library(tidyr)
install.packages("reshape2")
library(reshape2)
### For details on dcast see https://stackoverflow.com/questions/33051386/dcast-warning-aggregation-function-missing-defaulting-to-length

wide <- dcast(scents, ID + file ~ Peak, fun.aggregate = sum, value.var="Area")
write.csv(wide, "scents_wide.csv")


big <- left_join(master, wide, by=("ID"), all=TRUE)
big <- big[1:21,]
big <- subset(big,select=-c(11))

write.csv(big, "scents_master.csv")


