setwd("C:/Users/hanna/Documents/Universitetet/PhD/Me/GitHub/Cirsium/")

```{r}
library(tidyverse)
library(readxl)
library(dplyr)
```

```{r}



##getting data loaded in long format
path <- file.choose("Integration")
a <- path %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path = path, skip=4, sheet = .x,), .id = "sheet") %>% 
  rename(ID=sheet, Peak=`...1`) %>% select(Apex.RT, ID, Peak, Area) %>% 
  mutate(file="201103")
#View(a)

write_csv(a, "cirsium_data.csv")

#gets occurrence frequencies
counts<-a %>% group_by(Peak) %>% summarize(N=length(ID))
view(counts)

write_csv(counts, "cirsium_counts.csv")

#cleaning up names

new_names <- read_csv("Compound_dictionary.csv")
new_names_short <- new_names %>% filter(`New_name`!="NA")

#need to check all of these names for consistency

compound_list<-c("Eucalyptol", "Cis.Beta.Ocimene", "3.Hexenyl.acetate", "3.Hexen.1.ol", "Copaene", "Benzaldehyde", "Linalool",
                 "Elemene derivative.1", "Caryophyllene", "A.farnesene", "A.caryophyllene", "4.oxoisophorone", "A.cubebene",
                 "Germacrene.D", "Lilac.alcohol.", "Elemene derivative.2", "Lilac.alcohol.2", "B-farnesene", "Methyl.salicylate", 
                 "Lilac.alcohol.3", "Phenylethyl.acetate", "Lilac.alcohol.4", "Benzyl alcohol", "Cinnamaldehyde", "Phenylethyl.alcohol",
                 "Isoamyl benzoate", "Cis.jasmone", "Caryophyllene.oxide", "4.Methoxybenzaldehyde", "Benzenepropanol", "Benzyl.tiglate",
                 "3.Hexen.1.ol.benzoate", "Cis.3.Hexenyl.benzoate", "3.Hexen.1.ol.benzoate", "Eugenol", "Jasmine.lactone", "Cinnamyl.alcohol",
                 "Benzophenone", "Benzyl.benzoate", "Phenethyl.benzoate")


                 
scents_clean <- left_join(a, new_names_short, by="Peak") %>% 
  relocate(`New_name`, .after=ID) %>% select(ID, Apex.RT, `New_name`, Area, file) %>% rename(Peak=`New_name`)


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

scents34 <- read.csv("scents_34_wide.csv")

big <- left_join(master, wide, by=("ID"), all=TRUE)
big <- big[1:21,]
big <- subset(big,select=-c(11))

write.csv(big, "scents_master.csv")


