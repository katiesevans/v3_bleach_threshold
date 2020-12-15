# Date is NOV 15th, 2019


# L1-L4 assay V3

library(easysorter)
library(tidyverse)
library(ggplot2)
library(broom)
library(dplyr)
library(RColorBrewer)
library(ggridges)
library(cowplot)

source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_genotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_phenotype_plots.R")


#Define a vector of your experiement directories
dirs <- "~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/20191115_v3_mian/"


# Read in the data
raw_mian <- easysorter::read_data(dirs)

# Remove all data from the contaminated wells
raw_nocontam_mian<- easysorter::remove_contamination(raw_mian)

# look at strain distributions
raw_nocontam_mian %>%
  dplyr::filter(condition == "MIAN-0") %>%
  dplyr::select(strain, condition, EXT, TOF) %>%
  tidyr::gather(trait, value, EXT:TOF) %>%
  ggplot(.) +
  aes(x = value) +
  geom_histogram() +
  facet_grid(strain~trait, scales = "free") +
  theme_bw() +
  labs(title= "Mianserin")

# Summarize the data
# directories = FALSE because supplied data comes from one directory
summedraw_mian <- easysorter::sumplate(raw_nocontam_mian, directories = FALSE, quantiles = TRUE, v3_assay = TRUE) 

#Prune based on biological impossibilities
pruned_mian<- easysorter::bioprune(summedraw_mian) %>%
  tidyr::gather(trait, phenotype, -(date:col))%>%
  dplyr::filter(!grepl("red|green|yellow|f.|iqr", trait)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(control = ifelse(condition == "MIAN-0", "None", "MIAN-0")) %>% #mutate changes what you say; this is an if else statement 
  #and the usage is "ifelse(test, yes, no)" ,
  #so this is saying if condition has cis-0 then change the control to "none", if not mian-0, (=",") change the control to mian-0 %>%
  tidyr::drop_na(condition)

#look at bleaches
pruned_mian <- pruned_mian %>%
  mutate(assay = case_when(plate %in% c("1","2","3","4") ~
                             "A", 
                           plate %in% c("5","6","7","8") ~
                             "B", 
                           plate %in% c("9","10","11","12") ~
                             "C"))

#mutated assay column bc first 4 plates were bleach A, 5-8 were bleach B, and 9-12 were bleach C

a <- raw_nocontam_flx %>%
  mutate(assay = case_when(plate %in% c("1","2","3","4") ~
                             "A", 
                           plate %in% c("5","6","7","8") ~
                             "B", 
                           plate %in% c("9","10","11","12") ~
                             "C"))

ggplot(subset(a, !strain == "NA" & condition == "MIAN-0"), 
       aes(x=norm.EXT, y=norm.green)) + geom_jitter(size=0.5, alpha=0.5) + facet_wrap(~strain)

ggplot(subset(a, !strain == "NA" & condition == "MIAN-0"), 
       aes(x=log(TOF), y=log(EXT))) + 
  geom_jitter(aes(colour=assay), size=0.5, alpha=0.5) +
  facet_wrap(~strain) 

ggplot(subset(a, !strain == "NA" & condition == "MIAN-0"), 
       aes(x=log(TOF), y=log(EXT))) +
  geom_jitter(aes(colour=assay), size=0.5, alpha=0.5) + 
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~strain)

pruned_mian %>% 
  # dplyr::mutate(drug = stringr::str_split_fixed(condition, "_", 2)[,1],
  #background = ifelse(strain %in% c("N2", "ECA232", "EC377", "ECA386","ECA1113"), "N2", "CB"),
  #dose = stringr::str_split_fixed(condition, "_", 2)[,2]) %>%
  tidyr::drop_na(condition) %>%
  dplyr::filter(trait == "mean.EXT") %>%
  ggplot(.) +
  aes(x = factor(condition, levels = c("MIAN-25", "MIAN-50", "MIAN-100", "MIAN-200")), y = phenotype, fill = strain) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 0.5) +
  geom_boxplot(alpha = 0.5, outlier.color = NA) +
  theme_bw() +
  labs(x = "Dose", y = "mean.EXT") +
  scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "ECA377"="yellow", "ECA232"="red", "ECA386"="green", "ECA1113"="cyan" )) 


pruned_mian %>%
  dplyr::filter(condition == "MIAN-0") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.)+
  aes(x= mean.EXT , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") +
  #scale_color_manual(values=c("B"="black", "C"="turquoise")) +
  coord_flip()

# variability of titering: n / bleach / strain !!!!
ggplot(spread_mian) +
  aes(x=strain, y=n) +
  geom_boxplot() + geom_jitter(aes(fill=assay), width = 0.1, shape=21) +
  scale_color_brewer(palette = "Set1") + facet_grid(.~assay)

pruned_flx <- prune_outliers(pruned_flx)

assreg <- regress(pruned_mian, assay = TRUE)
