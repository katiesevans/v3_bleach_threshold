# 2020 12 14 --> V3 thresholds with Katie Evans

# all my assays are with lysate 

library(easysorter)
library(COPASutils)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dplyr)
library(cowplot)

source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_genotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_phenotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Loraina/Scripts/Base_theme.R")

#Error in read_template(strainsfile, type = "strains") : 
# The strains template file at ~/Dropbox/v3_bleach_threshold/data_ls/20191115_v3_flx/strains/V5.csv --> removed for now


#Define a vector of your experiement directories
dirs <-("~/Dropbox/v3_bleach_threshold/data_ls/20191003_mianserin_v3/")

dirs2 <- ("~/Dropbox/v3_bleach_threshold/data_ls/20191114_v3_flx/")
dirs4 <-  ("~/Dropbox/v3_bleach_threshold/data_ls/20191115_v3_mian/")  ### this folder had a processed folder i had to get rid of 
dirs5 <-("~/Dropbox/v3_bleach_threshold/data_ls/20200214_flx/") 
dirs6 <-("~/Dropbox/v3_bleach_threshold/data_ls/20200214_mian/")
dirs7 <-("~/Dropbox/v3_bleach_threshold/data_ls/20201105_flx/")
dirs8 <-("~/Dropbox/v3_bleach_threshold/data_ls/20201106_mianserin/")
dirs9 <-("~/Dropbox/v3_bleach_threshold/data_ls/20201112_flxmian_A1/") 
dirs10 <-("~/Dropbox/v3_bleach_threshold/data_ls/20201113_flxmian_A2/")
# reads in as a list ; extract out each element of list and rbind into a df orrrr do 9 separate easy sorters 


# Read in the data
raw <- easysorter::read_data(dirs)

# Remove all data from the contaminated wells
raw_nocontam1 <- easysorter::remove_contamination(raw)

raw_nocontam1 <- raw_nocontam1 %>%
  mutate(assay = case_when(plate %in% c("1","2","3", "4", "5") ~
                             "A",
                           plate %in% c("6","7","8", "9", "10", "11") ~
                             "B", 
                           plate %in% c("12", "13", "14", "15", "16", "17") ~
                             "C"))

raw2 <- easysorter::read_data(dirs2)
raw3 <- easysorter::read_data(dirs5)
raw4 <- easysorter::read_data(dirs4)
raw4 <- easysorter::read_data(dirs6)
raw5 <- easysorter::read_data(dirs7)
raw6 <- easysorter::read_data(dirs8)
raw7 <- easysorter::read_data(dirs9)
raw8 <- easysorter::read_data(dirs10)

raw_nocontam2 <- raw2 %>%
  mutate(assay = case_when(plate %in% c("1","2","3", "4") ~
                             "A",
                           plate %in% c("5","6","7", "8") ~
                             "B", 
                           plate %in% c("9", "10", "11","12") ~
                             "C"))

raw_nocontam3 <- raw3 %>%
  mutate(assay = case_when(plate %in% c("7","8","9") ~
                             "A"))

raw_nocontam4 <- raw4 %>%
  mutate(assay = case_when(plate %in% c("1","2","3") ~
                             "A",
                           plate %in% c("4","5","6") ~
                             "B"))


raw_nocontam5 <- raw5 %>%
  mutate(assay = case_when(plate %in% c("1","2","3") ~
                             "A",
                           plate %in% c("4","5","6") ~
                             "B",
                           plate %in% c("7", "8", "9")~
                             "C"))

raw_nocontam6 <- raw6 %>%
  mutate(assay = case_when(plate %in% c("1","2","3","4") ~
                             "A",
                           plate %in% c("5","6","7", "8") ~
                             "B",
                           plate %in% c("9", "10", "11")~
                             "C"))

raw_nocontam7 <- raw7 %>%
  mutate(assay = case_when(plate %in% c("1","2","3","4", "5","6") ~
                             "A",
                           plate %in% c("7", "8", "9", "10", "11", "12") ~
                             "B",
                           plate %in% c("13", "14", "15", "16", "17")~
                             "C"))

raw_nocontam8 <- raw8 %>%
  mutate(assay = case_when(plate %in% c("6","7", "8", "9", "10") ~
                             "B",
                           plate %in% c( "11", "12","13", "14", "15")~
                             "C"))


all_nocontam <- rbind(raw_nocontam1,raw_nocontam2, raw_nocontam3, raw_nocontam4, raw_nocontam5, raw_nocontam6, raw_nocontam7,raw_nocontam8)



# look at strain distributions
all_nocontam %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::filter(strain %in% c("N2", "CB4856")) %>%
  dplyr::select(date,strain, condition, EXT, TOF) %>%
  tidyr::gather(trait, value, EXT:TOF) %>%
  ggplot(.) +
  aes(x = value) +
  geom_histogram() +
  facet_grid(date~trait, scales = "free") +
  theme_bw() +
  labs(title= "DMSO")

mean(all_nocontam$TOF) # 351

mean(all_nocontam$EXT) #919


# less plates going into a bleach --> harder to tell when to stop? during bleach? 
all_nocontam <- all_nocontam %>%
  mutate(person="LS")


# use staging data to tell what percent of data is L1 , L2, L3, L4 