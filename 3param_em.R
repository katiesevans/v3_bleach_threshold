library(easysorter)
library(COPASutils)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mclust)



# read in data, filter to keep only controls
wd <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_kse/"
filelist <- list.files(wd)

kse_data <- NULL
for(file in filelist) {
  data <- easysorter::read_data(glue::glue("{wd}/{file}")) %>%
    dplyr::filter(contamination == FALSE,
                  call50 == "object", 
                  condition %in% c("DMSO", "water", "lyaste")) %>%
    dplyr::select(date:col, condition, strain, TOF, EXT, green,stage) %>%
    dplyr::mutate(person = "KSE")
  kse_data <- rbind(kse_data, data)
}
# mark unique assays
kse_data <- kse_data %>%
  na.omit() %>%
  dplyr::mutate(bleach = paste0(date, "_", assay)) %>%
  dplyr::filter(strain == "N2") %>%
  dplyr::filter(assay != "del")

kse_l1 <- kse_data %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) 

### # from 20210128_ls.R code --> added green to the df 
LS_data <- LS_data %>%
  dplyr::filter(strain == "N2")
ls_l1 <- LS_data %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) 

######

ls_log <- LS_data %>%
  dplyr::mutate(logGreen = log(green),logTOF = log(TOF), logEXT = log(EXT)) %>% # sorter doesnt cut off stages using EXT 
  dplyr::select(logTOF, logEXT, logGreen)  %>% # taking log separates better 
  mclust::Mclust(., G=3)

summary(ls_log)
# Clustering table:
#    1    2    3 
# 8076 1427 1130  
plot(ls_log)


kse_log <- kse_data %>%
  dplyr::mutate(green,logTOF = log(TOF), logEXT = log(EXT)) %>% # sorter doesnt cut off stages using EXT 
  dplyr::select(logTOF, logEXT, green)  %>% # taking log separates better 
  mclust::Mclust(., G=3)

summary(kse_log)
# Clustering table:
#  1    2    3 
#  5321 1220 4516 
plot(kse_log)


 












