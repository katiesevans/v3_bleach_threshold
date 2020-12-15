# Minimal directory name
# Date is October 7st, 2019
# Experiment name is "20191003_mianV3"

# 20191003_mianV3/
  
  #####################################


# L1-L4 assay V3
  
library(easysorter)
library(tidyverse)
library(ggplot2)
library(broom)
  
#source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts/NIL_genotype_plots.R")
#source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts/NIL_phenotype_plots.R")
  
  
#Define a vector of your experiement directories
dirs <- "~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/20191003_mianserin_v3"
  

# Read in the data
raw_mian <- easysorter::read_data(dirs)
  
# Remove all data from the contaminated wells
raw_nocontam_mian<- easysorter::remove_contamination(raw_mian)
  
# Summarize the data
# directories = FALSE because supplied data comes from one directory
summedraw_mian <- easysorter::sumplate(raw_nocontam_mian, directories = FALSE, quantiles = TRUE, v3_assay == TRUE) 
  
#Prune based on biological impossibilities
pruned_mian<- easysorter::bioprune(summedraw_mian) %>%
  tidyr::gather(trait, phenotype, -(date:col))%>%
  dplyr::filter(!grepl("red|green|yellow|f.|iqr", trait)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(control = ifelse(condition == "MIAN-0", "None", "MIAN-0")) %>% #mutate changes what you say; this is an if else statement 
  #and the usage is "ifelse(test, yes, no)" ,
  #so this is saying if condition has cis-0 then change the control to "none", if not mian-0, (=",") change the control to mian-0 %>%
  tidyr::drop_na(condition)
  
  
  
# save pruned data
save(pruned_mian, file = "~/Desktop/20190909_Mianserin/processed/pruned_mian.Rda")
  
# look at dose
pruned_mian %>% 
  # dplyr::mutate(drug = stringr::str_split_fixed(condition, "_", 2)[,1],
  #background = ifelse(strain %in% c("N2", "ECA232", "EC377", "ECA386","ECA1113"), "N2", "CB"),
  #dose = stringr::str_split_fixed(condition, "_", 2)[,2]) %>%
  tidyr::drop_na(condition) %>%
  dplyr::filter(trait == "mean.EXT") %>%
    ggplot(.) +
    aes(x = factor(condition, levels = c("MIAN-0", "MIAN-25", "MIAN-50", "MIAN-100", "MIAN-200")), y = phenotype, fill = strain) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 0.5) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) +
    theme_bw() +
    labs(x = "Dose", y = "mean.EXT") +
    scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "ECA377"="yellow", "ECA232"="red", "ECA386"="green", "ECA1113"="cyan" )) 
  
  
  
# regress for each dose
  regressed <- NULL
  conditiondf <- pruned_mian %>%
    dplyr::filter(!grepl("-0", condition)) #number doesnt start with zero
  for(c in unique(conditiondf$condition)) {
    newdf <- pruned_mian %>%
      dplyr::filter(condition %in% c(c, "MIAN-0") )#filters for c= mian doses and mian-0
# dplyr::mutate(condition = drug) %>%
    reg <- easysorter::regress(newdf)
    regressed <- rbind(regressed, reg)
  }
  
  
# plot each condition regression
trt <- "mean.EXT"
  mian_mean.EXT <- regressed %>%
    dplyr::filter(trait == trt) %>%
    ggplot(.) +
    aes(x = strain, 
        y = phenotype, fill = strain) +
    scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "ECA377"="yellow", "ECA232"="red", "ECA386"="green", "ECA1113"="cyan" )) +
    geom_jitter(width = 0.1, size = 0.5) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) +
    theme_bw() +
    facet_grid(~factor(condition, levels = c("MIAN-0", "MIAN-25", "MIAN-50", "MIAN-100", "MIAN-200")), scales = "free") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90)) +
    labs(x = "", y = "mean.EXT") 
  mian_mean.EXT
  
