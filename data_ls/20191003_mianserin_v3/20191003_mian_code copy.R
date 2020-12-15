
# Date is October 7st, 2019; assay happened on Thurs 10/03/219
# Experiment name is "20191003_mianV3"

# this was day before finding mites 
  
# L1-L4 assay V3
  
library(easysorter)
library(tidyverse)
library(ggplot2)
library(broom)
library(dplyr)
library(RColorBrewer)
library(ggridges)
library(cowplot)
  
#source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts/NIL_genotype_plots.R")
#source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts/NIL_phenotype_plots.R")
  
  
#Define a vector of your experiement directories
dirs <- "~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/20191003_mianserin_v3"
  

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

a <- raw_nocontam_mian %>%
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

ggplot(subset(a, !strain == "NA" & condition == "MIAN-50"), 
        aes(x=log(TOF), y=log(EXT))) +
        geom_jitter(aes(colour=assay), size=0.5, alpha=0.5) + 
        scale_color_brewer(palette = "Set2") +
        facet_wrap(~strain)


b <- summedraw_mian %>%
  ungroup(assay) %>%
  mutate(assay = case_when(plate %in% c("1","2","3","4") ~
                             "A", 
                           plate %in% c("5","6","7","8") ~
                             "B", 
                           plate %in% c("9","10","11","12") ~
                             "C"))

#looking at variability of titering bc n = # of worms / bleach/ strain 
ggplot(subset(b, !assay == "NA" & !strain == "NA"), aes(x=strain, y=n)) + 
    geom_boxplot() + geom_jitter(aes(fill=assay), width = 0.1, shape=21) +
    scale_color_brewer(palette = "Set1") + facet_grid(.~assay)


ext<- ggplot(subset(a, !assay == "NA" & !strain == "NA"), aes(x=EXT, y=strain)) + stat_density_ridges(alpha = 0.8, fill = "blue", scale = 0.98) +
  facet_grid(.~assay) 
tof <-ggplot(subset(a, !assay == "NA" & !strain == "NA"), aes(x=TOF, y=strain)) + stat_density_ridges(alpha = 0.8, fill = "blue", scale = 0.98) +
  facet_grid(.~assay)

green<- ggplot(subset(a, !assay == "NA" & !strain == "NA"), aes(x=green, y=strain)) + stat_density_ridges(alpha = 0.8, fill = "blue", scale = 0.98) +
  facet_grid(.~assay) 
green

plots <- list("ext"=ext, "tof" = tof)

plot_grid(plotlist=plots, nrow=1,ncol=2)

library(cowplot)

        

a %>%
  dplyr::filter(condition == "MIAN-0") %>%
  dplyr::select(strain, condition, EXT, TOF) %>%
  tidyr::gather(trait, value, EXT:TOF) %>%
  ggplot(.) +
  aes(x = value) +
  geom_histogram() +
  facet_grid(strain~trait, scales = "free") +
  theme_bw() +
  labs(title= "Mianserin")

raw_nocontam_flx %>%    #change this df for assay to plot  
  dplyr::filter(condition == "MIAN-0") %>%
  dplyr::select(strain, condition, EXT, TOF) %>%
  tidyr::gather(trait, value, EXT:TOF) #%>%
  #ggplot(.) +
  #aes(x = value) +
  #geom_histogram() +
  #facet_grid(strain~trait, scales = "free") +
  #theme_bw() +
  #labs(title= "Mianserin")

save(raw_nocontam_mian, file="~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/20191003_mianserin_v3/processed/raw_nocontam_mian")
  
  
# save pruned data
save(pruned_mian, file = "~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/20191003_mianserin_v3/processed/pruned_mian.Rda")

pruned_mian <- pruned_mian %>%
  dplyr::filter(! (assay=="A"))

pruned_mian <- prune_outliers(pruned_mian)

assreg <- regress(pruned_mian, assay = TRUE)

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
  conditiondf <- assreg %>%
    dplyr::filter(!grepl("-0", condition)) #number doesnt start with zero
  for(c in unique(conditiondf$condition)) {
    newdf <- assreg %>%
      dplyr::filter(condition %in% c(c, "MIAN-0") )#filters for c= mian doses and mian-0
# dplyr::mutate(condition = drug) %>%
    reg <- easysorter::regress(newdf)
    regressed <- rbind(regressed, reg)
  }
  
  
# plot each condition regression
  
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/pca_analysis.R")  



  
regressed25 <- regressed %>%
  dplyr::filter(trait != "n" & condition == "MIAN-25")

df <- pruned_mian %>%
  dplyr::select(date:phenotype) %>%
  dplyr::filter(condition=="MIAN-25") %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.)+
  aes(x= mean.EXT , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(~assay) +
  #scale_color_manual(values=c("B"="black", "C"="turquoise")) +
  coord_flip()
df

#regress by n 
n <- pruned_mian %>%
  tidyr::spread(trait, phenotype) %>%
  tidyr::gather( key = trait , value = phenotype, cv.EXT:var.TOF, -n)
  
regress
function (dataframe, assay = FALSE) 
{
  dataframe <- easysorter:::ensure_long(dataframe)
  dataframe <- dplyr::filter(dataframe, is.finite(phenotype), 
                             !is.na(strain))
  dataframe <- dataframe %>% dplyr::filter(trait != "n.sorted")
  
    # data <- dataframe %>% dplyr::filter(!is.na(control))
    # controls <- dataframe %>% dplyr::filter(is.na(control) | control == "None")
    # 
    # controls$control <- controls$condition
    # moltendata <- data
    # moltencontrols <- controls %>% 
    #   dplyr::group_by(strain, control, trait, assay) %>% 
    #   dplyr::summarize(controlphenotype = mean(phenotype, na.rm = TRUE))
    # 
    # fusedmoltendata <- dplyr::left_join(moltendata, moltencontrols, 
    #                                     by = c("strain", "control", "trait", "assay")) %>% 
    #   dplyr::filter(!is.na(phenotype), !is.na(controlphenotype))
  fusedmoltendata <- dataframe %>%
    dplyr::group_by(strain, trait, assay, condition) %>%
    dplyr::mutate(avg_n = mean(n, na.rm = T))
    
    regressed <- fusedmoltendata %>% 
      dplyr::group_by(condition, trait) %>% 
      dplyr::do(fit = lm(phenotype ~ avg_n -  1, .))
    
    regressedframe <- regressed %>% 
      broom::augment(fit) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(fusedmoltendata, ., by = c("condition", "trait", "phenotype", "avg_n")) %>% 
      dplyr::distinct(condition, trait, phenotype, avg_n, strain, row, col, plate, .keep_all = T) %>%
      dplyr::rename(resid = .resid) %>% 
      dplyr::mutate(phenotype = resid) %>% 
      dplyr::select(-resid)
  return(regressedframe)
}


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


  
#make PC
#calc_pc_reps is for NILS, noreps is for RILS
# only keeps PCs that capture 90% of phenotype 
pc <- calc_pc_reps(regressed25)[[2]]

trt <- "PC1" 
plot <- pc %>%
  dplyr::filter(trait == trt) %>%
  mutate(phenotype= phenotype* -1) %>%
  ggplot(.) +
  aes(x = factor(strain, levels= c("ECA386", "ECA1113", "ECA377", "ECA232", "N2", "CB4856")), 
      y = phenotype) + #add fill= strain later 
  #scale_fill_manual(values = c("N2" = "grey", "CB4856" = "blue", "ECA377"="yellow", "ECA232"="red", "ECA386"="green", "ECA1113"="cyan" )) +
  geom_jitter(aes(color=factor(assay)), width = 0.1, size = 2) +
  scale_color_manual(values=c("B"="black", "C"="turquoise")) +
  geom_boxplot(alpha = 0.5, outlier.color = NA) +  
  theme_bw() +
  facet_grid(~factor(condition, levels = c("MIAN-0", "MIAN-25", "MIAN-50", "MIAN-100", "MIAN-200")), scales = "free") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = trt) +
  coord_flip()
plot

cowplot::plot_grid(chroms, plot)



trt <- "q90.TOF" 
plot <- regressed %>%
    dplyr::filter(trait == trt & condition== "MIAN-25") %>%
    ggplot(.) +
    aes(x = factor(strain, levels= c("ECA386", "ECA1113", "ECA377", "ECA232", "N2", "CB4856")), 
        y = phenotype) + #add fill= strain later 
    #scale_fill_manual(values = c("N2" = "grey", "CB4856" = "blue", "ECA377"="yellow", "ECA232"="red", "ECA386"="green", "ECA1113"="cyan" )) +
    geom_jitter(aes(color=factor(assay)), width = 0.1, size = 2) +
  scale_color_manual(values=c( "B"="black", "C"="turquoise")) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) + #changed to violin to see distribution 
    theme_bw() +
    facet_grid(~factor(condition, levels = c("MIAN-0", "MIAN-25", "MIAN-50", "MIAN-100", "MIAN-200")), scales = "free") +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "", y = trt) +
  coord_flip()
plot



##
chroms <- nil_plot(c("N2" , "CB4856", "ECA377", "ECA232", "ECA386", "ECA230"), chr = "V", all.chr = TRUE) [[1]]

cowplot::plot_grid(chroms, plot)


trt <- "q90.TOF" 
plot <- assreg %>%
  dplyr::filter(trait == trt & condition== "MIAN-0") %>%
  ggplot(.) +
  aes(x = factor(strain, levels= c("ECA386", "ECA1113", "ECA377", "ECA232", "N2", "CB4856")), 
      y = phenotype) + #add fill= strain later 
  #scale_fill_manual(values = c("N2" = "grey", "CB4856" = "blue", "ECA377"="yellow", "ECA232"="red", "ECA386"="green", "ECA1113"="cyan" )) +
  geom_jitter(aes(color=factor(assay)), width = 0.1, size = 2) +
  scale_color_manual(values=c("A"="red", "B"="black", "C"="turquoise")) +
  geom_boxplot(alpha = 0.5, outlier.color = NA) + #changed to violin to see distribution 
  theme_bw() +
  facet_grid(~factor(condition, levels = c("MIAN-0", "MIAN-25", "MIAN-50", "MIAN-100", "MIAN-200")), scales = "free") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = trt) +
  coord_flip()
plot

cowplot::plot_grid(chroms, plot)
  
