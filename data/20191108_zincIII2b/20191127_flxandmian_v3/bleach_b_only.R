#assay from 20191127, 6 dmso plates, 6 flx and 6 mian 
#### look at just bleach b , 2 plates of each condition; ~30 reps each strain/ condition 


source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_genotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_phenotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Loraina/Scripts/Base_theme.R")
# L1-L4 assay V3

library(easysorter)
library(tidyverse)
library(ggplot2)
library(broom)
library(dplyr)
library(RColorBrewer)
library(ggridges)
library(cowplot)

#Define a vector of your experiement directories
dirs <- "~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/HTA/20191127_flxandmian_v3/"

# Read in the data
raw <- easysorter::read_data(dirs)

# Remove all data from the contaminated wells
raw_nocontam<- easysorter::remove_contamination(raw)

# look at strain distributions
raw_nocontam %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::select(strain, condition, EXT, TOF) %>%
  tidyr::gather(trait, value, EXT:TOF) %>%
  ggplot(.) +
  aes(x = value) +
  geom_histogram() +
  facet_grid(strain~trait, scales = "free") +
  theme_bw() +
  labs(title= "DMSO")


# Summarize the data
# directories = FALSE because supplied data comes from one directory
summedraw <- easysorter::sumplate(raw_nocontam, directories = FALSE, quantiles = TRUE, v3_assay = TRUE) 

#Prune based on biological impossibilities
pruned <- easysorter::bioprune(summedraw) %>%
  tidyr::gather(trait, phenotype, -(date:col))%>%
  dplyr::filter(!grepl("red|green|yellow|f.|iqr", trait)) %>%
  dplyr::ungroup()
#  dplyr::mutate(control = ifelse(condition == "FLX-0", "None", "FLX-0")) %>% #mutate changes what you say; this is an if else statement 
#and the usage is "ifelse(test, yes, no)" ,
#so this is saying if condition has cis-0 then change the control to "none", if not mian-0, (=",") change the control to mian-0 %>%
#  tidyr::drop_na(condition)



#look at bleaches
pruned<- pruned %>%
  mutate(assay = case_when(plate %in% c("1","2","3","4","5","6") ~
                             "A", 
                           plate %in% c("7","8","9","10","11","12") ~
                             "B", 
                           plate %in% c("13", "14", "15", "16", "17", "18") ~
                             "C")) 
pruned<- pruned %>%
  dplyr::filter(!(plate =="13")) # it was clogged during  plate #13 , pressure was off for DMSO 

pruned <- pruned %>%
  dplyr::filter(assay == "B") # only going to look at bleach b

#save pruned data
#save("pruned", file =  "~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/HTA/20191127_flxandmian_v3/processed/pruned")

#looking at number and trait
pruned %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.)+
  aes(x= mean.TOF , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") +
  coord_flip()

#change this x axis is wierd 
pruned %>% 
  tidyr::drop_na(condition) %>%
  dplyr::filter(trait == "mean.EXT") %>%
  ggplot(.) +
  aes(x = strain, y = phenotype, fill = strain) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 0.5) +
  geom_boxplot(alpha = 0.5, outlier.color = NA) +
  theme_bw() +
  facet_grid(~condition)
labs(x = "Dose", y = "mean.EXT") +
  scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "ECA377"="grey", "ECA232"="grey", "ECA386"="grey", "ECA1113"="grey" )) 




spread<-pruned %>%
  tidyr::spread(trait, phenotype)
# variability of titering: n / bleach / strain !!!!
ggplot(spread ) +
  aes(x=strain, y=n) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill=assay), width = 0.1, shape=21) +
  scale_color_brewer(palette = "Set1") + facet_grid(.~assay)

prunedout <- prune_outliers(pruned)

traita = "median.EXT"
prunedout%>%
  filter(trait==traita,
         condition=="DMSO",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "CB4856", "ECA377", "ECA232", "ECA386", "ECA1113"))) %>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = assay),width = 0.1)+
  theme_bw(24)+
  ylab(traita)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "right") +
  scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "ECA377"="grey", "ECA232"="grey", "ECA386"="grey", "ECA1113"="grey"))

# only one assay
# assreg <- regress(prunedout, assay = TRUE) #assay reg for bleaches effects 

#regress from DMSO using loop bc i have two conditions
regressed <- NULL
conditionloop <- prunedout %>%
  dplyr::filter(!grepl("DMSO", condition))
conditiondf <- prunedout
for(i in unique(conditionloop$condition)) {
  newdf <- conditiondf %>%
    dplyr::filter(condition %in% c(i,"DMSO") ) ##This section was pulling out both conditions at once so I changed it to only pull the one in the loop and DMSO 
  reg <- easysorter::regress(newdf)
  regressed <- rbind(regressed, reg)
}

chroms <- nil_plot(c("CB4856", "N2" ,"ECA232", "ECA377","ECA1113","ECA386"), chr = "V", all.chr = TRUE) [[1]]

DMSO <- plot_genopheno(prunedout, "DMSO", "median.EXT", "V", back=T, conf=3.9)

flxdose <- plot_genopheno(regressed, "FLX-31.25", "median.EXT", "V", back=T, conf=3.9)

miandose <- plot_genopheno(regressed, "MIAN-25", "median.EXT", "V", back = T, conf = 4.0)

DMSO
flxdose
miandose

traita <- "median.TOF"
medianEXT <- prunedout %>%
  filter(trait == traita & condition == "DMSO")
fit <- aov(phenotype ~ strain, data = medianEXT, projections = TRUE)
x <- as.data.frame(TukeyHSD(fit)$strain)
x$pair <- rownames(x)
dmso <-x %>%
  ggplot() +
  aes(colour=cut(`p adj`, c(0, 0.01, 0.05, 1), 
                 label=c("p<0.01","p<0.05","Non-Sig")))+
  geom_hline(yintercept = 0, lty="11",colour="grey30")+
  geom_errorbar(aes(pair,ymin=lwr,ymax=upr),width=0.3,size=2)+
  geom_point(aes(pair,diff))+
  labs(colour="", title = "DMSO") + 
  theme_classic(18) + 
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1)) +
  base_theme
dmso

medianEXT <- prunedout %>%
  filter(trait == traita & condition == "FLX-31.25")
fit <- aov(phenotype ~ strain, data = medianEXT, projections = TRUE)
x <- as.data.frame(TukeyHSD(fit)$strain)
x$pair <- rownames(x)

flx <-x %>%
  ggplot() +
  aes(colour=cut(`p adj`, c(0, 0.01, 0.05, 1), 
                 label=c("p<0.01","p<0.05","Non-Sig")))+
  geom_hline(yintercept = 0, lty="11",colour="grey30")+
  geom_errorbar(aes(pair,ymin=lwr,ymax=upr),width=0.3,size=2)+
  geom_point(aes(pair,diff))+
  labs(colour="", title= "Fluoxetine 31.25um") + 
  theme_classic(18) + 
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1)) +
  base_theme
flx

medianEXT <- prunedout %>%
  filter(trait == traita & condition == "MIAN-25")
fit <- aov(phenotype ~ strain, data = medianEXT, projections = TRUE)
x <- as.data.frame(TukeyHSD(fit)$strain)
x$pair <- rownames(x)

mian <-x %>%
  ggplot() +
  aes(colour=cut(`p adj`, c(0, 0.01, 0.05, 1), 
                 label=c("p<0.01","p<0.05","Non-Sig")))+
  geom_hline(yintercept = 0, lty="11",colour="grey30")+
  geom_errorbar(aes(pair,ymin=lwr,ymax=upr),width=0.3,size=2)+
  geom_point(aes(pair,diff))+
  labs(colour="", title= "Mianserin 25um") + 
  theme_classic(18) + 
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1)) +
  base_theme
mian

cowplot::plot_grid(DMSO, dmso, flxdose, flx, labels = c('A', 'B', 'C', 'D'), label_size = 10)

cowplot::plot_grid(DMSO, dmso, miandose, mian, labels = c('A', 'B', 'C', 'D'), label_size = 10)


#looking at numbers by plate
traita = "median.EXT"
prunedout%>%   #to look at drugs change df to regressed 
  filter(trait==traita,
         condition=="DMSO",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "CB4856", "ECA377", "ECA232", "ECA386", "ECA1113"))) %>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color =factor(plate)),width = 0.1)+
  theme_bw(24)+
  ylab(traita)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "right") 

#number by trait
D <- pruned %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.)+
  aes(x= median.EXT , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") +
  coord_flip()

FL <- regressed %>%
  dplyr::filter(condition == "FLX-31.25") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.)+
  aes(x= median.EXT , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") +
  coord_flip()

M <-regressed %>%
  dplyr::filter(condition == "MIAN-25") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.)+
  aes(x= median.EXT , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") +
  coord_flip()




cowplot::plot_grid(D, FL, M, labels = c('DMSO', 'FLX', 'MIAN'), label_size = 10)


# try to look at number of 1113 in flx-31.25 by well...


plate <-prunedout %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  dplyr::filter(strain == "ECA1113")

ggplot(plate)+
  aes(x= col , y= n) +
  geom_jitter(aes(color=factor(plate)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") 






