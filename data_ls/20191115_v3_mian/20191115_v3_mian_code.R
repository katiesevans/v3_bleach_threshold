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
source("~/Dropbox/AndersenLab/LabFolders/Loraina/Scripts/Base_theme.R")



#Define a vector of your experiement directories
dirs <- "~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/HTA/20191115_v3_mian/"

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
                             "B"))

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

ggplot(subset(a, !strain == "NA" & condition == "MIAN-0"), 
       aes(x=log(TOF), y=log(EXT))) +
  geom_jitter(aes(colour=assay), size=0.5, alpha=0.5) + 
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~strain)

pruned_mian %>% 
  tidyr::drop_na(condition) %>%
  dplyr::filter(trait == "mean.EXT") %>%
  ggplot(.) +
  aes(x = factor(condition, levels = c("MIAN-0", "MIAN-25", "MIAN-50", "MIAN-100", "MIAN-200")), y = phenotype, fill = strain) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 0.5) +
  geom_boxplot(alpha = 0.5, outlier.color = NA) +
  theme_bw() +
  labs(x = "Dose", y = "mean.EXT") +
  scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "ECA377"="grey", "ECA232"="grey", "ECA386"="grey", "ECA1113"="grey" )) 

pruned_mian %>%
  dplyr::filter(condition == "MIAN-0") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.) +
  aes(x= mean.EXT , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") +
  #scale_color_manual(values=c("B"="black", "C"="turquoise")) +
  coord_flip()

prunedout_mian <- prune_outliers(pruned_mian)

assreg_mian <- regress(prunedout_mian, assay = TRUE) #bleach assay regress

# variability of titering: n / bleach / strain !!!!
spread_mian <-pruned_mian %>%
  tidyr::spread(trait, phenotype)
ggplot(spread_mian) +
  aes(x=strain, y=n) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill=assay), width = 0.1, shape=21) +
  scale_color_brewer(palette = "Set1") + facet_grid(.~assay)

#this looks at trait by numbers 
prunedout_mian %>%
  dplyr::filter(condition == "MIAN-0") %>%
  dplyr::select(date:phenotype) %>%
  tidyr::spread (trait, phenotype) %>%
  ggplot(.)+
  aes(x= mean.EXT , y= n) +
  geom_jitter(aes(color=factor(strain)), width = 0.1, size = 2) +
  facet_grid(condition~assay, scales = "free") +
  coord_flip()

assreg_mian %>%
  dplyr::filter(condition == "MIAN-0" & trait == "mean.EXT") %>%
  ggplot(.)+
  aes(x= strain , y = phenotype) +
  geom_jitter(aes(color=factor(assay)), width = 0.1, size = 2) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  labs(title= "DMSO") +
  coord_flip() 

# aonly<- assreg_mian %>% #cb is wierd bc only has one point in assay b  
  dplyr::filter(! (strain=="CB4856" & assay=="B")) %>%
  dplyr::filter(condition == "MIAN-0" & trait == "mean.EXT") %>%
  ggplot(.)+
  aes(x= strain , y = phenotype) +
  geom_jitter(aes(color=factor(assay)), width = 0.1, size = 2) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  labs(title= "DMSO", ylab= "mean.EXT") +
  coord_flip() 
# aonly
#this is just messing around, kinda ignore?? ECA1113 is sick in dmso but other nils are fine 

# regress for each dose do this one
regressed <- NULL
conditionloop <- assreg_mian %>%
  dplyr::filter(!grepl("-0", condition))
conditiondf <- assreg_mian  #number doesnt start with zero
for(c in unique(conditionloop$condition)) {
  newdf <- conditiondf %>%
    dplyr::filter(condition %in% c(c, "MIAN-0") )#filters for c= flx doses and flx-0
  # dplyr::mutate(condition = drug) %>%
  reg <- easysorter::regress(newdf)
  regressed <- rbind(regressed, reg)
}


chroms <- nil_plot(c("N2" , "CB4856", "ECA377", "ECA232", "ECA386", "ECA230"), chr = "V", all.chr = TRUE) [[1]]


control <- plot_genopheno(assreg_mian, "MIAN-0", "median.EXT", "V", back=T, conf=4.0)
drug <- plot_genopheno(regressed, "MIAN-25", "median.EXT", "V", back=T, conf=4.0)




#looking at differences in dmso ! 
traita <- "median.EXT"
out_pruned_mian <- prune_outliers(assreg_mian)
medianEXT <- out_pruned_mian %>%
  filter(trait == traita & condition == "MIAN-0")
fit <- aov(phenotype ~ strain, data = medianEXT, projections = TRUE)
x <- as.data.frame(TukeyHSD(fit)$strain)
x$pair <- rownames(x)
dmso_stats <-x %>%
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




#ggsave(filename = "___regressed.png", plot = last_plot(), device = "png")


#diff trait stats on dose 
plot_genopheno(regressed, "MIAN-25", "median.TOF", "V", back=T, conf=4.0)

traita <- "median.TOF"
out_pruned_mian <- prune_outliers(regressed)
trait <- out_pruned_mian %>%
  filter(trait == traita & condition == "MIAN-25")
fit <- aov(phenotype ~ strain, data = trait, projections = TRUE)
x <- as.data.frame(TukeyHSD(fit)$strain)
x$pair <- rownames(x)
drugstats <-x %>%
  ggplot()+
  aes(colour=cut(`p adj`, c(0, 0.01, 0.05, 1), 
                 label=c("p<0.01","p<0.05","Non-Sig")))+
  geom_hline(yintercept = 0, lty="11",colour="grey30")+
  geom_errorbar(aes(pair,ymin=lwr,ymax=upr),width=0.3,size=2)+
  geom_point(aes(pair,diff))+
  labs(colour="", title= "Mianserin 25um")+
  theme_classic(18) +
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1)) +
  base_theme


cowplot::plot_grid(control, dmso_stats, labels = c('A', 'B'), label_size = 12)
cowplot::plot_grid(control, dmso_stats, drug, drugstats, labels = c('A', 'B', 'C', 'D' ), label_size = 12)




####
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/20180829/data/mianserin-GWER.chromosomal.annotated.Rda")
mianserin <-annotatedmap %>%
  filter(chr== "V" & marker=="V_4035606") 
library(cegwas)

mian_snps <- cegwas::snpeff("V:3742283-4074054",severity = "ALL",elements = "ALL") #for the marker of V_4035606
#CB is resistant in mianserin 
#for whole region
mian_high <- mian_snps %>%
  filter(impact == "HIGH", GT == "ALT", strain == "CB4856") 
highgenes <- as.tibble(unique(mian_high$gene_name))
mian_mod <- mian_snps %>%
  filter(impact == "MODERATE", GT == "ALT", strain == "CB4856") 
modgenes <- as.tibble(unique(mian_mod$gene_name))

####
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/findNILs.R")
# find nils 
findNILs("V:3839036-4034158")
nilgeno <- nilgeno %>%
  dplyr::filter(sample %in% c("ECA386", "ECA377") & chrom == "V")



#for break points between 377 and 386 
nil_region <- cegwas::snpeff("V:3120168-3955396	",severity = "ALL",elements = "ALL")  

high <- nil_region %>%
  filter(impact == "HIGH", GT == "ALT", strain == "CB4856") 
highgenes <- as.tibble(unique(high$gene_name))
mod <- nil_region %>%
  filter(impact == "MODERATE", GT == "ALT", strain == "CB4856") 
modgenes <- as.tibble(unique(mod$gene_name))

