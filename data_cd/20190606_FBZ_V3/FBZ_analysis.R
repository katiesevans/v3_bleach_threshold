library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)

###Read in day 1 and analyze
dir <- c("~/Desktop/20190606_FBZ_V3/")

raw <- read_data(dir)

raw_noncontam <- remove_contamination(raw)

summedplate <- sumplate(raw_noncontam, directories = FALSE, quantiles = TRUE)

summedplate_1 <- summedplate %>%
  mutate(norm.n = 3)
biopruned <- bioprune(summedplate_1)

out_pruned <- prune_outliers(biopruned)

out_1 <- out_pruned %>%
  filter(plate %in% x)
out_2 <- out_pruned %>%
  filter(plate %in% y)

full <- regress(out_pruned)
full_1 <- regress(out_1)
full_2 <- regress(out_2)
cols <- c("N2" = "orange", "920" = "red", "1081" = "blue", "1075" = "green", "1325" = "brown", "1327" = "yellow","N2" = "orange", "919" = "red", "1082" = "blue", "1076" = "green", "1326" = "brown", "1328" = "yellow","1139"="black","1137"="pink","1317"="gray","1097"="maroon")
x <- c(1,2,3,4,5,6,7,8,9,10)
y <- c(11,12,13,14,15,16,17,18,19,20)
z <- c("N2","1325","1327")
full_ele %>%
  filter(trait=="median.EXT",
         plate %in% y,
         #strain %in% z,
         !is.na(strain)) %>%
  #dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "919", "1082","1327", "1325","1075","1137","882","1139", "1317","1097")))%>%
  ggplot()+
  aes(x=strain, y = phenotype)+
  geom_boxplot(aes(fill=strain),outlier.shape = NA)+
  #scale_x_discrete(labels=c("DMSO" = "DMSO", "7_5" = "7.5µM",
  #"15" = "15µM", "30" = "30µM", "60" = "60µM", "120" = "120µM"))+
  theme_bw(24)+
  geom_jitter(width = 0.1)+
  #scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size=25),
        axis.title.y=element_text(size=30))
        #legend.position = "none")

