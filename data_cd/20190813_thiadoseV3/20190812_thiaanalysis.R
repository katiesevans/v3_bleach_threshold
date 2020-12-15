library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)

# colors
axis_color <- "#000F08"
highlight_color <- "#D7263D"
background_color <- "white"

# font
number_font <- "Helvetica"
axes_text_size <- 20
axes_title_font <- "Helvetica"
axes_title_size <- 20
title_size <- 20

base_theme <- theme(
  line = element_line(colour = axis_color, size = 0.5, linetype = 1), 
  rect = element_rect(fill = background_color, colour = axis_color, size = 0.5, linetype = 1), 
  text = element_text(family = axes_title_font, size = axes_text_size), 
  
  axis.text = element_text(family = number_font,size = rel(0.8), colour = "grey30", margin = unit(0.1, "cm")),
  axis.text.x = element_text(vjust = 1), 
  axis.text.y = element_text(hjust = 1), 
  axis.ticks = element_line(colour = "gray90"), 
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90), 
  axis.ticks.length = unit(0.15, "cm"),
  
  strip.text = element_text(size = rel(0.8)), 
  strip.background = element_rect(fill = background_color, colour = NA, size = 0.5, linetype = 1),
  
  plot.background = element_rect(fill = background_color, color = NA),
  
  legend.background = element_rect(fill=background_color,colour = background_color), 
  legend.spacing = unit(0.2, "cm"), 
  legend.key = element_rect(fill = background_color, colour = NA), 
  legend.key.size = unit(1, "lines"), 
  legend.key.height = NULL, 
  legend.key.width = NULL, 
  legend.text = element_text(size = rel(0.6)), 
  legend.text.align = NULL, 
  legend.title = element_text(size = rel(0.6), hjust = 0), 
  legend.title.align = NULL, 
  legend.position = "right", 
  legend.direction = NULL, 
  legend.justification = "center", 
  legend.box = NULL, 
  
  panel.background = element_rect(fill = background_color, colour = NA),  
  panel.grid.major = element_line(colour = "gray90"), 
  panel.grid.minor = element_blank(), 
  panel.spacing = unit(1, "lines"), 
  panel.margin.x = NULL, 
  panel.margin.y = NULL)


###Read in day 1 and analyze
dir <- c("~/20190812_thiadoseV3/")

raw <- read_data(dir)

raw_noncontam <- remove_contamination(raw)

raw_nongfp <- raw_noncontam[[1]]%>%
  mutate(row.col = paste(row,col))

summedraw <- sumplate(raw_nongfp,picked_plates = TRUE, directories = FALSE, quantiles = TRUE)

summedraw1 <- summedraw %>%
 mutate(norm.n = 3)
biopruned <- bioprune(summedraw1)

out_pruned <- prune_outliers(biopruned)


mylist <-NULL
for (x in unique(raw_nongfp$strain)){
  raw_nongfp_filt1<-raw_nongfp %>%
    filter(strain == x)
  for (con in unique(raw_nongfp_filt1$condition)){
    print(con)
    if (is.na(con)){
      next}else{
        TOF_DMSO_1 <- raw_nongfp_filt1 %>%
          filter(condition==con)
        TOF_DMSO <- TOF_DMSO_1%>%
          filter(EXT < (mean(TOF_DMSO_1$EXT) + 2*sd(TOF_DMSO_1$EXT))) %>%
          filter(EXT > (mean(TOF_DMSO_1$EXT) - 2*sd(TOF_DMSO_1$EXT))) %>%
          filter(green < 125) %>%
          filter(TOF > 100)
        mylist <- rbind(mylist,TOF_DMSO)
      }
  }}
N2_DMSO <- raw_nongfp %>%
  filter(condition == "DMSO")%>%
  filter(strain == "238")
N2_DMSO_filtered <- mylist %>%
  filter(condition == "DMSO")%>%
  filter(strain == "238")
hist(N2_DMSO$TOF,probability = TRUE,breaks = 50)
d <- density(N2_DMSO_filtered$TOF)
t <- density(N2_DMSO$TOF)
lines(t,col="green")
lines(d, col="blue")
mean_noF <- mean(N2_DMSO$TOF)
med_noF <- median(N2_DMSO$TOF)
mean_F <- mean(N2_DMSO_filtered$TOF)
med_F <- median(N2_DMSO_filtered$TOF)
abline(v=mean_noF, col="green")
abline(v=med_noF,col="green")
abline(v=med_F,col="red")
abline(v=mean_F,col="red")
summedplate <- sumplate(mylist,picked_plates = TRUE, directories = FALSE, quantiles = TRUE)

summedplate_1 <- summedplate %>%
  mutate(norm.n=3)
biopruned <- bioprune(summedplate_1)

biopruned_pr <- biopruned %>%
  dplyr::ungroup() %>%
  dplyr::select(-(date:assay)) %>%
  dplyr::filter(!is.na(condition)) %>%
  tidyr::gather(trait, value, -(plate:col))

control_dr <- biopruned_pr %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::group_by(strain, trait) %>%
  dplyr::summarise(control_value = mean(value, na.rm = T))

# RUN PCA ON ENTIRE EXPERIMENT
subtract_dr_to_pc <- biopruned_pr %>%
  dplyr::left_join(., control_dr, by = c("strain", "trait")) %>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(strain, condition, trait, delta_control, plate, row, col)%>%
  dplyr::mutate(u_strain = paste(strain, condition, plate,row,col,sep="_"))%>%
  tidyr::spread(trait, delta_control) %>%
  dplyr::select(-contains("green"),-contains("yellow"), -contains("red"),-contains("cv"),-contains("iqr"), -contains("var"))


subtract_dr_to_pc%>%
  filter(#trait == "median.EXT",
         #condition == "DMSO",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "CB4856", "238", "239")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","31_25","62_5","125","","250")))%>%
  ggplot()+
  aes(x=con_fix, y = q25.EXT)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  theme_bw(24)+
  #geom_jitter(aes(fill=fancy_strain), width =0.1)+
  ylab("Thiabendazole Response")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "31_25" = "31.25uM","62_5" = "62.5uM","125"="125uM","250"="250uM"))+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"))
