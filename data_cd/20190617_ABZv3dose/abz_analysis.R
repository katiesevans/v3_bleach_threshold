library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)

###Read in day 1 and analyze
dir <- c("~/Desktop/benzimidazole_assays/20190617_ABZv3dose/")

raw <- read_data(dir)

raw_noncontam <- remove_contamination(raw)

raw_nongfp <- raw_noncontam[[1]]%>%
  mutate(row.col = paste(row,col))
mylist <-NULL

for (x in unique(raw_nongfp$strain)){
  raw_nongfp_filt1<-raw_nongfp %>%
    filter(strain == x)
  for (con in unique(raw_nongfp_filt1$condition)){
    print(con)
    if (is.na(con)){
      next}else{
        TOF_DMSO_1 <- raw_nongfp_filt1 %>%
          filter(condition==con)%>%
          filter(norm.green < 0.2)%>%
          filter(TOF > 85)
        #if(x == "1326" | x == "1137" | x =="1097" | x == "1314"){
        # TOF_DMSO_1 <- TOF_DMSO_1 %>%
        #   filter(EXT > 300)
        #  }
        TOF_DMSO <- TOF_DMSO_1%>%
          #filter(TOF < (mean(TOF_DMSO_1$TOF) + 2*sd(TOF_DMSO_1$TOF))) %>%
          #filter(TOF > (mean(TOF_DMSO_1$TOF) - 2*sd(TOF_DMSO_1$TOF)))
          filter(EXT < (mean(TOF_DMSO_1$EXT) + 2*sd(TOF_DMSO_1$EXT))) %>%
          filter(EXT > (mean(TOF_DMSO_1$EXT) - 2*sd(TOF_DMSO_1$EXT)))
        #fit <- Mclust(TOF_DMSO[,c(10:11,13:15)],G=2)
        #if (!is.null(fit)){
        # TOF_FULL <- cbind(TOF_DMSO, fit$classification)
        #TOF_1 <- TOF_FULL %>%
        # filter(`fit$classification`==1)
        #TOF_2 <- TOF_FULL %>%
        # filter(`fit$classification`==2)
        #if(nrow(TOF_2) > nrow(TOF_1)){
        mylist <- rbind(mylist,TOF_DMSO)
        #}else{
        # mylist <- rbind(mylist,TOF_1)
      }
  }}
N2_DMSO <- raw_nongfp %>%
  filter(condition == "DMSO")%>%
  filter(strain == "N2")
N2_DMSO_filtered <- mylist %>%
  filter(condition == "DMSO")%>%
  filter(strain == "N2")
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

out_pruned <- prune_outliers(biopruned)

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

full <- out_pruned %>%
  filter(condition == "DMSO" | condition == "50")%>%
  mutate(control= case_when(condition == "50" ~ "DMSO", condition=="DMSO" ~ "None"))
save(full, file="ABZ_50uM.RData")
out_para_100 <- out_pruned %>%
  filter(strain %in% p)%>%
  filter(condition == "100" | condition =="DMSO")%>%
  mutate(control= case_when(condition == "100" ~ "DMSO", condition=="DMSO" ~ "None"))
full_100 <- regress(out_para_100)
f <- c("N2","882","919","920")
cols <- c("N2" = "orange","882"="grey","919"="red", "920" = "red", "1081" = "blue", "1075" = "green", "1325" = "brown","1326"="brown", "1327" = "yellow","N2" = "orange", "919" = "red", "1082" = "blue", "1076" = "green", "1326" = "brown", "1328" = "yellow","1139"="black","1140"="black","1137"="pink","1138"="pink","1314"="tan","1317"="tan","1097"="maroon","1098"="maroon")
x <- c("N2","882","919","1325","1327","1081","1097","1075","1137","1314","1139")
x_p <- c("N2","882","919","1325","1327","1081","1075")
y <-c("N2", "882","920","1326","1328","1082","1098","1076","1138","1317","1140")
y_p <- c("N2", "882","920", "1326","1328","1082","1076")
p <- c("N2", "882","919","920", "1325","1326","1327","1328","1082","1081","1075","1076")
cut <- c("1314","1139","1137","1075")
sub_01_mut <- subtract_dr_to_pc %>%
  #dplyr::filter(condition != "200")%>%
  dplyr::select(strain, condition, median.EXT,plate)%>%
  dplyr::mutate( med_norm = (median.EXT - min(median.EXT)) / (max(median.EXT) - min(median.EXT)))


subtract_dr_to_pc%>%
  filter(#condition!="200",
         #condition =="DMSO",
         #trait == "median.EXT",
         !(strain %in% cut),
         #strain %in% x,
         !is.na(strain))%>%
  #strain %in% parasite) %>%
  #condition=="DMSO") %>%
  #strain==c("N2","882","1325","1327","1082","1097","919")) %>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = median.TOF)+
  #geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  #geom_jitter(aes(color=strain),width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  ylab("Animal Length")+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'loess')+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25uM","12_5" = "12.5uM","25"="25uM","50"="50uM","100"="100uM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=25),
        axis.title.y=element_text(size=30),
        legend.position = "None")


save("subtract_dr_to_pc", )

ggsave("parasites_ABZmedTOF.png", device = "png", plot = last_plot(), width = 13, height = 8, units = "in") 

