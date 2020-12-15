library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)
library(mclust)

###Read in day 1 and analyze
dir <- c("~/Desktop/benzimidazole_assays/20190606_FBZ_V3/")

raw <- read_data(dir)

raw_noncontam <- remove_contamination(raw)
raw_nongfp <- raw_noncontam[[1]]%>%
  mutate(row.col = paste(row,col))

mylist <-matrix(0,ncol=27,nrow = 0)

mylist <-data.frame(mylist)

colnames(mylist) <- c("date","experiment","round","assay","plate","row","col","sort","TOF","EXT","time","green","yellow","red","norm.EXT","norm.green","norm.yellow","norm.red","object","call50","stage","strain","condition","control","contamination","row.col","fit$classification")

for (x in unique(raw_nongfp$strain)){
  raw_nongfp_filt1<-raw_nongfp %>%
    filter(strain == x)
  for (con in unique(raw_nongfp_filt1$condition)){
    print(con)
    if (is.na(con)){
      next}else{
        TOF_DMSO_1 <- raw_nongfp_filt1 %>%
          filter(condition==con)#%>%
         # filter(EXT > 300)%>%
         # filter(TOF < 550)%>%
          #filter(TOF > 85)
        #if(x == "1326" | x == "1137" | x =="1097" | x == "1314"){
         # TOF_DMSO_1 <- TOF_DMSO_1 %>%
         #   filter(EXT > 300)
      #  }
        TOF_DMSO <- TOF_DMSO_1%>%
          filter(TOF < (mean(TOF_DMSO_1$TOF) + 2*sd(TOF_DMSO_1$TOF))) %>%
          filter(TOF > (mean(TOF_DMSO_1$TOF) - 2*sd(TOF_DMSO_1$TOF)))
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
  filter(strain == "1325")
N2_DMSO_filtered <- mylist %>%
  filter(condition == "DMSO")%>%
  filter(strain == "1325")
hist(N2_DMSO$TOF,probability = TRUE,breaks = 100)
d <- density(N2_DMSO_filtered$TOF)
t <- density(N2_DMSO$TOF)
lines(t,col="black")
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
  mutate(norm.n = 3)
biopruned <- bioprune(summedplate_1)

out_pruned <- prune_outliers(biopruned)
out_para <- out_pruned %>%
  filter(strain %in% parasite)
out_full <- out_pruned %>%
  filter(strain %in% parasite | strain %in% elegans)
save(out_full, file= "FBZ_V3_fig.RData")
out_ele <- out_pruned %>%
  filter(strain %in% elegans)
save(out_pruned, file = "FBZ_V3_full.RData")
out_1 <- biopruned %>%
  filter(plate %in% plates_1)

save(out_11, file="FBZ_V3_plates1.RData")
out_11 <- prune_outliers(out_1)

out_2 <- biopruned %>%
  filter(plate %in% plates_2)
out_22 <- prune_outliers(out_2)
save(out_22, file="FBZ_V3_plates2.RData")

full <- regress(out_pruned)
full_fig <- regress(out_full)
full_para <- regress(out_para)
full_ele <- regress(out_ele)
full_1 <- regress(out_11)
full_2 <- regress(out_22)
cols <- c("N2" = "orange", "920" = "red", "1081" = "blue", "1075" = "green","1314"="grey", "1325" = "brown","1138"="pink", "1098" = "maroon", "1327" = "yellow","N2" = "orange", "919" = "red", "1082" = "blue", "1076" = "green","1075" = "green", "1326" = "brown", "1328" = "yellow","1139"="black","1140"="black","1137"="pink","1317"="gray","1097"="maroon")
x <- c("N2", "882", "1137", "1140", "1097","1317")
y <- c("N2", "882", "920", "1325","1327","1081","1075")
z <- c("N2","1325","1327")
plates_1 <- c(1,2,3,4,5,6,7,8,9,10)
plates_2 <- c(11,12,13,14,15,16,17,18,19,20)
full_22 <- full %>%
  mutate(plates = case_when())
elegans <- c("N2", "882","1137", "1138","1139", "1140", "1097","1098","1314","1317")
parasite <- c("N2","882","1325","1326","1327","1328","1081","1082","920","1075","1076")

full_ele %>%
  filter(trait=="median.TOF",
         #condition=="DMSO",
         #strain %in% parasite,
         #plate %in% plates_1,
         !is.na(strain))%>%
         #strain %in% parasite) %>%
  #condition=="DMSO") %>%
  #strain==c("N2","882","1325","1327","1082","1097","919")) %>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328","1326", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(aes(color=factor(plate)),width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  ylab("Optical Density")+
  #ylim(-100,50)+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1137"="Q131L","1076"="F167Y", "1326"="E198L","1328"="E198V","1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("full_DMSOmeanEXT_fbz.pdf", device = "pdf", plot = last_plot(), width = 13, height = 8, units = "in")
