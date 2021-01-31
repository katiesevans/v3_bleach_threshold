
#####################################################
#     trying to figure out V3 thresholds for subgroup 
#     2021 01 28 
#####################################################
library(easysorter)
library(COPASutils)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mclust)


source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_genotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_phenotype_plots.R")
source("~/Dropbox/AndersenLab/LabFolders/Loraina/Scripts/Base_theme.R")

load("~/Dropbox/AndersenLab/LabFolders/Loraina/v3_bleach_threshold/all_data.Rda")
# clean up all data to only include data with 3 bleaches, and all DMSO control

all_dmso <- all_data %>%
  dplyr::filter(!person == "DL") %>%
  dplyr::filter(!date %in% c("20190308", "20190315", "20190613", "20190722", "20190723", "20190822", "20190823", "20190606", "20190607",
                             "20190617")) %>%
  dplyr::filter(!is.na(assay)) %>%
  dplyr::filter(!assay %in% c("del", "dose", "tags"))



ls_data <- all_dmso %>% 
  dplyr::filter(person== "LS")

ls_l1 <- ls_data %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) 
  

ls_data_l1 <- ls_data %>% 
  dplyr::count(stage == "L1")
nrow(ls_data_l1)
# i have total of 7.7% L1's in my data 
5141 / ( 61045 + 5141) 
# =  0.07767504 == 7.76 %

# what is the percent of L1's in the data? 
percent_alldata <- all_dmso %>%
  dplyr::group_by(person , bleach) %>%
  dplyr::summarise(percent = (sum(stage == "L1")/n())*100) 


# if i take the dates of my own data with highest and lowest % of L1 -- what does clustering show ?
# bleach 20191003_B	has L1 % of 3.693114
# strain var : 1.8804096283072	
# bleach var: 0.37678515633937	
load("~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/HTA/20191003_mianserin_v3/raw_nocontam_mian.Rda")
raw_nocontam_mian <- raw_nocontam_mian %>%
  dplyr::filter(condition == "MIAN-0")

low_log <- raw_nocontam_mian %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT)) %>% # sorter doesnt cut off stages using EXT 
  dplyr::select(logTOF, logEXT)  %>% # taking log separates better 
  mclust::Mclust(., G=2)
summary(low_log)
#Clustering table:
# 1    2 
# 5191  782 
plot(low_log)
782 / (5191+782) = 0.1309225 # 13% of data in "smaller objects class"

low_class <- raw_nocontam_mian %>% # the classes fit the log(EXT) and log(TOF) into two classifications
  bind_cols(.,as_tibble(low_log$classification))  # from the model fit, bind back to df

l <- low_class %>% 
  ggplot(.) +
  aes(x= log(TOF), fill=factor(value)) +
  geom_histogram(bins = 50)  + # value is the num if in cluster 1 or 2
  theme_bw(15) +
  scale_fill_manual(values = c("red", "grey")) +
  theme(legend.position = "none ") 
  #facet_grid(~value)  

# high is high % of L1s and a bad bleach... 
# bleach 20201105_C	has L1 % of 16.322089 
# strain var: 6.62160617811502	
# bleach var: 1.14510518298996
load("~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/HTA/20201105_flx/raw_nocontam_flx.Rda")
raw_nocontam <- raw_nocontam %>%
  dplyr::filter(condition == "DMSO")

high_log <-  raw_nocontam  %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT)) %>% # sorter doesnt cut off stages using EXT 
  dplyr::select(logTOF, logEXT)  %>% # taking log separates better 
  mclust::Mclust(., G=2)

summary(high_log) # should be more in "true worm" 
#   Clustering table:
#   1     2 
#   12602  3281

plot(high_log) # choose #2 for classification 

3281 / (12602 + 3281)  0.2065731 # 20% of data in "smaller objects class"


high_class <- raw_nocontam %>% # the classes fit the log(EXT) and log(TOF) into two classifications
  bind_cols(.,as_tibble(high_log$classification))  # from the model fit, bind back to df

h <- high_class %>% 
  ggplot(.) +
  aes(x= log(TOF), fill=factor(value)) +
  geom_histogram(bins = 50)  + # value is the num if in cluster 1 or 2
  theme_bw(15) +
  scale_fill_manual(name = "Cluster", labels = c("small", "large"), values = c("red", "grey")) # check by eye how the data separates out 

cowplot::plot_grid(l,h, scale = 1)



#####

# variance explained by bleach / assay (assay bc grouped by date)
# in  all_dmso
dmso_date <- all_dmso %>%
  group_by(date)

all_expvar_date <- NULL 
for(d in unique(dmso_date$date)){ # date
  fit <- lm(TOF ~ strain + bleach, data = dplyr::filter(dmso_date, date == d))
  full_var <- sum(anova(fit)[[2]])
  strain_var <-anova(fit)[[2]][2]/full_var*100 # variables cant be same as column name which i tried 
  bleach_var <- anova(fit)[[2]][1]/full_var*100
  df <- data.frame( date = d, strain_variation = strain_var, bleach_variation=bleach_var) #naming col names by variable name
  all_expvar_date <- rbind(all_expvar_date, df)
}

################################################################################################################
# trying to use cluster to tell me something 
## using log bc thats what joy does, spreads/normalizes data (?)

all_dmso %>%
  tidyr::gather(trait, phenotype, TOF:EXT) %>%
  ggplot(.) +
  aes(x = phenotype, fill=stage) + # color by life stage
  geom_histogram(bins = 150) +
  theme_bw()+
  theme(axis.text = element_text(size = 10)) +
  facet_grid(person~trait, scales = "free")

all_dmso %>%
  tidyr::gather(trait, phenotype, TOF:EXT) %>%
  dplyr::filter(trait == "TOF") %>%
  ggplot(.) +
  aes(x = phenotype, fill=stage) + # color by life stage
  geom_histogram(bins = 150) +
  theme_bw()+
  theme(axis.text = element_text(size = 10)) # all data

all_dmso %>%
  ggplot(.) +
  aes(x = TOF, y= EXT, color=factor(stage)) + 
  geom_point() +
  theme_bw()+
  theme(axis.text = element_text(size = 10))
  



log <- all_dmso %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT)) %>% # sorter doesnt cut off stages using EXT 
  dplyr::select(logTOF, logEXT)  %>% # taking log separates better 
  mclust::Mclust(., G=2) # im telling the number of clusers, make sure it makes sense 

summary(log) # should be more in "true worm" 


plot(log) # choose #2 for classification 

class <- all_dmso %>% # the classes fit the log(EXT) and log(TOF) into two classifications
  bind_cols(.,as_tibble(log$classification))  # from the model fit, bind back to df

class %>% 
  ggplot(.) +
  aes(x= TOF) +
  geom_histogram() + # value is the num if in cluster 1 or 2
  facet_grid(~value)  # check by eye how the data separates out 

class %>%
  ggplot(.)+
  aes(x=log(TOF), y=log(EXT), color = stage)+
  geom_jitter() +
  facet_grid(~value)

two <- class %>%
  dplyr::filter(value == "2") 
unique(two$stage) # all stages in cluster 2, using elipses, sorter only has TOF cutoffs 
####
# try mclust, select log tof or ext, will have density plot
log2 <- all_dmso %>%
  dplyr::mutate(logTOF = log(TOF)) %>% # only looking at TOF, which is what sorter does, sorter has TOF cutoffs
  dplyr::select(logTOF)  %>% # taking log separates better 
  mclust::Mclust(., G=2) # JOY says 2 variables is better , the more the better lol
summary(log2)
plot(log2)
# bottom line is all data , blue is one part, red is other classified cluster
# plot density using 4
# 

class2 <- all_dmso %>%
  bind_cols(.,as_tibble(log2$classification))
class2 %>%
  ggplot(.)+
  aes(x=log(TOF), y=log(EXT), color = stage)+
  geom_jitter() +
  facet_grid(~value) 



#################################################################################################################
# found note from past subgroup about trying to cluster EM with L1's
# this tells me there isnt enough of a difference in what the soter calls "L1s" at 48 hours

l1_dmso <- all_dmso %>%
  dplyr::filter(stage == "L1") 

l1_m <- l1_dmso %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT)) %>%
  dplyr::select(logTOF, logEXT)  %>% # taking log separates better 
  mclust::Mclust(., G=2) # im telling the number of clusers, make sure it makes sense 
# played with using TOF/EXT and logTOF/logEXT

summary(l1_m)

plot(l1_m) # choose #2 for classification 
# not much diff in length for mclust to separate, but theres a diff in EXT... 

# when choosing 3 uncertainity -- bigger the dot, the more uncertain 

l1_dmso <- l1_dmso %>% 
  bind_cols(.,as_tibble(l1_m$classification))

l1_dmso %>% 
  ggplot(.) +
  aes(x= TOF) +
  geom_histogram() + # value is the num if in cluster 1 or 2
  facet_grid(~value)

# is looking at L1s at 48 houirs usefulll??? 
# NO  theyre shouldnt be L1s at 48 hours

# arrest L1's amount shouldnt cahnge from 24 hours to 48

# what traits at 48 hours is useful -- ?! 
# fraction of L1s is not correlated with expected results OR percent variance explained by bleach. 
# goal: define a bad bleach and determine data for assay should be thrown out
# clustering to get a better idea of fraction of "small objects "
# the clusters we see are L4's and smaller things
# --> how to use this to set thresh.... look at ratio? 
# what  about bad bleaches do we expect ? 


## just logEXT?
log3 <- all_dmso %>%
  dplyr::mutate(logEXT = log(EXT)) %>% # only looking at TOF, which is what sorter does, sorter has TOF cutoffs
  dplyr::select(logEXT)  %>% # taking log separates better 
  mclust::Mclust(., G=2) # JOY says 2 variables is better , the more the better lol
summary(log3)
plot(log3)

###############################################################################################################################################################
#####################################################

# try mclus analysis on bleaches that katies assay that the parent trait switched 
# in control -- is she seeing more L1's ? if not ... then how do we tell ? 
# use green or yellow?
# 20190926_b
# 20190927_b are the dates where katies parents flipped


flip_6 <- all_dmso %>% 
  dplyr::filter(bleach == "20190926_b") %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT))

f6_m <- flip_6 %>%
  dplyr::select(logTOF, logEXT) %>%
  mclust::Mclust(., G=2)
summary(f6_m) # 5777 vs 175 # of objects

175/ (5777 + 175)

f6class <- flip_6 %>%
  bind_cols(.,as_tibble(f6_m$classification)) 

# not what we expected, not a lot of animals in clust 2 (smaller objs)
six <- ggplot(f6class) +
  aes(x= log(TOF), fill=factor(value)) +
  geom_histogram(bins = 50)  + # value is the num if in cluster 1 or 2
  theme_bw(15) +
  scale_fill_manual(values = c("purple", "grey")) +
  theme(legend.position = "none ") 

plot(f6_m)


flip_7 <- all_dmso %>%
  dplyr::ungroup() %>%
  dplyr::filter(bleach =="20190927_b") %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT))

f7_m  <- flip_7 %>%
  dplyr::select(logTOF, logEXT) %>%
  mclust::Mclust(., G=2)

f7_class <-  flip_7 %>%
  bind_cols(.,as_tibble(f7_m$classification))

ggplot(f7_class) +
  aes(x=logTOF) +
  geom_histogram() +
  facet_grid(~value)
seven <- ggplot(f7_class) +
  aes(x= log(TOF), fill=factor(value)) +
  geom_histogram(bins = 50)  + # value is the num if in cluster 1 or 2
  theme_bw(15) +
  scale_fill_manual(values = c("purple", "grey")) +
  theme(legend.position = "none ") 


summary(f7_m) #4696 vs 388 
388 / (4700 + 384)
plot(f7_m) 


cowplot::plot_grid(six, seven)


## cleaned up a little i think from above, looking at my own data 
ls_data


# assay is in a b c
df3 <- ls_data %>%
  dplyr::filter(date == 20191003) %>%
  dplyr::ungroup()

# var explained 
fitA <- lm(TOF ~ strain + bleach, data = df3 ) # strain and bleach
anova(fitA)
summary(fitA)
full_varA <- sum(anova(fitA)[[2]]) #takes anova fit and adds sum sq -- which is full variance sum of element 2, full phenotypic variation

# first variable
# strain variance = 0.3767852 %
anova(fitA)[[2]][1]/full_varA*100   # sum sq of strain anova(fit)[[2]][1] = [[2]] is total strain variatin * strain var / total var * 100 
# bleach variance (second variable) :
anova(fitA)[[2]][2]/full_varA*100 # 1.88041





dfB <- ls_data %>%
  dplyr::filter(date == 20201105) %>%
  dplyr::ungroup()


# var explained 
fit <- lm(TOF ~ strain + bleach, data = dfB ) # strain and bleach
anova(fit)
summary(fit)
full_var <- sum(anova(fit)[[2]]) #takes anova fit and adds sum sq -- which is full variance sum of element 2, full phenotypic variation

# first variable
# strain variance = 1.145105 %
anova(fit)[[2]][1]/full_var*100   # sum sq of strain anova(fit)[[2]][1] = [[2]] is total strain variatin * strain var / total var * 100 
# bleach variance (second variable) :
anova(fit)[[2]][2]/full_var*100 # 6.621606







fit2 <- lm(TOF ~ strain + bleach, data = df ) #strain and bleach 
anova(fit2)
summary(fit2)
full_var2 <- sum(anova(fit2)[[2]])

# strain variance : 
anova(fit2)[[2]][1]/full_var2*100

# bleach variance (second variable) :
anova(fit2)[[2]][2]/full_var*100
# is 0.029 * remember i chose this assay bc all percentages of L1's were near 7% # consistent 







