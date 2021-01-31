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

# all V3 data from Katie, Clay, Deahan, myself 
# plot
all_data %>%
  tidyr::gather(trait, phenotype, TOF:EXT) %>%
  ggplot(.) +
  aes(x = phenotype) +
  geom_histogram(bins = 150) +
  theme_bw() +
  facet_grid(~trait, scales = "free")

all_data %>%
  tidyr::gather(trait, phenotype, TOF:EXT) %>%
  ggplot(.) +
  aes(x = phenotype, fill=stage) + # color by life stage
  geom_histogram(bins = 150) +
  theme_bw()+
  theme(axis.text = element_text(size = 10)) +
  facet_grid(person~trait, scales = "free")


mean <-all_data %>%
  dplyr::group_by(strain, person) %>%
  summarise(mean.TOF = mean(TOF))

# use life stages from sorter to figure out percentages of L1's in data 
# find percents 

# TOF: 230 units is late L4

# easysorter stages:
### L1: 60 > TOF < 90
### L2/L3: 90 > TOF < 200
### L4: 200 > TOF < 300
### Adult: TOF > 300


# stage by joy's data - grown with live liquid bacteria
### L1: 60 > TOF < 150
### L2/L3: 150 > TOF < 300
### L4: 300 > TOF < 500
### Adult: TOF > 500

ls_data <- all_data %>% dplyr::filter(person== "LS")
ls_data %>% dplyr::count(stage == "L1")

nrow(ls_data)
# i have total of 7.7% L1's in my data 
  
# testing for mine
percent <- ls_data %>%
  dplyr::filter(person == "LS") %>%
  dplyr::select(condition:bleach) %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(percent = (sum(stage == "L1")/n())) # this makes it a new df with one column with percents
              # per bleach 

# what is the percent of L1's in the data? 
percent_alldata <- all_data %>%
  dplyr::group_by(person,bleach) %>%
  dplyr::summarise(percent = (sum(stage == "L1")/n())*100) # full numbs on plot

percent_alldata %>%
  ggplot(.) +
  aes(x=percent) +
  geom_histogram() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
         axis.text.y = element_text(size=12) ) +
  labs(title = "Distribution of percent of L1's") +
  facet_wrap(~person)


unique(all_data$date)
unique(all_data$bleach)

############### me going through code and thinking about stuff 
# test, starting with assay from 20191114 -- most consistent % of L1's in 3 bleaches

df <- all_data %>%
  dplyr::filter(date == 20191114)

# var explained just strain differences
fit <- lm(TOF ~ strain, data = df )
anova(fit)
summary(fit)
full_var <- sum(anova(fit)[[2]]) #takes anova fit and adds sum sq -- which is full variance sum of element 2, full phenotypic variation

# first variable
# strain variance = 1.3 %
anova(fit)[[2]][1]/full_var*100   # sum sq of strain anova(fit)[[2]][1] = [[2]] is total strain variatin * strain var / total var * 100 


fit2 <- lm(TOF ~ strain + bleach, data = df ) #strain and bleach 
anova(fit2)
summary(fit2)
full_var2 <- sum(anova(fit2)[[2]])

# strain variance : 
anova(fit2)[[2]][1]/full_var2*100

# bleach variance (second variable) :
anova(fit2)[[2]][2]/full_var*100
# is 0.029 * remember i chose this assay bc all percentages of L1's were near 7% # consistent 

# now choose a more variable bleach 


# assay from 20201105
b <- all_data %>%
  dplyr::filter(date == 20201105)

fit3 <- lm(TOF ~ strain + bleach, data = b)
anova(fit3)
summary(fit3)
full_var3 <- sum(anova(fit3)[[2]])

#strain variance 
anova(fit3)[[2]][2]/full_var3*100 
# bleach var :
anova(fit3)[[2]][1]/full_var3*100

##

ls_grouped <- ls_data  %>%
  dplyr::group_by(date) 

# make a for-loop to see differences in strain var to bleach var for all assays

exp_var <- NULL # make an empty 
# d is my variable --> first value in the vector, the 2nd, 3rd etc
for (d in unique(ls_grouped$date)){ # tell it where "d in dates"   
  fit <-lm(TOF ~ strain + bleach, data = dplyr::filter(ls_grouped, date == d))
  full_var <- sum(anova(fit)[[2]])
  strain_var <-anova(fit)[[2]][2]/full_var*100 # variables cant be same as column name which i tried 
  bleach_var <- anova(fit)[[2]][1]/full_var*100
  df <- data.frame(strain_variation = strain_var, bleach_variation=bleach_var) #naming col names by variable name
  exp_var <- rbind(exp_var, df)
}

# exp_var is for all of my experiments 

# need to remove assays that dont have multiple bleaches
# all data 

# 20181122  20181123 20181220 adtes from DL are an issue
# 20190308


all_grouped <- all_data %>%
  dplyr::filter(!person == "DL") %>%
  dplyr::filter(!date %in% c("20190308", "20190315", "20190613", "20190722", "20190723", "20190822", "20190823", "20190606", "20190607",
                             "20190617")) %>%
  dplyr::filter(!is.na(assay)) %>%
  dplyr::filter(!assay %in% c("del", "dose", "tags")) %>%
   dplyr::group_by(date) 

all_expvar <- NULL 
for(d in unique(all_grouped$date)){
  fit <- lm(TOF ~ strain + bleach, data = dplyr::filter(all_grouped, date == d))
  full_var <- sum(anova(fit)[[2]])
  strain_var <-anova(fit)[[2]][2]/full_var*100 # variables cant be same as column name which i tried 
  bleach_var <- anova(fit)[[2]][1]/full_var*100
  df <- data.frame(strain_variation = strain_var, bleach_variation=bleach_var, date = d) #naming col names by variable name
  all_expvar <- rbind( d, all_expvar, df)
}

# cant figure out how to not have it doubled oh well 
all_expvar <- all_expvar %>%
  dplyr::slice(23:44)

# 20190308 20190315 20190613 20190722 20190723 20190822 20190823 20190606
[1] "20190308" "20190315" "20190613" "20190722" "20190723" "20190822" "20190823" "20190926"
[9] "20190927" "20190930" "20191001" "20191031" "20191101" "20191107" "20191108" "20200116"
[17] "20200117" "20200227" "20200228" "20190606" "20190607" "20190617" "20190618" "20190812"
[25] "20190813" "20200130" "20200131" "20200930" "20191003" "20191114" "20200214" "20201105"
[33] "20201106" "20201112" "20201113"
# had to get rid of SO MANY DATES that dont have multiple bleaches 

d <- "20190926"

fit <- lm(TOF ~ strain + bleach, data = dplyr::filter(all_grouped, date == d))
anova(fit)
summary(fit)
full_var <- sum(anova(fit)[[2]])
#strain variance 
anova(fit)[[2]][2]/full_var*100 
# bleach var :
anova(fit)[[2]][1]/full_var*100

 # from joy makes a df with % var explained 
aov(TOF~strain + bleach, data = dplyr::filter(all_grouped, date == d)) %>%
      summary() %>% .[[1]]  %>%
      dplyr::mutate(Terms = rownames(.), `% Var Explained` = round((`Sum Sq`/sum(`Sum Sq`))*100,2)) %>% # ' keeps % var ex as a phrase 
      dplyr::select(Terms, everything())


date <- data.frame(unique(all_grouped$date))

all <- cbind(date,all_expvar) %>%
  rename(date = unique.all_grouped.date.)

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

### 2021 01 15 subgroup notes 
# 2d 2 parameter and 
# 3D plots --> cluster of arrested L1's maybe separate 
# EM 

dmso <- all_grouped %>%
  dplyr::ungroup() %>%
  dplyr::filter(condition == "DMSO")

dmso %>%
  ggplot(.) +
  aes(x=TOF, y=EXT, color = factor(stage)) + 
  geom_point()

LS <- dmso %>% # just my data 
  dplyr::filter(person == "LS") 

ggplot(LS) +
  aes(x=TOF, y=EXT, color = factor(stage)) + 
  geom_point()

a <- LS %>%
  dplyr::filter(bleach == "20201112_A") %>%
  ggplot(.) +
  aes(x=TOF, y=EXT, color = factor(stage)) + 
  geom_point()

l1 <- all_grouped %>%
  dplyr::ungroup() %>%
  dplyr::filter(stage == "L1")

ggplot(l1) + 
  aes(x=TOF, y=EXT, color = factor(person)) + 
  geom_point() 


##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
### 2021 01 26 talking to Joy about mclust 
# i want to clsuter by 

dmso <- all_grouped %>%
  dplyr::ungroup() %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::filter(person == "LS")

dmso %>%
  ggplot(.) +
  aes(x = TOF, y= EXT) +
  geom_point(alpha = .1, size = .1) # start to see groups inthe data 

# using mclust package , 

m <- dmso %>%
  dplyr::select(TOF, EXT) %>% # the columns I want as axis
  mclust::Mclust(., G=2) # G=groups, opens prompts thru consol

summary(m)
plot(m) 


log <- dmso %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT)) %>%
  dplyr::select(logTOF, logEXT)  %>% # taking log separates better 
  mclust::Mclust(., G=2) # im telling the number of clusers, make sure it makes sense 
  
summary(log)


plot(log) # choose 2 for classification, clustering in elipses 


class <- dmso %>%
  bind_cols(.,as_tibble(log$classification))  # from the model fit, bind back to df

class %>% 
  ggplot(.) +
  aes(x= TOF) +
  geom_histogram() + # value is the num if in cluster 1 or 2
  facet_grid(~value) # check by eye how the data separates out 


#### whatever gets separated
# cant "remove" out from drug data


# just looking at Katie's data 
ks <- all_grouped %>%
  dplyr::ungroup() %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::filter(person == "KSE") %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT))  # using log gives higher percent 

ks_m <- ks %>%
  dplyr::select(logTOF, logEXT) %>%
  mclust::Mclust(., G=2)

ks_class <- ks %>%
  bind_cols(.,as_tibble(ks_m$classification)) 

ks_class %>%
  ggplot(.) +
  aes(x=logTOF) +
  geom_histogram() +
  facet_grid(~value)


km <- ks %>% # change to log?
   dplyr::select(TOF, EXT) %>% # the columns I want as axis
   mclust::Mclust(., G=2)
   
summary(km) # Katie's data 5% is in cluster 2 , at 48 hours we should not have L1's 



cd <- all_grouped %>%
  dplyr::ungroup() %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::filter(person == "CD")
cdm <- cd %>%
  dplyr::select(TOF, EXT) %>%
  mclust::Mclust(., G=2)
summary(cdm) # 4% in cluster 2 


ls <- dmso %>%
  dplyr::filter(bleach == "20201112_A") %>%
  
  dplyr::select(TOF, EXT) %>% 
  mclust::Mclust(., G=2)

summary(ls)


####
# from Joy, she makes a df with % var explained 
aov(TOF~strain + bleach, data = dplyr::filter(all_grouped, date == d)) %>%
  summary() %>% .[[1]]  %>%
  dplyr::mutate(Terms = rownames(.), `% Var Explained` = round((`Sum Sq`/sum(`Sum Sq`))*100,2)) %>% # ' keeps % var ex as a phrase 
  dplyr::select(Terms, everything())

##########################################################################################################
##########################################################################################################
##########################################################################################################
# try mclus analysis on bleaches that katies assay that the parent trait switched 
# in control -- is she seeing more L1's ? if not ... then how do we tell ? 
# use green or yellow?
# 20190926_b
# 20190927_b are the dates where katies parents flipped
  
  
flip_6 <- all_grouped %>% 
  dplyr::ungroup() %>%
  dplyr::filter(bleach == "20190926_b") %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT))

f6_m <- flip_6 %>%
  dplyr::select(logTOF, logEXT) %>%
  mclust::Mclust(., G=2)
summary(f6_m)

f6class <- flip_6 %>%
  bind_cols(.,as_tibble(f6_m$classification)) 
  
ggplot(f6class) +
  aes(x=logTOF) +
  geom_histogram() +
  facet_grid(~value) # not what we expected, not a lot of animals in clust 2 (smaller objs)

plot(f6_m)


flip_7 <- all_grouped %>%
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

summary(f7_m)
plot(f7_m)
