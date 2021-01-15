library(easysorter)
library(COPASutils)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dplyr)
library(cowplot)

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
  df <- data.frame(strain_variation = strain_var, bleach_variation=bleach_var) #naming col names by variable name
  all_expvar <- rbind(all_expvar, df)
}

# 20190308 20190315 20190613 20190722 20190723 20190822 20190823 20190606
[1] "20190308" "20190315" "20190613" "20190722" "20190723" "20190822" "20190823" "20190926"
[9] "20190927" "20190930" "20191001" "20191031" "20191101" "20191107" "20191108" "20200116"
[17] "20200117" "20200227" "20200228" "20190606" "20190607" "20190617" "20190618" "20190812"
[25] "20190813" "20200130" "20200131" "20200930" "20191003" "20191114" "20200214" "20201105"
[33] "20201106" "20201112" "20201113"
# had to get rid of SO MANY DATES that dont have multiple bleaches 

d <- "20190618"
fit <- lm(TOF ~ strain + bleach, data = dplyr::filter(all_grouped, date == d))


date <- data.frame(unique(all_grouped$date))

all <- cbind(date,all_expvar) %>%
  rename(date = unique.all_grouped.date.)

### 2021 01 15 subgroup notes 
# 2d 2 parameter and 
# 3D plots --> cluster of arrested L1's maybe separate 
# EM 
