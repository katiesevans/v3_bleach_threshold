library(easysorter)
library(COPASutils)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mclust)
library(factoextra)



# read in data, filter to keep only controls
wd <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_kse/"
filelist <- list.files(wd)

kse_data <- NULL
for(file in filelist) {
  data <- easysorter::read_data(glue::glue("{wd}/{file}")) %>%
    dplyr::filter(contamination == FALSE,
                  call50 == "object", 
                  condition %in% c("DMSO", "water", "lyaste")) %>%
    dplyr::select(date:col, condition, strain, TOF, EXT, green,stage) %>%
    dplyr::mutate(person = "KSE")
  kse_data <- rbind(kse_data, data)
}
# mark unique assays
kse_data <- kse_data %>%
  na.omit() %>%
  dplyr::mutate(bleach = paste0(date, "_", assay)) %>%
  dplyr::filter(strain == "N2") %>%
  dplyr::filter(assay != "del")
save(kse_data, "~/")

kse_l1 <- kse_data %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) 

### # from 20210128_ls.R code --> added green to the df 
LS_data <- LS_data %>%
  dplyr::filter(strain == "N2")
ls_l1 <- LS_data %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) 

clusterdata <- rbind(data_kse, data_ls)
save(clusterdata, file="/Users/loraina/Dropbox/AndersenLab/LabFolders/Loraina/v3_bleach_threshold/clusterdata.Rda")


##############################################################################

# G = 3
# --> maybe three clusters captures the L4-adult variability
# (thinking about how mean.TOF shifts between bleaches from katies lab meeting slide), 
#  this would only work if like the 2 clusters that are not "smaller objs: "overlap" are more similar than the 3rd
# thinking about my time my data separated out 45%, 40%, 15%.....

data_ls <- LS_data %>%
  dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT)) 

three_ls <- data %>% 
  dplyr::select(logTOF, logEXT, green)  %>% # taking log separates better 
  mclust::Mclust(., G=3)


summary(three_ls)
# Clustering table: 
#    1    2    3 
# 4915 1467 4251  
4251/(1467 + 4915 + 4251) 
# 13.8% 46.2.% 39.9 % ....# kind of like this 3d plot to see 
# majority in 86% 
plot(three_ls)

three_class <- data %>% 
  bind_cols(.,as_tibble(three_ls$classification)) # add back in classification number
# figure out how to group the data in a way that everytime "group 1" is the "small objects group" by
# attaching the class value to the smallest mean(tof) 
# thats some group, order, mutate shit 
# it'll help with visualization - consistent colors and "small" obejcts group 

# 3D, color by group, axes are traits 
three_class <- three_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)
                
test<- three_class %>%
  dplyr::summarise(mean=mean(TOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to
  
three_plot <- three_class  %>%
  merge(test, by = "value")  #  yes!  triple checked lol

colors <- as.factor(three_plot$new_g) # new grouping class, smallest object will always be group 1
# clean up code from above ???

plotly::plot_ly(three_plot, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")



#G=2
two_ls <- data %>%
  dplyr::select(logTOF, logEXT, green)  %>%  
  mclust::Mclust(., G=2)

summary(two_ls)
#   Clustering table:
#         1    2 
#       8631 2002
#   8631/(8631 + 2002)
#   18.8 %      81.2%

plot(two_ls, "classification")

two_class <- data %>% 
  bind_cols(.,as_tibble(two_ls$classification)) 

two_test <- two_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value) %>%
  dplyr::summarise(mean=mean(TOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

two_plot <- two_class  %>%
  merge(two_test, by = "value")  #  yes!  triple checked lol

colors <- as.factor(two_plot$new_g) # new grouping class, smallest object will always be group 1
# clean up code from above ???

plotly::plot_ly(two_plot, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("black","cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")



#####
data_kse <- kse_data %>%
  dplyr::mutate(green,logTOF = log(TOF), logEXT = log(EXT))


kse_log <- kse_data %>%
  dplyr::mutate(green,logTOF = log(TOF), logEXT = log(EXT)) %>% # sorter doesnt cut off stages using EXT 
  dplyr::select(logTOF, logEXT, green)  %>% # taking log separates better 
  mclust::Mclust(., G=3)

summary(kse_log)


plot(kse_log)

5321/(5321+ 1220 +4516)

kse_log2 <- kse_data %>%
  dplyr::mutate(green,logTOF = log(TOF), logEXT = log(EXT))# sorter doesnt cut off stages using EXT 
kse_log2m<- kse_log2 %>%
  dplyr::select(logTOF, logEXT, green)  %>% # taking log separates better 
  mclust::Mclust(., G=2)

summary(kse_log2m)
# 1    2 
# 1456 9601 
# 1456/(1456+96201) 15%

plot(kse_log2m, "classification")

kse_log2_class <- kse_log2 %>% 
  bind_cols(.,as_tibble(kse_log2m$classification))
# 3D, color by group, axes are traits 
kse_log2_class <- kse_log2_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

t <- kse_log2_class %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

kse_log2_plot <- kse_log2_class  %>%
  merge(t, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
kse_log2_plot %>% group_by(new_g) %>% summarise(n = n())

plotly::plot_ly(kse_log2_plot, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = ~new_g, 
                alpha = 0.5, size = .5, type="scatter3d")
# 85% of kse all data in L4-young adult 







###
# 3d plots, do the thing where I pull out classification value, rbind back in 
# color by class value, axis are parameters 


# more parameters the better _ L4's should cluster together across traits 
# i could - tof, ext, green, yellow
# 2 groups --> define which one is 

