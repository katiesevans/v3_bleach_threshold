library(tidyverse)
library(dplyr)
library(purrr)
library(mclust)
library(plotly)
load("~/Dropbox/AndersenLab/LabFolders/Loraina/v3_bleach_threshold/clusterdata.Rda")

percentl1 <- clusterdata %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) %>%
  dplyr::right_join(clusterdata, by = "bleach") %>%
  dplyr::select(bleach, L1_percent, person) %>%
  dplyr::distinct(L1_percent, bleach, person)

# this is the structure I want 
nesttest <- clusterdata %>%
  dplyr::select(date, assay, condition, strain, bleach, logTOF, logEXT, green, person) %>%
  dplyr::group_by(bleach) %>%
  tidyr::nest()

# PERFECT 
new_df <- NULL

for (b in 1:length(unique(nesttest$bleach[]))) {
  df <-nesttest$data[[b]] 
  m <- df %>%
    dplyr::select(logTOF, logEXT, green) %>%
    mclust::Mclust(., G=2)  # each nested df will have its own m clust obj 
  df <- df %>%
    bind_cols(.,as_tibble(m$classification)) %>% # add cluster class num
    dplyr::mutate(old_class = value) 
  dumb <- df %>% 
    dplyr::group_by(old_class) %>% 
    dplyr::summarise(mean=mean(logTOF)) %>% # find the mean of TOF per group 
    dplyr::arrange(mean) %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(new_g=rowname) %>%
    dplyr::select(!rowname)
  df <- df %>%
    merge(dumb) %>%
    dplyr::select(!c(old_class, value))
  badgroup_percent <- df %>% 
    dplyr::group_by(new_g) %>% 
    dplyr::summarise(n = n()) %>%
    dplyr::summarise(g1_percent = (sum(df$new_g == "1") /nrow(df)*100))%>%
    as.numeric()
  df2 <- data.frame(badgroup_percent=badgroup_percent, bleach=nesttest$bleach[[b]]) 
  test <- percentl1 %>%
    dplyr::filter(bleach == nesttest$bleach[[b]])
  final <-  merge(df2,percentl1)
  new_df <- rbind(new_df, final) 
}

view(new_df)
# SAVE THIS DATA  


plot <- ggplot2::ggplot(new_df) +
  aes(x=L1_percent, y=badgroup_percent, color=person, bleach = bleach) +
  ggplot2::geom_point() +
  ggplot2::theme_bw()
  
  
ggplotly(plot, tooltip = c("bleach","badgroup_percent","L1_percent","person"))
  
colnames(new_df)
  
cluster_cluster <- new_df %>%
  dplyr::select(L1_percent, badgroup_percent) %>%
  mclust::Mclust(., G=3) # pulls out high L1 percent ? are these bad assays in other ways than the high bad group percent ?
summary(cluster_cluster) 
plot(cluster_cluster, "classification")

cluster_cluster_df <- new_df %>%
  bind_cols(.,as_tibble(cluster_cluster$classification)) # add cluster class num
   
  
cluster_cluster_person <- new_df %>%
  dplyr::select(L1_percent, badgroup_percent, person) %>% # mm not helpful?
  mclust::Mclust(., G=2) 

summary(cluster_cluster_person)  
plot(cluster_cluster_person, "classification")
plot(cluster_cluster_person, "uncertainty")

cluster_cluster_person_df <- new_df %>%
  bind_cols(.,as_tibble(cluster_cluster_person$classification)) #

