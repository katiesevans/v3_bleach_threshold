# try purr
# ???
library(tidyverse)
library(dplyr)
library(purrr)
library(mclust)
load("~/Dropbox/AndersenLab/LabFolders/Loraina/v3_bleach_threshold/clusterdata.Rda")

percentl1 <- clusterdata %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) %>%
  right_join(clusterdata, by = "bleach") %>%
  dplyr::select(bleach, L1_percent, person) %>%
  dplyr::distinct(L1_percent, bleach, person)

#

b <- "20190926_a" # individual bleach
cl_data <- clusterdata %>%
  dplyr::filter(bleach== b) 
# EMcluster
m <- cl_data %>%
  dplyr::select(logTOF, logEXT, green) %>%
  mclust::Mclust(., G=2)  
plot(m, "classification")
summary(m)  # 40/(40+561) =0.06655574 # check

class <- cl_data %>%
  bind_cols(.,as_tibble(m$classification)) %>% # add back cluster grouping from mclust
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

t <-class %>% # change cluster grouping so group 1 is always group with smallest logTOF
  group_by(value) %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname)
cluster_order <- class  %>% 
  merge(t, by = "value") #add back

badgroup_percent <- cluster_order %>%
  dplyr::group_by(new_g) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::summarise(g1_percent = (sum(cluster_order$new_g == "1") /nrow(cluster_order)*100))%>%
  as.numeric()

df <- data.frame(badgroup_percent=badgroup_percent, bleach=b)


  
#######################################
# this is the structure I want 
nesttest <- clusterdata %>%
  dplyr::select(date, assay, condition, strain, bleach, logTOF, logEXT, green, person) %>%
  dplyr::group_by(bleach) %>%
  tidyr::nest()

# ok now for loop this shit 

# nesttest[[2]][[1]]
# nesttest$data[[1]] same results




nesttest$data[[b]] %>%
  select(green)
b <- 3

#for each bleach 

  
df <-nesttest$data[[b]] 
m <- df %>%
    dplyr::select(logTOF, logEXT, green) %>%
    mclust::Mclust(., G=2)  # each nested df will have its own m clust obj 
df <- df %>%
    bind_cols(.,as_tibble(m$classification)) %>% # add cluster class num
    dplyr::mutate(old_class = value) #%>%

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
    dplyr::summarise(g1_percent = (sum(cluster_order$new_g == "1") /nrow(cluster_order)*100))%>%
    as.numeric()

df2 <- data.frame(badgroup_percent=badgroup_percent, bleach=nesttest$bleach[[b]])
  
test <- percentl1 %>%
    dplyr::filter(bleach == nesttest$bleach[[b]])
final <-  merge(df2,percentl1) # YAS
b<-2





# PERFECT 
new_df <- NULL
  # ok now loop that shit # ok its not actually looping thru
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
    dplyr::summarise(g1_percent = (sum(cluster_order$new_g == "1") /nrow(cluster_order)*100))%>%
    as.numeric()
  df2 <- data.frame(badgroup_percent=badgroup_percent, bleach=nesttest$bleach[[b]])
  test <- percentl1 %>%
    dplyr::filter(bleach == nesttest$bleach[[b]])
  final <-  merge(df2,percentl1)
  new_df <- rbind(new_df, final) # almost 
}
















per_bleach <- function(x){
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
    dplyr::summarise(g1_percent = (sum(cluster_order$new_g == "1") /nrow(cluster_order)*100))%>%
    as.numeric()
  df2 <- data.frame(badgroup_percent=badgroup_percent, bleach=nesttest$bleach[[b]])
  test <- percentl1 %>%
    dplyr::filter(bleach == nesttest$bleach[[b]])
  final <-  merge(df2,percentl1)
  new_df <- rbind(new_df, final)
}




per_bleach(nesttest$data[[6]])


help <- purrr::map(.x = unique(nesttest$bleach[])), .f = per_bleach))




help <- purrr::pmap(.l = nesttest, .f = per_bleach(.))

for (i in 1:length(nesttest$bleach[[6]])) {
  per_bleach()
}



