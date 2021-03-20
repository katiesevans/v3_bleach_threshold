# this is looking at individual bleaches and clustering  PER SELECTED BLEACH
load("/Users/loraina/Dropbox/AndersenLab/LabFolders/Loraina/v3_bleach_threshold/clusterdata.Rda")

percentl1 <- clusterdata %>%
  dplyr::group_by(bleach) %>%
  dplyr::summarise(L1_percent = (sum(stage == "L1")/n())*100) %>%
  right_join(clusterdata, by = "bleach") %>%
  dplyr::select(bleach, L1_percent, person) %>%
  dplyr::distinct(L1_percent, bleach, person)


# 20191003_C	1.0638298	LS
# G=3
low_ls <- clusterdata %>%
  dplyr::filter(bleach == "20191003_C" & person == "LS") %>%
  dplyr::select(green, logTOF, logEXT) 
low_ls_m <- mclust::Mclust(low_ls, G=3)
summary(low_ls_m)
# Clustering table:
# 1  2  3 
# 38  5 51 
plot(low_ls_m, "classification")

low_ls_class <- low_ls %>% 
  bind_cols(.,as_tibble(low_ls_m$classification))
# 3D, color by group, axes are traits 
low_ls_class <- low_ls_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test<- low_ls_class %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

low_lss_plot <- low_ls_class  %>%
  merge(test, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
low_lss_plot %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         5 = 5%
# 2        38 = 40%
# 3        51 = 54%

51/(38 + 5 + 51)

colors <- as.factor(low_lss_plot$new_g) # new grouping class, smallest object will always be group 1
# clean up code from above ???

plotly::plot_ly(low_lss_plot, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = 1, type="scatter3d")



##############################################
# diff date high percent L1's' 14% 
# G=3
high_ls2 <- clusterdata %>%
  dplyr::filter(bleach == "20201105_C" & person == "LS") %>% # high percent of L1's
  dplyr::select(green, logTOF, logEXT) 
high_ls_m2 <- mclust::Mclust(high_ls2, G=3)
summary(high_ls_m2)
# Clustering table:
# 1   2   3 
# 55 129 424  # did it twice, got same results....

plot(high_ls_m2, "classification")

high_ls_class2 <- high_ls2 %>% 
  bind_cols(.,as_tibble(high_ls_m2$classification))
# 3D, color by group, axes are traits 
high_ls_class2 <- high_ls_class2 %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test2 <- high_ls_class2 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

high_ls_plot2 <- high_ls_class2  %>%
  merge(test2, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
high_ls_plot2 %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         129 = 21% %
# 2        55 = 9%
# 3        424 = 70%%

424/(55 + 129 + 424 )

colors <- as.factor(high_ls_plot2$new_g) # new grouping class, smallest object will always be group 1
# clean up code from above ???

plotly::plot_ly(high_ls_plot2, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")


##############################################
# low percent L1 date , different date bc first date didnt have a lot of data
# do diff day G= 2 and 3
# 20191003_B 2% L1's
low_ls3 <- clusterdata %>%
  dplyr::filter(bleach == "20191003_B" & person == "LS") %>% # high percent of L1's
  dplyr::select(green, logTOF, logEXT) 
low_ls_m3 <- mclust::Mclust(low_ls3, G=3)
summary(low_ls_m3)
# Clustering table:
#   1   2   3 
#   193  67 395 

plot(low_ls_m3, "classification")

low_ls_class3 <- low_ls3 %>% 
  bind_cols(.,as_tibble(low_ls_m3$classification))
# 3D, color by group, axes are traits 
low_ls_class3 <- low_ls_class3 %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test2 <- low_ls_class3 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

low_lss_plot3 <- low_ls_class3  %>%
  merge(test2, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
low_lss_plot3 %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         67 = 10% %
# 2        193 = 30%
# 3        395 = 60%%

395/(193 +  67  + 395 )

colors <- as.factor(low_lss_plot3$new_g) # new grouping class, smallest object will always be group 1
# clean up code from above ???

plotly::plot_ly(low_lss_plot3, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")
# in this the red group .... could be arrested L1's and ... carcas??

# do same bleach but G = 2  this bleach has 2% L1's
low_ls_m4 <- mclust::Mclust(low_ls3, G=2)
summary(low_ls_m4)
# Clustering table:
#   1   2   
#   589 66 

plot(low_ls_m4, "classification")

low_ls_class4 <- low_ls3 %>% 
  bind_cols(.,as_tibble(low_ls_m4$classification))
# 3D, color by group, axes are traits 
low_ls_class4 <- low_ls_class4 %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test3 <- low_ls_class4 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

low_lss_plot4 <- low_ls_class4  %>%
  merge(test3, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
low_lss_plot4 %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         66 = 10% %
# 2        589 = 90%

# 10 % of data in small objects --> this data must be "consistant" bc this happens in 2 or 3 cluster.....

589/(66+  + 589 )

colors <- as.factor(low_lss_plot4$new_g) # new grouping class, smallest object will always be group 1


plotly::plot_ly(low_lss_plot4, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")
#######

# diff date high percent L1's' 14% 
# G=2
high_ls_m3 <- mclust::Mclust(high_ls2, G=2)
summary(high_ls_m3)
# Clustering table:
# 1   2    
# 445 163  

plot(high_ls_m3, "classification")

high_ls_class3 <- high_ls2 %>% 
  bind_cols(.,as_tibble(high_ls_m3$classification))
# 3D, color by group, axes are traits 
high_ls_class3 <- high_ls_class3 %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test4 <- high_ls_class3 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

high_ls_plot3 <- high_ls_class3  %>%
  merge(test4, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
high_ls_plot3 %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         163 = 26% %
# 2        445 = 73%%

445/(163 + 445)

colors <- as.factor(high_ls_plot3$new_g) # new grouping class, smallest object will always be group 1
# clean up code from above ???

plotly::plot_ly(high_ls_plot3, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")


##############
#### NOW TRY KATIES DATA by bleach and then later on try cluster by assay... ???

############## 20200227_a	 11% L1's
# G = 3
high_kse <- clusterdata %>%
  dplyr::filter(bleach == "20200227_a" & person == "KSE") %>%
  dplyr::select(green, logTOF, logEXT) 
high_kse_m <- mclust::Mclust(high_kse, G = 3)
summary(high_kse_m)
# Clustering table:
# 1  2  3 
# 91  51 73 


plot(high_kse_m, "classification")

high_kse_class <- high_kse %>% 
  bind_cols(.,as_tibble(high_kse_m$classification))
# 3D, color by group, axes are traits 
high_kse_class <- high_kse_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test5 <- high_kse_class %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

high_kse_plot <- high_kse_class  %>%
  merge(test, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
high_kse_plot %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         51 = 23%
# 2        92 = 42%
# 3        73 = 33%

73/(51 + 92 + 73)

colors <- as.factor(high_kse_plot$new_g) # new grouping class, smallest object will always be group 1

plotly::plot_ly(high_kse_plot, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")
######
### same bleach , 11 % L1, G-2
high_kse_m2 <- mclust::Mclust(high_kse, G = 2)
summary(high_kse_m2)
# Clustering table:
# 1  2  
# 168 48

plot(high_kse_m2, "classification")

high_kse_class2 <- high_kse %>% 
  bind_cols(.,as_tibble(high_kse_m2$classification))
# 3D, color by group, axes are traits 
high_kse_class2 <- high_kse_class2 %>%
  dplyr::group_by(value) %>% # p
  dplyr::mutate(old_class = value)

test6 <- high_kse_class2 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) 

high_kse_plot2 <- high_kse_class2  %>%
  merge(test6, by = "value")  # 

# "updated" cluster table 
high_kse_plot2 %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         48 = 22%
# 2       168 = 78%


168/(48 + 168)

colors <- as.factor(high_kse_plot2$new_g) 


plotly::plot_ly(high_kse_plot2, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")

########################################################
# now katies data low L1 percent... 
# 20191031_b 1% L1 G=3
low_kse <- clusterdata %>%
  dplyr::filter(bleach == "20191031_b" & person == "KSE") %>%
  dplyr::select(green, logTOF, logEXT) 
low_kse_m <- mclust::Mclust(low_kse, G = 3)
summary(low_kse_m)
# Clustering table:
# 1  2  3 
# 137  122 25 


plot(low_kse_m, "classification")

low_kse_class <- low_kse %>% 
  bind_cols(.,as_tibble(low_kse_m$classification))
# 3D, color by group, axes are traits 
low_kse_class <- low_kse_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test7 <- low_kse_class %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

low_kse_plot <- low_kse_class  %>%
  merge(test7, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
low_kse_plot %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         25 = 0.08%
# 2        122 = 43%
# 3        137 = 48%

137/(25 + 122 + 137)

colors <- as.factor(low_kse_plot$new_g) # new grouping class, smallest object will always be group 1

plotly::plot_ly(low_kse_plot, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")
##############3
# G= 2, 1% L1

low_kse_m2 <- mclust::Mclust(low_kse, G = 2)
summary(low_kse_m2)
# Clustering table:
# 1  2   
# 159  125
plot(low_kse_m2, "classification")

low_kse_class2 <- low_kse %>% 
  bind_cols(.,as_tibble(low_kse_m2$classification))
# 3D, color by group, axes are traits 
low_kse_class2 <- low_kse_class2 %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

test8 <- low_kse_class2 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

low_kse_plot2 <- low_kse_class2  %>%
  merge(test8, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
low_kse_plot2 %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         125 = 44%
# 2        159 = 56%
159/(125 + 159)

colors <- as.factor(low_kse_plot2$new_g) # new grouping class, smallest object will always be group 1

plotly::plot_ly(low_kse_plot2, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")
# tbh this doesnt make sense that almost half qualifies as "small objects" 
# ok so new conspiracy theory : with some combo of L1 % and seeing if data fits 2 or 3 cluster better....
# now do it on an assay, color 3d plot by bleach 

## my assay from 20201105:
#   20201105_A	2.068558	
#   20201105_B	6.875000	
#   20201105_C	13.980263

ls_20201105 <- clusterdata %>%
  dplyr::filter(date == "20201105" & person == "LS") %>% # high percent of L1's
  dplyr::select(green, logTOF, logEXT) # if i add assay here it will use assay to find cluster.... not what i want, already structured data that will mess with cluster?
ls_20201105_m3 <- mclust::Mclust(ls_20201105, G=3)
summary(ls_20201105_m3)
# Clustering table:
#   1   2   3 
#   626  370 1784 

plot(ls_20201105_m3, "classification")

ls_20201105_class <- ls_20201105 %>% 
  bind_cols(.,as_tibble(ls_20201105_m3$classification))
# 3D, color by group, axes are traits 
ls_20201105_class <- ls_20201105_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

t <- ls_20201105_class %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

ls_20201105_plot <- ls_20201105_class  %>%
  merge(t, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
ls_20201105_plot %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         261 = 9%
# 2        355 = 13%
# 3        2164 = 78

2164/(261 +  355  + 2164 )

nov <- clusterdata %>%
  dplyr::filter(date == "20201105" & person == "LS") %>%
  dplyr::select(assay,green, logTOF, logEXT )

poo <- merge(nov, ls_20201105_plot) %>%
  group_by(assay)           # number of objs in each pair group
poo %>% filter(assay=="A" & new_g=="1") #116 
poo %>% filter(assay=="A" & new_g=="2") # 414
poo %>% filter(assay=="A" & new_g=="3") # 1,127
poo %>% filter(assay=="B" & new_g=="1") # 70
poo %>% filter(assay=="B" & new_g=="2") # 77
poo %>% filter(assay=="B" & new_g=="3") # 303
poo %>% filter(assay=="C" & new_g=="1") #156
poo %>% filter(assay=="C" & new_g=="2") # 100
poo %>% filter(assay=="C" & new_g=="3")# 324


poo$new_g <- as.factor(poo$new_g)
assay <- factor(poo$assay, levels = c("A", "B", "C"))

plotly::plot_ly(poo, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = ~new_g, 
                alpha = 0.5, size = .5, type="scatter3d")
# cant figure out how to change color and shape based on assay and cluster group
#ugh 
### g=2
ls_20201105_m2 <- mclust::Mclust(ls_20201105, G=2)
summary(ls_20201105_m2)
# Clustering table:
#   1   2   
# 529   2251

plot(ls_20201105_m2, "classification")

ls_20201105_class2 <- ls_20201105 %>% 
  bind_cols(.,as_tibble(ls_20201105_m2$classification))
# 3D, color by group, axes are traits 
ls_20201105_class2 <- ls_20201105_class2 %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

t <- ls_20201105_class2 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

ls_20201105_plot2 <- ls_20201105_class2  %>%
  merge(t, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
ls_20201105_plot2 %>% group_by(new_g) %>% summarise(n = n())
plotly::plot_ly(ls_20201105_plot2, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = c("red", "blue", "cyan"), 
                alpha = 0.5, size = .5, type="scatter3d")

# 531/(531+2249) 20% is in small objects -- or not L4/young adults 

########3
# kse assay where L1% is variable 20200116 #condition is water 
# 20200116_a	3.4965035	KSE
#	20200116_b	1.6483516	KSE
#	20200116_c	10.4060914	KSE

kse_20200116 <- clusterdata %>%
  dplyr::filter(date == "20200116" & person == "KSE") %>% # high percent of L1's
  dplyr::select(green, logTOF, logEXT) 
kse_20200116_m3 <- mclust::Mclust(kse_20200116, G=3)
summary(kse_20200116_m3)
# Clustering table:
#   1   2   3 
#   673  527 169 

plot(kse_20200116_m3, "classification")

kse_20200116_class <- kse_20200116 %>% 
  bind_cols(.,as_tibble(kse_20200116_m3$classification))
# 3D, color by group, axes are traits 
kse_20200116_class <- kse_20200116_class %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

t <- kse_20200116_class %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

kse_20200116_plot <- kse_20200116_class  %>%
  merge(t, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
kse_20200116_plot %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         169 = 12% %
# 2        527 = 38%
# 3        673 = 49%% 

169/(673  +527 +169 )

plotly::plot_ly(kse_20200116_plot, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = ~new_g, 
                alpha = 0.5, size = .5, type="scatter3d")

 ##############
# 2 clusters
kse_20200116_m2 <- mclust::Mclust(kse_20200116, G=2)
summary(kse_20200116_m2)
# Clustering table:
#   1   2  
# 1204  165

plot(kse_20200116_m2, "classification")

kse_20200116_class2 <- kse_20200116 %>% 
  bind_cols(.,as_tibble(kse_20200116_m2$classification))
# 3D, color by group, axes are traits 
kse_20200116_class2 <- kse_20200116_class2 %>%
  dplyr::group_by(value) %>% # find the mean of TOF per group
  dplyr::mutate(old_class = value)

t <- kse_20200116_class2 %>%
  dplyr::summarise(mean=mean(logTOF)) %>%
  dplyr::arrange(mean) %>%
  tibble::rownames_to_column() %>%
  dplyr::mutate(new_g=rowname) %>%
  dplyr::select(!rowname) # pretty sure this does what I want it to

kse_20200116_plot2 <- kse_20200116_class2  %>%
  merge(t, by = "value")  #  yes!  triple checked lol

# "updated" cluster table 
kse_20200116_plot2 %>% group_by(new_g) %>% summarise(n = n())
#new_g     n
# 1         165 = 12% %
# 2        1204 = 88%

165/(165  +1204)

plotly::plot_ly(kse_20200116_plot2, x = ~logTOF, y = ~logEXT, z = ~green,
                color = ~new_g, colors = ~new_g, 
                alpha = 0.5, size = .5, type="scatter3d")

############
# have plot code ready and organized below to go along with keynote slides
# look at correlation between L1% in data and small objects group #1 

