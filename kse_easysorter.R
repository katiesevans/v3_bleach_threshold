library(tidyverse)
library(easysorter)

# read in data, filter to keep only controls
wd <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_kse/"
filelist <- list.files(wd)

kse_data <- NULL
for(file in filelist) {
    data <- easysorter::read_data(glue::glue("{wd}/{file}")) %>%
        dplyr::filter(contamination == FALSE,
                      call50 == "object", 
                      condition %in% c("DMSO", "water", "lyaste")) %>%
        dplyr::select(date:col, condition, strain, TOF, EXT, stage) %>%
        dplyr::mutate(person = "KSE")
    kse_data <- rbind(kse_data, data)
}

# mark unique assays
kse_data <- kse_data %>%
    dplyr::mutate(bleach = paste0(date, "_", assay))
unique(kse_data$bleach)


# plot all strains together
kse_data %>%
    dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    # dplyr::group_by(trait) %>%
    # dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(date~trait, scales = "free")

# only N2
kse_data %>%
    dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    # dplyr::group_by(trait) %>%
    # dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(date~trait, scales = "free") 
    # geom_vline(aes(xintercept = quant))

# weird assay?
assays <- all_data %>%
    dplyr::filter(strain == "N2") %>%
    dplyr::group_by(date, experiment,round, assay) %>%
    dplyr::summarize(meanTOF = mean(TOF, na.rm = T))

##### Daehan srg crispr data

# read in data, filter to keep only controls
wd <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_dl/"
filelist <- list.files(wd)

DL_data <- NULL
for(file in filelist) {
    data <- easysorter::read_data(glue::glue("{wd}/{file}")) %>%
        dplyr::filter(contamination == FALSE,
                      call50 == "object", 
                      condition %in% c("DMSO", "water", "lyaste", "EtOH")) %>%
        dplyr::select(date:col, condition, strain, TOF, EXT, stage) %>%
        dplyr::mutate(person = "DL")
    DL_data <- rbind(DL_data, data)
}

# mark unique assays
DL_data <- DL_data %>%
    dplyr::mutate(bleach = paste0(date, "_", assay))
unique(DL_data$bleach)

# plot
DL_data %>%
    # dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    # dplyr::group_by(trait) %>%
    # dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(~trait, scales = "free")


##### Clay data

# read in data, filter to keep only controls
wd <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_cd/"
filelist <- list.files(wd)

CD_data <- NULL
for(file in filelist) {
    data <- easysorter::read_data(glue::glue("{wd}/{file}")) %>%
        dplyr::filter(contamination == FALSE,
                      call50 == "object", 
                      condition %in% c("DMSO", "water", "lyaste", "EtOH")) %>%
        dplyr::select(date:col, condition, strain, TOF, EXT, stage) %>%
        dplyr::mutate(person = "CD")
    
    # fix bleaches on 20200130
    if(data$date[1] == "20200130") {
        data <- data %>%
            dplyr::mutate(assay = ifelse(plate < 5, "a", "b"))
    }
    CD_data <- rbind(CD_data, data)
}

# mark unique assays
CD_data <- CD_data %>%
    dplyr::mutate(bleach = paste0(date, "_", assay))
unique(CD_data$bleach)

# plot
CD_data %>%
    # dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    # dplyr::group_by(trait) %>%
    # dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(~trait, scales = "free")

CD_data %>%
    # dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    # dplyr::group_by(trait) %>%
    # dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(date~trait, scales = "free")

############ loraina data
#Define a vector of your experiement directories
dirs <-("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20191003_mianserin_v3/")

dirs2 <- ("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20191114_v3_flx/")
dirs4 <-  ("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20191115_v3_mian/")  ### this folder had a processed folder i had to get rid of 
dirs5 <-("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20200214_flx/") 
dirs6 <-("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20200214_mian/")
dirs7 <-("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20201105_flx/")
dirs8 <-("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20201106_mianserin/")
dirs9 <-("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20201112_flxmian_A1/") 
dirs10 <-("~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data_ls/20201113_flxmian_A2/")
# reads in as a list ; extract out each element of list and rbind into a df orrrr do 9 separate easy sorters 


# Read in the data
raw <- easysorter::read_data(dirs)

# Remove all data from the contaminated wells
raw_nocontam1 <- easysorter::remove_contamination(raw)

raw_nocontam1 <- raw_nocontam1 %>%
    mutate(assay = case_when(plate %in% c("1","2","3", "4", "5") ~
                                 "A",
                             plate %in% c("6","7","8", "9", "10", "11") ~
                                 "B", 
                             plate %in% c("12", "13", "14", "15", "16", "17") ~
                                 "C"))

raw2 <- easysorter::read_data(dirs2)
raw3 <- easysorter::read_data(dirs5)
raw4 <- easysorter::read_data(dirs4)
raw4 <- easysorter::read_data(dirs6)
raw5 <- easysorter::read_data(dirs7)
raw6 <- easysorter::read_data(dirs8)
raw7 <- easysorter::read_data(dirs9)
raw8 <- easysorter::read_data(dirs10)

raw_nocontam2 <- raw2 %>%
    mutate(assay = case_when(plate %in% c("1","2","3", "4") ~
                                 "A",
                             plate %in% c("5","6","7", "8") ~
                                 "B", 
                             plate %in% c("9", "10", "11","12") ~
                                 "C"))

raw_nocontam3 <- raw3 %>%
    mutate(assay = case_when(plate %in% c("7","8","9") ~
                                 "A"))

raw_nocontam4 <- raw4 %>%
    mutate(assay = case_when(plate %in% c("1","2","3") ~
                                 "A",
                             plate %in% c("4","5","6") ~
                                 "B"))


raw_nocontam5 <- raw5 %>%
    mutate(assay = case_when(plate %in% c("1","2","3") ~
                                 "A",
                             plate %in% c("4","5","6") ~
                                 "B",
                             plate %in% c("7", "8", "9")~
                                 "C"))

raw_nocontam6 <- raw6 %>%
    mutate(assay = case_when(plate %in% c("1","2","3","4") ~
                                 "A",
                             plate %in% c("5","6","7", "8") ~
                                 "B",
                             plate %in% c("9", "10", "11")~
                                 "C"))

raw_nocontam7 <- raw7 %>%
    mutate(assay = case_when(plate %in% c("1","2","3","4", "5","6") ~
                                 "A",
                             plate %in% c("7", "8", "9", "10", "11", "12") ~
                                 "B",
                             plate %in% c("13", "14", "15", "16", "17")~
                                 "C"))

raw_nocontam8 <- raw8 %>%
    mutate(assay = case_when(plate %in% c("6","7", "8", "9", "10") ~
                                 "B",
                             plate %in% c( "11", "12","13", "14", "15")~
                                 "C"))


all_nocontam <- rbind(raw_nocontam1,raw_nocontam2, raw_nocontam3, raw_nocontam4, raw_nocontam5, raw_nocontam6, raw_nocontam7,raw_nocontam8)

# mark unique assays
LS_data <- all_nocontam %>%
    dplyr::mutate(condition = ifelse(condition %in% c("MIAN-0", "FLX-0"), "DMSO", condition)) %>%
    dplyr::filter(contamination == FALSE,
                  call50 == "object", 
                  condition %in% c("DMSO", "water", "lyaste", "EtOH")) %>%
    dplyr::select(date:col, condition, strain, TOF, EXT, stage) %>%
    dplyr::mutate(person = "LS") %>%
    dplyr::mutate(bleach = paste0(date, "_", assay))
unique(LS_data$bleach)

# plot
LS_data %>%
    # dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    # dplyr::group_by(trait) %>%
    # dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(~trait, scales = "free")

LS_data %>%
    # dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    # dplyr::group_by(trait) %>%
    # dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(bleach~trait, scales = "free")

#################################################
# put it all together
#################################################

all_data <- kse_data %>%
    dplyr::bind_rows(DL_data, CD_data, LS_data)

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
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(person~trait, scales = "free")
