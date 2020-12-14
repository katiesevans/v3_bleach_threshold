library(tidyverse)
library(easysorter)

# read in data, filter to keep only controls
wd <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/v3_bleach_threshold/data/"
filelist <- list.files(wd)

all_data <- NULL
for(file in filelist) {
    data <- easysorter::read_data(glue::glue("{wd}/{file}")) %>%
        dplyr::filter(contamination == FALSE,
                      call50 == "object", 
                      condition %in% c("DMSO", "water", "lyaste"),
                      strain %in% c("N2", "CB4856")) %>%
        dplyr::select(date:col, condition, strain, TOF, EXT, stage) %>%
        dplyr::mutate(person = "KSE")
    all_data <- rbind(all_data, data)
}

all_data %>%
    dplyr::filter(strain == "N2") %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    dplyr::group_by(trait) %>%
    dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(date~trait, scales = "free") +
    geom_vline(aes(xintercept = quant))

# weird assay?
assays <- all_data %>%
    dplyr::filter(strain == "N2") %>%
    dplyr::group_by(date, experiment,round, assay) %>%
    dplyr::summarize(meanTOF = mean(TOF, na.rm = T))

all_data %>%
    dplyr::filter(strain == "N2",
                  date) %>%
    tidyr::gather(trait, phenotype, TOF:EXT) %>%
    dplyr::group_by(trait) %>%
    dplyr::mutate(quant = quantile(phenotype, probs = 0.05)) %>%
    ggplot(.) +
    aes(x = phenotype) +
    geom_histogram(bins = 150) +
    theme_bw() +
    facet_grid(~trait, scales = "free") +
    geom_vline(aes(xintercept = quant))
