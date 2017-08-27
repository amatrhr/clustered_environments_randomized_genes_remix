library(tidyverse)
library(fdrtool)
library(dplyr)
library(magrittr)
library(TwoSampleMR)
library(MRInstruments)
library(knitr)
library(readr)

## Lists the Exposures and SNPs in MR-Base to which we want access.
ao <- available_outcomes()

snps_list <- c(
    "rs4988235",
    "rs708272",
    "rs662799",
    "rs1800588",
    "rs328",
    "rs4075131",
    "rs1007",
    "rs1800629",
    "rs1041981",
    "rs7291467",
    "rs1800795",
    "rs20417",
    "rs2234693",
    "rs1801133",
    "rs325400",
    "rs1501299",
    "rs854541",
    "rs1799884",
    "rs4343",
    "rs2071307",
    "rs662",
    "rs696217",
    "rs713598"
)




## Accesses the MR base resource, gets all outcomes for the 23 SNPs

traits_by_snps <- extract_outcome_data(snps = snps_list, outcomes = unique(ao$id))

traits_by_snps$SNP <- stringr::str_to_lower(traits_by_snps$SNP)
## Within each SNP, compute q-values using fdrtool
### Then, compute the median q-value within each category for each SNP
traits_fdrq <- traits_by_snps %>% filter(., !is.na(pval.outcome)) %>% group_by(., SNP) %>% mutate(., fdrq = fdrtool(pval.outcome, statistic = "pvalue")$qval)

traits_fdr_medians_by_subcat <- traits_fdrq %>% do(data.frame(aggregate(data = ., fdrq ~ subcategory.outcome, FUN = "median") ))
 

## Heatmap display of median q-values?
traitq_hm_data <- reshape2::melt(traits_fdr_medians_by_subcat)
medianq_heatmap <- ggplot( data = traitq_hm_data, aes( x = SNP, y  = subcategory.outcome, fill = log(value, 10)))+ geom_tile()  + theme_minimal() + theme(axis.text.x = element_text(angle = 90 )) + scale_fill_gradient2()

