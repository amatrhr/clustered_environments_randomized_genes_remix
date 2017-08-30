library(tidyverse)
library(fdrtool)
library(dplyr)
library(magrittr)
library(TwoSampleMR)
library(MRInstruments)
library(knitr)
library(readr)
library(stringr)

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

## clean rsids 
traits_by_snps$SNP <- stringr::str_to_lower(traits_by_snps$SNP)

## remove traits that aren't equally measured by _all_ SNPs
tbs_filter <- traits_by_snps %>% group_by(., outcome) %>% summarize(., non_na = sum(!is.na(pval.outcome) ) > 22) 

## Within each SNP, compute q-values using fdrtool
### Then, compute the median q-value within each category for each SNP
traits_fdrq <- traits_by_snps  %>% filter(., outcome %in% tbs_filter$outcome[tbs_filter$non_na], !is.na(pval.outcome)) %>% group_by(., SNP) %>% mutate(., fdrq = fdrtool(pval.outcome, statistic = "pvalue")$qval)

traits_fdr_medians <- traits_fdrq %>% summarize(., mfdrq = median(fdrq)) %>% arrange(., desc(mfdrq))

traits_fdr_medians_by_subcat <- traits_fdrq %>% do(data.frame(aggregate(data = ., fdrq ~ subcategory.outcome, FUN = "median")))

traits_fdr_medians_by_subcat_top7 <- traits_fdrq %>% do(data.frame(aggregate(data = ., fdrq ~ subcategory.outcome, FUN = "median"))) %>% filter(., SNP %in%  c('rs4988235', 'rs708272', 'rs3135506', 'rs1800588', 'rs328', 'rs1800629', 'rs1801133') )

## ($rs4988235$), CETP ($rs708272$), APO_AV ($rs3135506$), HL ($rs1800588$), LPL ($rs328$), TNF$\alpha$ ($rs1800629$), MTHFR ($rs1801133$) 

## Heatmap display of median q-values?
traitq_hm_data <- reshape2::melt(traits_fdr_medians_by_subcat_top7)

traitq_hm_data <- traitq_hm_data %>% mutate(gene_tagged = str_replace_all(SNP, c("rs4988235" = "MCM6" , "rs708272$" = "CETP",  "rs3135506" = "APO_AV",  "rs1800588$" = "HL",   "rs328" = "LPL", "rs1800629" = "TNFa",  "rs1801133" =  "MTHFR"))) %>% arrange(., desc(value))

medianq_heatmap <- ggplot( data = traitq_hm_data, aes( x = gene_tagged, y  = subcategory.outcome, fill = log(value, 10)))+ geom_tile()  + theme_minimal() + theme(axis.text.x = element_text(angle = 90 )) + scale_fill_gradient2()


#### Look at the individual categories to check number of associations and median p-value
categories_fdrq <- traits_by_snps %>% filter(., !is.na(pval.outcome)) %>% group_by(., subcategory.outcome) %>% mutate(., fdrq = fdrtool(pval.outcome, statistic = "pvalue")$qval)

categories_median_fdrq <- categories_fdrq %>% summarize(., medfdr = median(fdrq), n1 = n())


