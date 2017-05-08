                                        # This script
## Lists the Exposures and SNPs in MR-Base to which we want access.
## Accesses the MR base resource and limits the phenotypes to a sensible subset
## Performs the single-SNP analyses
### Counts, for each, the number of significant associations
#### First, from the 24 available phenotypes
#### And then a count of all the phenotypes which exceed a value?
##### Does doing that require us to know how many phenotypes there are in general? 
## Performs all 276 pairwise analyses (actually, is it 552?)
### Gets the p-values, SEs, and z statistics from each 2SMR

library(dplyr)
library(magrittr)

library(TwoSampleMR)
library(MRInstruments)
library(knitr)
library(readr)

## Lists the Exposures and SNPs in MR-Base to which we want access.
ao <- available_outcomes()

### For each outcome, need to specify a _single_ row.
### This means a years field might need to be added 

trait_list <- c(
    "Body mass index",
    "Height",
    "Waist-to-hip ratio",
    "Mean cell volume",
    "Haemoglobin concentration",
    "Platelet count",
    "Urea",## not enough GWS instruments
    "Urate", ## not enough GWS instruments
    "Phosphate", ## not enough GWS instruments
    "Albumin",
    "Bilirubin (E,Z or Z,E)*",
    "Total cholesterol",
    "Triglycerides",
    "Ascorbate (Vitamin C)",## not enough GWS instruments
    "Alpha-tocopherol",## not enough GWS instruments
    "Fasting insulin",
    "Fasting glucose",
    "Subjective well being",
    "Years of schooling",
    "Cigarettes smoked per day",
    "Weight",
    "Alcohol dependence", ## not enough GWS instruments
    "Birth weight",
    "Age at menopause",
    "Age at menarche",
    "Father"
    
)

year_list <- c(
    "2015",
    "2014",
    "2015",
    "2012",
    "2012",
    "2011",
    "2014",
    "2013",
    "2014",## not enough GWS instruments
    "2016",
    "2014",
    "2013",
    "2013",
    "2014",## not enough GWS instruments
    "2014",## not enough GWS instruments
    "2013",
    "2012",
    "2016",
    "2016",
    "2010",
    "2013",
    "2012",
    "2016",
    "2015",
    "2014",
    "2016"
    
)


grabtrait <- function(index, traitlist, yearlist,  exttraitframe){

    traitrows <-  intersect(grep(traitlist[index],exttraitframe$trait, fixed = TRUE),grep(yearlist[index], exttraitframe$year, fixed = TRUE))
    outrow <- which.max(exttraitframe$priority[traitrows])
    traitrows[outrow]
    
}

selected_traits <- as.vector(sapply(1:length(year_list), FUN = grabtrait, traitlist = trait_list, yearlist = year_list, exttraitframe = ao))
##
selected_outcomes <- ao[selected_traits,]

#### 
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

## Accesses the MR base resource and limits the phenotypes to a sensible subset

## Performs the single-SNP analyses
traits_by_snps <- extract_outcome_data(snps = snps_list, outcomes = selected_outcomes$id)

### Counts, for each, the number of significant associations
#### First, from the 24 available phenotypes
### Under a bonferroni correction and also an FDR for more 'bite


tbs_padj <- traits_by_snps %>% mutate(., bonpv = p.adjust(pval.outcome, method = "holm"), fdrqv = p.adjust(pval.outcome, method = "fdr"))
### Which columns from this to write out?
write_csv(tbs_padj, path = "./singlesnp_results_top25.csv")
capture.output(kable(tbs_padj), file  = "singlesnp_results_top25.md", split = TRUE)
#### And then a count of all the phenotypes which exceed a value?
##### Not at this time
##### Does doing that require us to know how many phenotypes there are in general? 
## Performs all 276 pairwise analyses (actually, is it 552?)
### Gets the p-values, SEs, and z statistics from each 2SMR

## First try to harmonize


for (idvar in setdiff(1:length(selected_outcomes$id), c(7,8,9,14,15,22)) ) {
    uuselsnps <- extract_instruments(selected_outcomes$id[idvar])
    vvselsnps <- extract_outcome_data(snps = uuselsnps$SNP, outcomes = selected_outcomes$id)
    snps_Only <- harmonise_data(exposure_dat = uuselsnps, outcome_dat = vvselsnps )
    mr_out<- mr(snps_Only, method_list=c("mr_ivw"))

    mr_out <- mr_out %>% mutate(.,  bonpv = p.adjust(pval, method = "holm"), fdrqv = p.adjust(pval, method = "fdr"))

    write_csv(mr_out, path = paste0("./genetic_correlations_2smr",selected_outcomes$id[idvar],".csv"))

    capture.output(kable(mr_out), file = paste0("./genetic_correlations_2smr",selected_outcomes$id[idvar],".md"),split = TRUE)
}



