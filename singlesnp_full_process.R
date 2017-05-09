library(tidyverse)
library(fdrtool)

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

filenames_list <- paste0("mr_base_", snps_list, ".csv")

### Read each into the filenames list 
### They all have the fdr and bonferroni-corrected 
for (idx in 1:length(filenames_list)) {
  temporary <- read_csv(filenames_list[idx])
  temporary <- temporary %>% mutate(., z = beta/se)
  
  jpeg(paste0("locfdr",snps_list[idx], ".jpg"), width = 800, height = 1300)
  lfdr <- fdrtool(temporary$z[complete.cases(temporary$z)], statistic = "normal")
  dev.off()
  
  print("###~####")
  print(snps_list[idx])
  print(lfdr$param)
  print("!!!!~!!!")
  write_csv(temporary, path = paste0("./locfdr", snps_list[idx],".csv"))
  
}

planned_analyses <- read_csv("./singlesnp_results_top25.csv")
strict_sigs <- planned_analyses %>% filter(., bonpv < 0.05) %>% select(., SNP, originalname.outcome, beta.outcome, se.outcome, bonpv) %>% arrange(., bonpv)



plannedplot <- planned_analyses %>% select(., SNP, originalname.outcome, bonpv, fdrqv)

bonfiplot <- ggplot(data = plannedplot, aes(x = SNP, y = originalname.outcome, fill = bonpv)) + geom_tile() + theme_minimal() + scale_fill_gradient2(midpoint = 0.025, low = "black", trans = "sqrt", mid = "blue", high = "white" ) + labs(ylab = "Outcome", fill = "p-value", title = "Bonferroni-Corrected p-values") + theme(axis.text.x = element_text(angle = 90, vjust = 0.7))

ggsave("bonfiplot.jpg")

fdriplot <- ggplot(data = plannedplot, aes(x = SNP, y = originalname.outcome, fill = fdrqv)) + geom_tile() + theme_minimal() + scale_fill_gradient2(midpoint = 0.025, low = "black", trans = "sqrt", mid = "blue", high = "white") + labs(ylab = "Outcome", fill = "Q-values", title = "FDR q-values") + theme(axis.text.x = element_text(angle = 90, vjust = 0.7))

ggsave("fdriplot.jpg")
topseven <- c("rs4988235",
              "rs708272",
              "rs662799",
              "rs1800588",
              "rs328",
              "rs1800629",
              "rs1801133"
              )


gene_corr <- read_csv("combined_results_2smr_all.csv")

gene_sigs <- gene_corr %>% filter(., bonpv < 0.05)  %>% select(., outcome, exposure, nsnp, b, se, bonpv) %>% arrange(., bonpv, outcome, exposure )


gcbonbonplot <- ggplot(data = gene_corr, aes(x = exposure, y = outcome, fill = bonpv)) + geom_tile() + theme_minimal() + scale_fill_gradient2(midpoint = 0.025, low = "black", trans = "sqrt", mid = "blue", high = "white") + labs(xlab = "Exposure", ylab = "Outcome", fill = "p-value", title = "Bonferroni-Adjusted p-values") + theme(axis.text.x = element_text(angle = 90, vjust = 0.6))

ggsave("gcbonplot.jpg")



gcfdrplot <- ggplot(data = gene_corr, aes(x = exposure, y = outcome, fill = fdrqv)) + geom_tile() + theme_minimal() + scale_fill_gradient2(midpoint = 0.025, low = "black", trans = "sqrt", mid = "blue", high = "white") + labs(xlab = "Exposure", ylab = "Outcome", fill = "q-value", title = "FDR q-values") + theme(axis.text.x = element_text(angle = 90, vjust = 0.6))

ggsave("gcfdrplot.jpg")

outcome_names <- unique(planned_analyses$id.outcome)
outfile_list <- paste0("./genetic_correlations_2smr",outcome_names,".csv")


for (idx in 1:length(outfile_list)) {
  if( file.exists(outfile_list[idx])) {
  temporary <- read_csv(outfile_list[idx])
  temporary <- temporary %>% mutate(., z = b/se)
  
  jpeg(paste0("locfdr",outcome_names[idx], ".jpg"), width = 800, height = 1300)
  lfdr <- fdrtool(temporary$z[complete.cases(temporary$z)], statistic = "normal")
  dev.off()
  
  print("###~####")
  print(outfile_list[idx])
  print(lfdr$param)
  print("!!!!~!!!")
  write_csv(temporary, path = paste0("./locfdr", outcome_names[idx],".csv"))
  }
}