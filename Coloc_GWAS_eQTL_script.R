# This code runs a colocalisation analysis with coloc package #



library(coloc)
library(data.table)
library(tidyverse)

#### Load the data
GWAS <- fread("/Users/jrogers3/Documents/Computational Genetics/20231121 Colocalization/Tcell_GWAS_Colocalization/GWAS_eur_maf0pt05_cleaned.tsv")

eQTL <- fread("/Users/jrogers3/Documents/Computational Genetics/20231121 Colocalization/Tcell_GWAS_Colocalization/output_formatted copy.txt")

Luo_TBprogression <- fread("/Users/jrogers3/Documents/Computational Genetics/20231121 Colocalization/Tcell_GWAS_Colocalization/Luo_TBprogression.txt")



#### Match variant ID of Luo_TBprogression with variant IDs of GWAS and eQTL by making new column and combining chr, position, ref and variant allele. Data frame will later be merged with eQTL and GWAS data frames. 
Luo_TBprogression$chr = "chr"

ID_column_Luo <- Luo_TBprogression %>%
  mutate(chr_num = paste(chr, sep = "", chromosome))

chr_num_loc_Luo <- ID_column_Luo %>%
  mutate(chr_num_loc = paste(chr_num, sep = ":", base_pair_location_grch37))

chr_num_loc_other_Luo <- chr_num_loc_Luo %>%
  mutate(chr_num_loc_other = paste(chr_num_loc, sep = "_", other_allele))

chr_num_loc_other_effect_Luo <- chr_num_loc_other_Luo %>%
  mutate(SNP = paste(chr_num_loc_other, sep = "_", effect_allele))



#### Format the SNPs in the eQTL data set by chromosome, position, first allele and second allele to match the SNP format from the GWAS so they can be merged. 
CHRPOSeQTL <- eQTL %>%
  mutate(Merge = paste(V2, sep = ":", V3))
         
CHR_POS_Allele1 <- CHRPOSeQTL %>%
  mutate(Merge2 = paste(Merge, sep = "_", V4))

CHR_POS_ALLele12 <- CHR_POS_Allele1 %>%
  mutate(SNP = paste(Merge2, sep = "_", V5))
  


#### First merge eQTL and GWAS SNPs, then merge new that data frame with Luo_TBprogression data frame  
eQTL_gwas_merge <- merge(CHR_POS_ALLele12, GWAS, by = "SNP")

eQTL_GWAS_Luo_merge <- merge(eQTL_gwas_merge, chr_num_loc_other_effect_Luo, by = "SNP")




#### Separate the eQTL and GWAS data now that they have been merged by SNPs 
eQTL_merged <- eQTL_GWAS_Luo_merge %>% 
  select(SNP, V1, V2, V3, V4, V5, V6, V7, V8, effect_allele_frequency)

GWAS_merged <- eQTL_GWAS_Luo_merge %>% 
  select(SNP, CHR, BP, A1, A2, AF_A1, BETA, SE, z.score, P, N_case, N_ctrl, N, effect_allele_frequency)



#### Square the standard error of beta for the eQTL and GWAS data set
eQTL_merged_SE2 <- transform(eQTL_merged, varbeta = V7^2)

GWAS_merged_SE2 <- transform(GWAS_merged, varbeta = SE^2)



# Add sample size and case control column to eQTL data set 
eQTL_N <- eQTL_merged_SE2 %>%
  mutate(N = 259)



# Rename the columns in eQTL and GWAS data sets to match coloc names. 
GWAS_new_name <- GWAS_merged_SE2 %>%
  rename(snp = SNP, beta = BETA, MAF = effect_allele_frequency)

eQTL_new_name <- eQTL_N %>%
  rename(snp = SNP, beta = V8, MAF = effect_allele_frequency)



# Remove duplicated SNPs from GWAS and eQTL data frames
unique_GWAS <- unique(GWAS_new_name)

eQTL_duplicated_rows <- eQTL_new_name[duplicated(eQTL_new_name$snp) | duplicated(eQTL_new_name$snp, fromLast = TRUE), ]

eQTL_no_duplicates <- eQTL_new_name[!duplicated(eQTL_new_name$snp), ]



# Make the dataframes lists
GWASlist <- as.list(unique_GWAS)

GWASlist$CHR <- as.character(GWASlist$CHR)

eQTLlist <- as.list(eQTL_no_duplicates)



# Add variables as described by coloc package 
type_add <- list(type = "cc")
quant_add <- list(type = "quant")

GWASlist <- append(GWASlist, type_add)
eQTLlist <- append(eQTLlist, quant_add)



#### Coloc loop script

# Note sex chromosomes not listed
chr_list_mi <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

output_coloc <- data.frame(matrix(ncol = 22, nrow = 0))
colnames(output_coloc) <- c("snp", "pvalues.df1", "MAF.df1", "N.df1","V.df1","z.df1", "r.df1","lABF.df1", "pvalues.df2","MAF.df2","N.df2","V.df2","z.df2","r.df2","lABF.df2","internal.sum.lABF","SNP.PP.H4","ID", "CHR", "GENPOS","gene_id")

#Completing the colocalisation in a for loop for each chromosome
for (chr_l in chr_list_mi){
  print(paste("ANALYSING CHROMOSOME: ", chr_l))
  D1 <- GWASlist[GWASlist$CHR == chr_l,]
  D1$P <- 10**(D1$LOG10P)
  D1 <- D1[D1$LOG10P >= 5,]
  D2 <- eQTLlist[eQTLlist$CHR == chr_l,]
  
  DC <- merge(x=D1, y=D2, by=c("ID"), all=FALSE)
  
  maf <- DC$maf
  
  my.res <- coloc.abf(dataset1=list(pvalues=DC$P.val, beta=DC$T_geffect_comp, N=407046,type="cc", s=0.047),dataset2=list(pvalues=DC$pval_nominal, beta=DC$slope, N=670,type="quant"),p1 = 1e-04, p2 = 1e-03, p12 = 1e-05, MAF=maf)
  print(my.res) 
  
  h4 <-subset(my.res$results,SNP.PP.H4>0.01)
  print(h4)
  
}

write.csv(output_coloc, "Tcell_eQTL_Crit_COVID_GWAS_colocalization.csv", row.names = FALSE)
















