# Aim: Perform a two sample summary-data MR analysis to find SNPs causing critical COVID-19 phenotype using pseudobulk memory T-cell eQTL summary stats from this paper https://www.nature.com/articles/s41586-022-04713-1 and critical COVID-19 GWAS summary statistics (Release 3, March 2022) https://genomicc.org/data/r3/ 
# The package above will need to be ammended.  

# Used MRAPSS_Rpackage for MR.
# R Version: 

# Written by: JR 17 Aug 2023
# Last updated: JR 17 Sep 2023


## PART 1: Install Packages ######

install.packages("devtools")
install.packages("readr")
devtools::install_github("YangLabHKUST/MR-APSS")

## Load packages
# Tidyverse
library(tidyverse)
# Readr
library(readr)
# Load MR_APSS packages
library(MRAPSS)




## PART 2: Install data ######
# Dataset containing rs numbers from COVID19 HGI B2 data set
COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv <- read.delim("~/Downloads/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz", header=FALSE, comment.char="#")

# T-cell eQTL summary stats data. 
output_formatted <- read.table("~/Documents/MR/scRNAseq_data/output_formatted.txt", quote="\"", comment.char="")

# GWAS critical COVID data.
meta_critical_allcohorts_cleaned_tsv <- read.delim("~/Documents/MR/meta.critical.allcohorts.cleaned.tsv", header=FALSE)




## PART 3: Reduce data set to rs number and MarkerName ######
# Take rs reference data frame and make new "MarkerName" column
all_COVID_merge <- COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv %>%
  mutate(MarkerName = paste(V1, sep = "_", V2))

# Create new data frame with only "MarkerName and rs number 
twoC_all_COVID_merge <- all_COVID_merge %>%
  select(V18, MarkerName)




## PART 4: Merge rs data and crit gwas data to align chr position and rs number. I only want the rs number for this analysis. 
# Add sample size to GWAS dataframe and rename V3
crit_gwas_data_renamed <- meta_critical_allcohorts_cleaned_tsv %>%
  rename("MarkerName" = V1) %>%
  rename("samplesize_col" = V16) %>%
  mutate(Ncontrol = 1755569, Ncase = 24202) 

# Merge "new_crit_gwas_data" and "twoC_all_COVID_merge" data by chromosome position "MarkerName"
rs_crit_gwas_merged <- merge(crit_gwas_data_renamed, twoC_all_COVID_merge, by = "MarkerName")




## PART 5: Reformat rs dataframe for T-cell data ####
# Add chr to beginning of "MarkerName" in COVID19_HGI_... data set. It will be hard to delete chr from the entire column in the sc data set.
chr_twoC_all_COVID_merge <- twoC_all_COVID_merge %>%
  mutate(Chr = "chr")

# Join columns
joined_twoC <- chr_twoC_all_COVID_merge %>%
  mutate(MarkerName = paste(Chr, sep = "", MarkerName))

# Rid extra columns
rs_chr_position <- joined_twoC %>% 
  select(-Chr)




## PART 6: Merge rs data frame and T-cell data frame ####
# Join chromosome number and position
join_col_output_formatted <- output_formatted %>%
  mutate(MarkerName = paste(V2, sep = "_", V3))

# Remove columns
Tcell_data <- join_col_output_formatted %>%
  select(-V1, -V2, -V3, -V6)

# Add N column
Tcell_sample_data <- Tcell_data %>%
  mutate(N = 259) 

# Merge rs data and Tcell data. Takes a long time to compute!
merged_rs_Tcell <- merge(Tcell_sample_data, rs_chr_position, by = "MarkerName")





## PART 7: Reformat GWAS and T-cell dataset for MR_APSS ####
# Format gwas data with MR_APSS code.
MR_APSS_formatted_crit_COVID_GWAS <- format_data(rs_crit_gwas_merged, snp_col = "V18", b_col = "V8", se_col = "V9", A1 = "V6", A2 = "V7", p_col = "V10", ncase_col = "Ncase", ncontrol_col = "Ncontrol")

# Format T-cell data with MR_APPS code
MR_APSS_formatted_Tcell_eQTL_summary <- format_data(merged_rs_Tcell, snp_col = "V18", b_col = "V8", A1 = "V4", A2 = "V5", p_col = "V7", n_col = "N")




## PART 8: Harmonize the formatted data sets and estimate nuisance parameters ####
# Function is based on European LD scores. 
paras = est_paras(dat1 = MR_APSS_formatted_Tcell_eQTL_summary, dat2 = MR_APSS_formatted_crit_COVID_GWAS, trait1.name = "Tcell_eQTL", trait2.name = "Critical_COVID", ldscore.dir = "./eur_w_ld_chr")

# Check the estimates with the following commands
paras$Omega
paras$C

# Harmonized data set will be used for LD clumping 
head(paras$dat)




## PART 9: LD Clumping ####
MRdat = clump(paras$dat, IV.Threshold = 5e-05, SNP_col = "SNP", pval_col = "pval.exp", clump_kb = 1000, clump_r2 = 0.001)
head(MRdat)

## PART 10: Fit MRAPSS
MRres = MRAPSS(MRdat, exposure = "Tcell_eQTL", outcome = "Critical_COVID", C = paras$C, Omega = paras$Omega, Cor.SelectionBias = T)

# Visualize MR-APSS analysis 
MRplot(MRres, exposure = "Tcell_eQTL", outcome = "Critical_COVID")

# Sensivity analysis 
sensitivity(MRdat, Omega = paras$Omega, C = paras$C, exposure = "Tcell_eQTL", outcome = "Critical_COVID")