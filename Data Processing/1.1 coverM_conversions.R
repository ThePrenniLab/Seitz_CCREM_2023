
## ---------------------------
## Script merge_coverM_outputs
## 
## Purpose merge MetaG coverM outputs
##
## By Mikayla Borton 
## Date adapted: May 1, 2023
##
## ---------------------------
## Script adapted from place_coassemblies_to_treatment script
## Author: Ikaia Leleiwi
##
## Date Created: May 16, 2022
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
## ---------------------------

##Libraries##
library(tidyverse)
library(edgeR)

#mapping files
#default coverM output with min_covered_fraction set
covered_fraction <- read.csv(paste0( "coverm_min75_reformattedForR.csv"))
#coverM output with -m reads_per_base
reads_per_base <- read.csv(paste0( "coverm_reads_per_base_reformattedForR.csv"))
#coverM output with -m trimmed_mean
trimmed_mean <- read.csv(paste0( "coverm_trimmed_mean_reformattedForR.csv"))

##Functions##
relabund <- function(df, columns = c(NA)) 
  #takes feature table and calculates relative abundance of columns, omitting NA's
  #considers only columns listed in columns argument (character vector of column names). 
  #columns (default) = all columns
{
  if (NA %in% columns){
    df <- sweep(df, 2, colSums(df, na.rm = TRUE), '/')
  }
  else {
    df_relabund <- df %>% select(all_of(columns)) %>%
      sweep(., 2, colSums(., na.rm = TRUE), '/')
    df[,columns] <- df_relabund
  }
  
  return(df)
}

#rename function for reads_per_base
rename_cols <- function(x){
  x <- str_remove(x, "97peridmapped")
  x <- str_remove(x, "_99perMAGs_POSSORT Reads per base")
  x <- str_replace_all(x, " ", "_")
  return(x)
}

#rename function for default coverM
rename_cols_mean <- function(x){
  x <- str_remove(x, "97peridmapped")
  x <- str_remove(x, "_99perMAGs_POSSORT")
  x <- str_remove(x, paste0(" ", " Relative Abundance (%)"))
  return(x)
}

#rename function for trimmed mean
rename_cols_tm <- function(x){
  x <- str_remove(x, "97peridmapped")
  x <- str_remove(x, "_99perMAGs_POSSORT")
  x <- str_remove(x, " Trimmed Mean")
  return(x)
}

#modify longer function
#note your sample names should not have underscores, this is looking for the first underscore and then will split into sample id and map_type
mod_long <- function(df){
  df <- df %>%
    rename_with(rename_cols) %>%
    pivot_longer(cols = -Genome,
                 values_to = "count",
                 names_to = "sample") %>%
    separate(col = "sample",
             into = c("sample_id", "map_type"),
             extra = "merge")
}


##Combine data##
df_list <- list(covered_fraction,
                trimmed_mean,
                reads_per_base)

combined_df <- lapply(df_list, mod_long) %>%
  reduce(rbind)

combined_df_wide <- combined_df %>%
  pivot_wider(names_from = "map_type",
              values_from = "count")


##Filter data##
#changed names here to match what i used above 
bin_counts_matrix <- combined_df_wide %>%
  mutate(trimmed_mean = ifelse( covered_fraction > 0, trimmed_mean, 0), #calculate genomes with 75% MAG coverage
         trimmed_mean = ifelse(reads_per_base*151 >= 3, trimmed_mean, 0)) %>% #calculate 3x coverage per base
  select(Genome, sample_id, trimmed_mean) %>%
  pivot_wider(values_from = trimmed_mean,
              names_from = sample_id) %>%
  arrange(Genome) %>%
  column_to_rownames(var = "Genome") %>%
  as.matrix()

bin_counts_out <- bin_counts_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "bin.id")

write_csv(bin_counts_out, "strict_mapping_table_95id.csv")

