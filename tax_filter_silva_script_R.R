#loading libraries

library(ggplot2)
library(tidyverse)
library(rafalib)
library(devtools)
library(RColorBrewer)
library(gplots)
library(downloader)
library(genefilter)
library(SpikeInSubset)
library(RCy3)
library(readr)
library(dplyr)


# Taxonomic filtering

# List of files from local directory (SILVA-ADHD, SILVA-Control)
files = list.files(path = "/Users/dhyeyparikh/MyData/6th_Semester/BME_4550/Project/silva/control", pattern = "*.csv", full.names = T)

# Taxonomy Column Names
into_var = c("domain", "phylum", "class", "order", "family", "genus")

# Aggregate lists of genus and phyla for all files
genus_list = list() 
phylum_list = list()

# Loop to generate genus and phyla columns for aggregation
for (i in files) {
  tbl = sapply(i, read_csv, simplify=FALSE) %>% 
    bind_rows(.id = "id")
  
  tax_col = tbl[,9]
  tax_col = separate(tax_col, 1, into_var, sep = ";")
  tax_col[tax_col == ""] <- NA
  
  genus_count_table = as.data.frame(table(tax_col$genus))
  genus_count_table = genus_count_table[order(genus_count_table$Freq, decreasing = TRUE),]
  genus_list[[i]] = genus_count_table
  
  phylum_count_table = as.data.frame(table(tax_col$phylum))
  phylum_count_table = phylum_count_table[order(phylum_count_table$Freq, decreasing = TRUE),]
  phylum_list[[i]] = phylum_count_table
}

# Aggregation of phyla and genus values
genus_big_data = do.call(rbind, genus_list)
genus_big_data_ag = aggregate(x = genus_big_data$Freq, by = list(genus_big_data$Var1), FUN = sum)
genus_big_data_ag_sort = genus_big_data_ag[order(genus_big_data_ag$x, decreasing = TRUE),]

phylum_big_data = do.call(rbind, phylum_list)
phylum_big_data_ag = aggregate(x = phylum_big_data$Freq, by = list(phylum_big_data$Var1), FUN = sum)
phylum_big_data_ag_sort = phylum_big_data_ag[order(phylum_big_data_ag$x, decreasing = TRUE),]

# Top 5 taxonomic breakdown of all files
genus_big_data_ag_sort = genus_big_data_ag_sort[1:5, ]
phylum_big_data_ag_sort = phylum_big_data_ag_sort[1:5, ]

# Subset of 2 highest genus taxonomy 
genus_big_data_tax_genus1 = subset(genus_big_data, genus_big_data$Var1 == genus_big_data_ag_sort[1,1] )
genus_big_data_tax_genus2 = subset(genus_big_data, genus_big_data$Var1 == genus_big_data_ag_sort[2,1] )

# Subset of 2 highest phylum taxonomy
phylum_big_data_tax_phylum1 = subset(phylum_big_data, phylum_big_data$Var1 == phylum_big_data_ag_sort[1,1] )
phylum_big_data_tax_phylum2 = subset(phylum_big_data, phylum_big_data$Var1 == phylum_big_data_ag_sort[2,1] )


