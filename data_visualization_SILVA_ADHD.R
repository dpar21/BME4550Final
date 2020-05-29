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
library(stringr)
library(ape)
library(igraph)
library(phytools)
library(phyclust)
library(metacoder)
library(taxize)
library(HMPTrees)
library(plyr)
library(phyloseq)

# Taxonomic filtering

# List of files from local directory (SILVA-ADHD)
files = list.files(path = "/Users/dhyeyparikh/MyData/6th_Semester/BME_4550/Project/silva/ADHD", pattern = "*.csv", full.names = T)

# Taxonomy Column Names
into_var = c("domain", "phylum", "class", "order", "family", "genus")

# Aggregate lists of OTUs
otu_list = list() 

# Loop to generate OTUs column for aggregation
for (i in files) {
  tbl = sapply(i, read_csv, simplify=FALSE) %>% 
    bind_rows(.id = "id")
  
  tax_col = tbl[,9]
  
  otu_count_table = as.data.frame(table(tax_col))
  otu_count_table = otu_count_table[order(otu_count_table$Freq, decreasing = TRUE),]
  otu_list[[i]] = otu_count_table

}

# Aggregation of OTUs from all samples
otu_big_data = reduce(otu_list, full_join, by = "tax_col")
otu_big_data[is.na(otu_big_data)] = 0
colnames(otu_big_data)[1] = c("OTUs")
for (i in 2:ncol(otu_big_data)) {
  col_name_otu = paste("ADHD: Sample", i-1)
  colnames(otu_big_data)[i] = col_name_otu
}

# Taxonomic table breakdown from OTUs --> K, P, C, O, F, G
tax_table = data.frame(matrix(ncol = 6, nrow = 0))
colnames(tax_table) = into_var

for (i in 1:nrow(otu_big_data)) {
  count = str_count(otu_big_data[i,1] , ';')
  if (count == 5) {
    new_row_tax = separate(as.data.frame(otu_big_data[i,1]), 1, into_var[c(1, 2, 3, 4, 5, 6)], sep = ";")
    new_row_tax[new_row_tax == ""] = paste("Unclassified genus:", new_row_tax$family)
  } 
  else if (count == 6) {
    new_row_tax = separate(as.data.frame(otu_big_data[i, 1]), 1, into_var, sep = ";")
  } 
  
  
  tax_table = rbind.fill(tax_table, new_row_tax)
}

tax_table = cbind(as.data.frame(otu_big_data[,1]), tax_table)
colnames(tax_table)[1] = "OTUs"

# Reformatting tables for phyloseq
otus_vec = vector()

for (i in 1:nrow(otu_big_data)) {
  val = paste("OTU", i)
  otus_vec = c(otus_vec, val)
}

otus_vec = as.data.frame(otus_vec)
otu_big_data[,1] = otus_vec
tax_table[,1] = otus_vec

# Creation of sample table for phyloseq
samples_mat = data.frame(matrix(ncol = 4, nrow = 19))
colnames(samples_mat) = c("Sample", "Type", "Instrument Model", "Scientific Name")
for (i in 1:nrow(samples_mat)) {
  samples_mat[i,1] = paste("ADHD: Sample", i)
}
samples_mat[,2] = "case"
samples_mat[,3] = "454 GS FLX"
samples_mat[,4] = "Homo Sapiens"

# 1st column binding to row names
row.names(otu_big_data) = otu_big_data$OTUs
otu_big_data = otu_big_data %>% select (-OTUs) 
row.names(tax_table) = tax_table$OTUs
tax_table = tax_table %>% select (-OTUs) 
row.names(samples_mat) = samples_mat$Sample
samples_mat = samples_mat %>% select (-Sample) 

# Convert tables to matricies
otu_big_data = as.matrix(otu_big_data)
tax_table = as.matrix(tax_table)

# Converting matricies to phyloseq objects
OTU = otu_table(otu_big_data, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
samples = sample_data(samples_mat)
carbom = phyloseq(OTU, TAX, samples)

# Normalize number of reads in each sample using median sequencing depth and relative abundance phyloseq object
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
carbom_clean = prune_taxa(taxa_sums(carbom)>150, carbom)
carbom_ra = transform_sample_counts(carbom_clean, function(x) x/sum(x))

# Genera Analysis 
genus_taxa = tax_glom(carbom_clean, taxrank = "genus")
genus_taxa_sort = as.data.frame(names(sort(taxa_sums(genus_taxa), TRUE)[1:5]))
genus_top_sums = data.frame(matrix(ncol = 1, nrow = 5))
for (i in 1:nrow(genus_taxa_sort)) {
  genus_taxa_sort_aval = sum(as.data.frame(get_sample(genus_taxa, toString(genus_taxa_sort[i,1]))))
  genus_top_sums[i,1] = genus_taxa_sort_aval
  rownames(genus_top_sums)[i] = carbom_clean@tax_table[row.names(carbom_clean@tax_table) == toString(genus_taxa_sort[i,1]), colnames(carbom_clean@tax_table@.Data) == "genus"]
}
OTU_table_ra_genus = as.data.frame(otu_table(carbom_ra))
OTU_table_ra_means_genus = data.frame(matrix(ncol = 1, nrow = 5))
for (i in 1:nrow(genus_taxa_sort)) {
  taxa_ra_mean = rowMeans(OTU_table_ra_genus[row.names(OTU_table_ra_genus) == toString(genus_taxa_sort[i,1]), ])
  OTU_table_ra_means_genus[i,1] = taxa_ra_mean
  rownames(OTU_table_ra_means_genus)[i] = carbom_clean@tax_table[row.names(carbom_clean@tax_table) == toString(genus_taxa_sort[i,1]), colnames(carbom_clean@tax_table@.Data) == "genus"]
}

# Phylum analysis
phylum_taxa = tax_glom(carbom_clean, taxrank = "phylum")
phylum_taxa_sort = as.data.frame(names(sort(taxa_sums(phylum_taxa), TRUE)[1:5]))
phylum_top_sums = data.frame(matrix(ncol = 1, nrow = 5))
for (i in 1:nrow(phylum_taxa_sort)) {
  phylum_taxa_sort_aval = sum(as.data.frame(get_sample(phylum_taxa, toString(phylum_taxa_sort[i,1]))))
  phylum_top_sums[i,1] = phylum_taxa_sort_aval
  rownames(phylum_top_sums)[i] = carbom_clean@tax_table[row.names(carbom_clean@tax_table) == toString(phylum_taxa_sort[i,1]), colnames(carbom_clean@tax_table@.Data) == "phylum"]
}
OTU_table_ra_phylum = as.data.frame(otu_table(carbom_ra))
OTU_table_ra_means_phylum = data.frame(matrix(ncol = 1, nrow = 5))
for (i in 1:nrow(phylum_taxa_sort)) {
  taxa_ra_meanP = rowMeans(OTU_table_ra_phylum[row.names(OTU_table_ra_phylum) == toString(phylum_taxa_sort[i,1]), ])
  OTU_table_ra_means_phylum[i,1] = taxa_ra_meanP
  rownames(OTU_table_ra_means_phylum)[i] = carbom_clean@tax_table[row.names(carbom_clean@tax_table) == toString(phylum_taxa_sort[i,1]), colnames(carbom_clean@tax_table@.Data) == "phylum"]
}

# Heatmap of OTUs that represent at least 20% of reads in at least one sample
carbom_abund = filter_taxa(carbom_clean, function(x) sum(x > total*0.20) > 0, TRUE)
plot_heatmap(carbom_abund, method = "MDS", distance = "bray", title = "Top Genera Taxa in ADHD Case Samples",taxa.label = "genus", taxa.order = "genus", trans=NULL, low="lightblue", high="blue", na.value="lightblue") 

# Alpha diversity plot of OTUs
plot_richness(carbom_abund, measures=c("Chao1", "Shannon"), title = "Alpha Diversity Plots for Top Genera Taxa in ADHD Case Samples")
