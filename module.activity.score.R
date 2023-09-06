# METAMODULE ACTIVITY SCORING

#########
# SET UP
#########

library(stringr)
library(dplyr)
library(Seurat)
library(tidyr)
library(readxl)
library(tidyverse)

### RECOMMENDATIONS FOR METAATLAS ###
# Save a folder for each dataset int your metaatlas, in which all individual Seurat objects and all outputs of individual processing (correlation matrix, cluster assignment dataframe, gene score table) will be saved
# The module gene list used in this script will be generated in the last step of metamodule.generation.R (metaatlas_metamodules, saved with the suffix "_metamodules.rds")

# Arguments: 
# 1. directory containing all folders associated with each dataset to be added in metaatlas
# 2. module gene list table, with each row containing a module gene (column named: labels) and its assigned module (column named: metamodules) 
# 3. filepath to a normalized and scaled seurat object for which module activities are to be calculated 

# Execute as:
## Rscript module.activity.score.R [directory] [module gene list table] [filepath to query dataset] 

# Parse arguments
args = commandArgs(trailingOnly=TRUE)

if (is.na(args[1])| is.na(args[2]) | is.na(args[3])) {
  stop("argument missing.n", call.=FALSE)
  } 

print(
	paste("directory of datasets to be included in metaatlas:", args[1]),
	paste("module gene list table:", args[2]),
	paste("filepath to query dataset", args[3])	
	)

#########
# SET UP
#########

setwd(args[1])
module.gene.list.table <- args[2]
query.dataset <- args[3]

#########
# Activity of gene lists in dataset
#########

## Load list of module genes

modules <- group_by(readRDS(module.gene.list.table), metamodules)
head(modules)
modules$metamodules <- as.character(modules$metamodules)

## Establish a function to calculate avg log-normd CPMs for each module gene list

module.act <- 
  function(file) 
    {
  
  dataset <- readRDS(file)
  dataset.name <- strsplit(file, ".rds")[1]
  lognormd.dataset <- dataset@assays$RNA@data
  rownames(lognormd.dataset)<-toupper(rownames(lognormd.dataset)) # genes need to be uppercase in order to match with the human-derived meta-atlas gene list
  lognormd.dataset[1:5,1:5]
  
  cell.names <- rownames(dataset@meta.data)
  modules$metamodules <- as.character(modules$metamodules)
  
  all.modules_avg.lognormd.dataset <- as.data.frame(matrix(nrow = length(cell.names), ncol = n_distinct(modules$metamodules)))
  rownames(all.modules_avg.lognormd.dataset) <- cell.names
  colnames(all.modules_avg.lognormd.dataset) <- unique(modules$metamodules)
  
  for (i in unique(modules$metamodules)){ 
  print(paste0("-------------------------------------------", i))
  
  # Extract module gene list genes
  modules.genes <- subset(modules, metamodules == i)

  # Extract module gene list genes that are present in the dataset
  shared_modules.genes <- intersect(rownames(lognormd.dataset), modules.genes$genes)
  print("module genes not in dataset:")
  print(c((100-(length(shared_modules.genes)/length(modules.genes$genes))*100), "% of module genes"))
  print(modules.genes[!(modules.genes %in% shared_modules.genes)]$genes)
  
  # Extract lognormd CPMs for each module gene in each cell
  modules.lognormd.dataset <- lognormd.dataset[shared_modules.genes,]
  
  # Calculate the total lognormd CPMs for the entire gene.list
  total_modules.lognormd.dataset <- as.data.frame(colSums(modules.lognormd.dataset))
  
  # Divide lognormd CPMs by number of genes in module
  avg_modules.lognormd.dataset <-  total_modules.lognormd.dataset
  avg_modules.lognormd.dataset$`colSums(modules.lognormd.dataset)` <-  total_modules.lognormd.dataset$`colSums(modules.lognormd.dataset)` / length(shared_modules.genes)
  colnames(avg_modules.lognormd.dataset) <- "avg_modules.lognormd.dataset"
  
  # If cell names retained, append module gene list-specific columns into an aggregated avg gene.list lognormd CPM table 
  if(!(all(colnames(lognormd.dataset) == rownames(avg_modules.lognormd.dataset)))){
    print(paste("error", i))
    } else {
    all.modules_avg.lognormd.dataset[,i] <- avg_modules.lognormd.dataset$avg_modules.lognormd.dataset
    }
  print(all.modules_avg.lognormd.dataset[1:5,])
  rm(i, modules.genes, shared_modules.genes, modules.lognormd.dataset, total_modules.lognormd.dataset, avg_modules.lognormd.dataset)
  }
  
  # Sort module column names in numerical order
  column.order <- as.character(sort(as.numeric(colnames(all.modules_avg.lognormd.dataset))))
  all.modules_avg.lognormd.dataset <- all.modules_avg.lognormd.dataset[,column.order]
  all.modules_avg.lognormd.dataset[1:5,1:5]
  
  # Save table as csv
  write.csv(
    all.modules_avg.lognormd.dataset,
    file = paste0("all.modules_avg.lognormd_", dataset.name,".csv")
    )

  }

## Run function for query dataset

module.act(query.dataset)
