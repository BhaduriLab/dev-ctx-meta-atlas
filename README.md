# dev-ctx-meta-atlas

## System Requirements

### Dependencies
R version 4 or later and Seurat version 4 or later. 
Packages required for each step are indicated in each script, and include:
    *Seurat, tidyverse, WGCNA, factoextra, dynamicTreeCut, data.table, ggplot2, stringr, dplyr, readxl*

### Versions tested on:
  R v4.0.2, 4.1.0, 4.3.1
  Seurat v4.4.0, v5.0.1

## Instructions for use

## indvd.seurat_hierarchical.clustering_cluster.markers.R
Identification of gene signatures within each individual via hierarchical clustering.

For Bhaduri lab default, execute script as:
*Rscript indvd.seurat_hierarchical.clustering_cluster.markers.R [filename of individual Seurat object, saved as rds file]* 

Else execute script as:
*Rscript indvd.seurat_hierarchical.clustering_cluster.markers.R [filename of individual Seurat object, saved as rds file] [working directory] [deepSplit (integer from 0 – 4)]*

Arguments: 
 1. filename of individual Seurat object, saved as rds file 
 2. the filepath of the working directory, containing all inputfiles and in which outputs will be stored 
 3. deepSplit argument for adaptive branch pruning of hierarchical clustering dendrograms – controls sensitive of cluster splitting, with 0 producing the least number of clusters and 4 producing the most (Bhaduri default: 4)

## metamodule.generation.R
Aggregation of individual cluster markers and gene scores, filtration of cluster markers based on gene score, and hierarchical clustering of cluster markers to generate metamodules

For Bhaduri lab default, execute script as:
*Rscript metamodule.generation.R [filepath to .lis file] [directory] *

Else, execute as:
*Rscript metamodule.generation.R [filepath to .lis file] [directory] [file prefix for metaatlas files] [gene score percentile threshold for cluster marker filtration]*

Arguments: 
 1. filepath to a .lis file with the names of all the folders associated with each dataset to be added in metaatlas (one line for each filepath)
 2. directory containing all folders associated with each dataset to be added in metaatlas
 3. (optional) file prefix for metaatlas files
 4. (optional) gene score percentile threshold for cluster marker filtration (Bhaduri lab default: 0.90)

## module.activity.score.R
Score cells in datasets for meta-module activity (avg log-normd CPMs)

Arguments: 
  1. filepath to module gene list table, with each row containing a module gene (column named: labels) and its assigned module (column named: metamodules) 
  2. directory containing the query dataset (argument #3)
  3. filename of query dataset: a normalized and scaled seurat object for which module activities are to be calculated 

Execute as:
*Rscript module.activity.score.R [filepath to module gene list table] [directory] [filename of query dataset]**

## module.specificity.score.R
Score the ability of modules to represent specific groupings of cells (ex: cell types). 

Arguments: 
  1. working directory containing CSVs of query dataset metadata and module activity scores (arguments #2 and #3)
  2. filename of a CSV of metadata for each cell in query dataset
  3. filename a CSV of module activity scores for each cell in query dataset (cells in rows; modules in columns)
  4. column name of desired metadata group across which to calculate specificity scores
  5. cutoff value for module activity binarization (default: 0.9)
  6. output filename (needs to end in .csv)

Execute as:
*Rscript module.specificity.score.R [directory] [csv of metadata for each cell] [csv of module activity scores] [desired column of metadata across which to calculate specificity scores] [cutoff value for module activity binarization (default: 0.9)] [output filename]* 

