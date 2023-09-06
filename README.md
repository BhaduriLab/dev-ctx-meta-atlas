# dev-ctx-meta-atlas

Tools for meta-module generation and activity scoring, as described in "A Meta-Atlas of the Developing Human Cortex Identifies Modules Driving Cell Subtype Specification", Nano PR et al., bioRxiv, 2023.

This repository contains three scripts that can be used in sequence to generate meta-modules from a meta-atlas, and score module activity in a variety of query datasets.

*Recommended for meta-atlas generation*: Save a folder for each dataset in your meta-atlas, in which all individual Seurat objects and all outputs of individual processing (correlation matrix, cluster assignment dataframe, gene score table) will be saved.

## indvd.seurat_hiearchical.clustering_cluster.markers.R
Identification of gene signatures within each individual via hierarchical clustering.

For Bhaduri lab default, execute script as:
*Rscript indvd.seurat_hierarchical.clustering_cluster.markers.R [filepath of individual Seurat object, saved as rds file]*

Else, execute as: 
*Rscript indvd.seurat_hierarchical.clustering_cluster.markers.R [filepath of individual Seurat object, saved as rds file] [working directory] [deepSplit (integer from 0 – 4)]*

Arguments: 
1. filepath of individual Seurat object, saved as rds file 
2. the filepath of the working directory, in which all outputs will be stored (Bhaduri default: present working directory) 
3. deepSplit argument for adaptive branch pruning of hierarchical clustering dendrograms – controls sensitive of cluster splitting, with 0 producing the least number of clusters and 4 producing the most (Bhaduri default: 4)

## metamodule.generation.R
Aggregation of individual cluster markers and gene scores, filtration of cluster markers based on gene score, and hierarchical clustering of cluster markers to generate metamodules

For Bhaduri lab default, execute script as:
*Rscript metamodule.generation.R [directory] [filepath to .lis file]*

Else, execute as:
*Rscript metamodule.generation.R [directory] [filepath to .lis file] [file prefix for metaatlas files] [gene score percentile threshold for cluster marker filtration]*

Arguments: 
1. directory containing all folders associated with each dataset to be added in meta-atlas
2. filepath to a .lis file with the names of all the folders associated with each dataset to be added in meta-atlas (one line for each filepath)
3. (optional) file prefix for meta-atlas files
4. (optional) gene score percentile threshold for cluster marker filtration (Bhaduri lab default: 0.90)

## module activity.score.R
Score cells in datasets for meta-module activity (avg log-normd CPMs)

Execute as:
*Rscript module.activity.score.R [directory] [module gene list table] [filepath to query dataset]* 



