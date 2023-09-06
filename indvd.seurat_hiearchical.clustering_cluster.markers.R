# IDENTIFICATION OF INDIVIDUAL GENE SIGNATURES via HIERARCHICAL CLUSTERING

# For Bhaduri lab default, execute script as:
## Rscript indvd.seurat_hierarchical.clustering_cluster.markers.R [filepath of individual Seurat object, saved as rds file]

# Else, execute as: 
## Rscript indvd.seurat_hierarchical.clustering_cluster.markers.R [filepath of individual Seurat object, saved as rds file] [working directory] [deepSplit (integer from 0 – 4)]

# Arguments: 
# 1. filepath of individual Seurat object, saved as rds file 
# 2. the filepath of the working directory, in which all outputs will be stored (Bhaduri default: present working directory) 
# 3. deepSplit argument for adaptive branch pruning of hierarchical clustering dendrograms – controls sensitive of cluster splitting, with 0 producing the least number of clusters and 4 producing the most (Bhaduri default: 4)

### RECOMMENDATIONS FOR METAATLAS ###
# Save a folder for each dataset int your metaatlas, in which all individual Seurat objects and all outputs of individual processing (correlation matrix, cluster assignment dataframe, gene score table) will be saved

# Parse arguments
args = commandArgs(trailingOnly=TRUE)

if (is.na(args[1])) {
  stop("input Seurat object must be supplied.n", call.=FALSE)
} else if (is.na(args[2]) | is.na(args[3])) {
  print("using default parameters")
  args[2] = getwd()
  args[3] = 4
}

print(
paste("the filepath of input Seurat object:", args[1]),
paste("the destination of output files:", args[2]),
paste("deepSplit:", args[3])
)

# SET UP

inputfile <- args[1]
setwd(args[2])

library(Seurat)
library(tidyverse)
library(WGCNA)
library(factoextra)
library(dynamicTreeCut)
library(data.table)

# HIERARCHICAL CLUSTERING WITHIN EACH INDIVIDUAL

# This step will perform hierarchical clustering on the correlation matrices of expression
# values. The dynamicTreeCut package will be used to define the cluster membership of each
# cell. Note: cuttreeDynamic() is modulable based on the depth of cutting that you wish to 
# achieve. The deepSplit argument ranges from 0 to 4, with 4 being the most granular 
# (smallest clusters). See dynamicTreeCut package guide for more information.

# Load Seurat object
Indvd_SeuratObject <- readRDS(inputfile) 
print(paste0(inputfile, " Seurat object before Cluster ID"))
Indvd_SeuratObject@meta.data[1:5,]

# Convert normalized expression data into dense matrix for correlation
Indvd_SeuratObject_matrix <- as.matrix(Indvd_SeuratObject@assays$RNA@data)
print(paste0(inputfile, " matrix"))
Indvd_SeuratObject_matrix[1:5, 1:5]
saveRDS(Indvd_SeuratObject_matrix,file=paste0("matrix_", inputfile))
Indvd_SeuratObject_matrix[1:5, 1:5]

# Perform correlation
Indvd_SeuratObject_cor <- cor(Indvd_SeuratObject_matrix, Indvd_SeuratObject_matrix)
print(paste0(inputfile, " self-correlated"))

# Generate correlation-based distance matrix
class(Indvd_SeuratObject_cor) <- "numeric"
print(paste0(inputfile, " converted to numeric"))
dist_Indvd_SeuratObject_cor <- get_dist(Indvd_SeuratObject_cor, method = "pearson", stand = FALSE)
print(paste0(inputfile, " distance matrix computed"))
Indvd_SeuratObject_cor_mat <- as.matrix(dist_Indvd_SeuratObject_cor)
print(paste0(inputfile, " correlation-based distance matrix generated"))
Indvd_SeuratObject_cor_mat[1:5, 1:5]

# Perform hierarchical clustering on the distance matrix, use average linkage
hclust_Indvd_SeuratObject_cor <- hclust(dist_Indvd_SeuratObject_cor, method = 'average')

# Generate clusters
Indvd_SeuratObject_cor_ds <- cutreeDynamic(hclust_Indvd_SeuratObject_cor, cutHeight= NULL, minClusterSize= 10, method = "hybrid", pamStage = TRUE, distM = Indvd_SeuratObject_cor_mat, deepSplit = as.integer(args[3]))
print(paste0(inputfile, " max # clusters"))
max(Indvd_SeuratObject_cor_ds) #shows you number of clusters

# Make dataframe with cell name and cluster ID, saving this is recommended
labels <- hclust_Indvd_SeuratObject_cor$labels # extracts cell names
Indvd_SeuratObject_cuts <- data.frame(labels, Indvd_SeuratObject_cor_ds) # generates dataframe linking cells (rows) to their cluster assignment
print(paste0(inputfile, " cuts"))

# Add cluster ID to original Seurat object as metadata
rownames(Indvd_SeuratObject_cuts) <- Indvd_SeuratObject_cuts$labels
Indvd_SeuratObject <- AddMetaData(Indvd_SeuratObject, metadata = Indvd_SeuratObject_cuts[colnames(Indvd_SeuratObject@assays$RNA@data),])
print(paste0(inputfile, " Seurat object with Cluster ID"))
Indvd_SeuratObject@meta.data[1:5,]

# Saving this Seurat object with cluster ID is recommended
saveRDS(Indvd_SeuratObject_cuts,file=paste0("cuts_", inputfile))
saveRDS(Indvd_SeuratObject, file=inputfile)
print(paste0(inputfile, " hierarchical clustering saved and completed"))

# FIND CLUSTER MARKERS FOR EACH CLUSTER 

# This step will use Seurat to find the cluster marker genes using the cluster ID’s
# derived from the previous step. Based on cluster markers, a GeneScore for each cluster
# marker will be calculated. The GeneScore is calculated as follows: (% of cells in one
# cluster that express one gene / % of all other cells outside of that cluster that
# express the same gene) * avglogFC of that same gene. The GeneScore thereby serves as a
# measure of specificity and expression value of a cluster marker.

#feature selection
Indvd_SeuratObject <- FindVariableFeatures(Indvd_SeuratObject, selection.method = "vst", nfeatures = 2000)

#scale the data
all.genes <- rownames(Indvd_SeuratObject)
Indvd_SeuratObject <- ScaleData(Indvd_SeuratObject, features = all.genes)

#PCA
Indvd_SeuratObject <- RunPCA(Indvd_SeuratObject, features = VariableFeatures(object = Indvd_SeuratObject))
print(paste0(inputfile, " scaling and PCA complete"))

# set the identity to Cluster_ID from the previous step
Idents(Indvd_SeuratObject) <- Indvd_SeuratObject@meta.data$Indvd_SeuratObject_cor_ds
Indvd_SeuratObject.markers <- FindAllMarkers(Indvd_SeuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print(paste0(inputfile, " cluster markers complete"))
head(Indvd_SeuratObject.markers, n=5)

#Calculate GeneScores
Indvd_SeuratObject_genescores <- Indvd_SeuratObject.markers
Indvd_SeuratObject_genescores$GeneScore <- ((Indvd_SeuratObject_genescores $pct.1)/ (Indvd_SeuratObject_genescores$pct.2))* Indvd_SeuratObject_genescores$avg_log2FC

#Create a column that shows “Indvd_Cluster_ID” to use for downstream analysis
Indvd_SeuratObject_genescores$Indvd_Cluster_ID <- paste0(inputfile, sep="_", Indvd_SeuratObject_genescores$Indvd_SeuratObject_cor_ds)
print(paste0(inputfile, " gene scores complete and prepared for downstream analysis"))
head(Indvd_SeuratObject_genescores, n=5)

# Save the Seurat object now scaled with PCA, cluster markers, and gene scores and the gene scores table
saveRDS(Indvd_SeuratObject, file= inputfile)
saveRDS(Indvd_SeuratObject_genescores, file=paste0("genescores_", inputfile))