# METAMODULE GENERATION

### RECOMMENDATIONS FOR METAATLAS ###
# Save a folder for each dataset int your metaatlas, in which all individual Seurat objects and all outputs of individual processing (correlation matrix, cluster assignment dataframe, gene score table) will be saved

# Arguments: 
# 1. directory containing all folders associated with each dataset to be added in metaatlas
# 2. filename of .lis file with the names of all the folders associated with each dataset to be added in metaatlas (one line for each filepath)
# 3. (optional) file prefix for metaatlas files
# 4. (optional) gene score percentile threshold for cluster marker filtration (Bhaduri lab default: 0.90)

# For Bhaduri lab default, execute script as:
## Rscript metamodule.generation.R [directory] [filename of .lis file]

# Else, execute as:
## Rscript metamodule.generation.R [directory] [filename of .lis file] [file prefix for metaatlas files] [gene score percentile threshold for cluster marker filtration]

# Parse arguments
args = commandArgs(trailingOnly=TRUE)

if (is.na(args[1])| is.na(args[2])) {
  stop("directory and filepath to a .lis file need to be supplied.n", call.=FALSE)
  } else {
  if (is.na(args[3])){
  	args[3] = "metaatlas"
  	}
  }

if (is.na(args[4])){
	args[4] = 0.90
	}

print(paste("directory of datasets to be included in metaatlas:", args[1]))
print(paste(".lis file listing datasets:", args[2]))
print(paste("prefix for metaatlas output files:", args[3]))
print(paste("gene score pctl minimum for cluster marker filtration:", args[4]))

#########
# SET UP
#########

setwd(args[1])
dataset.lis <- read.table(args[2])

library(ggplot2)
library(tidyverse)
library(WGCNA)
library(factoextra)
library(dynamicTreeCut)
library(data.table)

#########
# AGGREGATION OF INDIVIDUAL CLUSTER MARKERS AND GENE SCORES
#########

## Aggregate the cluster markers and gene scores of individuals from the same dataset

metaatlas.genescores <- NULL
for (i in dataset.lis$V1){
	dataset.folder <- i
	for_aggregation <- list.files(dataset.folder, pattern = "^genescores_") # generate a list of all the gene score tables for a given dataset
	for_aggregation <- paste0(dataset.folder, "/", for_aggregation) # extract full file path of these gene score tables
	print(for_aggregation)
	
	# aggregate the gene scores of all the individuals from the same dataset
	dataset.genescores <- NULL
	for (j in 1:length(for_aggregation)){
		indvd <- readRDS(for_aggregation[j])
		print(for_aggregation[j])
		dataset.genescores <- rbind(dataset.genescores, indvd)
		rm(indvd)
		}

    	#change all infinite gene score values to 5000. Infinite gene scores arise when one gene is highly specific to one cluster, such that expression in other clusters is not detectable. 
    	dataset.genescores$GeneScore[is.infinite(dataset.genescores$GeneScore)] <- 5000
	
	# create a column that shows "Dataset" to use for downstream analysis
	dataset.genescores$Dataset <- dataset.folder
	
	# optional but highly recommended – validate the succesful transfer of individual gene scores to aggregated gene score table	
	
	print(paste("validating aggregation of", dataset.folder, "gene scores"))
	counter = 0
	for (k in for_aggregation){
		print(k)
		indvd.genescore <- readRDS(k)
		indvd_in_dataset.genescores <- subset(dataset.genescores, startsWith(Indvd_Cluster_ID, gsub(pattern = paste0(dataset.folder, "/genescores_"), replacement = "", k)))
		if(
			nrow(indvd.genescore) == nrow(indvd_in_dataset.genescores) &
			sum(indvd.genescore$GeneScore == indvd_in_dataset.genescores$GeneScore)
			){
				print(paste(k, "succesfully transferred to aggregate gene score table"))
				print("# clusters:")
				print(max(as.numeric(indvd.genescore$cluster)))
				} else {
					counter = counter + 1 # to tally how many times an individual gene score table was not successfully integrated
					if (
						sum(sort(unique(as.character(indvd.genescore$cluster))) == sort(unique(as.character(indvd_in_dataset.genescores$cluster)))) == length(unique(as.character(indvd.genescore$cluster))) 
						){
							print(paste(k, "not succesfully transferred to aggregate gene score table"))
							} else {
							 	print(paste("Not all clusters in", k, "transferred to aggregate gene score table"))}
							 	}
  		rm(indvd.genescore, indvd_in_dataset.genescores)	
  	}
  	
  	if (counter != 0) {
		print("at least more than one error – not saved")
	} else {
		## Aggregate the cluster markers and gene scores of all datasets
		metaatlas.genescores <- rbind(metaatlas.genescores, dataset.genescores)
	}
	
	rm(for_aggregation, dataset.genescores, counter)
}

#########
# FILTRATION OF CLUSTER MARKERS BASED ON GENE SCORE
#########

# NOTE: The following assumes cluster markers will be filtered using a pan-metaatlas gene score cutoff. Your metaatlas may instead require filtration of cluster markers using dataset-specific gene score cutoffs. 

# Analyses to guide filtration cutoff

print("clusters per dataset")
print(summarise(group_by(metaatlas.genescores,Dataset), n_distinct(Indvd_Cluster_ID)))

print("cluster markers per dataset")
print(tally(group_by(metaatlas.genescores,Dataset)))
		
# distribution of gene scores in aggregated table and in each dataset
all.gene.scores <- as.data.frame(metaatlas.genescores$GeneScore)
colnames(all.gene.scores) <- "gene.scores"
print(ggplot(all.gene.scores, aes(x=gene.scores)) 
		+ geom_density() 
		+ xlim(0,15) 
		+ geom_vline(aes(xintercept = median(gene.scores)), color = "blue", linetype = "dashed", size = 0.5) 
		+ labs(title = "Distribution of Gene Scores in metaatlas")
		)

print(ggplot(metaatlas.genescores[,c("GeneScore", "Dataset")], aes(x=GeneScore)) 
		+ geom_density() 
		+ facet_wrap(~Dataset) 
		+ xlim(0,15) 
		+ geom_vline(aes(xintercept = median(GeneScore)), color = "blue", linetype = "dashed", size = 0.5) 
		+ labs(title = "Distribution of Gene Scores in each dataset")
		)		 

# Filtering of aggregated cluster marker X gene scores table
filtered_metaatlas.genescores <- filter(metaatlas.genescores, GeneScore >= quantile(metaatlas.genescores$GeneScore, as.numeric(args[4])))

# Effects of filtering on cluster/cluster marker distribution

print("clusters per dataset")
print(summarise(group_by(filtered_metaatlas.genescores,Dataset), n_distinct(Indvd_Cluster_ID)))

print("cluster markers per dataset")
print(tally(group_by(filtered_metaatlas.genescores,Dataset)))
		
# distribution of gene scores in aggregated table and in each dataset
all.gene.scores <- as.data.frame(filtered_metaatlas.genescores$GeneScore)
colnames(all.gene.scores) <- "gene.scores"
print(ggplot(all.gene.scores, aes(x=gene.scores)) 
		+ geom_density() 
		+ geom_vline(aes(xintercept = median(gene.scores)), color = "blue", linetype = "dashed", size = 0.5) 
		+ labs(title = "Distribution of Gene Scores in metaatlas post-filter")
		)

print(ggplot(filtered_metaatlas.genescores[,c("GeneScore", "Dataset")], aes(x=GeneScore)) 
		+ geom_density() 
		+ facet_wrap(~Dataset) 
		+ geom_vline(aes(xintercept = median(GeneScore)), color = "blue", linetype = "dashed", size = 0.5) 
		+ labs(title = "Distribution of Gene Scores in each dataset post-filter")
		)		

saveRDS(filtered_metaatlas.genescores, file = paste0(args[3], "_filtered.genescores.rds"))

#########
# HIERARCHICAL CLUSTERING OF CLUSTER MARKERS TO GENERATE METAMODULES
#########

# Reformat the filtered gene score table

## Trim the filtered aggregate gene score table to only contain the following columns: gene, GeneScore, Indvd_Cluster_ID
trimmed_filtered_metaatlas.genescores <- filtered_metaatlas.genescores[,c("gene", "GeneScore", "Indvd_Cluster_ID")]
head(trimmed_filtered_metaatlas.genescores)

## Use pivot wider to create a table with the following dimensions: gene x Indvd_Cluster_ID.
## After this step, you should have a table of thefiltered cluster marker genes in all individuals from all datasets, with the gene scores for those genes in each cluster for which it is a marker. 
trimmed_filtered_metaatlas.genescores <- pivot_wider(trimmed_filtered_metaatlas.genescores, names_from = Indvd_Cluster_ID, values_from = GeneScore)

## If a gene does not pass the gene score threshold in a particular cluster (column), designated that gene score to be a 0.
trimmed_filtered_metaatlas.genescores[is.na(trimmed_filtered_metaatlas.genescores)] <- 0
head(trimmed_filtered_metaatlas.genescores)

saveRDS(trimmed_filtered_metaatlas.genescores, file = paste0(args[3], "_trimmed.filtered.genescores.rds"))

# Convert the gene score table into a distance matrix
trimmed_filtered_metaatlas.genescores_matrix <- as.matrix(trimmed_filtered_metaatlas.genescores)
print("trimmed filtered gene score table converted into matrix")
trimmed_filtered_metaatlas.genescores_matrix[1:5, 1:5]

rownames(trimmed_filtered_metaatlas.genescores_matrix) <- trimmed_filtered_metaatlas.genescores_matrix[,"gene"]
trimmed_filtered_metaatlas.genescores_matrix <- trimmed_filtered_metaatlas.genescores_matrix[,-1]
class(trimmed_filtered_metaatlas.genescores_matrix) <- "numeric"
dist <- get_dist(trimmed_filtered_metaatlas.genescores_matrix, method = "pearson", stand = FALSE)
dist_matrix <- as.matrix(dist)
print("distance matrix generated")
dist_matrix[1:5, 1:5]

# Perform hierarchical clustering on the distance matrix, using average linkage
hclust <- hclust(dist, method = 'average')

# Generate metamodules
metamodules <- cutreeDynamic(hclust, cutHeight= NULL, minClusterSize= 10, method = "hybrid", pamStage = TRUE, distM = dist_matrix, deepSplit = 1)
print(" max # clusters")
print(max(metamodules)) #shows you number of metamodules
genes <- hclust$labels

# Save the results as a table with each cluster of cluster marker genes (metamodules) and the cluster marker genes within each of these (labels) 
metaatlas_metamodules <- data.frame(genes, metamodules)
print(head(metaatlas_metamodules))

# optional but highly recommended – validate that each cluster marker has been assigned to a metamodule 
validation.table <- tally(group_by(metaatlas_metamodules, metamodules))
print(validation.table)
if (sum(validation.table$n) == n_distinct(trimmed_filtered_metaatlas.genescores$gene)){
	saveRDS(metaatlas_metamodules, file=paste0(args[3], "_metamodules.rds"))
	} else {
		print("error – did not save")
		}

