# METAMODULE SPECIFICITY SCORING

# Arguments: 
# 1. working directory
# 2. filename of a CSV of metadata for each cell in query dataset
# 3. filename a CSV of module activity scores for each cell in query dataset (cells in rows; modules in columns)
# 4. column name of desired metadata group across which to calculate specificity scores
# 5. cutoff value for module activity binarization (default: 0.9)
# 6. output filename (needs to end in .csv)

# Execute as:
## Rscript module.specificity.score.R [directory] [csv of metadata for each cell] [csv of module activity scores] [desired column of metadata across which to calculate specificity scores] [cutoff value for module activity binarization (default: 0.9)] [output filename] 

#####
# SET UP 
#####

library(tidyverse)
library(dplyr)

# Parse arguments
args = commandArgs(trailingOnly=TRUE)

if (is.na(args[1]) |
    is.na(args[2]) |
    is.na(args[3]) |
    is.na(args[4]) |
    is.na(args[5]) |
    is.na(args[6])) {
  stop("argument missing.n", call. = FALSE)
} 

print(paste("working directory:", args[1]))
print(paste("metadata CSV", args[2]))
print(paste("module activity table CSV", args[3]))
print(paste("grouping", args[4]))
print(paste("cutoff value for binarization: ", args[5]))

#########
# SET UP
#########

setwd(args[1])
meta <- read.csv(args[2], row.names = 1, check.names = FALSE)
mod.act <- read.csv(args[3], row.names = 1, check.names = FALSE)
meta.column <- args[4]
cutoff <- as.numeric(args[5])
output <- args[6]

#####
# CLEAN UP MOD. ACT. TABLE
#####

# remove NA's
mod.act[is.na(mod.act)] <- 0

# make sure there's a dataset metadata column
if(!("Dataset" %in% colnames(meta))) {
  meta$Dataset <- "dataset"
}

mod.act <- cbind(meta[rownames(mod.act),c("Dataset", meta.column)],
                 mod.act)
colnames(mod.act)[colnames(mod.act) == meta.column] <- "meta.column"

#####
# CALCULATE AVERAGE MODULE ACTIVITY PER METADATA GROUPING
#####

summary_mod.act <- as.data.frame(summarise_at(group_by(mod.act, meta.column), colnames(mod.act)[-c(1:2)], mean))
rownames(summary_mod.act) <- summary_mod.act$meta.column

#####
# BINARIZE MODULE ACTIVITY per CELL
#####

## establish cutoff value by dataset
long_mod.act.table <- as.data.frame(
  pivot_longer(
    mod.act,
    cols = colnames(mod.act)[-c(1:2)],
    names_to = "module",
    values_to = "mod.activity"
  )
)

cutoff.by.dataset <- as.data.frame(summarise(
  group_by(long_mod.act.table, Dataset),
  "pctile" = quantile(mod.activity, cutoff)
))

for (i in unique(mod.act$Dataset)){
  print(i)
  binarized.dataset <- subset(mod.act, Dataset == i)
  binarized.dataset <-
    ifelse(binarized.dataset[, colnames(mod.act)[-c(1:2)]] >= (as.numeric(
      subset(cutoff.by.dataset, Dataset == i)$pctile
    )), 1, 0)
  
  if(max(binarized.dataset) == 1){
    if(i == unique(long_mod.act.table$Dataset)[1]){
      mod.bin <- binarized.dataset
    } else {
      mod.bin <- rbind(mod.bin, binarized.dataset)
    }
  }else{
    print("error")
    print(min(mod.bin, Dataset == i))
  }
  rm(i)
  rm(binarized.dataset)
}

# add metadata
mod.bin <- cbind(meta[rownames(mod.bin),c("Dataset", meta.column)],
                 mod.bin)


#####
# CALCULATE MODULE SPECIFICITY
#####

# Module specificity score of Module A in Cell Type 1 = [(% of cells in Cell Type 1 with > (default 90.pctl) of Module A activity) / (% of cells not in Cell Type 1 with > (default 90.pctl) of Module A activity)] * log2-fold-change of average Module A expression between Cell Type 1 cell vs all other cells in dataset

# calculate active cells per target cell group
active.cells <-
  as.data.frame(summarise_at(
    group_by(mod.bin, mod.bin[[meta.column]]),
    vars(colnames(mod.bin)[-c(1:2)]),
    sum
  ))
colnames(active.cells)[colnames(active.cells) == "mod.bin[[meta.column]]"] <- meta.column
rownames(active.cells) <- active.cells[[meta.column]]

# initiate module specificity table
mod.spec.table <- data.frame(matrix(nrow = nrow(active.cells), ncol = (ncol(mod.bin)-2)))
rownames(mod.spec.table) = rownames(active.cells)
colnames(mod.spec.table) = colnames(mod.bin)[-c(1:2)]

# calculate necessary metrics and populate
for(r in rownames(mod.spec.table)){
  
  print("######## NEW GROUP")
  print(r)
  
  for(c in colnames(mod.spec.table)){
    
    print(c)
    
    pct.active.IN.cell.type <- subset(active.cells, active.cells[[meta.column]] ==r)[[c]] / sum(active.cells[[c]])
    pct.active.OUT.cell.type <- sum(subset(active.cells, active.cells[[meta.column]] !=r)[[c]]) / sum(active.cells[[c]])
    
    avg.activity.IN.cell.type <- summary_mod.act[r,c]
    avg.activity.OUT.cell.type <- mean(subset(mod.act, mod.act$meta.column != r)[[c]])
    
    print(c(pct.active.IN.cell.type, 
            pct.active.OUT.cell.type,
            avg.activity.IN.cell.type, 
            avg.activity.OUT.cell.type))
    
    mod.spec.table[r,c] <- (pct.active.IN.cell.type / pct.active.OUT.cell.type) * log(avg.activity.IN.cell.type / avg.activity.OUT.cell.type, 2)
    
  }  
  
}

write.csv(mod.spec.table,
          file = output,
          row.names = TRUE)

# for visualization, retain only positive values
mod.spec.table[is.na(mod.spec.table)] <- 0
mod.spec.table[mod.spec.table < 0] <- 0

# resolve infinite values
max.value <- max(as.numeric(unlist(mod.spec.table[mod.spec.table != Inf])))
mod.spec.table[mod.spec.table == Inf] <- max.value

write.csv(mod.spec.table,
          file = paste0("for.visualization_", output),
          row.names = TRUE)


