#This project is the final AI and Bioinformatics project whereby I chose to lead a group of 5 GSE117613
#for Transcription Profiling of Malaria Biomarkers

#clear memory
gc()

#Set working directory
getwd()
setwd("C:/Users/Dorothy/Documents/AI_Omics_Internship_2025/Malaria/")

# Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery","limma","oligo","lumi","illuminaio","annotate",
  "AnnotationDbi","org.Hs.eg.db","clusterProfiler","enrichplot","DOSE", 
  "EnhancedVolcano", "AnnotationDbi","org.Hs.eg.db","clusterProfiler","enrichplot", 
))
BiocManager::install("hgu133a.db")

# CRAN helpers
install.packages(c("tidyverse","pheatmap","EnhancedVolcano","ggpubr","janitor","stringr","caret","pROC","RColorBrewer"))

#Import libraries and necessary packages
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(pheatmap) 
library(EnhancedVolcano)
library(arrayQualityMetrics)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(hgu133a.db)
library(dplyr)
library(tidyr)


## ##################################################################
#                                                                   #
##            1) Quality Control and Preprocessing                 ##
#                                                                   #
## ##################################################################

# Load data from GEO NCBI into working space
gse_1124 <- getGEO("GSE1124", GSEMatrix = TRUE)[[1]]

#Extract expression, feature and phenotype data
expression_1124 <- exprs(gse_1124)

feature_1124 <- fData(gse_1124)

phenotype_1124 <- pData(gse_1124)


#Inspect column names in phenotype data
names(phenotype_1124)

#Inspect column names in phenotype data
head(phenotype_1124$source_name_ch1)

# Check missing values in sample annotation
sum(is.na(phenotype_1124$source_name_ch1)) 


## ----------------------------------------------------------------
## 1.b) Data quality metrics BEFORE normalization 
## ----------------------------------------------------------------

  quality_metrics_raw <- arrayQualityMetrics(
    expressionset = gse_1124,
    outdir = "Results/GC_Metrics_Raw",
    force = TRUE
)

#Normalize data set using quality metrics libraries
expression_1124_norm  <- normalizeBetweenArrays(expression_1124, method = "quantile")

# Wrap back into an ExpressionSet so downstream code looks identical
normalized_data <- ExpressionSet(
  assayData = expression_1124_norm,
  phenoData = AnnotatedDataFrame(phenotype_1124),
  featureData = AnnotatedDataFrame(feature_1124)
)


## ----------------------------------------------------------------
## 1.c.) Data quality metrics AFTER normalization 
## ----------------------------------------------------------------

quality_metrics_norm <- arrayQualityMetrics(
  expressionset = normalized_data,
  outdir = "Results/GC_Metrics_Normalized",
  force = TRUE
)

#Extract normalized data and convert from matrix to dataframe
preprocessed_data <- as.data.frame(exprs(normalized_data))

dim(preprocessed_data)

## ----------------------------------------------------------------
## 1.c.) Filtration of Low transcript reads 
## ----------------------------------------------------------------

row_median <- rowMedians(as.matrix(preprocessed_data))

Median_intensity_distribution <- hist(row_median,
                                      breaks=100,
                                      freq=FALSE,
                                      main="Mean Intensity Distribution",
                                      xlim = c(0, max(2000, na.rm = TRUE)))
                                  
                                    
pdf(file = "Plots/Median Intensity Distribution.pdf", width = 8, height = 6)


dev.off()  


# Setting threshold for filtration for probes to be 3.5
threshold <- 300

# Indicating Threshold on plot
abline(v=threshold, col="red", lw=2)

indx <- row_median > threshold

filtered_data <- preprocessed_data[indx, ]

dim(filtered_data)

#Renaming columns
# Verifying if length of rows and column names match
length(colnames(filtered_data))
length(rownames(phenotype_1124))

#Overwriting column names with filtered_data
colnames(filtered_data) <- rownames(phenotype_1124)

# Viewing newly formatted filtered_data df
View(filtered_data)

# Passing filtered_data df into "Processed_data"
Processed_data <- filtered_data

dim(Processed_data)

save(Processed_data, file="GSE1124.RData")


