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


data <- load("GSE1124.RData")


##################### End of Pre-processing and Normalization #################


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#Check annotation slot of the data set
annotation(gse_1124)
ls("package:hgu133a.db")

columns(hgu133a.db)
keytypes(hgu133a.db)

colnames(feature_1124)[1:20]

probe_ids <- rownames(Processed_data)

colnames(feature_1124)[1] <- "Probe_ID"

head(feature_1124$Probe_ID)

#Get ALL probe-to-gene mappings (22,283 rows) in one data frame
gene_symbols <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multivals = "first"
)


# Converting mapping to data frame and renaming columns
gene_map_df <- gene_symbols %>%
               as.data.frame() %>%
               tibble::rownames_to_column("PROBE_ID") %>%
               dplyr::rename(SYMBOL = 2)


#Address many-to-one gene by summarize probe signals

#Investigating how many probes map to more than one gene -Summarize number of probes per gene symbol
duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

# Identify gene associated with many probes
duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene >1)

sum(duplicate_genes$probes_per_gene)

#Verify gene probe IDs are found in processed data
all(gene_map_df$PROBE_ID == row.names(Processed_data))

#Merge annotation (SYMBOL)
Processed_data_df <- Processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)


#Remove probes without valid gene symbol annotation
Processed_data_df <- Processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

#select only numerical data
expr_only <- Processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

############################################################################

###Collapse multiple probes per gene using average expression###

############################################################################

# limma::avereps() computes the average for probes representing the same gene
averaged_data <- limma::avereps(expr_only, ID = Processed_data_df$SYMBOL)

dim(averaged_data)

# Convert averaged expression data to matrix format
data <- as.data.frame(averaged_data)
data <- data.matrix(data)
str(data)        # Structure check
is.numeric(data) # Confirm numeric matrix

nrow(data)
##################### DIFFERENTIAL GENE EXPRESSION ANALYSIS ###################

## ------------------------------------------------------------------
## 7) Phenotype grouping (diagnosis: Normal / Severe Malaria / Cerebral Malaria)
## ------------------------------------------------------------------
# Find the characteristics column that contains 'diagnosis:'
phenotype_1124$characteristics_ch1.2

groups <- factor(phenotype_1124$characteristics_ch1.2,
                 levels = c("diesease status: healthy", 
                            "diesease status: with asymptomatic Plasmodium falciparum infection",
                            "diesease status: with uncomplicated malaria",
                            "diesease status: with malaria associated with severe anemia",
                            "diesease status: with cerebral malaria"),
                 labels = c("Healthy", 
                            "Asymptomatic", 
                            "Uncomplicated", 
                            "SevereAnemia",
                            "Cerebral"))


class(groups)


############################################################################

###Create design matrix for linear modeling###
 

############################################################################ 
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Define contrast to compare disease states vs healthy samples
contrast_matrix <- makeContrasts(
  Asymptomatic_vs_Healthy = Asymptomatic - Healthy,
  Uncomplicated_vs_Healthy = Uncomplicated - Healthy,
  SevereAnemia_vs_Healthy = SevereAnemia - Healthy,
  Cerebral_vs_Healthy = Cerebral - Healthy,
  levels = design)

## Extract DEGs for each contrast separately
contrast_names <- colnames(contrast_matrix)

# Apply contrasts and compute moderated statistics
# Fit linear model to expression data
fit_1 <- lmFit(data, design)
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)
fit_2 <- eBayes(fit_contrast, trend =TRUE)

############################################################################
### Extract list of differential expressed genes (DEGs) ###
############################################################################  
## left here
# Define function for your DEG processing pipeline
process_contrast <- function(fit_object, contrast_name) {
  # Extract DEG results
  deg_results <- topTable(fit_object, 
                          coef = contrast_name,
                          number = Inf,
                          adjust.method = "BH",
                          sort.by = "p")
  
  # Your cleaning and thresholding pipeline
  deg_results <- deg_results %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::separate(Gene, into = c("PROBEID", NA), sep = "\\.", extra = "drop", fill = "right") %>%
    tibble::column_to_rownames("PROBEID")
  
  deg_results$threshold <- as.factor(ifelse(
    deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1,  "Upregulated",
    ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
           "No")
  ))
  
  # Add contrast name as a column
  deg_results$Contrast <- contrast_name
  
  return(deg_results)
}


# Process all contrasts and combine into one dataframe
combined_results <- do.call(rbind, lapply(contrast_names, function(contrast) {
  process_contrast(fit_2, contrast)
}))

# Add row names as a proper column for the combined data
combined_results <- combined_results %>%
  tibble::rownames_to_column("PROBEID")

# -------------------------------------------------------------
# SAVE SINGLE COMBINED FILE
# -------------------------------------------------------------
write.csv(combined_results, file = "Results/All_Contrasts_DEGs_Combined.csv", row.names = FALSE)

# Optional: Also create separate files for upregulated/downregulated across all contrasts
upregulated_all <- combined_results %>% filter(threshold == "Upregulated")
downregulated_all <- combined_results %>% filter(threshold == "Downregulated")
deg_updown_all <- combined_results %>% filter(threshold %in% c("Upregulated", "Downregulated"))

write.csv(upregulated_all, file = "Results/All_Upregulated_DEGs.csv", row.names = FALSE)
write.csv(downregulated_all, file = "Results/All_Downregulated_DEGs.csv", row.names = FALSE)
write.csv(deg_updown_all, file = "Results/All_Updown_DEGs.csv", row.names = FALSE)

# Print summary
cat("Combined results summary:\n")
cat("Total rows:", nrow(combined_results), "\n")
cat("Unique contrasts:", unique(combined_results$Contrast), "\n")
cat("Upregulated genes:", nrow(upregulated_all), "\n")
cat("Downregulated genes:", nrow(downregulated_all), "\n")
cat("Total DEGs:", nrow(deg_updown_all), "\n")


png("Plots/volcano_plot.png", width = 2000, height = 1500, res = 300)
volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(P-value)",
       color = "Regulation")

print(volcano_plot)


dev.off()

# -------------------------------------------------------------
# Heatmap of Top Differentially Expressed Genes
# -------------------------------------------------------------

# Select top genes with smallest adjusted p-values
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 10)

# Subset averaged expression matrix for selected genes
heatmap_data <- data[top_genes, ]

# Generate unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

# Assign formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names

# Save heatmap as PNG
png("Plots/heatmap_top50_DEGs.png", width = 2000, height = 1500, res = 300)


# Generate heatmap without additional scaling
pheatmap(
  heatmap_data,
  scale = "none", # for already normalized data
  cluster_rows = FALSE,              # Disable row clustering
  cluster_cols = TRUE,               # Cluster samples
  show_rownames = TRUE,              # Display gene names
  show_colnames = TRUE,              # Display sample labels
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 10 Differentially Expressed Genes"
)

dev.off()



