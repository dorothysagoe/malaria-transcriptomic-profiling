##################### DIFFERENTIAL GENE EXPRESSION ANALYSIS ###################

## ------------------------------------------------------------------
## Phenotype grouping (diagnosis: Healthy / Asymptomatic / SevereAnemia / Cerebral)
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



