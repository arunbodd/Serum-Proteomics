# Load required libraries
library(readxl)          # For reading Excel files
library(dplyr)           # For data manipulation
library(limma)           # For linear modeling and differential expression
library(ggplot2)         # For data visualization
library(tidyr)           # For reshaping data
library(pheatmap)        # For heatmap visualization
library(RColorBrewer)    # For color palettes
library(bit64)           # For handling large integer data
library(ggrepel)         # For improved text label placement in plots
library(plotly)          # For interactive plots
library(data.table)      # For data manipulation
library(htmlwidgets)     # For saving interactive widgets
library(tidyverse)       # Collection of data manipulation and visualization tools
library(edgeR)           # For RNA-Seq and proteomics analysis
library(ggvenn)          # For creating Venn diagrams
library(ComplexHeatmap)  # For complex heatmap visualization
library(circlize)        # For circular heatmaps and annotations

# Read the proteomics data
serum_prot = read.csv("../Raw_Counts_serum_proteomics.csv", header = T, sep = ",", row.names = 1)

# Transform the data: replace NA with 0, pivot longer, and round values to integers
serum_prot = serum_prot %>% 
  rownames_to_column(var = "tmp") %>%          # Convert row names to a column for pivoting
  pivot_longer(cols = -"tmp") %>%              # Transform wide data to long format
  mutate(new_value = if_else(is.na(value), 0, value)) %>%  # Replace NA with 0
  mutate(val_integers = ceiling(new_value))    # Round values to the nearest integer

# Reshape the data back to wide format, using integer values
serum_prot_integerdat = serum_prot %>% 
  pivot_wider(names_from = "name", values_from = "val_integers", id_cols = "tmp") %>%
  column_to_rownames(var = "tmp")  # Restore row names

# Rename columns to remove unnecessary prefixes
colnames(serum_prot_integerdat) = gsub(".*Sample..", "", colnames(serum_prot_integerdat))

# Load metadata for the samples
metadata = read.csv("../Metadata_TimePoints.csv", header = T, sep = ",")

# Update specific column names in the dataset using metadata
colnames(serum_prot_integerdat)[13:18] = metadata$Samples[34:39]

# Ensure metadata order matches the sample order in the dataset
metadata <- metadata %>% arrange(match(Samples, colnames(serum_prot_integerdat)))
rownames(metadata) = metadata[,1]  # Set the sample names as row names
all(colnames(serum_prot_integerdat) == rownames(metadata))  # Verify alignment between metadata and dataset

# ---- Before Outlier Removal Plots and Normalization ----

# Total Ion Current (TIC) normalization function
tic_normalize <- function(data) {
  totals <- colSums(data)  # Calculate column sums (Total Ion Current)
  normalized_data <- sweep(data, 2, totals, FUN = "/") * 1e6  # Normalize and scale to 1e6
  return(normalized_data)
}

# Apply TIC normalization to the original data
original_filtered_norm_data <- tic_normalize(serum_prot_integerdat)

# Filter out proteins with a high proportion of zeros
threshold <- 0.8  # Set threshold for maximum allowable proportion of zeros
zero_proportion <- rowSums(original_filtered_norm_data == 0) / ncol(original_filtered_norm_data)
original_filtered_data <- original_filtered_norm_data[zero_proportion <= threshold, ]

# Apply log2 normalization to the filtered data
original_log2_normalized_data <- log2(original_filtered_data + 1)

# Perform PCA on log2-normalized data
pca <- prcomp(t(original_log2_normalized_data), scale. = FALSE)
pca_data <- as.data.frame(pca$x)  # Extract PCA results
pca_data$Samples <- rownames(pca_data)  # Add sample names for plotting
pca_data <- merge(pca_data, metadata, by.x = "Samples", by.y = "Samples")  # Merge with metadata

# Calculate variance explained by each principal component
explained_variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)

# Plot PCA results
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Samples)) +
  geom_point(size = 3) +                # Plot data points
  geom_text_repel(size = 3) +           # Add text labels to points
  labs(
    title = "Before Outlier Removal and Before Voom Normalization",
    x = paste0("PC1 (", explained_variance[1], "% variance)"),
    y = paste0("PC2 (", explained_variance[2], "% variance)")
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# ---- After Normalization Before Outlier Removal ----

# Create the design matrix for the linear model
design_org <- model.matrix(~ 0 + metadata$Group)  # Create dummy variables for groups
colnames(design_org) <- c("Control", "Patient_12", "Patient_24", "Patient_Base")  # Rename columns

# Create DGEList object for differential expression analysis
dge_org <- DGEList(counts = original_filtered_data)

# Apply `voom` normalization to transform data for linear modeling
v_org <- voom(dge_org, design_org, plot = TRUE, normalize.method = "quantile")

# Perform PCA on voom-normalized data
pca <- prcomp(t(v_org$E), scale. = FALSE)
pca_data <- as.data.frame(pca$x)  # Extract PCA results
pca_data$Samples <- rownames(pca_data)  # Add sample names
pca_data <- merge(pca_data, metadata, by.x = "Samples", by.y = "Samples")  # Merge with metadata

# Calculate variance explained by each principal component
explained_variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)

# Plot PCA results
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Samples)) +
  geom_point(size = 3) +                # Plot data points
  geom_text_repel(size = 3) +           # Add text labels to points
  labs(
    title = "Before Outlier Removal and After Voom Normalization",
    x = paste0("PC1 (", explained_variance[1], "% variance)"),
    y = paste0("PC2 (", explained_variance[2], "% variance)")
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# ---- After Outlier Removal ----

# Identify outlier samples to remove
samples_to_remove <- c("E1", "E2", "E3", "HS6")

# Filter out these samples from the gene matrix and metadata
filtered_gene_matrix <- serum_prot_integerdat[, !(colnames(serum_prot_integerdat) %in% samples_to_remove)]
filtered_metadata <- metadata[!(rownames(metadata) %in% samples_to_remove), ]

# Ensure alignment between the filtered data and metadata
stopifnot(all(colnames(filtered_gene_matrix) == rownames(filtered_metadata)))

# Apply TIC normalization to the filtered data
filtered_norm_data <- tic_normalize(filtered_gene_matrix)

# Generate a boxplot to visualize scaled data distribution
boxplot(filtered_norm_data, las = 2, col = "lightblue", main = "Scaled Data Distribution")

# Compute pairwise distances between samples and plot a heatmap
dist_matrix <- dist(t(filtered_norm_data))
dist_matrix <- as.matrix(dist_matrix)
pheatmap(dist_matrix, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         main = "Sample Distance Heatmap")

# Filter out proteins with a high proportion of zeros from normalized data
zero_proportion <- rowSums(filtered_norm_data == 0) / ncol(filtered_norm_data)
filtered_data <- filtered_norm_data[zero_proportion <= threshold, ]

# Apply log2 normalization to the filtered data
log2_normalized_data <- log2(filtered_data + 1)

# ---- Perform PCA After Outlier Removal ----

# Perform PCA on log2-normalized data after outlier removal
pca <- prcomp(t(log2_normalized_data), scale. = FALSE)

# Extract PCA results and add sample names
pca_data <- as.data.frame(pca$x)
pca_data$Samples <- rownames(pca_data)

# Merge PCA results with filtered metadata
pca_data <- merge(pca_data, filtered_metadata, by.x = "Samples", by.y = "Samples")

# Calculate the percentage of variance explained by each principal component
explained_variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)

# Plot PCA results after outlier removal
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, shape = Group, label = Samples)) +
  geom_point(size = 3) +  # Plot data points
  geom_text_repel(size = 3) +  # Add sample labels
  labs(
    title = "PCA of TIC Normalized Data After Outlier Removal",
    x = paste0("PC1 (", explained_variance[1], "% variance)"),
    y = paste0("PC2 (", explained_variance[2], "% variance)")
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_shape_manual(values = c(15, 16, 17, 18))  # Custom shapes for groups

# ---- Create Interactive 3D PCA Plot ----

# Use Plotly to create an interactive 3D PCA plot
fig <- plot_ly(
  pca_data,
  x = ~PC1,  # PC1 values
  y = ~PC2,  # PC2 values
  z = ~PC3,  # PC3 values
  color = ~Group,  # Group information for coloring
  symbol = ~Group,  # Group information for shapes
  text = ~paste("Sample:", Samples, "<br>Group:", Group),  # Hover text
  hoverinfo = "text",  # Information to display on hover
  type = "scatter3d",  # 3D scatter plot
  mode = "markers"  # Use markers for points
) %>%
  layout(
    title = "Interactive 3D PCA of TIC Normalized Data After Outlier Removal",
    scene = list(
      xaxis = list(title = paste0("PC1 (", explained_variance[1], "% variance)")),
      yaxis = list(title = paste0("PC2 (", explained_variance[2], "% variance)")),
      zaxis = list(title = paste0("PC3 (", explained_variance[3], "% variance)"))
    )
  )

# Save the interactive 3D plot to an HTML file
#saveWidget(fig, "../Plots/3D_PCA_After_Normalized_AfterOutlierRemoval.html", selfcontained = TRUE)

# ---- Create Design Matrix and Contrasts ----

# Create the design matrix for linear modeling
design <- model.matrix(~ 0 + filtered_metadata$Group)  # Dummy variables for group comparisons
levels(filtered_metadata$Group) <- factor(unique(filtered_metadata$Group))  # Ensure correct factor levels
colnames(design) <- c("Control", "Patient_12", "Patient_24", "Patient_Base")  # Rename columns

# Define contrasts for group comparisons
contrast_matrix <- makeContrasts(
  "Control_vs_Base" = Patient_Base - Control,
  "Control_vs_12" = Patient_12 - Control,
  "Base_vs_12" = Patient_Base - Patient_12,
  "Base_vs_24" = Patient_Base - Patient_24,
  "Control_vs_24" = Patient_24 - Control,
  levels = design
)

# Print the contrast matrix for verification
print(contrast_matrix)

# ---- Differential Expression Analysis ----

# Create DGEList object for differential expression analysis
dge <- DGEList(counts = filtered_data)

# Apply voom normalization to prepare data for linear modeling
v <- voom(dge, design, plot = TRUE, normalize.method = "quantile")

# Fit a linear model using the voom-normalized data
fit <- lmFit(v, design)

# Apply contrasts to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)

# Perform empirical Bayes moderation to improve variance estimation
fit2 <- eBayes(fit2)

# Extract results for each comparison
results_base <- topTable(fit2, coef = "Control_vs_Base", number = Inf)
results_12 <- topTable(fit2, coef = "Control_vs_12", number = Inf)
results_24 <- topTable(fit2, coef = "Control_vs_24", number = Inf)

# ---- Generate Plots for Significant Genes ----

generate_plots <- function(results_data, metadata, log2_data, group_order, output_prefix, threshold = 0.05) {
  # Filter significant genes based on adjusted p-value
  significant_genes <- results_data %>%
    filter(adj.P.Val < threshold) %>%
    arrange(adj.P.Val) %>%
    pull(Gene)
  
  # Subset expression data for significant genes
  significant_expression_data <- log2_data[rownames(log2_data) %in% significant_genes, ]
  
  # Prepare data for heatmap visualization
  metadata$Group <- factor(metadata$Group, levels = group_order)
  metadata <- metadata %>% arrange(Group)  # Order metadata by group
  time_ordered_data <- significant_expression_data[, rownames(metadata)]  # Align data with metadata
  
  # Define annotations for the heatmap
  annotation_col <- data.frame(Group = metadata$Group)
  rownames(annotation_col) <- rownames(metadata)
  set2_palette <- brewer.pal(length(group_order), "Set2")
  annotation_colors <- list(Group = setNames(set2_palette, group_order))
  
  # Create and save heatmap
  pdf(paste0(output_prefix, "_Heatmap.pdf"), height = 14, width = 12)
  pheatmap(
    as.matrix(time_ordered_data),
    name = "Expression",
    cluster_rows = TRUE,  # Cluster genes
    cluster_cols = FALSE,  # Do not cluster samples
    scale = "row",  # Row-wise z-score scaling
    color = colorRampPalette(c("blue", "white", "red"))(100),
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    main = paste0("Heatmap Showing Significant Genes - ", output_prefix)
  )
  dev.off()
  
  # Prepare data for boxplots
  long_data <- significant_expression_data %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
    left_join(metadata, by = c("Sample" = "Samples"))
  
  long_data$Group <- factor(long_data$Group, levels = group_order)
  
  # Create and save boxplots
  ggsave(
    paste0(output_prefix, "_BoxPlots.pdf"),
    ggplot(long_data, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6, size = 1.5) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
      scale_fill_brewer(palette = "Set2") +
      labs(
        title = paste0("Boxplot Showing Expression Levels - ", output_prefix),
        x = "Group (Time-Point)",
        y = "Expression (Row-wise Z-Score)"
      ) +
      theme_minimal() +
      theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
  )
}

# Generate heatmaps and boxplots for comparisons
group_order <- c("Control", "Patient_Base", "Patient_12", "Patient_24")
generate_plots(results_base, filtered_metadata, log2_normalized_data, group_order, "../Plots/Base_vs_Control")
generate_plots(results_12, filtered_metadata, log2_normalized_data, group_order, "../Plots/Patient_12_vs_Control")
generate_plots(results_24, filtered_metadata, log2_normalized_data, group_order, "../Plots/Patient_24_vs_Control")

# ---- Generate Volcano Plots ----

generate_volcano_plot <- function(results, output_prefix, pvalue_threshold = 0.05, logfc_threshold = 1) {
  # Annotate significance levels
  results <- results %>%
    mutate(
      Significance = case_when(
        adj.P.Val < pvalue_threshold & logFC > logfc_threshold ~ "Significant & Upregulated",
        adj.P.Val < pvalue_threshold & logFC < -logfc_threshold ~ "Significant & Downregulated",
        TRUE ~ "Not Significant"
      ),
      Gene_Label = ifelse(adj.P.Val < pvalue_threshold & abs(logFC) > logfc_threshold, Gene, NA)
    )
  
  # Create volcano plot
  ggsave(
    paste0("../Plots/Volcano_Plot_", output_prefix, ".pdf"),
    ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
      geom_point(alpha = 0.8, size = 3) +
      scale_color_manual(
        values = c(
          "Significant & Upregulated" = "red",
          "Significant & Downregulated" = "blue",
          "Not Significant" = "gray"
        )
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      geom_text_repel(aes(label = Gene_Label), size = 3) +
      labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value"
      ) +
      theme_minimal()
  )
}

# Generate volcano plots
generate_volcano_plot(results_base, "Base_vs_Control")
generate_volcano_plot(results_12, "12_vs_Control")
generate_volcano_plot(results_24, "24_vs_Control")

# ---- Generate Venn Diagrams ----

# Define significant genes for Venn diagram
top_genes_base <- results_base %>%
  filter(adj.P.Val < 0.05) %>%
  pull(Gene)

top_genes_12 <- results_12 %>%
  filter(adj.P.Val < 0.05) %>%
  pull(Gene)

top_genes_24 <- results_24 %>%
  filter(adj.P.Val < 0.05) %>%
  pull(Gene)

# Combine gene sets into a list
gene_list <- list(
  Base_vs_Control = top_genes_base,
  Patient_12_vs_Control = top_genes_12,
  Patient_24_vs_Control = top_genes_24
)

# Create and save Venn diagram
ggsave(
  "../Plots/Venn_TopGenes_ggvenn.pdf",
  ggvenn(gene_list, fill_color = brewer.pal(3, "Set2"))
)

# ---- Save Results to Files ----

# Save unique genes and differential expression results for each comparison
write.csv(top_genes_base, "../Results_Files/Unique_Proteins_Base_vs_Control.csv", row.names = FALSE)
write.csv(top_genes_12, "../Results_Files/Unique_Proteins_12_vs_Control.csv", row.names = FALSE)
write.csv(top_genes_24, "../Results_Files/Unique_Proteins_24_vs_Control.csv", row.names = FALSE)
write.csv(results_base, "../Results_Files/Differential_Proteins_Base_vs_Control.csv", row.names = FALSE)
write.csv(results_12, "../Results_Files/Differential_Proteins_12_vs_Control.csv", row.names = FALSE)
write.csv(results_24, "../Results_Files/Differential_Proteins_24_vs_Control.csv", row.names = FALSE)
