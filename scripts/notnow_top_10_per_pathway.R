#need to run the 'transcriptome' code first

# Load required libraries
library(biomaRt)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(tidyr)

# Connect to Ensembl (initialize `hsmart`)
hsmart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Ensure gene_name column exists in results_limma
results_limma$gene_name <- rownames(results_limma)

# Function to extract top 5 and bottom 5 genes for a specific pathway
get_top_bottom_genes <- function(results, pathway_id) {
  pathway_genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "reactome",
    values = pathway_id,  # Reactome Pathway ID
    mart = hsmart
  )
  
  # Filter results to keep only genes from the pathway
  pathway_results <- results %>%
    filter(gene_name %in% pathway_genes$hgnc_symbol | rownames(results) %in% pathway_genes$ensembl_gene_id)
  
  # Extract top 5 and bottom 5 genes
  top_genes <- pathway_results %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(5)
  bottom_genes <- pathway_results %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>% head(5)
  
  return(c(rownames(top_genes), rownames(bottom_genes)))
}

# Define Reactome Pathway IDs
syndecan_pathway_id <- "R-HSA-3000170"
extracellular_pathway_id <- "R-HSA-1474244"
hemostasis_pathway_id <- "R-HSA-109582"
platelet_degranulation_pathway_id <- "R-HSA-76005"

# Extract the correct genes per pathway (Fixed issue: Now `gene_name` is available)
genes_syndecan <- get_top_bottom_genes(results_limma, syndecan_pathway_id)
genes_extracellular <- get_top_bottom_genes(results_limma, extracellular_pathway_id)
genes_hemostasis <- get_top_bottom_genes(results_limma, hemostasis_pathway_id)
genes_platelet <- get_top_bottom_genes(results_limma, platelet_degranulation_pathway_id)


# Function to generate heatmaps
generate_heatmap <- function(gene_list, title, file_name) {
  
  # Extract relevant gene expression data
  syndecan_reactome <- BSI_transcriptome_2[, c("ICU_ID_from_datasource", "Group", gene_list)]
  
  # Extract Group column (use dplyr::select() explicitly to avoid namespace conflicts)
  group <- syndecan_reactome %>% dplyr::select(ICU_ID_from_datasource, Group)
  rownames(group) <- group$ICU_ID_from_datasource
  
  # Remove non-gene columns
  gene_data <- syndecan_reactome[, -c(1:2)]
  rownames(gene_data) <- syndecan_reactome$ICU_ID_from_datasource
  
  # Z-score normalization
  gene_data_normalized <- scale(gene_data)
  
  # Convert to data frame
  gene_data_normalized <- as.data.frame(gene_data_normalized)
  
  # Convert to matrix and transpose
  gene_matrix <- t(as.matrix(gene_data_normalized))
  
  # Create annotation for groups
  annotation <- data.frame(Group = factor(group$Group))
  rownames(annotation) <- rownames(gene_data_normalized)
  
  # Define annotation colors
  ann_colors <- list(Group = c("1" = "#17becf", "4" = "#ff7f0e"))  # Cyan for Group1, Orange for Group4
  
  # Define color scale
  colors <- colorRampPalette(c("blue", "white", "red"))(200)
  
  # Function to reorder columns within each group
  reorder_within_group <- function(gene_matrix, annotation, group_label) {
    group_indices <- which(annotation$Group == group_label)
    group_matrix <- gene_matrix[, group_indices]
    group_means <- colMeans(group_matrix)
    ordered_indices <- group_indices[order(-group_means)]
    return(ordered_indices)
  }
  
  # Reorder columns for each group
  ordered_indices_group_1 <- reorder_within_group(gene_matrix, annotation, "1")
  ordered_indices_group_4 <- reorder_within_group(gene_matrix, annotation, "4")
  
  # Combine ordered indices
  final_ordered_indices <- c(ordered_indices_group_1, ordered_indices_group_4)
  
  # Reorder the gene matrix and annotation dataframe
  gene_matrix_ordered <- gene_matrix[, final_ordered_indices]
  annotation_ordered <- annotation[final_ordered_indices, , drop = FALSE]
  
  # Define custom breaks to center at 0
  breaks <- c(seq(-2, 0, length.out = 100), seq(0, 3, length.out = 100)[-1])
  
  # Create and save heatmap
  p <- pheatmap(
    gene_matrix_ordered, 
    cluster_cols = FALSE, 
    cluster_rows = FALSE, 
    annotation_col = annotation_ordered, 
    annotation_colors = ann_colors, 
    color = colors,
    breaks = breaks, # Ensure breaks cover range with 0 as midpoint
    annotation_names_col = TRUE,
    annotation_names_row = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = title
  )
  
  ggsave(file_name, plot = p, width = 8, height = 4)
}

# Generate heatmaps for each pathway
generate_heatmap(genes_syndecan, "Top 5 and Bottom 5 Genes - Syndecan Interactions", "heatmap_syndecan.svg")
generate_heatmap(genes_extracellular, "Top 5 and Bottom 5 Genes - Extracellular Matrix", "heatmap_extracellular.svg")
generate_heatmap(genes_hemostasis, "Top 5 and Bottom 5 Genes - Hemostasis", "heatmap_hemostasis.svg")
generate_heatmap(genes_platelet, "Top 5 and Bottom 5 Genes - Response to elevated platelet cytosolic Ca2+", "heatmap_platelet_degranulation.svg")  # NEW PATHWAY
