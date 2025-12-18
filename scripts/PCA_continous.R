# ============================================================
# Figure 3 PCA (continuous Syndecan-1, unified color scale)
# ============================================================

rm(list = ls())
pacman::p_load(tidyverse, ggfortify, ggrepel, readr, rio, naniar, grDevices, scales, grid, patchwork)

# ============================================================
# --- 1) Load data and define global color palette ----------
# ============================================================
merged_data <- read.csv("original_data/final_data.csv")

## ---------- Global color scale for Syndecan-1 ----------
# Create diverging palette (blue → white → red)
sdc1_all  <- log10(merged_data$Syndecan1)
sdc1_all  <- sdc1_all[is.finite(sdc1_all)]
q_sdc1    <- quantile(sdc1_all, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
vals_sdc1 <- scales::rescale(q_sdc1, from = range(q_sdc1))

pca_div_cols <- c(
  "#2D5F9A",  # deep blue (decrease)
  "#B4C5E4",  # light blue
  "#FFFFFF",  # white (negligible)
  "#F4B7A9",  # light red
  "#B22222"   # deep red (increase)
)
# ============================================================


# ============================================================
# --- 2) Coagulation Activation PCA --------------------------
# ============================================================
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

markers <- c("Syndecan1", "Antitrombin", "PROC", "Tissue factor", "Ddimer", "Platelets", "PT")
markers <- c("Syndecan1", "Tissue factor", "Ddimer", "Platelets", "PT")
markers_data <- merged_data[, markers] %>%
  mutate(across(all_of(markers), ~ log10(.))) %>%
  na.omit()

endo_pca <- prcomp(markers_data %>% select(-Syndecan1), center = TRUE, scale. = TRUE)

arrow_df <- as.data.frame(endo_pca$rotation[, 1:2])
arrow_df$label <- rownames(arrow_df)
arrow_df <- arrow_df %>%
  mutate(PC1 = PC1 * 4, PC2 = PC2 * 4) %>%
  mutate(label = recode(label,
                        "Tissue factor" = "Tissue Factor",
                        "Ddimer" = "D-dimer",
                        "Platelets" = "Platelet count",
                        "PT" = "PT",
                        "Antitrombin" = "Antithrombin",
                        "PROC" = "Protein-C"))

SDC1 <- markers_data$Syndecan1
PC1 <- endo_pca$x[, 1]; PC2 <- endo_pca$x[, 2]
cor_PC1 <- cor.test(SDC1, PC1); cor_PC2 <- cor.test(SDC1, PC2)

label_text <- paste0(
  "PC1 vs Syndecan-1: r = ", round(cor_PC1$estimate, 2),
  ", p = ", format(cor_PC1$p.value, scientific = TRUE, digits = 2), "\n",
  "PC2 vs Syndecan-1: r = ", round(cor_PC2$estimate, 2),
  ", p = ", format(cor_PC2$p.value, scientific = TRUE, digits = 2)
)

p_coa <- ggplot(as.data.frame(endo_pca$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = markers_data$Syndecan1), shape = 21, size = 4, colour = "black") +
  scale_fill_gradientn(
    colours = pca_div_cols,
    values  = vals_sdc1,
    limits  = range(q_sdc1),
    oob     = scales::squish,
    name    = "Syndecan-1\n(low → high)"
  ) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = arrow_df,
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = label),
            color = "black", size = 4, family = "sans") +
  annotate("text", x = min(PC1), y = min(PC2), label = label_text,
           hjust = 0, vjust = 0, size = 4, family = "sans") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "Coagulation Activation",
       x = "PC1 (37.4%)", y = "PC2 (24.3%)")



# ============================================================
# --- 3) Systemic Inflammation and Organ Damage PCA ----------
# ============================================================
inflam_markers <- c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1",
                    "IL10", "IL18", "IL1ra", "IL6", "IL8", "Syndecan1")

inflam_data <- merged_data[, inflam_markers] %>%
  mutate(across(all_of(inflam_markers), ~ log10(.))) %>%
  na.omit()

inflam_pca <- prcomp(inflam_data %>% select(-Syndecan1), center = TRUE, scale. = TRUE)

arrow_df <- as.data.frame(inflam_pca$rotation[, 1:2])
arrow_df$label <- rownames(arrow_df)
arrow_df <- arrow_df %>%
  mutate(PC1 = PC1 * 5, PC2 = PC2 * 5) %>%
  mutate(label = recode(label,
                        "MMP8" = "MMP-8",
                        "TenascinC" = "Tenascin C",
                        "TREM1" = "sTREM-1",
                        "IL10" = "IL-10",
                        "IL18" = "IL-18",
                        "IL1ra" = "IL-1RA",
                        "IL6" = "IL-6",
                        "IL8" = "IL-8"))

SDC1 <- inflam_data$Syndecan1
PC1 <- inflam_pca$x[, 1]; PC2 <- inflam_pca$x[, 2]
cor_PC1 <- cor.test(SDC1, PC1); cor_PC2 <- cor.test(SDC1, PC2)

label_text <- paste0(
  "PC1 vs Syndecan-1: r = ", round(cor_PC1$estimate, 2),
  ", p = ", format(cor_PC1$p.value, scientific = TRUE, digits = 2), "\n",
  "PC2 vs Syndecan-1: r = ", round(cor_PC2$estimate, 2),
  ", p = ", format(cor_PC2$p.value, scientific = TRUE, digits = 2)
)

arrow_labels <- arrow_df %>%
  mutate(label = ifelse(abs(PC1) > 0.4 | abs(PC2) > 0.4, label, NA))

p_inflam <- ggplot(as.data.frame(inflam_pca$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = SDC1), shape = 21, size = 4, colour = "black") +
  scale_fill_gradientn(
    colours = pca_div_cols,
    values  = vals_sdc1,
    limits  = range(q_sdc1),
    oob     = scales::squish,
    name    = "Syndecan-1\n(low → high)"
  ) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1) +
  ggrepel::geom_text_repel(data = arrow_labels %>% filter(!is.na(label)),
                           aes(x = PC1 * 1.05, y = PC2 * 1.05, label = label),
                           size = 4, family = "sans", color = "black",
                           box.padding = 0.5, segment.color = "grey30") +
  annotate("text", x = min(PC1), y = min(PC2), label = label_text,
           hjust = 0, vjust = 0, size = 4, family = "sans") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "Systemic Inflammation and Organ Damage",
       x = paste0("PC1 (", round(summary(inflam_pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(inflam_pca)$importance[2,2]*100,1), "%)"))



# ============================================================
# --- 4) Endothelial Cell Activation and Dysfunction PCA -----
# ============================================================
endo_markers <- c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1", "Syndecan1")

endo_data <- merged_data[, endo_markers] %>%
  mutate(across(all_of(endo_markers), ~ log10(.))) %>%
  na.omit()

endo_pca <- prcomp(endo_data %>% select(-Syndecan1), center = TRUE, scale. = TRUE)

arrow_df <- as.data.frame(endo_pca$rotation[, 1:2])
arrow_df$label <- rownames(arrow_df)
arrow_df <- arrow_df %>%
  mutate(PC1 = PC1 * 4, PC2 = PC2 * 4) %>%
  mutate(label = recode(label,
                        "Ang1" = "ANG1",
                        "Ang2" = "ANG2",
                        "CX3CL1" = "Fractalkine",
                        "Thrombomodulin" = "sThrombomodulin",
                        "ESM1" = "Endocan"))

SDC1 <- endo_data$Syndecan1
PC1 <- endo_pca$x[, 1]; PC2 <- endo_pca$x[, 2]
cor_PC1 <- cor.test(SDC1, PC1); cor_PC2 <- cor.test(SDC1, PC2)

label_text <- paste0(
  "PC1 vs Syndecan-1: r = ", round(cor_PC1$estimate, 2),
  ", p = ", format(cor_PC1$p.value, scientific = TRUE, digits = 2), "\n",
  "PC2 vs Syndecan-1: r = ", round(cor_PC2$estimate, 2),
  ", p = ", format(cor_PC2$p.value, scientific = TRUE, digits = 2)
)

p_endo <- ggplot(as.data.frame(endo_pca$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = SDC1), shape = 21, size = 4, colour = "black") +
  scale_fill_gradientn(
    colours = pca_div_cols,
    values  = vals_sdc1,
    limits  = range(q_sdc1),
    oob     = scales::squish,
    name    = "Syndecan-1\n(low → high)"
  ) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = arrow_df,
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = label),
            color = "black", size = 4, family = "sans") +
  annotate("text", x = min(PC1), y = min(PC2), label = label_text,
           hjust = 0, vjust = 0, size = 4, family = "sans") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "Endothelial Cell Activation and Dysfunction",
       x = paste0("PC1 (", round(summary(endo_pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(endo_pca)$importance[2,2]*100,1), "%)"))



# ============================================================
# --- 5) Combine and export ----------------------------------
# ============================================================
combined_plot <- (p_endo + p_inflam + p_coa) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "right")

combined_plot <- combined_plot + plot_annotation(tag_levels = 'A')
combined_plot

ggsave("figures/Figure3_pca_continuous_diverging.svg", combined_plot, width = 15.5, height = 6)






#the percentage of pca
rm(list = ls())
pacman::p_load(tidyverse, readr, rio)
merged_data <- read.csv("original_data/final_data.csv")

# ============================================================
# Function: calculate % contribution of each biomarker to PC1/PC2
# ============================================================
get_pca_contribution <- function(pca_obj) {
  loadings <- pca_obj$rotation[, 1:2]     # Raw PCA loadings (PC1, PC2)
  load_sq  <- loadings^2                  # Square loadings
  
  # Normalize within each PC so the contributions sum to 100%
  contrib  <- sweep(load_sq, 2, colSums(load_sq), FUN = "/") * 100
  
  contrib_df <- as.data.frame(round(contrib, 2))
  contrib_df$Marker <- rownames(contrib_df)
  contrib_df <- contrib_df[, c("Marker", "PC1", "PC2")]
  return(contrib_df)
}

endo_markers <- c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1", "Syndecan1")

endo_data <- merged_data[, endo_markers] %>%
  mutate(across(all_of(endo_markers), ~ log10(.))) %>%
  na.omit()

endo_pca <- prcomp(endo_data %>% select(-Syndecan1), center = TRUE, scale. = TRUE)

inflam_markers <- c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC",
                    "TREM1", "IL10", "IL18", "IL1ra", "IL6", "IL8", "Syndecan1")

inflam_data <- merged_data[, inflam_markers] %>%
  mutate(across(all_of(inflam_markers), ~ log10(.))) %>%
  na.omit()

inflam_pca <- prcomp(inflam_data %>% select(-Syndecan1), center = TRUE, scale. = TRUE)

colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

coa_markers <- c("Tissue factor", "Ddimer", "Platelets", "PT", "Antitrombin", "PROC", "Syndecan1")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

coa_data <- merged_data[, coa_markers] %>%
  mutate(across(all_of(coa_markers), ~ log10(.))) %>%
  na.omit()

coa_pca <- prcomp(coa_data %>% select(-Syndecan1), center = TRUE, scale. = TRUE)


# ============================================================
# Get percentage contribution tables
# ============================================================

endo_contrib   <- get_pca_contribution(endo_pca)
inflam_contrib <- get_pca_contribution(inflam_pca)
coa_contrib    <- get_pca_contribution(coa_pca)

endo_contrib
inflam_contrib
coa_contrib

export(
  list(
    Endothelial = endo_contrib,
    Inflammation = inflam_contrib,
    Coagulation = coa_contrib
  ),
  "figures/PCA_contribution_tables.xlsx"
)
