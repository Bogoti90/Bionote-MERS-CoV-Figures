#=============================================================================
# MERS-CoV Bionote Rapid Antigen Test Validation: Complete Figure Generation
#=============================================================================
# This code generates both main figures from the manuscript:
# "On-site Detection of MERS-CoV Infections in a Camel Slaughterhouse in Kenya 
# Using a Commercial Rapid Antigen Test"
#
# FIGURE 1: Clade comparison and validation analysis (Combined A+B panels)
# FIGURE 2: Four-panel validation analysis (Panels A, B, C, D)
#
# Author: Brian Ogoti
# Institution: University of Nairobi, Kenya
# Date: January 2025
# 
# Required data files:
# - clade_bionote_analysis.csv: Clade comparison data (A vs C strains)
# - Bionote_results.csv: Field validation study results
#
# Dependencies: tidyverse, ggplot2, cowplot, scales, pROC
#=============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(pROC)

# Set working directory
setwd("/Users/admin/Desktop/HpDiskE/LocalDisk/Manuscripts/MERS_COV Rapid Kit validation")

#=============================================================================
# UTILITY FUNCTIONS
#=============================================================================

# Safe log10 transformation function (handles zero and negative values)
safe_log10 <- function(x) {
  ifelse(is.na(x) | x <= 0, NA, log10(x))
}

# Statistical testing function for clade comparisons
perform_stat_test <- function(data, group_var, value_var) {
  # Perform Mann-Whitney U test (non-parametric)
  test_result <- wilcox.test(
    data[[value_var]][data[[group_var]] == levels(data[[group_var]])[1]],
    data[[value_var]][data[[group_var]] == levels(data[[group_var]])[2]]
  )
  
  p_value <- test_result$p.value
  
  # Format p-value for display
  if (p_value < 0.001) {
    p_text <- "p<0.001"
  } else if (p_value < 0.01) {
    p_text <- paste0("p=", sprintf("%.3f", p_value))
  } else {
    p_text <- paste0("p=", sprintf("%.3f", p_value))
  }
  
  return(list(p_value = p_value, p_text = p_text))
}

# Sensitivity analysis function for threshold-based performance evaluation
calc_sensitivity_by_threshold <- function(data, thresholds) {
  results <- data.frame()
  
  for (threshold in thresholds) {
    # Filter data above threshold and calculate detection rate
    subset_data <- data %>% 
      filter(MERS_RNA_copies_per_ml >= threshold) %>%
      summarize(total = n(), detected = sum(Bionote_positive))
    
    if (subset_data$total > 0) {
      # Calculate sensitivity and confidence intervals
      sens <- subset_data$detected / subset_data$total * 100
      binom_test <- binom.test(subset_data$detected, subset_data$total)
      ci_lower <- binom_test$conf.int[1] * 100
      ci_upper <- binom_test$conf.int[2] * 100
    } else {
      sens <- NA; ci_lower <- NA; ci_upper <- NA
    }
    
    results <- rbind(results, data.frame(
      threshold = threshold,
      log10_threshold = log10(threshold),
      sensitivity = sens,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    ))
  }
  
  return(results[!is.na(results$sensitivity), ])
}

#=============================================================================
# DATA LOADING AND PREPROCESSING
#=============================================================================
cat("Loading and preprocessing datasets...\n")

# Load clade comparison data (Figure 1A)
clade_data <- read.csv("clade_bionote_analysis.csv", header = TRUE, stringsAsFactors = FALSE)
clade_processed <- clade_data %>%
  rename(
    sample_id = `Sample.ID`,
    virus_dose = `Virus.dose..expected.TCID50.mL.`,
    test_result = `Rapid.Ag.Test.result..15.min.`,
    tcid50_measured = `TCID50.mL..VeroE6.T.titration.`,
    genome_copies = `genome.copies.mL`
  ) %>%
  mutate(
    clade = case_when(
      str_detect(virus_dose, "EMC") ~ "EMC (Clade A)",
      str_detect(virus_dose, "Kenya/9954") ~ "Kenya/9954 (Clade C)",
      TRUE ~ "Unknown"
    ),
    log10_genome_copies = safe_log10(genome_copies),
    log10_tcid50 = safe_log10(tcid50_measured),
    test_result_simplified = case_when(
      test_result %in% c("+", "faint", "very faint") ~ "Positive",
      test_result == "-" ~ "Negative",
      TRUE ~ "Negative"
    ),
    test_result_numeric = case_when(
      test_result_simplified == "Positive" ~ 0.6,
      test_result_simplified == "Negative" ~ 0.4,
      TRUE ~ NA_real_
    ),
    clade = factor(clade, levels = c("EMC (Clade A)", "Kenya/9954 (Clade C)")),
    test_result_simplified = factor(test_result_simplified, levels = c("Negative", "Positive"))
  )

# Load validation data (Figure 1B and Figure 2)
validation_data <- read.csv("Bionote_results.csv", header = TRUE, stringsAsFactors = FALSE)

# Check column structure
cat("Validation dataset columns:\n")
print(colnames(validation_data))

# Standardize column names
colnames(validation_data)[6] <- "Viral_isolation"
validation_processed <- validation_data %>%
  rename(
    No = No,
    Animal_ID = Animal_ID,
    final_PCR_result = Final_result,
    MERS_RNA_copies_per_ml = MERS_RNA_copiesPerml,
    Bionote_result = Bionote_result,
    Viral_isolation = Viral_isolation
  ) %>%
  filter(!is.na(final_PCR_result) & !is.na(Bionote_result)) %>%
  mutate(
    final_PCR_result = factor(final_PCR_result, levels = c("Negative", "Positive")),
    Bionote_result = factor(Bionote_result, levels = c("Negative", "Positive")),
    Viral_isolation = factor(Viral_isolation, levels = c("Negative", "Positive")),
    PCR_positive = ifelse(final_PCR_result == "Positive", 1, 0),
    Bionote_positive = ifelse(Bionote_result == "Positive", 1, 0),
    Viral_isolation_positive = ifelse(Viral_isolation == "Positive", 1, 0),
    MERS_RNA_copies_per_ml = as.numeric(MERS_RNA_copies_per_ml),
    log10_viral_load = safe_log10(MERS_RNA_copies_per_ml)
  )

# Prepare data for analysis
pcr_pos_data <- validation_processed %>% 
  filter(PCR_positive == 1) %>%
  arrange(MERS_RNA_copies_per_ml)

pcr_pos_with_isolation <- validation_processed %>%
  filter(PCR_positive == 1 & !is.na(Viral_isolation)) %>%
  arrange(MERS_RNA_copies_per_ml)

# Prepare combined clade data for box plots (Figure 1A)
clade_combined_data <- clade_processed %>%
  select(sample_id, clade, test_result_simplified, test_result_numeric, genome_copies, tcid50_measured) %>%
  pivot_longer(
    cols = c(genome_copies, tcid50_measured),
    names_to = "measurement_type",
    values_to = "viral_load"
  ) %>%
  mutate(
    measurement_type = case_when(
      measurement_type == "genome_copies" ~ "RNA Copies/mL",
      measurement_type == "tcid50_measured" ~ "TCID50/mL",
      TRUE ~ measurement_type
    ),
    measurement_type = factor(measurement_type, levels = c("RNA Copies/mL", "TCID50/mL")),
    log10_viral_load = safe_log10(viral_load)
  ) %>%
  filter(!is.na(log10_viral_load))

cat("Clade data loaded:", nrow(clade_combined_data), "measurements\n")
cat("Validation data loaded:", nrow(validation_processed), "samples\n")
cat("PCR-positive samples:", nrow(pcr_pos_data), "\n")

#=============================================================================
# STATISTICAL ANALYSES FOR FIGURE 1
#=============================================================================
cat("\nPerforming statistical analyses for Figure 1...\n")

# Calculate statistics for Panel A (Clade A vs Clade C within each group)
stat_results <- list()
for (mtype in c("RNA Copies/mL", "TCID50/mL")) {
  for (test_result in c("Negative", "Positive")) {
    subset_data <- clade_combined_data %>%
      filter(measurement_type == mtype, test_result_simplified == test_result)
    
    if (nrow(subset_data) > 0 && length(unique(subset_data$clade)) == 2) {
      stat_test <- perform_stat_test(subset_data, "clade", "log10_viral_load")
      stat_results[[paste(mtype, test_result, sep = "_")]] <- stat_test
    }
  }
}

# Statistical test for Panel B (Negative vs Positive)
pcr_pos_with_isolation_categorized <- pcr_pos_with_isolation %>%
  mutate(
    x_position = ifelse(Bionote_result == "Positive", 0.6, 0.4),
    result_category = case_when(
      Bionote_result == "Negative" ~ "Negative",
      Bionote_result == "Positive" & Viral_isolation == "Negative" ~ "Positive",
      Bionote_result == "Positive" & Viral_isolation == "Positive" ~ "Positive+Isolation",
      TRUE ~ "Other"
    ),
    result_category = factor(result_category,
                             levels = c("Negative", "Positive", "Positive+Isolation")),
    point_shape = case_when(
      result_category %in% c("Negative", "Positive") ~ "Regular",
      result_category == "Positive+Isolation" ~ "Isolation",
      TRUE ~ "Other"
    )
  )

stat_test_panel_b <- perform_stat_test(pcr_pos_with_isolation_categorized, "Bionote_result", "log10_viral_load")

# Print Figure 1 statistical results
cat("=== FIGURE 1 STATISTICAL RESULTS ===\n")
for (name in names(stat_results)) {
  result <- stat_results[[name]]
  cat(sprintf("%s: %s\n", name, result$p_text))
}
cat(sprintf("Panel B - Negative vs Positive: %s\n", stat_test_panel_b$p_text))

#=============================================================================
# STATISTICAL ANALYSES FOR FIGURE 2
#=============================================================================
cat("\nPerforming statistical analyses for Figure 2...\n")

# Define viral load thresholds for sensitivity analysis (Panel A)
thresholds <- c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7)
sensitivity_by_threshold <- calc_sensitivity_by_threshold(pcr_pos_data, thresholds)

# ROC curve analysis (Panel B)
roc_data <- pcr_pos_data %>% select(log10_viral_load, Bionote_positive) %>% na.omit()
roc_curve <- roc(roc_data$Bionote_positive, roc_data$log10_viral_load)
auc_value <- auc(roc_curve)

# Limit of Detection (LOD) analysis (Panel C)
fine_thresholds <- 10^seq(2, 9, by = 0.1)
detection_rates <- calc_sensitivity_by_threshold(pcr_pos_data, fine_thresholds)

# Determine optimal LOD (LOD90 - 90% detection rate)
if (max(detection_rates$sensitivity, na.rm = TRUE) >= 90) {
  detection_rates$diff_from_90 <- abs(detection_rates$sensitivity - 90)
  lod_row <- detection_rates[which.min(detection_rates$diff_from_90),]
  lod_threshold <- lod_row$threshold
  lod_detection_rate <- lod_row$sensitivity
} else {
  max_sens_row <- detection_rates[which.max(detection_rates$sensitivity),]
  lod_threshold <- max_sens_row$threshold
  lod_detection_rate <- max_sens_row$sensitivity
}

#=============================================================================
# FIGURE 1: CLADE COMPARISON AND VALIDATION ANALYSIS
#=============================================================================
cat("\nGenerating Figure 1: Clade Comparison and Validation Analysis...\n")

# Define color palettes and labels
clade_palette <- c("EMC (Clade A)" = "#373A40", "Kenya/9954 (Clade C)" = "#DC5F00")
clade_labels <- c("EMC (Clade A)" = "Clade A", "Kenya/9954 (Clade C)" = "Clade C")
result_colors <- c("Negative" = "forestgreen", "Positive" = "firebrick", "Positive+Isolation" = "firebrick")
shape_values <- c("Regular" = 16, "Isolation" = 21)

# Define LOD values
lod_clade_a_rna <- log10(2.0e7)  # Clade A RNA copies/mL
lod_clade_c_rna <- log10(9.0e7)  # Clade C RNA copies/mL
lod_clade_a_tcid50 <- log10(1.14e4)  # Clade A TCID50/mL
lod_clade_c_tcid50 <- log10(5.32e4)  # Clade C TCID50/mL
lod_validation_rna <- 6.19       # Validation study LOD (log10)

# Create annotation data for facet-specific p-values
annotation_data <- data.frame(
  measurement_type = c("RNA Copies/mL", "RNA Copies/mL", "TCID50/mL", "TCID50/mL"),
  group = c("Negative", "Positive", "Negative", "Positive"),
  x = c(1.0, 2.0, 1.0, 2.0),
  y = c(9.3, 10.0, 3.8, 7.0),
  segment_x1 = c(0.7, 1.7, 0.7, 1.7),
  segment_x2 = c(1.3, 2.3, 1.3, 2.3),
  segment_y = c(9.1, 9.8, 3.6, 6.8),
  p_text = c(
    ifelse("RNA Copies/mL_Negative" %in% names(stat_results), stat_results[["RNA Copies/mL_Negative"]]$p_text, ""),
    ifelse("RNA Copies/mL_Positive" %in% names(stat_results), stat_results[["RNA Copies/mL_Positive"]]$p_text, ""),
    ifelse("TCID50/mL_Negative" %in% names(stat_results), stat_results[["TCID50/mL_Negative"]]$p_text, ""),
    ifelse("TCID50/mL_Positive" %in% names(stat_results), stat_results[["TCID50/mL_Positive"]]$p_text, "")
  ))

# Panel A: Clade comparison with box plots, LOD lines and p-values
panel_a_figure1 <- ggplot(clade_combined_data, aes(x = test_result_simplified, y = log10_viral_load)) +
  # Box plots
  geom_boxplot(aes(group = interaction(test_result_simplified, clade)),
               fill = "white", color = "black", alpha = 1, 
               outlier.shape = NA, width = 0.6,
               position = position_dodge(0.8), size = 0.8) +
  # Error bars
  stat_boxplot(aes(group = interaction(test_result_simplified, clade)),
               geom = "errorbar", width = 0.3, size = 0.8,
               position = position_dodge(0.8), color = "black") +
  # Data points
  geom_point(aes(color = clade), size = 2, alpha = 0.8,
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  # LOD lines for RNA panel
  geom_hline(data = subset(clade_combined_data, measurement_type == "RNA Copies/mL"),
             aes(yintercept = lod_clade_a_rna), 
             linetype = "dashed", color = "#373A40", size = 1, alpha = 0.8) +
  geom_hline(data = subset(clade_combined_data, measurement_type == "RNA Copies/mL"),
             aes(yintercept = lod_clade_c_rna), 
             linetype = "dashed", color = "#DC5F00", size = 1, alpha = 0.8) +
  # LOD lines for TCID50 panel
  geom_hline(data = subset(clade_combined_data, measurement_type == "TCID50/mL"),
             aes(yintercept = lod_clade_a_tcid50), 
             linetype = "dashed", color = "#373A40", size = 1, alpha = 0.8) +
  geom_hline(data = subset(clade_combined_data, measurement_type == "TCID50/mL"),
             aes(yintercept = lod_clade_c_tcid50), 
             linetype = "dashed", color = "#DC5F00", size = 1, alpha = 0.8) +
  # Facet-specific p-value annotations
  geom_segment(data = annotation_data,
               aes(x = segment_x1, xend = segment_x2, y = segment_y, yend = segment_y),
               color = "black", size = 0.3, inherit.aes = FALSE) +
  geom_text(data = annotation_data,
            aes(x = x, y = y, label = p_text),
            size = 4, family = "Arial", fontface = "bold", 
            inherit.aes = FALSE) +
  # Aesthetics and theming
  scale_color_manual(values = clade_palette, name = "Clade", labels = clade_labels) +
  facet_wrap(~ measurement_type, scales = "free_y", ncol = 1) +
  labs(x = "Bionote Test Result", y = "Log10 Viral Load") +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = "black", size = 1),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "black")
  ) +
  scale_y_continuous(breaks = seq(1, 11, by = 1), labels = seq(1, 11, by = 1)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Panel B: Validation study with viral isolation data
panel_b_figure1 <- ggplot(pcr_pos_with_isolation_categorized,
                          aes(x = x_position, y = log10_viral_load)) +
  # Box plots
  geom_boxplot(aes(group = cut(x_position, breaks = c(0.35, 0.5, 0.65))),
               fill = "white", color = "black", alpha = 1, 
               outlier.shape = NA, width = 0.15, size = 0.8) +
  # Error bars
  stat_boxplot(aes(group = cut(x_position, breaks = c(0.35, 0.5, 0.65))),
               geom = "errorbar", width = 0.08, size = 0.8, color = "black") +
  # LOD line
  geom_hline(yintercept = lod_validation_rna, 
             linetype = "dashed", color = "#D55E00", size = 1, alpha = 0.8) +
  # P-value annotation
  annotate("segment", x = 0.35, xend = 0.65, y = 8.7, yend = 8.7,
           color = "black", size = 0.3) +
  annotate("text", x = 0.5, y = 8.9,
           label = stat_test_panel_b$p_text,
           size = 4, family = "Arial", fontface = "bold") +
  # Data points
  geom_point(aes(color = result_category, shape = point_shape,
                 fill = ifelse(result_category == "Positive+Isolation", "white",
                               ifelse(result_category == "Positive", "firebrick", "forestgreen"))),
             size = 3, alpha = 0.8, stroke = 1.5,
             position = position_jitter(width = 0.02, height = 0)) +
  scale_color_manual(values = result_colors, name = "Result") +
  scale_shape_manual(values = shape_values, guide = "none") +
  scale_fill_identity() +
  scale_x_continuous(breaks = c(0.4, 0.6), labels = c("Negative", "Positive"), limits = c(0.3, 0.7)) +
  labs(x = "Bionote Test Result", y = "Log10 RNA Copies/mL") +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = "black", size = 1),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black")
  ) +
  scale_y_continuous(breaks = seq(2, 9, by = 1), labels = seq(2, 9, by = 1)) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(16, 16, 21),
        fill = c("forestgreen", "firebrick", "white"),
        color = c("forestgreen", "firebrick", "firebrick"),
        size = 4
      ),
      ncol = 1, byrow = TRUE
    )
  )

# Combine Figure 1 panels
figure1_combined <- plot_grid(
  panel_a_figure1, panel_b_figure1,
  labels = c("A", "B"),
  label_size = 20, label_fontface = "bold", label_fontfamily = "Arial",
  label_x = c(0, -0.05), label_y = c(1, 1),
  ncol = 2, align = "h", rel_widths = c(1, 1))

figure1_final <- ggdraw(figure1_combined) +
  theme(plot.background = element_rect(fill = "white", color = NA))

#=============================================================================
# FIGURE 2: FOUR-PANEL VALIDATION ANALYSIS (CORRECTED VERSION)
#=============================================================================
cat("Generating Figure 2: Four-Panel Validation Analysis...\n")

# Define consistent theme (from clean code)
base_theme <- theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = "black", size = 1),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black")
  )

mers_palette <- c("Negative" = "forestgreen", "Positive" = "firebrick")

# Panel A: Sensitivity vs Viral Load Threshold (CORRECTED)
panel_a_figure2 <- ggplot(sensitivity_by_threshold, aes(x = log10_threshold, y = sensitivity)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#0072B2", alpha = 0.2) +
  geom_line(size = 1.2, color = "#0072B2") +
  geom_point(size = 3, color = "#0072B2") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "#D55E00", size = 1) +
  labs(x = "Log10 Viral Load Threshold", y = "Sensitivity (%)") +
  base_theme +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  scale_x_continuous(breaks = seq(2, 7, by = 1), labels = seq(2, 7, by = 1))

# Panel B: ROC Curve (CORRECTED)
panel_b_figure2 <- ggplot() +
  geom_path(
    data = data.frame(
      specificity = 1 - roc_curve$specificities,
      sensitivity = roc_curve$sensitivities
    ),
    aes(x = specificity, y = sensitivity), size = 1.2, color = "black"
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  base_theme +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  annotate("text", x = 0.6, y = 0.2,
           label = paste("AUC =", round(auc_value, 3)),
           size = 4, fontface = "bold", family = "Arial")

# Panel C: Limit of Detection Analysis (CORRECTED)
panel_c_figure2 <- ggplot(detection_rates, aes(x = log10_threshold, y = sensitivity)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#0072B2", alpha = 0.2) +
  geom_line(size = 1.2, color = "#0072B2") +
  geom_point(size = 2, alpha = 0.6, color = "#0072B2") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "#D55E00", size = 1) +
  geom_vline(xintercept = log10(lod_threshold), linetype = "dashed", color = "#D55E00", size = 1) +
  geom_point(aes(x = log10(lod_threshold), y = lod_detection_rate),
             color = "#D55E00", size = 4) +
  labs(x = "Log10 Viral Load Threshold", y = "Detection Rate (%)") +
  base_theme +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  scale_x_continuous(breaks = seq(2, 8, by = 1), labels = seq(2, 8, by = 1)) +
  annotate("text", x = log10(lod_threshold) - 1.2, y = 15,
           label = paste0("LOD90 = ", formatC(lod_threshold, format = "e", digits = 1)),
           color = "#D55E00", size = 4, fontface = "bold", family = "Arial")

# Panel D: Detection Probability Model (CORRECTED)
if (nrow(pcr_pos_data) > 5) {
  logistic_model <- glm(Bionote_positive ~ log10_viral_load,
                        family = binomial(link = "logit"), data = pcr_pos_data)
  
  pred_data <- data.frame(
    log10_viral_load = seq(min(pcr_pos_data$log10_viral_load, na.rm = TRUE),
                           max(pcr_pos_data$log10_viral_load, na.rm = TRUE),
                           length.out = 100)
  )
  pred_data$detection_prob <- predict(logistic_model, newdata = pred_data, type = "response")
  
  panel_d_figure2 <- ggplot(pred_data, aes(x = log10_viral_load, y = detection_prob)) +
    geom_line(size = 1.2, color = "#0072B2") +
    geom_point(data = pcr_pos_data,
               aes(x = log10_viral_load, y = Bionote_positive, color = Bionote_result),
               alpha = 0.7, size = 2) +
    scale_color_manual(values = mers_palette) +
    labs(x = "Log10 Viral Load", y = "Detection Probability") +
    base_theme +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = seq(2, 9, by = 1), labels = seq(2, 9, by = 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "#D55E00", size = 1)
} else {
  panel_d_figure2 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Insufficient data\nfor probability\nmodeling",
             size = 5, family = "Arial", hjust = 0.5, vjust = 0.5) +
    labs(x = "Log10 Viral Load", y = "Detection Probability") +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", color = "black", size = 1),
          axis.title = element_text(size = 13, color = "black"))
}

# Combine Figure 2 panels (CORRECTED)
figure2_combined <- plot_grid(
  panel_a_figure2, panel_b_figure2, panel_c_figure2, panel_d_figure2,
  labels = c("A", "B", "C", "D"),
  label_size = 20, label_fontface = "bold", label_fontfamily = "Arial",
  label_x = c(-0.025, -0.025, -0.025, -0.025), 
  label_y = c(1, 1, 1, 1),
  ncol = 2, nrow = 2, align = "hv")

figure2_final <- ggdraw(figure2_combined) +
  theme(plot.background = element_rect(fill = "white", color = NA))

#=============================================================================
# DISPLAY AND SAVE RESULTS
#=============================================================================
# Display both figures
cat("\n=== FIGURE GENERATION COMPLETE ===\n")
print("Figure 1: Clade Comparison and Validation Analysis")
print(figure1_final)
cat("\n")
print("Figure 2: Four-Panel Validation Analysis")
print(figure2_final)

# Print comprehensive summary
cat("\n=== COMPREHENSIVE ANALYSIS SUMMARY ===\n")
cat("Dataset Information:\n")
cat("- Clade comparison measurements:", nrow(clade_combined_data), "\n")
cat("- Validation samples analyzed:", nrow(validation_processed), "\n")
cat("- PCR-positive samples:", nrow(pcr_pos_data), "\n")
cat("- Samples with viral isolation data:", nrow(pcr_pos_with_isolation), "\n")

cat("\nFigure 1 Results:\n")
cat("- Clade A RNA LOD:", formatC(2.0e7, format = "e", digits = 1), "copies/mL\n")
cat("- Clade C RNA LOD:", formatC(9.0e7, format = "e", digits = 1), "copies/mL\n")
cat("- Clade A TCID50 LOD:", formatC(1.14e4, format = "e", digits = 2), "TCID50/mL\n")
cat("- Clade C TCID50 LOD:", formatC(5.32e4, format = "e", digits = 2), "TCID50/mL\n")

cat("\nFigure 2 Results:\n")
cat("- Overall test sensitivity:", round(mean(pcr_pos_data$Bionote_positive) * 100, 1), "%\n")
cat("- Area Under ROC Curve:", round(auc_value, 3), "\n")
cat("- Estimated Limit of Detection (LOD90):", formatC(lod_threshold, format = "e", digits = 2), "RNA copies/mL\n")
cat("- Detection rate at LOD90:", round(lod_detection_rate, 1), "%\n")

cat("\nStatistical Comparisons (Figure 1):\n")
for (name in names(stat_results)) {
  result <- stat_results[[name]]
  cat(sprintf("- %s: %s\n", name, result$p_text))
}

# Optional: Save figures (uncomment to use)
cat("\n=== SAVING FIGURES ===\n")
# Save Figure 1
ggsave("/Users/admin/Desktop/HpDiskE/LocalDisk/Manuscripts/MERS_COV Rapid Kit validation/Figure1_Clade_Comparison_Validation.png",
       figure1_final, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("/Users/admin/Desktop/HpDiskE/LocalDisk/Manuscripts/MERS_COV Rapid Kit validation/Figure1_Clade_Comparison_Validation.pdf",
       figure1_final, width = 12, height = 8, bg = "white")

# Save Figure 2
ggsave("/Users/admin/Desktop/HpDiskE/LocalDisk/Manuscripts/MERS_COV Rapid Kit validation/Figure2_Four_Panel_Validation.png",
       figure2_final, width = 12, height = 10, dpi = 300, bg = "white")
ggsave("/Users/admin/Desktop/HpDiskE/LocalDisk/Manuscripts/MERS_COV Rapid Kit validation/Figure2_Four_Panel_Validation.pdf",
       figure2_final, width = 12, height = 10, bg = "white")

cat("Figures saved to:", getwd(), "\n")

#=============================================================================
# SESSION INFO FOR REPRODUCIBILITY
#=============================================================================
cat("\n=== SESSION INFORMATION FOR REPRODUCIBILITY ===\n")
print(sessionInfo())

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Both Figure 1 and Figure 2 have been generated and saved successfully.\n")
cat("All statistical analyses are documented and reproducible.\n")
cat("Data files required: clade_bionote_analysis.csv, Bionote_results.csv\n")
