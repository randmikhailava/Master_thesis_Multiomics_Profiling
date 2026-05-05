library(GeomxTools)
library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(BiocParallel)
library(NanoStringNCTools)
library(lme4)
library(pbapply)
library(lmerTest)
library(tidyr)

# 1. Load your data (already in your script)
ssgsea_pathways <- fread("/Volumes/h30492/farkkilab2/9_EyeMT/Data_analysis/geomx_processing/batch123-2808/pathway_analysis/gsea/ssgsea_norm_harmony_batch_corr_q3_norm_all_msigdb.csv")
aoi_summary <- fread("/Volumes/h30492/farkkilab2/4_CellCycle/geomx/Data/proliferative_aois.csv")


wide_ssgsea <- ssgsea_pathways %>%
  pivot_wider(
    id_cols = c(dcc_filename, Sample),
    names_from = pathway,
    values_from = ssgsea_score
  )

wide_ssgsea <- wide_ssgsea %>%
  rename(dcc_filename_1 = dcc_filename) %>%
  mutate(dcc_filename_1 = as.character(dcc_filename_1))

merged_df_ssgsea <- wide_ssgsea %>%
  inner_join(aoi_summary, by = "dcc_filename_1")


merged_df_ssgsea <- merged_df_ssgsea %>%
  filter(AOI_category %in% c("Low P / Low C", "High P / High C"))

merged_df_ssgsea$AOI_category <- factor(
  merged_df_ssgsea$AOI_category,
  levels = c("Low P / Low C", "High P / High C")
)

vars_to_exclude <- c("dcc_filename_1", "Sample", "AOI_category", "Tumor_total_count", 
                     "Tumor_Ki67_count", "Tumor_Ki67_pRb_count", "Tumor_PH3_count", 
                     "Ki67_pRb_niche_count", "Ki67_pRb_niche_prop", "Ki67_prop", "Ki67_pRb_prop", "PH3_prop", "V1")

all_pathways <- setdiff(colnames(merged_df_ssgsea), vars_to_exclude)

############## function for the lmm running #################
run_gene_test_satterthwaite <- function(gene, df) {
  tryCatch({
    y_values <- as.numeric(trimws(as.character(df[[gene]])))
    print(y_values)
    # remove the genes where the number of available values is 0
    if (sum(!is.na(y_values)) < 3) {
      return(list(
        gene = gene,
        status = "error",
        message = "Not enough non-NA values",
        result = NA
      ))
    }
    #### removing the genes where the variance is 0
    if (var(y_values, na.rm = TRUE) == 0) {
      return(list(
        gene = gene,
        status = "error",
        message = "Zero variance",
        result = NA
      ))
    }
    
    # Fit model
    mod <- lmer(y_values ~ AOI_category + (1 | Sample), data = df)
    summary_mod <- summary(mod)$coefficients
    
    # Extract contrast
    target_row <- grep("AOI_categoryHigh P / High C", rownames(summary_mod), value = TRUE)
    
    if (length(target_row) == 0) {
      return(list(
        gene = gene,
        status = "error",
        message = "Contrast not found",
        result = NA
      ))
    }
    
    stats <- summary_mod[target_row, ]
    return(list(
      gene = gene,
      status = "success",
      message =NA ,
      result = c(
        Estimate = stats["Estimate"],
        Std_Error = stats["Std. Error"],
        df = stats["df"],
        t_value = stats["t value"],
        P_Value = stats["Pr(>|t|)"]
      )
    ))
    
  }, error = function(e) {
    return(list(
      gene = gene,
      status = "error",
      message = e$message,
      result = NA
    ))
  })
}

# --- 5. Run Analysis ---
cat("Starting LMM analysis for", length(all_pathways), "pathways...\n")


# Using pbapply for progress bar
results_list <- pblapply(all_pathways, run_gene_test_satterthwaite, df = merged_df_ssgsea)


# 1. Assign gene names to the list (ensure all_genes matches results_list length)
names(results_list) <- all_pathways


# 2. Remove NULLs and genes that returned an 'error' status
# This ensures we only try to bind successful models
clean_results <- results_list[sapply(results_list, function(x) !is.null(x) && x$status == "success")]

# 3. Extract the 'result' vector from each entry and stack them
matrix_data <- do.call(rbind, lapply(clean_results, `[[`, "result"))

# 4. Convert to a clean Data Frame
final_df <- as.data.frame(matrix_data)

# 5. Add the Gene names as a proper column (instead of just row names)
final_df$Pathway <- rownames(final_df)

# 6. Ensure all statistical columns are numeric (rbind sometimes makes them 'any')
final_df[, 1:5] <- lapply(final_df[, 1:5], as.numeric)

# 7. Add FDR correction (Essential for 14,000+ genes!)

colnames(final_df) <- c("Estimate", "Std_Error", "df", "t_value", "P_Value", "Pathway")

final_df$FDR <- p.adjust(final_df$P_Value, method = "fdr")
# 8. Reorder columns to put Gene first and sort by P-value
final_df <- final_df[, c("Pathway", "Estimate", "Std_Error", "t_value", "P_Value", "FDR")]
final_df <- final_df[order(final_df$P_Value), ]

# --- 7. Save ---
write.csv(final_df, "/Volumes/h30492/farkkilab2/4_CellCycle/geomx/GeoMx_DGE_LMM_Satterthwaite_Results_Raw_Pathways.csv", row.names = FALSE)
