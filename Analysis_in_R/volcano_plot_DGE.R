library(GeomxTools)
library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(BiocParallel)
library(NanoStringNCTools)
library(lme4)
library(pbapply)
library(lmerTest)
library(EnhancedVolcano)
library(dplyr)
library(plotly)


#metadata <- fread("/Volumes/h30492/farkkilab2/9_EyeMT/Data/geomx/batch123/metadata/dcc_metadata_batch123_cleaned.csv")

# Fix column name
#setnames(metadata, "dcc_filename", "dcc_filename_1")


dge_limma_genes <-fread("/Volumes/h30492/farkkilab2/4_CellCycle/geomx/GeoMx_DGE_LMM_Satterthwaite_Results_Bulk.csv")

# 1. Load Data
dge_limma_genes <- fread("/Volumes/h30492/farkkilab2/4_CellCycle/geomx/GeoMx_DGE_LMM_Satterthwaite_Results_Bulk.csv")

# Ensure numeric types
dge_limma_genes[, `:=`(Estimate = as.numeric(Estimate), P_Value = as.numeric(P_Value))]

# 2. Prepare Data Subsets
# Colored hits
sig_data_colored <- dge_limma_genes %>% 
  filter(P_Value < 0.01 & abs(Estimate) > 0.3)

# Labeling: Top 25 Up and Top 25 Down
top_up <- dge_limma_genes %>% filter(P_Value < 0.01) %>% arrange(desc(Estimate)) %>% head(15)
top_down <- dge_limma_genes %>% filter(P_Value < 0.01) %>% arrange(Estimate) %>% head(15)
top_genes <- bind_rows(top_up, top_down)

# 3. Build the Interactive Plot
# Start with Background points
p_interactive <- plot_ly(data = dge_limma_genes, 
                         x = ~Estimate, 
                         y = ~-log10(P_Value), 
                         text = ~Gene, 
                         type = 'scatter', 
                         mode = 'markers',
                         marker = list(size = 5, opacity = 0.2, color = 'grey'),
                         hoverinfo = 'text',
                         name = "Non-Significant") %>%
  
  # Add colored significant markers
  add_markers(data = sig_data_colored,
              x = ~Estimate,
              y = ~-log10(P_Value),
              text = ~Gene,
              marker = list(
                size = 7, 
                color = ~ifelse(Estimate > 0.3, 'firebrick', 'royalblue'), 
                opacity = 0.8
              ),
              name = "Significant DEGs") %>%
  
  # Add labels with optimized directional arrows
  add_annotations(data = top_genes,
                  x = ~Estimate,
                  y = ~-log10(P_Value),
                  text = ~Gene,
                  showarrow = TRUE,
                  arrowhead = 0.5,
                  arrowsize = 0.2,
                  arrowwidth=1,
                  # Spread arrows further out (ax) and staggered height (ay)
                  ax = ~ifelse(Estimate > 0, 10, -10),
                  ay = ~seq(-20, -30, length.out = nrow(top_genes)), # Staggers heights to prevent text overlap
                  font = list(size = 7, color = "black"),
                  xanchor = ~ifelse(Estimate > 0, "left", "right")) %>%
  
  # Add Layout and Threshold Lines
  layout(
    title = list(text = "<b>High Proliferation vs Low Proliferation Niche</b><br><sup>LMM Satterthwaite Results</sup>"),
    xaxis = list(title = "Log2 Fold Change (Estimate)", zeroline = FALSE),
    yaxis = list(title = "-log10(P-Value)"),
    shapes = list(
      # Vertical lines
      list(type = "line", x0 = -0.3, x1 = -0.3, y0 = 0, y1 = max(-log10(dge_limma_genes$P_Value), na.rm=T), 
           line = list(color = "black", dash = "dot", width = 1)),
      list(type = "line", x0 = 0.3, x1 = 0.3, y0 = 0, y1 = max(-log10(dge_limma_genes$P_Value), na.rm=T), 
           line = list(color = "black", dash = "dot", width = 1)),
      # Horizontal line (P=0.05)
      list(type = "line", x0 = min(dge_limma_genes$Estimate), x1 = max(dge_limma_genes$Estimate), 
           y0 = -log10(0.01), y1 = -log10(0.01), 
           line = list(color = "black", dash = "dot", width = 1))
    ),
    showlegend = TRUE
  )

# 4. Final Output
p_interactive

# 5. Save Results
significant_genes <- dge_limma_genes %>% 
  filter(P_Value < 0.05 & abs(Estimate) > 0.3) %>% 
  arrange(P_Value)

write.csv(significant_genes, "Significant_DEGs_HighP_vs_LowP_Bulk.csv", row.names = FALSE)
print(paste("Results saved in:", getwd()))