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


# 1. Load Pathway Results
dge_pathway_results <- fread("/Volumes/h30492/farkkilab2/4_CellCycle/geomx/GeoMx_DGE_LMM_Satterthwaite_Results_Raw_Pathways.csv")

# Ensure numeric types and rename 'Gene' to 'Pathway' if necessary for clarity
dge_pathway_results[, `:=`(Estimate = as.numeric(Estimate), P_Value = as.numeric(P_Value))]
if("Gene" %in% colnames(dge_pathway_results)) setnames(dge_pathway_results, "Gene", "Pathway")

# 2. Prepare Data Subsets (Using 0.01 and 0.03 thresholds)
sig_pathways_colored <- dge_pathway_results %>% 
  filter(P_Value < 0.01 & abs(Estimate) > 0.02)

# Labeling: Top 15 Up and Top 15 Down
top_up <- dge_pathway_results %>% filter(P_Value < 0.01) %>% arrange(desc(Estimate)) %>% head(15)
top_down <- dge_pathway_results %>% filter(P_Value < 0.01) %>% arrange(Estimate) %>% head(15)
top_pathways <- bind_rows(top_up, top_down)

# 3. Build the Interactive Plot
p_interactive <- plot_ly(data = dge_pathway_results, 
                         x = ~Estimate, 
                         y = ~-log10(P_Value), 
                         text = ~Pathway, 
                         type = 'scatter', 
                         mode = 'markers',
                         marker = list(size = 6, opacity = 0.2, color = 'grey'),
                         hoverinfo = 'text',
                         name = "Non-Significant") %>%
  
  # Add colored significant pathway markers
  add_markers(data = sig_pathways_colored,
              x = ~Estimate, 
              y = ~-log10(P_Value),
              text = ~Pathway,
              marker = list(
                size = 8, 
                color = ~ifelse(Estimate > 0.02, 'firebrick', 'royalblue'), 
                opacity = 0.9
              ),
              name = "Significant Pathways") %>%
  
  # Add labels
  add_annotations(data = top_pathways,
                  x = ~Estimate,
                  y = ~-log10(P_Value),
                  text = ~Pathway,
                  showarrow = TRUE,
                  arrowhead = 1,
                  arrowsize = 0.3,
                  arrowwidth = 0.5,
                  arrowcolor = "lightgrey",
                  ax = ~ifelse(Estimate > 0, 25, -25),
                  ay = ~seq(-10, -150, length.out = nrow(top_pathways)), # Spread labels vertically
                  font = list(size = 8, color = "black"),
                  xanchor = ~ifelse(Estimate > 0, "left", "right")) %>%
  
  # Use plotly::layout to fix "unused arguments" error
  plotly::layout(
    title = list(text = "<b>Pathway Enrichment: High vs Low Proliferation</b><br><sup>Thresholds: |ES| > 0.02 & P < 0.01</sup>"),
    xaxis = list(title = "Differential Enrichment Score (Estimate)", zeroline = FALSE),
    yaxis = list(title = "-log10(P-Value)"),
    shapes = list(
      # Vertical line at -0.03
      list(type = "line", x0 = -0.02, x1 = -0.02, y0 = 0, y1 = max(-log10(dge_pathway_results$P_Value), na.rm=T), 
           line = list(color = "black", dash = "dot", width = 1)),
      # Vertical line at 0.03
      list(type = "line", x0 = 0.02, x1 = 0.02, y0 = 0, y1 = max(-log10(dge_pathway_results$P_Value), na.rm=T), 
           line = list(color = "black", dash = "dot", width = 1)),
      # Horizontal line at P=0.01
      list(type = "line", x0 = min(dge_pathway_results$Estimate, na.rm=T), x1 = max(dge_pathway_results$Estimate, na.rm=T), 
           y0 = -log10(0.01), y1 = -log10(0.01), 
           line = list(color = "black", dash = "dot", width = 1))
    ),
    showlegend = TRUE
  )

# Display the plot
p_interactive

# 4. Save Results
write.csv(sig_pathways_colored %>% arrange(P_Value), 
          "Significant_Pathways_HighP_vs_LowP.csv", row.names = FALSE)

print(paste("Results saved in:", getwd()))