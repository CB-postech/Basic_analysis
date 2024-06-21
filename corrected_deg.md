There is a logical leap in identifying cell types with large transcriptome changes across conditions based on the number of differentially expressed genes.
To overcome this, we propose a function that obtains a linear model between the number of cells and the number of DEGs for each cell type and utilizes the residuals to identify cell types with large transcriptome changes.

```R
library(Seurat)
library(ggplot2)
library(dplyr)

### Input; dataframe including celltype, num_deg, num_cell -----
meta = seurat@meta.data
celltype_list = levels(meta$celltype) # cell type vector
DEG_df = data.frame(celltype = celltype_list,
                    num_deg = rep(0, length(celltype_list)),
                    num_cell = rep(0, length(celltype_list)))

for(i in 1:length(celltype_list)) {
  DEG = get(paste0('DEG_', celltype_list[i])) # Get your DEGs
  DEG_sig = DEG[DEG$p_val_adj < 0.05 & abs(DEG$avg_log2FC) > 0.5,] # DEG threshold (user-defined)
  DEG_df$num_deg[i] <- nrow(DEG_sig) # Number of DEGs for each cell type
  DEG_df$num_cell[i] <- nrow(meta[meta$celltype == celltype_list[i],]) # Number of cells for each cell type
}

### Function -----
# "celltype_color" is a vector containing color information for each cell type used in the annotation visualization.

corrected_deg <- function(DEG_df, celltype_color) {
  # modeling linear regression
  model <- lm(num_deg ~ num_cell, data = DEG_df)
  
  # calculating residuals for each cell type
  DEG_df$residuals <- residuals(model)
  
  # ranking cell type
  DEG_df_rank <- DEG_df %>%
    arrange(desc(residuals)) %>%
    mutate(rank = row_number())
  
  # plotting rankplot
  ggplot(DEG_df_rank, aes(x = rank, y = residuals, color = celltype)) +
    geom_point(size = 3) +
    geom_line(aes(group = celltype)) +
    geom_text(aes(label = celltype), vjust = -1, hjust = 0.2, size = 3, color = "black") +
    scale_color_manual(values = celltype_color) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x = "Rank",
         y = "Corrected Number of DEGs") +
    theme_classic()
}

corrected_deg(DEG_df, celltype_color)
```

![image](https://github.com/CB-postech/Basic_analysis/assets/98519284/491f23ff-c5db-44fc-be7e-336f960673de)
