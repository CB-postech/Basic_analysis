## Plot OddRatio with fdr-corrected fisher exact test p-value
### Contact
seongju6787@postech.ac.kr

If you want any toy data please send an email

### Motivation
If a binarized regulon matrix is given, how do we quantify their enrichment in a certain condition?

### Method
<img src="https://github.com/user-attachments/assets/cccb133e-3a72-443e-ac24-1a26f7dfbaeb" width="400">

1. Calculate OddRatio from above matrix for each regulon.
In above example, (324/34) / (79/782)

2. Do fisher-exact test for each regulon, then do fdr correction
### Example Result 
Condition : Epithelial Lineage or Not

Highlight : Epithelial-Enriched Regulon from hierarchical clustering

<img src="https://github.com/user-attachments/assets/a90d259c-b6a8-4973-a27d-46969b757a97" width="700">

### Requirements
Requirement 1. Binarized Regulon Matrix (row : cell, column : regulon). Output of pySCENIC

Requirement 2. Seurat object. It should have Condition (or Lineage, Lineage&Exercise, ... something like that) in metadata

Reguirement 3. Regulon Module. Output from heirchical clustering of pHeatmap. (row : regulon, column : module)

Requirement 4. Color vector to highlight the regulon module

### There are 4 functions in this script
1. get_OddRatio : Calculate OddRatio, p-value, then return them as DataFrame
2. plot_OddRatio_text_on_right : plot OddRatio DataFrame. Text on right side
3. plot_OddRatio_text_on_left : same as 2, but text on left side
4. get_OddRatio_plot: A utility function to generate plots. It returns top 5 odd ratio TFs

### Example Script
```R
library(Seurat)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(data.table)

# load custom function saved below
source('odd_ratio.R')

save_path <- '/home/sjcho/projects/lung_exercise/24_06_07_after_delete_con4_ex4/1.transcriptomic_difference/1.5.scenic/1.5.6.OddRatio_binarized_regulon/'

### load datas
# Requirement 1. Binarized Regulon Matrix
pyscenic_path = '/home/sjcho/projects/lung_exercise/24_06_07_after_delete_con4_ex4/1.transcriptomic_difference/1.5.scenic/20240828_wo_con4_ex4_wo_erythrocyte/pyscenic_output.csv'
bin_regulon <- as.data.frame(fread(pyscenic_path))
rownames(bin_regulon) <- bin_regulon[, 1] # convert 1st column to rownames
bin_regulon <- bin_regulon[, -1]

# Requirement 2. Seurat object
so <- readRDS('/home/sjcho/projects/lung_exercise/24_06_07_after_delete_con4_ex4/wo_erythrocyte/wo_erythrocyte.rds')

# Requirement 3. Regulon Module
load('/home/sjcho/projects/lung_exercise/24_06_07_after_delete_con4_ex4/1.transcriptomic_difference/1.5.scenic/1.5.5.scenic_clean_w_function/TF_annotations.Rdata')
# save as regulon_modules

# Requirement 4. Color vector
cols_highlight = c('Exercise-Immune' = '#00008B', 'Exercise-non-Immune' = '#8B0000', 'Control-Epithelial' = '#228B22')

### parameter setting
condition_name = 'Condition' # or 'lineage' ...
percentage_cutoff = 0 # only use when regulon is activated then 0% of total cells, so... by default, use all the regulons!

top5_genes <- get_OddRatio_plot(so, bin_regulon, condition_name, cols_highlight, percentage_cutoff, regulon_modules, save_path)
# it took almost 1 min for 50000 cells
```

### Functions
```R
get_OddRatio <- function(bin_regulon, control_cells, experimental_cells, percentage_cutoff, regulon_modules) {
    ### initialize
    total_cell_number = bin_regulon %>% nrow
    regulon_test_results <- data.frame(matrix(ncol = 5, nrow = ncol(bin_regulon)))
    colnames(regulon_test_results) <- c('p-value', 'adjusted_p-value', 'odd_ratio', 'sum', 'activated_percentage')
    rownames(regulon_test_results) <- colnames(bin_regulon)

    for (regulon in colnames(bin_regulon)) {            
        sum_control <- sum(bin_regulon[control_cells, regulon], na.rm = TRUE)
        # it save the number of regulon positive control cells
        sum_experimental <- sum(bin_regulon[experimental_cells, regulon], na.rm = TRUE)
        # it save the number of regulon positive experimental cells
        
        regulon_positive = rownames(bin_regulon)[bin_regulon[, regulon] > 0]
        regulon_negative = rownames(bin_regulon)[bin_regulon[, regulon] == 0]

        mat = matrix(c(
            length(intersect(regulon_positive, control_cells)), length(intersect(regulon_positive, experimental_cells)),
            length(intersect(regulon_negative, control_cells)), length(intersect(regulon_negative, experimental_cells))
        ), nrow = 2, byrow = TRUE)
        p_value = fisher.test(mat)$p.value
        
        # ignore case where regulon is not active in both conditions
        # or very sparsely activated
        # sum_control == 0 || sum_experimental == 0) || 
        if(((sum_control + sum_experimental) / total_cell_number < percentage_cutoff)) {
            odd_ratio = NA
        }
        else {
            odd_ratio <- (sum_experimental / length(experimental_cells)) / (sum_control / length(control_cells))
        }
        regulon_test_results[regulon, ] <- c(p_value, NA, odd_ratio, sum_control + sum_experimental, (sum_control+ sum_experimental) / total_cell_number)
    }
    regulon_test_results <- regulon_test_results[!is.na(regulon_test_results[, 'odd_ratio']), ]
    regulon_test_results[, 'odd_ratio'] <- as.numeric(regulon_test_results[, 'odd_ratio'])
    order <- regulon_test_results[order(regulon_test_results[, 'odd_ratio'], decreasing = TRUE), ]
    order[, 'rank'] = 1:nrow(order)
    order[, 'TFs'] = rownames(order)
    order[, 'module'] = regulon_modules[rownames(order), 1]
    
    order$'adjusted_p-value' =  p.adjust(order$'p-value', method='fdr') # fdr correction
    order$'-log10(adj.p-value)' = -log10(order$'adjusted_p-value')
    order[order$'-log10(adj.p-value)' > 300, '-log10(adj.p-value)'] = 300 # set maximum -log10(adj.p-value) as 300
    order$size = order$'-log10(adj.p-value)' # size of dot is -log10(adj.p-value)
    order[order$'odd_ratio' == Inf, 'odd_ratio'] = 10^8 # set maximum odd ratio as 10^8
    e = 10^-8 # to avoid log10(0)
    order$log_ratio = log10(order$odd_ratio + e) # log10 of odd ratio

    return(order)
}

plot_OddRatio_text_on_right <- function(order, highlihgt, highlight_cols, basal_ratio) {
    p <- ggplot(order, aes(x = rank, y = log_ratio)) +
        geom_point(aes(size = size), color = 'black') +
        geom_point(data = subset(order, module == highlihgt), 
                    aes(size = size),
                    color = highlight_cols[[highlihgt]]) +
        geom_text_repel(data = subset(order, module == highlihgt) %>% head(n = 5), 
                        aes(label = TFs),
                        size = 5, color = highlight_cols[[highlihgt]],
                        box.padding = 1.5,
                        point.padding = 1.5,
                        force = 10,
                        nudge_x = 0.5 * length(rownames(order)),
                        nudge_y = -0.1
                        ) + 
        labs(title = "",
            x = "rank",
            y = "log10(Odd Ratio)",
            size = 'adjusted p-value') +
        geom_hline(yintercept = basal_ratio, color = 'darkred', linetype = 'dashed') +
        theme_classic() +
        theme(legend.position = "right", 
            axis.text.x = element_text(size = 15),
            axis.title.x = element_text(size = 17), 
            axis.text.y = element_text(size = 15), 
            axis.title.y = element_text(size = 17),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20)) +
        guides(size = guide_legend(override.aes = list(color = "black")))
        # scale_size_continuous(breaks = c(1, 2, 3, 4, 5)) 
    p <- p + scale_size_binned(breaks = c(0, 50, 100, 150, 200, 250, 300), range = c(0.1, 5))
    p <- p + ylim(min(order$log_ratio) * 1.15, max(order$log_ratio) * 1.15) # to visualize better
    return(p)
}

plot_OddRatio_text_on_left <- function(order, highlihgt, highlight_cols, basal_ratio) {
    p <- ggplot(order, aes(x = rank, y = log_ratio)) +
        geom_point(aes(size = size), color = 'black') +
        geom_point(data = subset(order, module == highlihgt), 
                    aes(size = size),
                    color = highlight_cols[[highlihgt]]) +
        geom_text_repel(data = subset(order, module == highlihgt) %>% head(n = 5), 
                        aes(label = TFs),
                        size = 5, color = highlight_cols[[highlihgt]],
                        box.padding = 1.5,
                        point.padding = 1.5,
                        force = 10,
                        nudge_x = -0.5 * length(rownames(order)),
                        nudge_y = 1
                        ) + 
        labs(title = "",
            x = "rank",
            y = "log10(Odd Ratio)",
            size = 'adjusted p-value') +
        geom_hline(yintercept = basal_ratio, color = 'darkred', linetype = 'dashed') +
        theme_classic() +
        theme(legend.position = "right", 
            axis.text.x = element_text(size = 15),
            axis.title.x = element_text(size = 17), 
            axis.text.y = element_text(size = 15), 
            axis.title.y = element_text(size = 17),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20)) +
        guides(size = guide_legend(override.aes = list(color = "black")))
        # scale_size_continuous(breaks = c(1, 2, 3, 4, 5)) 
    p <- p + scale_size_binned(breaks = c(0, 50, 100, 150, 200, 250, 300), range = c(0.1, 5))
    p <- p + ylim(min(order$log_ratio) * 1.15, max(order$log_ratio) * 1.15) # to visualize better
    return(p)
}


get_OddRatio_plot <- function(so, bin_regulon, condition_name, cols_highlight, percentage_cutoff, regulon_modules, save_path) {
    Idents(so) = so[[condition_name]][[1]]
    top5_genes = list()
    
    for (condition in levels(Idents(so))) {
        experimental_cells = WhichCells(so, ident = condition)
        control_cells = setdiff(Cells(so), experimental_cells)

        order <- get_OddRatio(bin_regulon, control_cells, experimental_cells, percentage_cutoff, regulon_modules)

        basal_ratio = length(experimental_cells)/length(Cells(so)) # its means ratio of experimental cells over total cells
        ## plot them separately
        for (highlihgt in names(cols_highlight)) {
            p <- plot_OddRatio_text_on_right(order, highlihgt, cols_highlight, basal_ratio)
            ggsave(p, file = paste0(save_path, condition, '_RightText_', highlihgt, '.png'), width = 10, height = 5)

            p <- plot_OddRatio_text_on_left(order, highlihgt, cols_highlight, basal_ratio)
            ggsave(p, file = paste0(save_path, condition, '_LeftText_', highlihgt, '.png'), width = 10, height = 5)
        }
        
        top5_genes[[condition]] = rownames(order %>% head(n = 5))
    }
    return(top5_genes)
}

```
