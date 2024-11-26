

```R
library(Seurat)
library(magrittr)
library(ggplot2)
library(ComplexHeatmap)

save_path <- '/home/sjcho/projects/lung_exercise/24_06_07_after_delete_con4_ex4/1.transcriptomic_difference/1.5.scenic/1.5.6.OddRatio_binarized_regulon/'

### load datas
so <- readRDS('/home/sjcho/projects/lung_exercise/24_06_07_after_delete_con4_ex4/wo_erythrocyte/wo_erythrocyte.rds')

pyscenic_path= '/home/sjcho/projects/lung_exercise/24_06_07_after_delete_con4_ex4/1.transcriptomic_difference/1.5.scenic/20240828_wo_con4_ex4_wo_erythrocyte/pyscenic_output.csv'
bin_regulon <- as.data.frame(fread(pyscenic_path))
rownames(bin_regulon) <- bin_regulon[, 1] # convert 1st column to rownames
bin_regulon <- bin_regulon[, -1]

### parameter setting
# so[['lineage_exercise']] = paste0(so[['lineage']][[1]], '_', so[['Condition']][[1]])
condition_name = 'Condition' # or 'lineage'

# cols_condition = c('Lineage-Endothelial' = '#2980b9', 'Lineage-Stromal' = '#50c444', 'Lineage-Immune' = '#f1c40f', 'Lineage-Epithelial' = 'red')
cols_condition = c('Exercise-Immune' = '#00008B', 'Exercise-non-Immune' = '#8B0000', 'Control-Epithelial' = '#228B22')
# highlight other lineage

# only use when regulon is activated then 0% of total cells
# so... default : use all the cells!
percentage_cutoff = 0
top5_genes = list()

e = 10^-15 # to avoid log10(0)
for (condition in names(condition_cells_list)) {
    experimental_cells = condition_cells_list[[condition]]
    control_cells = setdiff(Cells(so), experimental_cells)

    # fisher_pvalue = get_fisher_test(control_cells, experimental_cells, bin_regulon)
    order <- get_odd_ratio(bin_regulon, control_cells, experimental_cells, percentage_cutoff)
    order$'adjusted_p-value' =  p.adjust(order$'p-value', method='fdr') # fdr correction
    order$'-log10(adj.p-value)' = -log10(order$'adjusted_p-value')
    order[order$'-log10(adj.p-value)' > 300, '-log10(adj.p-value)'] = 300 # set maximum -log10(adj.p-value) as 300
    order$size = order$'-log10(adj.p-value)' # size of dot is -log10(adj.p-value)

    order[order$'odd_ratio' == Inf, 'odd_ratio'] = 10^8 # set maximum odd ratio as 10^8
    order$log_ratio = log10(order$odd_ratio + e) # log10 of odd ratio

    basal_ratio = length(experimental_cells)/sum(length(control_cells), length(experimental_cells)) # its means ratio of experimental cells over total cells
    ## plot them separately
    for (highlihgt in names(cols_condition)) {
        p <- plot_order_text_on_right(order, highlihgt, cols_condition, basal_ratio)
        p <- p + ylim(min(order$log_ratio) * 1.15, max(order$log_ratio) * 1.15) # to visualize better
        ggsave(p, file = paste0(save_path, condition, '_RightText_', highlihgt, '.png'), width = 10, height = 5)

        p <- plot_order_text_on_left(order, highlihgt, cols_condition, basal_ratio)
        p <- p + ylim(min(order$log_ratio) * 1.15, max(order$log_ratio) * 1.15) # to visualize better
        ggsave(p, file = paste0(save_path, condition, '_LeftText_', highlihgt, '.png'), width = 10, height = 5)
        }
    top5_genes[[condition]] = rownames(order %>% tail(n = 5))
}
```

```R
get_odd_ratio <- function(bin_regulon, control_cells, experimental_cells, percentage_cutoff) {
    ### initialize
    total_cell_number = bin_regulon %>% nrow
    regulon_test_results <- data.frame(matrix(ncol = 5, nrow = ncol(bin_regulon)))
    colnames(regulon_test_results) <- c('p-value', 'adjusted_p-value', 'odd_ratio', 'sum', 'activated_percentage')
    rownames(regulon_test_results) <- colnames(bin_regulon)
    for (regulon in colnames(bin_regulon)) {            
        sum_control <- sum(bin_regulon[control_cells, regulon], na.rm = TRUE)
        sum_experimental <- sum(bin_regulon[experimental_cells, regulon], na.rm = TRUE)

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
    order[, 'module'] = col_annotation[rownames(order), 'TFs']
    return(order)
}

plot_order_text_on_right <- function(order, highlihgt, highlight_cols, basal_ratio) {
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
    return(p)
}

plot_order_text_on_left <- function(order, highlihgt, highlight_cols, basal_ratio) {
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
    return(p)
}
```

![image](https://github.com/CB-postech/Basic_analysis/assets/98519284/491f23ff-c5db-44fc-be7e-336f960673de)
