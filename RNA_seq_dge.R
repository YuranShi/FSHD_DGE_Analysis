# Reference: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

library(DESeq2)
library(tidyverse)
library(openxlsx)

## Read gene signatures
signatures <- read.csv( file = "gene_signatures.csv")
sig_names <- colnames(signatures)
sig_df_list <- c()
pval_df_list <- c()

for (name in sig_names) {
  sig <- signatures[[name]]
  sig <- sig[sig != ""]
  sig_df <- data.frame(sig, row.names = 1)
  pval_df <- data.frame(sig, row.names = 1)
  sig_df_list <- c(sig_df_list, list(sig_df))
  pval_df_list <- c(pval_df_list, list(pval_df))
}

## Iterate through GSE_list and use DESeq2 for differential expression
GSE_list <- c("GSE140261", "GSE115650", "GSE56787") # "GSE56787", "GSE115650", "GSE140261" 

for (GSE in GSE_list) {
  ## Design formula
  design <- ~ sampletype
  
  ## Read data and metadata
  data <- read.csv(file = sprintf("data/%s/%s_raw_counts.csv", GSE, GSE),
             header = T)
  data = data %>% distinct(Gene_Symbol, .keep_all = TRUE)  # remove duplicate gene names
  data = data %>% remove_rownames %>% column_to_rownames(var="Gene_Symbol")  # set Gene_Symbol as row.names
  meta <- read.csv(
    file = sprintf("meta/%s_meta.csv", GSE),
    header = T,
    row.names = 1
  )
  # Match meta and data orders if needed
  if (all(colnames(data) == rownames(meta)) == F) {
    idx <- match(colnames(data), rownames(meta))
    meta$index <- idx
    meta <- meta[meta$index, ]
  }
  
  ## Create DESeq object
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
  
  ## Run analysis
  dds <- DESeq(dds)
  
  ## Define contrasts, extract results table, and shrink the log2 fold changes
  contrast <- c("sampletype", "FSHD", "Control")
  res_table_unshrunken <- results(dds, contrast=contrast, alpha = 0.05)
  res_table <- lfcShrink(dds, coef = 2, res=res_table_unshrunken)
  
  summary(res_table)
  
  res_table_tb <- res_table %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  # Get top 500 most differentially expressed genes
  ordered_res <- res_table_tb[order(res_table_tb$padj),]
  write(ordered_res$gene[1:500], sprintf("%s_500DEG.csv", GSE))
  
  # Determine whether genes in gene signatures are differentially expressed
  for (i in 1:length(sig_df_list)) {
    table <- res_table_tb %>% filter(gene %in% rownames(sig_df_list[[i]]))  # filter for genes in signature
    idx <- match(rownames(sig_df_list[[i]]), table$gene)
    # Based on adjusted p-value cutoff = 0.05
    de <- as.integer(table$padj[idx] < 0.05) # DE genes as 1 and non-DE genes as 0
    pval <- table$padj[idx]  # record the p-value
    sig_df_list[[i]][[GSE]] <- de  # Add to dataframe
    pval_df_list[[i]][[GSE]] <- pval
  }
}

# Write to Excel
dataset_names <- list('DUX4' = sig_df_list[[1]], 'D4Z4' = sig_df_list[[2]], 'PAX7' = sig_df_list[[3]])
#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(dataset_names, file = 'DE_counts.xlsx', rowNames = TRUE)

# Write to Excel
dataset_names <- list('DUX4' = pval_df_list[[1]], 'D4Z4' = pval_df_list[[2]], 'PAX7' = pval_df_list[[3]])
#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(dataset_names, file = 'DE_pvalue.xlsx', rowNames = TRUE)

# 