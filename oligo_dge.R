library(limma)
library(tidyverse)
library(qdapTools)

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
GSE_list <- c("GSE15090") # "GSE26852", "GSE3307","GSE9397","GSE36398","GSE10760", "GSE15090" (3 groups)

for (GSE in GSE_list) {
  ## Read data and metadata
  data <-
    read.csv(
      file = sprintf("data/%s/%s.csv", GSE, GSE),
      header = T,
      row.names = 1
    )
  data <- subset(data, select = -c(Ensembl_Gene_ID, Probe_ID))  # remove ensembl gene id
  data <- data %>% distinct(Gene_Symbol, .keep_all = TRUE)  # remove duplicate gene names
  data <- data %>% remove_rownames %>% column_to_rownames(var = "Gene_Symbol")  # set Gene_Symbol as row.names
  meta <- read.csv(
    file = sprintf("meta/%s_meta.csv", GSE),
    header = T,
    row.names = 1
  )
  meta <- subset(meta, select = sampletype)  # keep only the sampletype column
  # Match meta and data orders if needed
  if (all(colnames(data) == rownames(meta)) == F) {
    idx <- match(colnames(data), rownames(meta))
    meta$index <- idx
    meta <- meta[meta$index,]
  }
  
  tmp <- data.matrix(meta)  # 2 -- FSHD, 1 -- Control // 1 -- Affected, 2 -- 	Asymptomatic, 3 -- Control
  design <- mtabulate(tmp)  # parameterization of sample type
  rownames(design) <- rownames(meta)
  if (GSE == "GSE15090") {  
    colnames(design) <- c("Affected", "Asymptomatic", "Control")  # the study with 3 groups
  } else {
    colnames(design) <- c("Control", "FSHD")  # all other studies 
  }
  
  # DEG Analysis by Limma
  fit <- lmFit(data, design)
  if (GSE == "GSE15090") {
    cont.matrix <- makeContrasts(Asymptomatic-Control, levels=design)
  } else {
    cont.matrix <- makeContrasts(FSHDvsControl=FSHD-Control, levels=design)
  }
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  top_table <- topTable(fit2, adjust="BH", sort="none", n=Inf)
  top_table <- tibble::rownames_to_column(top_table, "gene") # make Genes as a column
  
  # Determine whether genes in gene signatures are differentially expressed
  for (i in 1:length(sig_df_list)) {
    table <-  top_table %>% filter(gene %in% rownames(sig_df_list[[i]]))  # filter for genes in signature
    idx <- match(rownames(sig_df_list[[i]]), table$gene)
    de <- as.integer(table[, 'adj.P.Val'][idx] < 0.05) # DE genes as 1 and non-DE genes as 0
    pval <- table[, 'adj.P.Val'][idx]  # record the p-value
    sig_df_list[[i]][[GSE]] <- de  # Add to dataframe
    pval_df_list[[i]][[GSE]] <- pval
  }
}

# Write to Excel
dataset_names <- list('DUX4' = sig_df_list[[1]], 'D4Z4' = sig_df_list[[2]], 'PAX7' = sig_df_list[[3]])
#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(dataset_names, file = 'DE_counts_oligo.xlsx', rowNames = TRUE)

# Write to Excel
dataset_names <- list('DUX4' = pval_df_list[[1]], 'D4Z4' = pval_df_list[[2]], 'PAX7' = pval_df_list[[3]])
#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(dataset_names, file = 'DE_pvalue_oligo.xlsx', rowNames = TRUE)
