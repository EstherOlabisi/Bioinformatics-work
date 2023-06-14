#Esther: RNA-seq Count Data Analysis

BiocManager::install(c("AnnotationHub", "DESeq2", "ensembldb", "dblyr", "EnhancedVolcano"))
library(DESeq2)
library(ensembldb)
library(dplyr)
library(AnnotationHub)
library(EnhancedVolcano)
library(tidyverse)

########DATA DOWNLOAD########
#read RNA-seq count matrix downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450. The first column contains "EntrezGeneIDs" which will be our rownames. Remove "Length" column as it is not needed 
counts_f <- "GSE60450_Lactation-GenewiseCounts.txt"
raw_counts <- read.table(counts_f, header=T, row.names = 1)
raw_counts <- raw_counts[,-1]

#column names to rownames for new dataframe
ID <- colnames(raw_counts)
col.data <- as.data.frame(ID)
#cell type column
col.data$celltype <- c("basal", "basal", "basal", "basal", "basal", "basal", "luminal", "luminal", "luminal", "luminal", "luminal", "luminal")
#status of individual
col.data$status <- c("virgin", "virgin", "pregnant", "pregnant", "lactating", "lactating", "virgin", "virgin", "pregnant", "pregnant", "lactating", "lactating")

#########DOWNLOAD ANNOTATIONS#######
#connect database
annot_hub <- AnnotationHub::AnnotationHub()

#Query the db and view the latest(last on the list) version
mm <- AnnotationHub::query(annot_hub, c("Mus musculus", "EnsDb"))
mm

#select the latest
mm <- mm[["AH109655"]]

#view annotations and other info. Some columns of interest include entrezid and description
mm_annots <- genes(mm, return.type = "data.frame")

#select specific columns
mm_annots <- genes(mm, return.type = "data.frame")%>%
  select(gene_name, description, entrezid)

  
###########PCA PLOT#######
#create DESeqDataset first
de.dataset <- DESeqDataSetFromMatrix(countData=raw_counts, colData=col.data, design= ~celltype + status)

#Transform dataset 
?rlog
rlog_data <- rlog(de.dataset)

#PCA by cell type
plotPCA(rlog_data, intgroup="celltype")
#PCA by status
plotPCA(rlog_data, intgroup="status")


########DE ANALYSIS######
#run differential expression analysis with DESeq()
de.analysis <- DESeq(de.dataset)

#view results. Set tidy to true to view as a dataframe. EntrezIDs are in column 1.
de_results <- results(de.analysis, tidy=T)


########Add annotations########
#add annotations but first convert entrezid column from list to character for successful merging.  
mm_annots$entrezid <- as.character(mm_annots$entrezid)
annotated_df <- left_join(de_results, mm_annots, by=join_by(row==entrezid))

#More observations in annotated_df than the original de_results which suggests we have duplicates. Checking...; yes, we do.
annotated_df %>%
  group_by_all() %>%
  filter(n()>1)

#remove duplicates
annotated_df <- annotated_df[!duplicated(annotated_df$row), ]

#remove NA adjusted-pvalues 
annotated_df <- annotated_df %>%
  filter(!is.na(padj))


#########VOLCANO PLOT#######
#extract names of the most significant (both up and down-regulated genes)
most_sig_genes <- top_n(annotated_df, -20, padj)$gene_name
most_sig_genes

#do volcano plot highlighting some of the most significant genes 
EnhancedVolcano(annotated_df,
               lab = annotated_df$gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
               pCutoff = 0.05,
               FCcutoff = 1.5,
               selectLab = c("Pigr", "Wap", "Pik3r1", "Glycam1", "Dsc2", "Slc24a3","Lao1","Fgfr3", "Kit", "Csn1s1","Csn2", "Csn1s2a", "Csn1s2b", "Naaa", "Sh2b2", "Atp6v1b1", "Fcgbp", "Rab11fip1"),
               title = "Differential Expression of luminal and basal mammary gland cell",
               subtitle = ""
                )

#add table for description
t_most_sig_genes <- top_n(annotated_df, -20, padj)[c(8,6,9)]
write.table(t_most_sig_genes, file = "most_sig.csv", sep=",", row.names = F)


#Links that were helpful
#https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#https://hbctraining.github.io/DGE_workshop_salmon/lessons/genomic_annotation.html
