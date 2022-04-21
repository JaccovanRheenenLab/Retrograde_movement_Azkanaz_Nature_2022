#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.13")
##BiocManager::install(c("limma","Biobase","GenomeInfoDb","edgeR","DESeq2","convert","phyloseq","DOSE","enrichplot","EnhancedVolcano","NOISeq","clusterProfiler","pathview","org.Mm.eg.db","org.Hs.eg.db","sojourner"))
##BiocManager::install("DESeq2")
#devtools::install_github("ingewortel/celltrackR")
# Loading important libraries

library(ggplot2)
library(tidyverse)
library(reshape2)
library(pipeR)
library(devtools)
library(data.table)

library(limma)
library(Biobase)
library(GenomeInfoDb)
library(RCurl)
library(edgeR)
library(DESeq2)
library(convert)
library(devtools)
library(phyloseq)

library(RColorBrewer) # Load a package giving more colors
library(pheatmap) # load a package for making heatmaps
library(GOplot)
library(enrichR)
library(DOSE)
library(enrichplot)
library(EnhancedVolcano)
library(ggpubr)
library(rstatix)
library(NOISeq)

organism = "org.Mm.eg.db"
library(clusterProfiler)
library(wordcloud)
library(pathview)
library(organism, character.only = TRUE)
library(sojourner)
library(celltrackR)
library(flowcatchR)
library(ggbeeswarm)
library(lm.beta)

# Setting working directory
getwd()
setwd("Retrograde_movement_Azkanaz_Nature_2022/Data/")

# Reading input files
raw_counts <- fread("raw_counts.csv", header = T)

#Setting conditions
condition <- c("SI_high", "SI_high", "SI_high", 
               "SI_medium", "SI_medium", "SI_medium", 
               "SI_low", "SI_low", "SI_low",
               "LI_high", "LI_high", "LI_high", 
               "LI_medium", "LI_medium", "LI_medium", 
               "LI_low", "LI_low", "LI_low")

replicates <- c("SI_high_1", "SI_high_2", "SI_high_3",
                "SI_medium_1", "SI_medium_2", "SI_medium_3",
                "SI_low_1", "SI_low_2", "SI_low_3",
                "LI_high_1", "LI_high_2", "LI_high_3",
                "LI_medium_1", "LI_medium_2", "LI_medium_3",
                "LI_low_1", "LI_low_2", "LI_low_3")




#  _       _     _       _                   _                 _     _  #
# | |     | |   (_)     | |                 | |               (_)   | | #
# | |     | |__  _  __ _| |__   __   _____  | |      _ __ ___  _  __| | #
# | |     | '_ \| |/ _` | '_ \  \ \ / / __| | |     | '_ ` _ \| |/ _` | #
# | |_____| | | | | (_| | | | |  \ V /\__ \ | |_____| | | | | | | (_| | #
# \_____(_)_| |_|_|\__, |_| |_|   \_/ |___/ \_____(_)_| |_| |_|_|\__,_| #
#                   __/ |                                               #
#                  |___/                                                #
coldata_LLarge_LMid <- data.frame(id = colnames(raw_counts[,c(12:14,18:20)]), 
                                                           condition = condition[10:15], 
                                                           replicates = replicates[10:15])

countdata_LLarge_LMid <- raw_counts[,c(1,12:17)] %>%
  textshape::column_to_rownames("gene_id") %>%
  magrittr::set_colnames(NULL)

## Conditions must be in the header fo the columns
dds_LLarge_LMid <- DESeqDataSetFromMatrix(countData = countdata_LLarge_LMid, 
                                          colData = coldata_LLarge_LMid, 
                                          design = ~ condition)

dds_LLarge_LMid_df <- as.data.frame(counts(dds_LLarge_LMid))

# Filter step to remove genes with zero total counts 
LARGE_summary_data_samples_raw_counts <- function(dds_LLarge_LMid_df){
  
  # Here we don't filter out any values
  data_noneFilter <- dds_LLarge_LMid_df %>%
    dplyr::transmute(
      genes = rownames(.),
      LI_high = rowMeans(.[,1:3]),
      LI_medium = rowMeans(.[,4:6])
    ) %>%
    summarise(n = n(),
              raw_count_mean_LI_high = mean(LI_high),
              raw_count_mean_LI_medium = mean(LI_medium)
    )
  
  ## Here we exclude if any sample has a 0 read 
  data_none0 <- dds_LLarge_LMid_df %>% 
    dplyr::mutate(
      LI_high_0 = case_when(.[,1] == 0 ~ 0,
                            .[,2] == 0 ~ 0,
                            .[,3] == 0 ~ 0, TRUE ~ 1),
      LI_medium_0 = case_when(.[,4] == 0 ~ 0,
                              .[,5] == 0 ~ 0,
                              .[,6] == 0 ~ 0, TRUE ~ 1)
    ) %>% 
    dplyr::filter((LI_high_0 + LI_medium_0) != 0) %>%
    dplyr::select(-c(LI_high_0, LI_medium_0)) %>%
    dplyr::transmute(
      genes = rownames(.),
      LI_high = rowMeans(.[,1:3]),
      LI_medium = rowMeans(.[,4:6])
    ) %>%
    summarise(n = n(),
              raw_count_mean_LI_high = mean(LI_high),
              raw_count_mean_LI_medium = mean(LI_medium)
    )
  
  
  # Here we filter out genes that have mean raw values of 0
  data_noneMean0 <- dds_LLarge_LMid_df %>%
    dplyr::transmute(
      genes = rownames(.),
      LI_high = rowMeans(.[,1:3]),
      LI_medium = rowMeans(.[,4:6])
    ) %>%
    filter(LI_high != 0 & LI_medium != 0) %>%
    summarise(n = n(),
              raw_count_mean_LI_high = mean(LI_high),
              raw_count_mean_LI_medium = mean(LI_medium)
    )
  
  ##########################################################
  # Here we filter out genes that have mean raw values of 10
  data_noneMeanLess10 <- dds_LLarge_LMid_df %>%
    dplyr::transmute(
      genes = rownames(.),
      LI_high = rowMeans(.[,1:3]),
      LI_medium = rowMeans(.[,4:6])
    ) %>%
    filter(LI_high >= 10 & LI_medium >= 10) %>%
    summarise(n = n(),
              raw_count_mean_LI_high = mean(LI_high),
              raw_count_mean_LI_medium = mean(LI_medium)
    )
  ##########################################################
  
  dataframe_results <- rbind(
    cbind(data_noneFilter, data.frame(x = "Nofilter")),
    cbind(data_none0, data.frame(x = "Anyvalueis0")),
    cbind(data_noneMean0, data.frame(x = "Meanof0")),
    cbind(data_noneMeanLess10, data.frame(x = "MeanLess10"))
  )
  
  return(dataframe_results)
}

LARGE_summary_data_samples_raw_counts(dds_LLarge_LMid_df)

dds_LLarge_LMid <- dds_LLarge_LMid[rowMeans(counts(dds_LLarge_LMid)) >= 10,]

dds_LLarge_LMid <- DESeq(dds_LLarge_LMid)

dds_LLarge_LMid$condition <- factor(dds_LLarge_LMid$condition, levels = c("LI_high", "LI_medium"))
resultsNames(dds_LLarge_LMid)[2]
res_LLarge_LMid <- results(dds_LLarge_LMid)

table_creation <- function(res_object){
  p0.05_log2FC0.5 <- data.frame(res_object) %>% 
    dplyr::filter(pvalue < 0.05 & (log2FoldChange < -0.5 | log2FoldChange > 0.5)) %>%
    dplyr::summarise(
      n = n()
    )
  
  q0.05_log2FC0.5 <- data.frame(res_object) %>% 
    dplyr::filter(padj < 0.05 & (log2FoldChange < -0.5 | log2FoldChange > 0.5)) %>%
    dplyr::summarise(
      n = n()
    )
  
  p0.05_log2FC2 <- data.frame(res_object) %>% 
    dplyr::filter(pvalue < 0.05 & (log2FoldChange < -2 | log2FoldChange > 2)) %>%
    dplyr::summarise(
      n = n()
    )
  
  q0.05_log2FC2 <- data.frame(res_object) %>% 
    dplyr::filter(padj < 0.05 & (log2FoldChange < -2 | log2FoldChange > 2)) %>%
    dplyr::summarise(
      n = n()
    )
  
  table <- data.frame(x = c("pvalue0.05_log2FC0.5", "qvalue0.05_log2FC0.5", 
                            "pvalue0.05_log2FC2", "qvalue0.05_log2FC2"),
                      y = c(as.numeric(p0.05_log2FC0.5$n), as.numeric(q0.05_log2FC0.5$n),
                            as.numeric(p0.05_log2FC2$n), as.numeric(q0.05_log2FC2$n))
  )
  return(table)
}
table_creation(res_LLarge_LMid)

# VOLCANOPLOT
genes_LI.df_LLarge_LMid <- bitr(rownames(res_LLarge_LMid), fromType = "ENSEMBL",
                                toType = c("SYMBOL"),
                                OrgDb = org.Mm.eg.db)

res_modified_LLarge_LMid <- res_LLarge_LMid %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::rename(ENSEMBL = gene_id) %>%
  dplyr::left_join(., genes_LI.df_LLarge_LMid, by = "ENSEMBL")

EnhancedVolcano(res_modified_LLarge_LMid,
                lab = res_modified_LLarge_LMid$SYMBOL,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 0.001,
                FCcutoff = 2,
                labSize = 4.0,
                shape = c(6, 6, 19, 16),
                title = "LI Lgr5 mid vs LI Lgr5 high",
                subtitle = "Differential expression analysis",
                caption = "FC cutoff, 2; p-value cutoff, 0.001",
                legendPosition = "right",
                legendLabSize = 14,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 0.9,
                hline = c(10e-8)) +
  ggplot2::coord_cartesian(xlim=c(-6, 6)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1))
dev.off()



#  _____  _     _       _                   _____            _     _  #
# /  ___|| |   (_)     | |                 /  ___|          (_)   | | #
# \ `--. | |__  _  __ _| |__   __   _____  \ `--.  _ __ ___  _  __| | #
#  `--. \| '_ \| |/ _` | '_ \  \ \ / / __|  `--. \| '_ ` _ \| |/ _` | #
# /\__/ /| | | | | (_| | | | |  \ V /\__ \ /\__/ /| | | | | | | (_| | #
# \____(_)_| |_|_|\__, |_| |_|   \_/ |___/ \____(_)_| |_| |_|_|\__,_| #
#                  __/ |                                              #
#                 |___/                                               #
coldata_SLarge_SMid <- data.frame(id = colnames(raw_counts[,3:8]), 
                                  condition = condition[1:6], 
                                  replicates = replicates[1:6])

countdata_SLarge_SMid <- raw_counts[,c(1,3:8)] %>%
  textshape::column_to_rownames("gene_id") %>%
  magrittr::set_colnames(NULL)

## Conditions must be in the header fo the columns
dds_SLarge_SMid <- DESeqDataSetFromMatrix(countData = countdata_SLarge_SMid, 
                                          colData = coldata_SLarge_SMid, 
                                          design = ~ condition)


dds_SLarge_SMid_df <- as.data.frame(counts(dds_SLarge_SMid))

# Filter step to remove genes with zero total counts 
SMALL_summary_data_samples_raw_counts <- function(dds_SLarge_SMid_df){
  
  # Here we don't filter out any values
  data_noneFilter <- dds_SLarge_SMid_df %>%
    dplyr::transmute(
      genes = rownames(.),
      SI_high = rowMeans(.[,1:3]),
      SI_medium = rowMeans(.[,4:6])
    ) %>%
    summarise(n = n(),
              raw_count_mean_SI_high = mean(SI_high),
              raw_count_mean_SI_medium = mean(SI_medium)
    )
  
  ## Here we exclude if any sample has a 0 read 
  data_none0 <- dds_SLarge_SMid_df %>% 
    dplyr::mutate(
      SI_high_0 = case_when(.[,1] == 0 ~ 0,
                            .[,2] == 0 ~ 0,
                            .[,3] == 0 ~ 0, TRUE ~ 1),
      SI_medium_0 = case_when(.[,4] == 0 ~ 0,
                              .[,5] == 0 ~ 0,
                              .[,6] == 0 ~ 0, TRUE ~ 1)
    ) %>% 
    dplyr::filter((SI_high_0 + SI_medium_0) != 0) %>%
    dplyr::select(-c(SI_high_0, SI_medium_0)) %>%
    dplyr::transmute(
      genes = rownames(.),
      SI_high = rowMeans(.[,1:3]),
      SI_medium = rowMeans(.[,4:6])
    ) %>%
    summarise(n = n(),
              raw_count_mean_SI_high = mean(SI_high),
              raw_count_mean_SI_medium = mean(SI_medium)
    )
  
  
  # Here we filter out genes that have mean raw values of 0
  data_noneMean0 <- dds_SLarge_SMid_df %>%
    dplyr::transmute(
      genes = rownames(.),
      SI_high = rowMeans(.[,1:3]),
      SI_medium = rowMeans(.[,4:6])
    ) %>%
    filter(SI_high != 0 & SI_medium != 0) %>%
    summarise(n = n(),
              raw_count_mean_SI_high = mean(SI_high),
              raw_count_mean_SI_medium = mean(SI_medium)
    )
  
  ##########################################################
  # Here we filter out genes that have mean raw values of 10
  data_noneMeanLess10 <- dds_SLarge_SMid_df %>%
    dplyr::transmute(
      genes = rownames(.),
      SI_high = rowMeans(.[,1:3]),
      SI_medium = rowMeans(.[,4:6])
    ) %>%
    filter(SI_high >= 10 & SI_medium >= 10) %>%
    summarise(n = n(),
              raw_count_mean_SI_high = mean(SI_high),
              raw_count_mean_SI_medium = mean(SI_medium)
    )
  ##########################################################
  
  dataframe_results <- rbind(
    cbind(data_noneFilter, data.frame(x = "Nofilter")),
    cbind(data_none0, data.frame(x = "Anyvalueis0")),
    cbind(data_noneMean0, data.frame(x = "Meanof0")),
    cbind(data_noneMeanLess10, data.frame(x = "MeanLess10"))
  )
  
  return(dataframe_results)
}
SMALL_summary_data_samples_raw_counts(dds_SLarge_SMid_df)


dds_SLarge_SMid <- dds_SLarge_SMid[rowMeans(counts(dds_SLarge_SMid)) >= 10,]

dds_SLarge_SMid <- DESeq(dds_SLarge_SMid)
####IMPORTANT
dds_SLarge_SMid$condition <- factor(dds_SLarge_SMid$condition, levels = c("SI_high", "SI_medium"))
resultsNames(dds_SLarge_SMid)[2]
res_SLarge_SMid <- results(dds_SLarge_SMid)

table_creation <- function(res_object){
  p0.05_log2FC0.5 <- data.frame(res_object) %>% 
    dplyr::filter(pvalue < 0.05 & (log2FoldChange < -0.5 | log2FoldChange > 0.5)) %>%
    dplyr::summarise(
      n = n()
    )
  
  q0.05_log2FC0.5 <- data.frame(res_object) %>% 
    dplyr::filter(padj < 0.05 & (log2FoldChange < -0.5 | log2FoldChange > 0.5)) %>%
    dplyr::summarise(
      n = n()
    )
  
  p0.05_log2FC2 <- data.frame(res_object) %>% 
    dplyr::filter(pvalue < 0.05 & (log2FoldChange < -2 | log2FoldChange > 2)) %>%
    dplyr::summarise(
      n = n()
    )
  
  q0.05_log2FC2 <- data.frame(res_object) %>% 
    dplyr::filter(padj < 0.05 & (log2FoldChange < -2 | log2FoldChange > 2)) %>%
    dplyr::summarise(
      n = n()
    )
  
  table <- data.frame(x = c("pvalue0.05_log2FC0.5", "qvalue0.05_log2FC0.5", 
                            "pvalue0.05_log2FC2", "qvalue0.05_log2FC2"),
                      y = c(as.numeric(p0.05_log2FC0.5$n), as.numeric(q0.05_log2FC0.5$n),
                            as.numeric(p0.05_log2FC2$n), as.numeric(q0.05_log2FC2$n))
  )
  return(table)
}
table_creation(res_SLarge_SMid)


# Principal component analysis
vsd_SLarge_SMid <- vst(dds_SLarge_SMid, blind=FALSE)
rld_SLarge_SMid <- rlog(dds_SLarge_SMid, blind=FALSE)

#vsn::meanSdPlot(assay(dds_SLarge_SMid)) # raw counts
#vsn::meanSdPlot(assay(vsd_SLarge_SMid)) # variance stabilizing transformation
#vsn::meanSdPlot(assay(rld_SLarge_SMid)) # regularized log transformation

## 
#plotPCA(rld_SLarge_SMid, intgroup=c("condition", "id")) # #Principal componen analysis with reg log transformation
#plotMA(res_SLarge_SMid, ylim= c(-2,2), alpha = 0.05)


# VOLCANOPLOT
genes_LI.df_SLarge_SMid <- bitr(rownames(res_SLarge_SMid), fromType = "ENSEMBL",
                                toType = c("SYMBOL"),
                                OrgDb = org.Mm.eg.db)

res_modified_SLarge_SMid <- res_SLarge_SMid %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::rename(ENSEMBL = gene_id) %>%
  dplyr::left_join(., genes_LI.df_SLarge_SMid, by = "ENSEMBL")



#####Volcano plot
#png("VP_SI_Lgr5_high_vs_SI_Lgr5_mid.png", width = 800, height = 714)
EnhancedVolcano(res_modified_SLarge_SMid,
                lab = res_modified_SLarge_SMid$SYMBOL,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 0.001,
                FCcutoff = 2,
                #xlim = c(-5.5, 5.5),
                #ylim = c(0, -log10(10e-12)),
                # pointSize = c(ifelse(res$log2FoldChange>2, 8, 1)),
                labSize = 6.0,
                shape = c(6, 6, 19, 16),
                title = "SI Lgr5 mid vs SI Lgr5 high",
                subtitle = "Differential expression analysis",
                caption = "FC cutoff, 2; p-value cutoff, 0.001",
                legendPosition = "right",
                legendLabSize = 14,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 0.9,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                hline = c(10e-8)) +
  ggplot2::coord_cartesian(xlim=c(-6, 6)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1))

dev.off()






# Principal component analysis
countdata_PCA <- raw_counts[,c(1,3:20)] %>%
  textshape::column_to_rownames("gene_id") %>%
  magrittr::set_colnames(NULL)

coldata_PCA <- data.frame(id = colnames(raw_counts[,3:20]), 
                                  condition = condition, 
                                  replicates = replicates)

dds_PCA <- DESeqDataSetFromMatrix(countData = countdata_PCA, 
                                          colData = coldata_PCA, 
                                          design = ~ condition)


vsd_PCA <- vst(dds_PCA, blind=FALSE)
rld_PCA <- rlog(dds_PCA, blind=FALSE)
plotPCA(rld_PCA, intgroup=c("condition")) # #Principal componen analysis with reg log transformation
plotMA(res_PCA, ylim= c(-2,2), alpha = 0.05)

