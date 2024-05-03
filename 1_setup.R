# 0.0 install packages (If you don't have it)
# install.packages("BiocManager")
# install.packages("tidyverse")
# install.packages("readxl")
# install.packages("writexl")
# install.packages("pheatmap")
# install.packages("ggpubr")
# BiocManager::install("preprocessCore")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("DOSE")
# BiocManager::install("enrichplot")
# BiocManager::install("simplifyEnrichment")
# BiocManager::install("ggridges")
# BiocManager::install("biomaRt")

# 0.1 loda libraries
library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(biomaRt)
library(pheatmap)
library(ggpubr)
library(preprocessCore)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(simplifyEnrichment)
library(ggridges)

# 0.2 The species name of your data
org <- org.Mm.eg.db # change to your organism
KEGGorg <- "mmu" # search kegg organism for your species

# 0.3 prepare biomaRt
# access the ensembl database
ensembl <- biomaRt::useEnsembl(
  biomart = "genes", # retrieve gene data
  dataset = "mmusculus_gene_ensembl", # designate spiecies !!IMPORTANT!!
  mirror = "asia"
) # version ensures consistency
# 111 is updated in 2024
attributes <- listAttributes(ensembl)

# 0.4 mapID function
map <- function(input) {
  ensembl <- input$proteinID %>%
    mapIds(org.Mm.eg.db, keys = ., column = "ENSEMBL", keytype = "UNIPROT") %>% # change org.Mm.eg.db to your organism
    ifelse(is.na(.) | duplicated(.), names(.), .)
  input$Ensembl_id <- ensembl
  return(input)
}

# 0.5 prepare annotation function
annotate <- function(input) {
  input <- input %>% dplyr::left_join(
    getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "description"),
      filters = "ensembl_gene_id",
      values = .$Ensembl_id,
      mart = ensembl
    ),
    by = c("Ensembl_id" = "ensembl_gene_id")
  )
}

# function for MAplot
plotma <- function(input, fname) {
  label <- input %>%
    mutate(md = abs(baseMean) + abs(log2FoldChange)) %>%
    filter(abs(log2FoldChange) > 1) %>%
    arrange(desc(md)) %>%
    slice(1:15)
  label <- label$name
  ma <- ggmaplot(input,
    main = paste0(fname, " MA-plot"),
    fc = 2, size = 1,
    palette = c("#B31B21", "#1465AC", "darkgray"),
    genenames = as.vector(input$name),
    legend = "bottom", top = 0,
    font.label = c("bold", 12), label.rectangle = TRUE,
    font.legend = "bold",
    font.main = "bold",
    label.select = label,
    ggtheme = ggplot2::theme_minimal()
  )
  png(paste0(fname, "_MAplot.png"), height = 5000, width = 5000, res = 600)
  print(ma)
  dev.off()
}

# 0.7 GSEA analysis
# This include perform and visualise GSEA analysis
performenrichment <- function(genelist, fname) {
  # GO biological process
  bp <- gseGO(
    geneList = genelist,
    ont = "BP",
    OrgDb = org,
    keyType = "ENTREZID",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    by = "fgsea",
    nPermSimple = 10000
  )
  bp2 <- DOSE::setReadable(bp, OrgDb = org, keyType = "ENTREZID")

  # GO cecullar component
  cc <- gseGO(
    geneList = genelist,
    ont = "CC",
    OrgDb = org,
    keyType = "ENTREZID",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    by = "fgsea",
    nPermSimple = 100000
  )
  cc2 <- DOSE::setReadable(cc, OrgDb = org, keyType = "ENTREZID")

  # GO molecular function
  mf <- gseGO(
    geneList = genelist,
    ont = "MF",
    OrgDb = org,
    keyType = "ENTREZID",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    by = "fgsea",
    nPermSimple = 100000
  )
  mf2 <- DOSE::setReadable(mf, OrgDb = org, keyType = "ENTREZID")

  # KEGG pathway
  kk <- gseKEGG(
    geneList = genelist,
    organism = "mmu", # search KEGG organism name
    keyType = "kegg",
    minGSSize = 1,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = TRUE,
    by = "fgsea",
    nPermSimple = 100000
  )
  kk2 <- DOSE::setReadable(kk, OrgDb = org, keyType = "ENTREZID")

  # visualise the results as dotplot
  res <- c(bp, cc, mf, kk)
  ridgeplot_filenames <- c("bpridge.png", "ccridge.png", "mfridge.png", "keggridge.png")
  dot_filenames <- c("bpdot.png", "ccdot.png", "mfdot.png", "keggdot.png")
  for (i in seq_along(res)) {
    png(paste0(fname, "_", ridgeplot_filenames[i]),
      height = 5000, width = 5000, res = 600
    )
    print(enrichplot::ridgeplot(res[[i]], showCategory = 10))
    dev.off()

    png(paste0(fname, "_", dot_filenames[i]),
      height = 5000, width = 5000, res = 600
    )
    print(enrichplot::dotplot(res[[i]], showCategory = 5, split = ".sign") +
      facet_grid(. ~ .sign))
    dev.off()
  }

  # save the raw results in a .xlsx file
  write_xlsx(
    list(
      gsea_bp = as.data.frame(bp2),
      gsea_cc = as.data.frame(cc2),
      gsea_mf = as.data.frame(mf2),
      gsea_kegg = as.data.frame(kk2)
    ),
    path = paste0(fname, "_GSEA_Results.xlsx")
  )
}
