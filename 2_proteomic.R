# 0.1 loda libraries
# go to 0_setup.R to install and load libraries

# 0.2 preparing expression table
# prepare your data for input with following criteria: Unique peptide > 1
# data should contain: Accession and Abundance of each group

# 1.0 read data
# specify your file path
filepath <- "23122201-Protein Peptide Report(1)mass檔案.xlsx" # your file path
colnames(read_excel(filepath, sheet = "BD"))

# do the following for each group
df1 <- read_excel(filepath, sheet = "BD") %>% # change sheet name to your data
  filter(!is.na(Accession)) %>% # filter out Accession with NA
  tibble::column_to_rownames("Accession") %>% # designate Accession as rownames
  filter(`# Unique Peptides` > 1) %>% # filter out Accession with Unique Peptides < 1 # nolint: line_length_linter.
  select("Abundance: BD") # select the column with abundance data

df2 <- read_excel(filepath, sheet = "LG") %>%
  filter(!is.na(Accession)) %>%
  column_to_rownames("Accession") %>%
  filter(`# Unique Peptides` > 1) %>%
  select("Abundance: LG")

df3 <- read_excel(filepath, sheet = "HG") %>%
  filter(!is.na(Accession)) %>%
  column_to_rownames("Accession") %>%
  filter(`# Unique Peptides` > 1) %>%
  select("Abundance: HG")

# Perform the above for all groups

# prepare expression table
# merge all data into one table
df <- merge(df1, df2, by = "row.names", all = TRUE) %>%
  column_to_rownames("Row.names")
df <- merge(df, df3, by = "row.names", all = TRUE) %>%
  column_to_rownames("Row.names")

# designate your column names
# change the column names to your group names
col <- c("BD", "LG", "HG")
colnames(df) <- col

# 2.0 data preprocessing
# imputation and normalisation
# 2.1 imputation
df[is.na(df)] <- 0
df$BD <- as.numeric(df$BD) + 1
df$LG <- as.numeric(df$LG) + 1
df$HG <- as.numeric(df$HG) + 1

# 2.2 quantile noralisation
df_norm <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(df)))

# you can find the quantiles of your data before and after transformation
sapply(df, function(x) quantile(x, probs = seq(0, 1, 1 / 4)))
sapply(df_norm, function(x) quantile(x, probs = seq(0, 1, 1 / 4)))

# specify your column and row names, as your index and column names will lost during the process
names(df_norm) <- col
rownames(df_norm) <- rownames(df)

# 2.3 calculate log2 Fold Change
df_log2FC <- df_norm %>%
  mutate(
    LG_log2fc = log2(LG / BD), # change IL, ML to your group names
    HG_log2fc = log2(HG / BD),
    HGvsLG_log2fc = log2(HG / LG),
    padj = 0 # add padj column for MA plot
)

# save your expression table as .xlsx
# you can find DEG inside e.g. |log2FC| > 1
write_xlsx(df_log2FC %>%
  rownames_to_column("proteinID") %>%
  map() %>%
  annotate(),
  "data_normalised.xlsx"
)

# 2.4 map your rownames(now it's Uniprot accession id) to NCBI ENTREZ id
# ENTREZ id is required for KEGG pathway analysis
df_log2FC$ENTREZID <- rownames(df_log2FC) %>%
  mapIds(org, keys = ., column = "ENTREZID", keytype = "UNIPROT") %>%  # change keytype to your accession type
  ifelse(is.na(.) | duplicated(.), names(.), .)

# 3.0 start analyse
# 3.1 heatmap
pheat <- pheatmap::pheatmap(df_norm,
                            scale = "row",
                            show_rownames = FALSE
)
png("heatmap.png", res = 600, width = 6000, height = 8000)
pheat
dev.off()

# 3.2 MA plot
# prepare input
# specify your group names
LG <- df_log2FC %>% 
  rownames_to_column("name") %>% # designate rownames as column
  dplyr::select("name", "LG", "LG_log2fc", "padj") # select columns
colnames(LG) <- c("name", "baseMean", "log2FoldChange", "padj") # rename columns for ggmaplot

# do the same for other groups
HG <- df_log2FC %>%
  rownames_to_column("name") %>%
  dplyr::select("name", "HG", "HG_log2fc", "padj")
colnames(HG) <- c("name", "baseMean", "log2FoldChange", "padj")

HGvsLG <-df_log2FC %>%
  rownames_to_column("name") %>%
  dplyr::select("name", "HG", "HGvsLG_log2fc", "padj")
colnames(HGvsLG) <- c("name", "baseMean", "log2FoldChange", "padj")

# plot and save the MA plot
plotma(LG, "LGvsBD")
plotma(HG, "HGvsBD")
plotma(HGvsLG, "HGvsLG")

# 3.3 enrichment analysis
# construct geneList as input for ClusterProfiler 4.0
# set log2FC as values and index of df_norm as names
LG_ranked <- df_log2FC %>%
  arrange(desc(LG_log2fc)) # arrange by log2FC
LG_genelist <- LG_ranked$LG_log2fc # select the column
names(LG_genelist) <- LG_ranked$ENTREZID # designate the index as names
LG_genelist <- LG_genelist[LG_genelist != 0] # filter out 0 values

# do the same for other groups
HG_ranked <- df_log2FC %>%
  arrange(desc(HG_log2fc))
HG_genelist <- HG_ranked$HG_log2fc
names(HG_genelist) <- HG_ranked$ENTREZID
HG_genelist <- HG_genelist[HG_genelist != 0]

HGvsLG_ranked <- df_log2FC %>%
  arrange(desc(HGvsLG_log2fc))
HGvsLG_genelist <- HGvsLG_ranked$HGvsLG_log2fc
names(HGvsLG_genelist) <- HGvsLG_ranked$ENTREZID
HGvsLG_genelist <- HGvsLG_genelist[HGvsLG_genelist != 0]


# perform enrichment analyses
# every thing is wrapped as function so it's easier to read and execute
performenrichment(LG_genelist, "LGvsBD")
performenrichment(HG_genelist, "HGvsBD")
performenrichment(HGvsLG_genelist, "HGvsLG")

view(LG_genelist)
view(df_log2FC)