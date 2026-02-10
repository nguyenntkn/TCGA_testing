############################################################
## TCGA-CHOL RNA-seq Analysis
## Paired Tumor vs Normal using limma-voom
## Reference: https://cran.r-project.org/web/packages/easybio/vignettes/example_limma.html
############################################################




# ============ 0. Set up ============

wd = "/Users/tknguyenhollandcollege.com/Documents/R/TestingTCGA"

# install.packages(c('data.table', 'tidyverse'))
# BiocManager::install(c('TCGAbiolinks', 'easybio', 'SummarizedExperiment',
#                        'DESeq2', 'edgeR', 'EnhancedVolcano', 'clusterProfiler',
#                        'org.Hs.eg.db', 'enrichplot'))


library('TCGAbiolinks')
library('easybio')
library('data.table')
library('SummarizedExperiment')
library('DESeq2')
library('tidyverse')
library('edgeR')
library('EnhancedVolcano')
library('org.Hs.eg.db')
library('clusterProfiler')
library('enrichplot')



# ============ 1. Query TCGA metadata and download count data ============
query <- GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)

# Download associated data. Then convert into a RangedSummarizeExperiment object.
GDCdownload(query = query)
data <- GDCprepare(query = query)





# ============ 2. Exploring metadata to identify paired patients ============
# Extract sample-level metadata
metadata <- query[[1]][[1]]

# Check sample composition
unique(metadata$sample_type)
unique(metadata$experimental_strategy)
unique(metadata$analysis_workflow_type)
# There's only Primary Tumor and Solid Tissue Normal
# All RNA-Seq. All used STAR - Counts work flow.

# Total samples and patients
n_samples  <- length(unique(metadata$cases))
n_patients <- length(unique(metadata$cases.submitter_id))

# Get number of samples
length(unique(metadata$cases))
# 44 samples total from 36 patients

# Count tumor and normal samples per patient
sample_count <- metadata %>% 
  group_by(cases.submitter_id, sample_type) %>% 
  summarise(n_samples = n()) %>% 
  pivot_wider(
    names_from = sample_type,
    values_from = n_samples
  )

# Rename columns for clarity
colnames(sample_count) <- c("case.submitter_id", "tumor", "normal")

# Keep patients with at least one tumor and one normal sample
paired_patients <- sample_count %>%
  filter(tumor >= 1 & normal >= 1) %>%
  pull(case.submitter_id)





# ============= 3. Subset paired gene expression data  ==================
# Extract raw count matrix (genes x samples)
exprs_data <- assay(data, "unstranded")
dim(exprs_data) # 60660 genes and 44 samples

# Subset metadata to paired patients
paired_metadata <- metadata %>% 
  filter(cases.submitter_id %in% paired_patients) 

# Clean sample_type labels
paired_metadata$sample_type <- gsub(" ", "", paired_metadata$sample_type)

# Subset count matrix to paired samples
paired_exprs_data <- exprs_data[, paired_metadata$cases]





# ============ 4. Exploring expression data ============

# ----------- 4.1. PCA with VST transformation ---------

# Some notes before doing PCA: 
#   1. PCA needs expression data table to be transposed (row = samples, col = feature)
#   2. Gene expression data must be transformed/normalized. This is to prevent the 
#      highly expressed and/or highly variance genes to dominate the PCs, and correct
#      for sequencing depth. 
#         NOTE: This is not needed for DEA since some tools requires raw count input. 

# For visualization (PCA or heatmap), 3 normalization techniques are recommended:
#   1. Transcript per million (TPM)
#   2. Regularized Log (rlog) or Variance Stabilizing Transformation (VST)
#   3. Log2(x+1)
# Apparently it's best to use VST. BUT WHY? READ MORE INTO THIS!!

# Perform PCA
pca_data <- prcomp(paired_exprs_data %>% vst() %>% t())

# PCA data is stored in "x" within the prcomp object (col-PCs, row-samples)
pca_data$x %>% 
  ggplot(aes(x=PC1, y=PC2, colour = paired_metadata$sample_type)) + 
  geom_point(size = 4) +
  scale_color_discrete(name="Sample type") +
  theme_classic()

# -------------- 4.2. MDS plot ----------------
# WHY?
# But need to make sure the column names in count table has the same order as
# the samples in meta data table.
colnames(paired_exprs_data) == paired_metadata$cases
plotMDS(dge, col = as.numeric(as.factor(paired_metadata$sample_type)), pch = 16)





# ==================== 5. Filter lowly-expressed genes ===================

# Step 1: Store raw counts in a DGEList object
#    - Contains: raw counts, sample info (group, patient IDs), library sizes
#    - Needed for limma-voom normalization
dge <- DGEList(counts = paired_exprs_data)

# Step 2: Calculate normalization factors
#    - This step does NOT normalize the counts
#    - Provides info for voom to account for library size differences
dge <- calcNormFactors(dge)

# Step 3: Filter out lowly expressed genes
#    - Low-expression genes are usually uninformative and increase computational load
#    - Use counts-per-million (CPM) rather than raw counts to account for library size
#    - Here we keep genes with max CPM >= 1 (roughly equivalent to ~10 raw counts per 10M reads)
cpm_threshold <- 1
low_expr_genes <- which(apply(cpm(dge), 1, max) < cpm_threshold)
dge <- dge[-low_expr_genes, ]

# Step 4: Check remaining dimensions
dim(dge$counts)



# ============== 6. limma-voom differential analysis =============

# Define experimental group
group <- factor(paired_metadata$sample_type)
patient <- factor(gsub("-", "_", paired_metadata$cases.submitter_id))

# Notes for Voom transformation:
# 1. Counts are transformed to log2 counts per million reads (CPM), where “per million reads” 
# is defined based on the normalization factors we calculated earlier.
# 2. A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated.
# 3. A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression 
# (see red line in plot)
# 4. The smoothed curve is used to obtain weights for each gene and sample that are passed into 
# limma along with the log2 CPMs.


# ----------- 6.1. Unpaired analysis --------------
# This is just to see the difference between unpaired and paired analysis

# Design model matrix (no intercept) 
mm_unpaired <- model.matrix(~ 0 + group)

# Voom transformation
v_unpaired <- voom(dge, mm_unpaired, plot=T)

# Fit linear model
fit_unpaired <- lmFit(v_unpaired, mm_unpaired)

# Define contrast: Tumor vs Normal
contrast_unpaired <- makeContrasts(
  groupPrimaryTumor - groupSolidTissueNormal,
  levels = colnames(coef(fit_unpaired))
)

# Apply the specified contrast
tmp_unpaired <- contrasts.fit(fit_unpaired, contrast_unpaired)

# Performs empirical Bayes moderation
tmp_unpaired <- eBayes(tmp_unpaired)

# Extract results
top.table_unpaired <- topTable(tmp_unpaired, sort.by = "P", n = Inf)

# Export DGE results
write.csv(top.table_unpaired, file = file.path(wd, "Data", "DE_unpaired.csv"), row.names = TRUE)


# ----------- 6.2. Paired analysis --------------

# Design model matrix (no intercept) 
mm_paired <- model.matrix(~ 0 + group + patient)

# Voom transformation. 
v_paired <- voom(dge, mm_paired, plot=T)

# Fit linear model
fit_paired <- lmFit(v_paired, mm_paired)

# Define contrast: Tumor vs Normal
contrast_paired <- makeContrasts(
  groupPrimaryTumor - groupSolidTissueNormal,
  levels = colnames(coef(fit_paired))
)

# Apply the specified contrast
tmp_paired <- contrasts.fit(fit_paired, contrast_paired)

# Performs empirical Bayes moderation
tmp_paired <- eBayes(tmp_paired)

# Extract results
top.table_paired <- topTable(tmp_paired, coef = 1, sort.by = "P", n = Inf)

# Export DGE results
write.csv(top.table_paired, file = file.path(wd, "Data", "DE_paired.csv"), row.names = TRUE)


# ----------- 6.3. Volcano plots --------------
EnhancedVolcano(top.table_paired,
                lab = rownames(top.table_paired),
                x = 'logFC',
                y = 'P.Value',
                pointSize = 0.5)





# ============== 7. GO enrichment analysis ==================
paired_gene_list <- top.table_paired$logFC

# Clean ensembl ID to remove everything after the "."
names(paired_gene_list) <- gsub("\\..*", "", row.names(top.table_paired))

paired_gene_list <- na.omit(paired_gene_list)

sorted_paired_gene_list = sort(paired_gene_list, decreasing = TRUE)

head(sorted_paired_gene_list)
is.numeric(sorted_paired_gene_list)
!is.null(names(sorted_paired_gene_list))

set.seed(1234)

# gse <- gseGO(geneList=gene_list, 
#              ont ="ALL",              # Ontology: biological process (BP), cellular components (CC), 
#              # molecular function (MF) or all 3 (ALL)
#              keyType = "ENTREZID",    # Gene list is named by Entrez IDs
#              pvalueCutoff = 0.05,     
#              verbose = TRUE,          # Print progress in console
#              OrgDb = hs,              # Organism database: human
#              seed = TRUE,             # Makes the permutation process reproducible
#              pAdjustMethod = "none")  # "none" for raw p-value or 
# # "BH" (Benjamin-Hochberg) for FDR correction.

gse <- gseGO(geneList = sorted_paired_gene_list, 
             ont ="ALL",                # Ontology: biological process (BP), cellular components (CC), 
                                        # molecular function (MF) or all 3 (ALL)
             keyType = "ENSEMBL",       # Gene list is named by Ensembl IDs
             pvalueCutoff = 0.05,       # Only keep GO terms with p < 0.05
             verbose = TRUE,            # Print progress in console
             OrgDb = org.Hs.eg.db,      # Organism database: human
             seed = TRUE,               # Makes the permutation process reproducible
             pAdjustMethod = "BH")

dotplot(gse,
        x = "GeneRatio",
        color = "p.adjust",
        title = "Top 15 of GO Enrichment",
        showCategory = 10,
        label_format = 80,
        split=".sign") +
  facet_grid(.~.sign)




# ========== Supplement =============

# Check differences before and after filtering lowly expressed genes
dge_0 <- DGEList(counts = paired_exprs_data)
dge_0 <- calcNormFactors(dge_0)

cpm_threshold <- 1
low_expr_genes <- which(apply(cpm(dge_0), 1, max) < cpm_threshold)
dge <- dge_0[-low_expr_genes, ]

v_0 <- voom(dge_0, mm_paired, plot=T)    # Pre-filter
v <- voom(dge, mm_paired, plot=T)        # Post-filter
