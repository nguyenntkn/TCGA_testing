############################################################
## TCGA-CHOL RNA-seq Analysis
## Paired Tumor vs Normal using limma-voom
## Reference: https://cran.r-project.org/web/packages/easybio/vignettes/example_limma.html
############################################################




# ============ 0. Set up ============

wd = "/Users/nguyennguyen/Documents/TCGA_testing"

# install.packages(c('data.table', 'tidyverse'))
# BiocManager::install(c('TCGAbiolinks', 'easybio', 'SummarizedExperiment',
#                        'DESeq2', 'edgeR', 'EnhancedVolcano'))

library('TCGAbiolinks')
library('easybio')
library('data.table')
library('SummarizedExperiment')
library('DESeq2')
library('tidyverse')
library('edgeR')
library('EnhancedVolcano')





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
  geom_point()

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

tmp_unpaired <- contrasts.fit(fit_unpaired, contrast_unpaired)
tmp_unpaired <- eBayes(tmp_unpaired)


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

tmp_paired <- contrasts.fit(fit_paired, contrast_paired)
tmp_paired <- eBayes(tmp_paired)


top.table_paired <- topTable(tmp_paired, coef = 1, sort.by = "P", n = Inf)

# Export DGE results
write.csv(top.table_paired, file = file.path(wd, "Data", "DE_paired.csv"), row.names = TRUE)






# ============== Volcano plot ==================
EnhancedVolcano(top.table_paired,
                lab = rownames(top.table_paired),
                x = 'logFC',
                y = 'P.Value',
                pointSize = 0.5)


















# ========== Supplement =============

# Check differences before and after filtering lowly expressed genes
dge_0 <- DGEList(counts = paired_exprs_data)
dge_0 <- calcNormFactors(dge_0)

cpm_threshold <- 1
low_expr_genes <- which(apply(cpm(dge_0), 1, max) < cpm_threshold)
dge <- dge_0[-low_expr_genes, ]

v_0 <- voom(dge_0, mm_paired, plot=T)    # Pre-filter
v <- voom(dge, mm_paired, plot=T)        # Post-filter

