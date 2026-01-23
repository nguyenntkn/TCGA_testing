library('TCGAbiolinks')
library('easybio')
library('data.table')
library('SummarizedExperiment')
library('DESeq2')
library('tidyverse')

# Following this vignette: 
# https://cran.r-project.org/web/packages/easybio/vignettes/example_limma.html

# Get query (metadata) for downloading
query <- GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)
# ============= exploring metadata =============
metadata <- query[[1]][[1]]
# Just checking what sample types there are
unique(metadata$sample_type)
# There's only Primary Tumor and Solid Tissue Normal

# Experimental stratery (RNA-seq, microarray?...)
unique(metadata$experimental_strategy)
# All RNA-Seq. Wonder if people still use microarray here...

# Analysis workflow?
unique(metadata$analysis_workflow_type)
# They're all STAR - Counts. Makes sense that they use STAR aligner for RNA-seq
# Seemed like they don't use HISAT2 aligner or Kallisto/Salmon pseudoaligner...

# Get number of samples
length(unique(metadata$cases))
# 44 samples total

# Get number of patients
length(unique(metadata$cases.submitter_id))
# 36 patients

# Since there is some miss match between patient number and samples, might be good
# to check more closely
sample_count <- metadata %>% 
  group_by(cases.submitter_id, sample_type) %>% 
  summarise(n_samples = n()) %>% 
  pivot_wider(
    names_from = sample_type,
    values_from = n_samples
  )
print(sample_count)
# Each patient have different number of tumor vs normal samples:
#   - Tumor only (patient only has 1 tumor data)
#   - Normal only (patient only has 1 normal data)
#   - Paired (patient has 1 tumor and 1 normal data)


# Might want to do two of the following approach:
# 1. Paired analysis only: do DEA and GSEA on paired samples only 
# 2. Supplement: treat all samples as unpaired

# Get a list of patients with paired data
# Clean up col names
colnames(sample_count) <- c("cases.submitter_id", "tumor", "normal") 
paired_patient <- sample_count %>% 
  filter(tumor == 1 & normal == 1) %>% 
  select(cases.submitter_id) %>% 
  as.list()

paired_patient <- paired_patient[[1]]

# ============= Download data and some more exploration ==================
# Download associated data. 
GDCdownload(query = query)

# Convert the downloaded data into an R object suitable for analysis.
# Since these are RNA-seq data, it will be a RangedSumarizeExperiment object
data <- GDCprepare(query = query)

# Some exploration of count data
exprs_data <- assay(data, "unstranded")
dim(exprs_data)
# 60660 genes and 44 samples


# ================= PCA =================

# First prep data for PCA

# Some notes before doing PCA: 
#   1. PCA needs expression data table to be transposed (row = samples, col = feature)
#   2. Gene expression data must be transformed/normalized. This is to prevent the 
#      highly expressed and/or highly variance genes to dominate the PCs, and correct
#      for sequencing depth. 
#         NOTE: This is not needed for DEA since some tools requires raw count input. 

# For visualization (PCA or heatmap), 3 normalization techniques are recommended
#   1. Transcript per million (TPM)
#   2. Regularized Log (rlog) or Variance Stabilizing Transformation (VST)
#   3. Log2(x+1)
# Apparently it's best to use VST. BUT WHY? READ MORE INTO THIS!!


# PCA:
pca_data <- prcomp(
  exprs_data %>% vst() %>% t()        
)

# PCA data is stored in "x" within the prcomp object. This contains:
# col: PCs
# row: samples
pca_data$x %>% 
  ggplot(aes(x=PC1, y=PC2, colour = metadata$sample_type)) + 
  geom_point()


# ================= DEA =================

# First do paired analysis
# Extract paired data
paired_metadata <- metadata %>% 
  filter(cases.submitter_id %in% paired_patient)

# Subseting expression data to contain only paired samples
paired_exprs_data <- exprs_data[, paired_metadata$cases]

# Exploratory on PCA:
pca_paired <- prcomp(
  paired_exprs_data %>% vst() %>% t()        
)

# PCA data is stored in "x" within the prcomp object. This contains:
# col: PCs
# row: samples
pca_paired$x %>% 
  ggplot(aes(x=PC1, y=PC2, colour = paired_metadata$sample_type)) + 
  geom_point()




