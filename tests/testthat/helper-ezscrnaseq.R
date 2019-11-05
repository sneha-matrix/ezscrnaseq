library(BiocParallel)
library(knitr)
library(pheatmap)
library(rprojroot)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scater)
library(scran)
library(splatter)
library(testthat)

options(stringsAsFactors = FALSE)

suppressWarnings(set.seed(seed = 100, sample.kind = "Rounding"))
data("sc_example_counts")
sc_example_counts_shift <- sc_example_counts
sc_example_counts_shift[] <- sample(c(rep(0, 10), 1:10), size=nrow(sc_example_counts)*ncol(sc_example_counts), replace=TRUE)
# Estimate parameters from example data
params <- splatEstimate(cbind(sc_example_counts, sc_example_counts_shift+sc_example_counts))
# Simulate data using estimated parameters
sce <- splatSimulate(params, group.prob=c(0.5, 0.5), method="groups", verbose=FALSE)
rownames(sce) <- readLines(find_testthat_root_file("testing_data", "ensembl_ids.txt", path=".."))
rowData(sce)$Gene[1:10] <- paste0("MT-", rowData(sce)$Gene[1:10])
