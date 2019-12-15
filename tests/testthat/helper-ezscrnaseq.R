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
#data("sc_example_counts")
sc_example_counts <- readRDS(find_testthat_root_file("testing_data","sc_example_counts.rda", path=".."))
sc_example_counts_shift <- sc_example_counts
sc_example_counts_shift[] <- sample(c(rep(0, 10), 1:10), size=nrow(sc_example_counts)*ncol(sc_example_counts), replace=TRUE)
# Estimate parameters from example data
params <- splatEstimate(cbind(sc_example_counts, sc_example_counts_shift+sc_example_counts))
# Simulate data using estimated parameters
sce <- splatSimulate(params, group.prob=c(0.5, 0.5), method="groups", verbose=FALSE)
rownames(sce) <- readLines(find_testthat_root_file("testing_data", "ensembl_ids.txt", path=".."))
rowData(sce)$Gene[1:10] <- paste0("MT-", rowData(sce)$Gene[1:10])
#Define phase and pairs in data
genes.list <- head(rownames(sce),200)
G1 <- c(3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 17, 18, 19, 21, 23, 24, 26, 27, 28, 30, 31, 32, 34, 35, 36, 37, 39, 40, 42, 
43, 44, 45, 46, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61, 63, 64, 65, 66, 69, 73, 74, 75, 77, 79, 80)
S <- c(1, 4, 9, 14, 16, 20, 25, 29, 33, 38, 41, 47, 59, 62, 67, 68, 70, 71, 72, 76, 78)
G2M <- c(2, 22)
