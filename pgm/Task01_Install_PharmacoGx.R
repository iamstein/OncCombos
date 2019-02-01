#http://bioconductor.org/packages/release/bioc/html/PharmacoGx.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("PharmacoGx", version = "3.8")
