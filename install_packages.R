if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
remove.packages("Matrix")
install.packages("BiocManager")
BiocManager::install("Matrix") 
install.packages("C:/Users/Daniel Chang/Downloads/Matrix_1.6-4.tar.gz", repos = NULL, type="source")
install.packages("Matrix")
install.packages("Rcpp")
install.packages("glue")
devtools::install_github("saeyslab/nichenetr")
install.packages("htmltools")
install.packages('Seurat')
install.packages("nichenetr")