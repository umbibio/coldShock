library(Matrix)

# this scripts reads the filtered_feature_bc_matrix generated by Cellranger Count
# and generates the feature_count matrix and save it in RData folder.

input.dir <- "../Input/coldShock/filtered_feature_bc_matrix/"
matrix.dirs <- list.files(input.dir)
num.total.files <- length(matrix.dirs)

out.dir <- "../Input/coldShock/counts/"


for(i in 1:num.total.files){
  d <- matrix.dirs[i]
  cat(paste('processing file', d))
  cat('\n')
  matrix_dir <- paste(input.dir, d, sep = "")
  barcode.path <- paste(matrix_dir, "/barcodes.tsv.gz", sep = "")
  features.path <- paste(matrix_dir, "/features.tsv.gz", sep = "")
  matrix.path <- paste(matrix_dir, "/matrix.mtx.gz", sep = "")
  
  mat <- readMM(file = matrix.path)
  
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  expr <- as.data.frame(as.matrix(mat))
  
  expr.csv.name <- paste(d, ".expr.csv", sep = "")
  write.csv(expr,paste(out.dir, expr.csv.name, sep = ""))
}


