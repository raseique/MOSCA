# Normalization of gene expression data
# By João Sequeira
# Jan 2019

normalize_gene_expression <- function(df, output_file, norm_method, imput_method="LLS") {
  # TMM or RLE normalization -> for RNA-Seq
  if (norm_method == "TMM" || norm_method == "RLE") {
    print(paste("Performing", if (norm_method == "TMM") {"Trimmed Mean of M-values"} else {"Relative Log Expression"},
                "normalization.", sep=' '))
    library("edgeR")
    df[is.na(df)] <- 0
    factors <- calcNormFactors(df, method=norm_method)
    write.table(factors, file=paste0(dirname(output_file), "/norm_factors.txt"), sep='\n', row.names=FALSE, col.names=FALSE)
    df[, 1:ncol(df)] <- mapply("*", df[, 1:ncol(df)], factors)
  } else if (norm_method == "VSN") {
    df <- as.matrix(df)
    print("Performing Variance Stabilizing Normalization.")
    library("vsn")
    library("pcaMethods")
    df[df == 0] <- NA
    df <- df[rowSums(is.na(df)) < ncol(df), ]
    norm <- exprs(justvsn(ExpressionSet(df)))
    if (imput_method == "LLS") {
      print("Performing Local Least Squares Imputation.")
      allVariables <- FALSE
      if (sum(complete.cases(norm)) / nrow(norm) < 0.5) {
        allVariables <- TRUE
      }
      imputed <- llsImpute(t(norm), correlation="pearson", allVariables=allVariables)
      df <- t(completeObs(imputed))
    } else if (imput_method == "MIN") {
      print("Imputing missing values with the minimum value.")
      df[is.na(df)] <- min(norm, na.rm=TRUE)
    }
  } else {
    stop("Error: normalization method must be either TMM, RLE, or VSN")
  }
  print('Writing normalized results.')    # but first convert back to counts (normalized)
  write.table(2^df, file=output_file, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
}

print("Reading data to normalize.")
df1 <- read.table(snakemake@input[[1]], header=TRUE, sep="\t", row.names=1)
normalize_gene_expression(df1, snakemake@output[[1]], snakemake@params$norm_method)
if (length(snakemake@input) > 1) {
  print("Reading more data to normalize.")
  df2 <- read.table(snakemake@input[[2]], header=TRUE, sep="\t", row.names=1)
  # check if mt is in the name (snakemake@input[[2]]) and mp is not
  if (grepl("mt", snakemake@input[[2]]) && !grepl("mp", snakemake@input[[2]])) {
      normalize_gene_expression(df2, snakemake@output[[2]], snakemake@params$norm_method)
  } else {      # we're dealing with mp
      normalize_gene_expression(
        df2, snakemake@output[[2]], "VSN", imput_method=snakemake@params$imput_method)
  }
}
