#' @title Compute Correlation for Sparse Matrix
#' @description
#' Efficiently computes the correlation matrix for a sparse matrix using
#' memory-optimized calculations that avoid dense matrix conversions.
#'
#' @param x A sparse matrix (dgCMatrix) where rows represent features and
#' columns represent samples
#'
#' @return A correlation matrix with the same dimensions as the input matrix
#'
#' @details
#' This function is particularly useful for single-cell RNA-seq data where
#' expression matrices are typically sparse. It handles the sparsity efficiently
#' by only processing non-zero elements and using mathematical optimizations
#' to compute correlations without converting to dense format.
#'
#' @examples
#' \dontrun{
#' # Create a sparse matrix
#' library(Matrix)
#' x <- Matrix::rsparsematrix(1000, 500, density = 0.1)
#'
#' # Compute correlation matrix
#' cor_matrix <- sparse.cor(x)
#' dim(cor_matrix)
#' }
#'
#' @export
sparse.cor <- function(x) {
    n <- nrow(x)
    m <- ncol(x)
    ii <- unique(x@i) + 1 # rows with a non-zero element

    Ex <- Matrix::colMeans(x)
    nozero <- as.vector(x[ii, ]) - rep(Ex, each = length(ii)) # colmeans

    covmat <- (crossprod(matrix(nozero, ncol = m)) +
        crossprod(t(Ex)) * (n - length(ii))) /
        (n - 1)
    sdvec <- sqrt(Matrix::diag(covmat))
    covmat / crossprod(t(sdvec))
}
