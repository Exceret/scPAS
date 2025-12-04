#' @title The function of imputaion.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#' @param method The method for imputation. The default is 'RNA'.
#' @param verbose Logical, whether to print messages
#'
#' @return  A seurat object after imputaion.
#' @export
#' @family scPAS
#' @family imputation
#'
imputation2 <- function(
    obj,
    assay = 'RNA',
    method = c('KNN', 'ALRA'),
    verbose = SigBridgeRUtils::getFuncOption('verbose')
) {
    switch(
        method,
        'KNN' = {
            if (verbose) {
                ts_cli$cli_alert_info(
                    "Imputation of missing values in single cell RNA-sequencing data with {.val KNN}"
                )
            }
            imputation_KNN2(obj = obj, assay = assay, LogNormalized = TRUE)
        },
        'ALRA' = {
            if (verbose) {
                ts_cli$cli_alert_info(
                    "Imputation of missing values in single cell RNA-sequencing data with {.val ALRA}"
                )
            }
            imputation_ALRA2(obj = obj, assay = assay)
        },
        {
            cli::cli_warn(
                'The {.val {method}} method does not exist, so imputaion is invalid!',
                "Skip imputation"
            )
            obj
        }
    )
}

#' @title  A method for imputation of missing values in single cell RNA-sequencing data based on ALRA.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#'
#' @return  A seurat object after imputaion.
#'
#' @export
#' @family scPAS
#' @family imputation
#'
imputation_ALRA2 <- function(obj, assay = 'RNA') {
    # library(ALRA)
    # library(Matrix)
    # library(Seurat)
    rlang::check_installed("ALRA")
    # data <- GetAssayData(object = obj, assay = assay, slot = 'data')
    data <- SeuratObject::LayerData(obj, assay = assay)
    alra <- getExportedValue("ALRA", "alra")

    data_alra <- Matrix::t(alra(Matrix::t(as.matrix(data)))[[3]])
    # data_alra <- as.matrix(data) %>%
    #     Matrix::t() %>%
    #     alra() %>%
    #     .[[3]] %>%
    #     Matrix::t()
    colnames(data_alra) <- colnames(data)
    data_alra <- Matrix::Matrix(data_alra, sparse = T)

    obj[["imputation"]] <- SeuratObject::CreateAssayObject(data = data_alra)
    SeuratObject::DefaultAssay(obj) <- "imputation"
    obj
}

#' @title A method for imputation of missing values in single cell RNA-sequencing data based on the average expression value of nearest neighbor cells.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#' @param LogNormalized Whether the data is LogNormalized.
#'
#' @return A seurat object after imputaion.
#'
#' @export
#' @family scPAS
#' @family imputation
#'
imputation_KNN2 <- function(obj, assay = 'RNA', LogNormalized = T) {
    # library(Matrix)
    # exp_sc <- Seurat::GetAssayData(object = obj, assay = assay, slot = 'data')
    exp_sc <- SeuratObject::LayerData(obj, assay = assay)
    # nn_network <- obj@graphs[[paste0(assay, "_nn")]]
    # nn_network <- obj@graphs$RNA_nn
    nn_network <- SeuratObject::Graphs(obj, slot = paste0(assay, "_nn"))

    if (!methods::is(object = exp_sc, class2 = "sparseMatrix")) {
        exp_sc <- methods::as(exp_sc, "sparseMatrix")
    }
    if (!methods::is(object = nn_network, class2 = "sparseMatrix")) {
        nn_network <- methods::as(nn_network, "sparseMatrix")
    }
    if (LogNormalized) {
        exp_sc <- methods::as(exp(exp_sc) - 1, "sparseMatrix")
    }

    network_count <- methods::as(
        Matrix::Diagonal(x = 1 / Matrix::rowSums(nn_network)),
        "sparseMatrix"
    )
    exp_sc_mean <- tcrossprod(
        x = tcrossprod(x = exp_sc, y = nn_network),
        y = network_count
    )
    if (LogNormalized) {
        exp_sc_mean <- log1p(exp_sc_mean)
    }
    colnames(exp_sc_mean) <- colnames(exp_sc)
    obj[["imputation"]] <- SeuratObject::CreateAssayObject(data = exp_sc_mean)
    SeuratObject::DefaultAssay(obj) <- "imputation"

    obj
}
