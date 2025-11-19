#' @title Single-Cell Phenotype Association Score Prediction
#' @description
#' Uses a trained scPAS model to make predictions on independent single-cell data.
#' Calculates phenotype association scores and identifies significant cells based
#' on false discovery rate thresholding.
#'
#' @param model Seurat object containing the trained scPAS model
#' @param test.data Matrix or Seurat object containing test single-cell RNA-seq data
#' @param assay Name of assay to use for prediction (default: 'RNA')
#' @param FDR.threshold FDR threshold for identifying phenotype-associated cells (default: 0.05)
#' @param imputation Whether to perform data imputation (default: FALSE)
#' @param imputation_method Method for data imputation (default: 'KNN')
#' @param independent Whether test data is independent from training (default: TRUE)
#'
#' @return Seurat object with prediction results added as metadata, or data frame
#' containing prediction results if input was a matrix
#'
#' @examples
#' \dontrun{
#' # Predict using trained model on new data
#' predictions <- scPAS.prediction(
#' model = trained_model,
#' test.data = new_sc_data,
#' FDR.threshold = 0.05
#' )
#' }
#'
#' @export
scPAS.prediction <- function(
    model,
    test.data,
    assay = 'RNA',
    FDR.threshold = 0.05,
    imputation = F,
    imputation_method = 'KNN',
    independent = T
) {
    model <- SeuratObject::Misc(model, slot = 'scPAS_para')

    if (inherits(test.data, 'Seurat')) {
        if (imputation) {
            test.data <- imputation2(
                test.data,
                assay = assay,
                method = imputation_method
            )
            assay <- SeuratObject::DefaultAssay(test.data)
        }
        test.exp <- SeuratObject::GetAssayData(
            object = test.data,
            assay = assay,
            slot = 'data'
        )
        Expression_cell <- test.exp
        rownames(Expression_cell) <- rownames(test.exp)
        colnames(Expression_cell) <- colnames(test.exp)
    } else {
        test.exp <- as.matrix(test.data)
        Expression_cell <- methods::as(
            SigBridgeRUtils::normalize.quantiles(test.exp),
            'dgCMatrix'
        )
        rownames(Expression_cell) <- rownames(test.exp)
        colnames(Expression_cell) <- colnames(test.exp)
    }

    Coefs <- model$Coefs
    common <- SeuratObject::intersect(names(Coefs), rownames(Expression_cell))

    if (sum(Coefs != 0) < 20) {
        stop(
            "There are too few valid features and the test data may not be suitable for the model!"
        )
    }

    Coefs <- Coefs[common]
    Expression_cell <- Expression_cell[common, ]
    scaled_exp <- Seurat:::FastSparseRowScale(
        Expression_cell,
        display_progress = F
    )
    colnames(scaled_exp) <- colnames(Expression_cell)
    rownames(scaled_exp) <- rownames(Expression_cell)
    #scaled_exp <- methods::as(scaled_exp, "sparseMatrix")

    scaled_exp[which(is.na(scaled_exp))] <- 0
    risk_score <- crossprod(scaled_exp, Coefs)

    set.seed(12345)

    randomPermutation <- vapply(
        X = seq_len(2000),
        FUN = function(x) {
            set.seed(1234 + x)
            sample(Coefs, length(Coefs), replace = F)
        },
        FUN.VALUE = numeric(length(Coefs))
    )
    randomPermutation <- methods::as(randomPermutation, "sparseMatrix")
    risk_score.background <- crossprod(scaled_exp, randomPermutation)

    if (independent) {
        mean.background <- SeuratObject::rowMeans(risk_score.background)
        sd.background <- apply(risk_score.background, 1, sd)
    } else {
        mean.background <- mean(as.matrix(risk_score.background))
        sd.background <- stats::sd(as.matrix(risk_score.background))
    }

    Z <- (risk_score[, 1] - mean.background) / sd.background

    p.value <- stats::pnorm(q = abs(Z), mean = 0, sd = 1, lower.tail = F)
    q.value <- stats::p.adjust(p = p.value, method = 'BH')

    risk_score_data.frame <- data.frame(
        sample = colnames(Expression_cell),
        scPAS_RS = risk_score[, 1],
        scPAS_NRS = Z,
        scPAS_Pvalue = p.value,
        scPAS_FDR = q.value
    )
    risk_score_data.frame$scPAS <- ifelse(
        Z > 0 & q.value <= FDR.threshold,
        'Positive',
        ifelse(Z < 0 & q.value <= FDR.threshold, 'Negative', 'Neutral')
    )

    if (!inherits(test.data, 'Seurat')) {
        return(risk_score_data.frame)
    }

    test.data <- SeuratObject::AddMetaData(
        test.data,
        metadata = risk_score_data.frame[-"sample"]
    )

    test.data
}
