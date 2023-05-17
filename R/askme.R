#' Master data frame for all model objects
#'
#' This is an internal object which will be updated on every new release
#' 
#' @name askme_models
#' @docType data
#' @return master sheet of askme model objects
#' @export
NULL

#' Master data frame for all prediction labels
#'
#' This is an internal object which will be updated on every new release
#' 
#' @name prediction_labels
#' @docType data
#' @return master sheet of prediction labels
#' @export
NULL

#' Model for CNS cancer type classification
#'
#' @name m_cancertype_CNS66
#' @docType data
#' @return random forest model
#' @export
NULL

#' Model for TCGA cancer type classification
#'
#' @name m_cancertype_TCGA33
#' @docType data
#' @return random forest model
#' @export
NULL

#' classify cancer type
#'
#' @param betas DNA methylation betas
#' @param model cancer classification model
#' @return predicted cancer type label
#' @examples
#' library(sesameData)
#' betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
#' askme_cancertype(betas, model=m_cancertype_TCGA33)
#' ## expect PAAD
#' @import randomForest
#' @export
askme_cancertype <- function(betas, model=m_cancertype_TCGA33) {
    require("randomForest")
    betas <- betas[match(rownames(model$importance), names(betas))]
    predict(model, newdata=betas)
}
