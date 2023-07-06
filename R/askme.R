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
    require("dtplyr")
    betas <- betas[match(rownames(model$importance), names(betas))]
    res <- sort(predict(model, newdata=betas, type="prob")[1,], decreasing=TRUE)
    tibble(response = names(res)[1], prob = res[1])
}

#' classify sample.
#'
#' @param sample DNA methylation beta
#' @param model classification model
#' @param Probe_IDs methylation probes used
#' @return predicted cancer type label
#' @examples
#' library(sesameData)
#' betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
#' model <- readRDS("rfc_model.rds")
#' probes <- readRDS("HM450.rds")
#' askme_classify(betas[1], model, probes)
#' ## expect PAAD
#' @import randomForest
#' @import e1071
#' @export
askme_classify <- function(sample, model, Probe_IDs) {
    sample <- t(as.data.frame(sample))
    colnames(sample) <- Probe_IDs

    if (grepl("randomForest", model$call[1])) {
        require(randomForest)
        sample <- sample[,rownames(model$importance)]
        res <- sort(predict(model, newdata = sample, type = "prob")[1, ], decreasing = TRUE)
        tibble(response = names(res)[1], prob = res[1])
    }
    else if (grepl("svm", model$call[1])) {
        require(e1071)
        sample <- t(as.data.frame(sample[,attr(model$terms, "term.labels")]))
        res <- as.character(predict(model, newdata = sample))
        probs <- attr(predict(model, newdata = sample, probability = TRUE), "probabilities")
        prob_max <- apply(probs, MARGIN = 1, FUN = max)[1]
        tibble(response = as.character(res), prob = prob_max)
    }

}
