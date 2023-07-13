#' Master data frame for all model objects
#'
#' This is an internal object which will be updated on every new release
#' 
#' @name askme_models
#' @docType data
#' @return master sheet of askme model objects
#' @examples print(askme_models$ModelID)
#' @export
NULL

#' Master data frame for all prediction labels
#'
#' This is an internal object which will be updated on every new release
#' 
#' @name prediction_labels
#' @docType data
#' @return master sheet of prediction labels
#' @examples print(prediction_labels)
#' @export
NULL
#' Model for CNS cancer type classification
#'
#' @name m_cancertype_CNS66
#' @docType data
#' @return random forest model
#' @examples print(m_cancertype_CNS66$ntree)
#' @export
NULL

#' Model for TCGA cancer type classification
#'
#' @name m_cancertype_TCGA33
#' @docType data
#' @return random forest model
#' @examples print(m_cancertype_TCGA33$ntree)
#' @export
NULL


#' classify sample.
#'
#' @param betas DNA methylation beta
#' @param model classification model
#' @return predicted cancer type label
#' @examples
#' library(sesameData)
#' betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
#' askme_classify(betas, m_cancertype_TCGA33)
#' ## expect PAAD
#' @import randomForest
#' @import e1071
#' @export
askme_classify <- function(betas, model) {
    betas <- t(as.data.frame(betas))
    if (grepl("randomForest", class(model)[1])) {
        require(randomForest)
        feature <- rownames(model$importance)
        betas <- betas[,feature]
        res <- sort(predict(model, newdata = betas, type = "prob")[1, ], decreasing = TRUE)
        tibble::tibble(response = names(res)[1], prob = res[1])
    }
    else if (grepl("svm", class(model)[1])) {
        require(e1071)
        betas <- t(as.data.frame(betas[,attr(model$terms, "term.labels")]))
        res <- as.character(predict(model, newdata = betas))
        probs <- attr(predict(model, newdata = betas, probability = TRUE), "probabilities")
        prob_max <- apply(probs, MARGIN = 1, FUN = max)[1]
        tibble::tibble(response = as.character(res), prob = prob_max)
    }
    else if(grepl("xgb", class(model)[1])) {
        require(xgboost)
        feature <- as.data.frame(xgboost::xgb.importance(model=model))$Feature
        betas <- betas[, feature]
        betas <- xgb.DMatrix(t(as.matrix(betas)))
        pred_probabilities <- predict(xgbModel, betas)

    }

}

