#' Master data frame for all model objects
#'
#' This is an internal object which will be updated on every new release
#'
#' @name cmi_models
#' @docType data
#' @return master sheet of CytoMethIC model objects
#' @examples print(cmi_models$ModelID)
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

#' Impute Missing Values with Mean
#' This function replaces missing values (NA) in a data frame or matrix, default is col means.
#'
#' @param df A dataframe or matrix
#' @param axis A single integer. Use 1 to impute column means (default), and 2 to impute row means.
#' @return A data frame or matrix with missing values imputed.
#' @examples
#' df <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 2, 3, 4))
#' df <- impute_mean(df)
#' df <- impute_mean(df, axis = 2)
#' @export
impute_mean <- function(df, axis = 1) {
  if (axis == 1) {
    for (i in 1:ncol(df)) {
      df[is.na(df[, i]), i] <- mean(df[, i], na.rm = TRUE)
    }
  } else if (axis == 2) {
    for (i in 1:nrow(df)) {
      df[i, is.na(df[i, ])] <- mean(df[i, ], na.rm = TRUE)
    }
  } else {
    stop("Invalid axis. Use 1 for columns or 2 for rows.")
  }
  return(df)
}



#' classify sample.
#'
#' @param betas DNA methylation beta
#' @param model classification model
#' @param feature list of features if not stored within model
#' @param label_levels factor-label mapping if not stored within model
#' @return predicted cancer type label
#' @examples
#' library(sesameData)
#' betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
#' cmi_classify(betas, m_cancertype_CNS66)
#' ## expect PAAD
#' @import randomForest
#' @import stats
#' @import tools
#' @export
cmi_classify <- function(betas, model, feature = NULL, label_levels = NULL) {
  betas <- t(as.data.frame(betas))
  if (grepl("randomForest", class(model)[1])) {
    if (!requireNamespace("randomForest", quietly = TRUE)) stop("randomForest not installed")
    feature <- rownames(model$importance)
    betas <- betas[, feature]
    res <- sort(predict(model, newdata = betas, type = "prob")[1, ], decreasing = TRUE)
    tibble::tibble(response = names(res)[1], prob = res[1])
  } else if (grepl("svm", class(model)[1])) {
    if (!requireNamespace("e1071", quietly = TRUE)) stop("e1071 not installed")
    betas <- t(as.data.frame(betas[, attr(model$terms, "term.labels")]))
    res <- as.character(predict(model, newdata = betas))
    probs <- attr(predict(model, newdata = betas, probability = TRUE), "probabilities")
    prob_max <- apply(probs, MARGIN = 1, FUN = max)[1]
    tibble::tibble(response = as.character(res), prob = prob_max)
  } else if (grepl("xgb", class(model)[1])) {
    if (!requireNamespace("xgboost", quietly = TRUE)) stop("xgboost not installed")
    if (setequal(feature, NULL)) stop("Must provide feature parameter with xgboost model")
    if (setequal(label_levels, NULL)) stop("Must provide label_levels parameter with xgboost model")
    betas <- betas[, feature]
    betas <- xgboost::xgb.DMatrix(t(as.matrix(betas)))
    pred_probabilities <- predict(model, betas)
    num_classes <- length(pred_probabilities)
    pred_prob_matrix <- matrix(pred_probabilities, nrow = 1, ncol = num_classes, byrow = TRUE)
    max_probability <- apply(pred_prob_matrix, 1, max)
    pred_label <- label_levels[apply(pred_prob_matrix, 1, which.max)]
    tibble::tibble(response = pred_label, prob = max_probability)
  } else if (grepl("keras", class(model)[1])) {
    if (!requireNamespace("keras", quietly = TRUE)) stop("keras not installed")
    if (!requireNamespace("tensorflow", quietly = TRUE)) stop("tensorflow not installed")
    if (setequal(feature, NULL)) stop("Must provide feature parameter with Keras model")
    if (setequal(label_levels, NULL)) stop("Must provide label_levels parameter with Keras model")
    betas <- betas[, feature]
    betas <- t(as.matrix(betas))
    pred_prob_matrix <- model %>% predict(betas)
    max_probability <- apply(pred_prob_matrix, 1, max)
    highest_prob_prediction <- apply(pred_prob_matrix, 1, function(x) which.max(x))
    pred_label <- label_levels[highest_prob_prediction]
    tibble::tibble(response = pred_label, prob = max_probability)
  } else {
    stop("Package not supported")
  }
}


