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
        df <- data.frame(lapply(df, function(x) {
            x[is.na(x)] <- mean(x, na.rm = TRUE)
            return(x)
        }))
    } else if (axis == 2) {
        df <- t(apply(df, 1, function(x) {
            x[is.na(x)] <- mean(x, na.rm = TRUE)
            return(x)
        }))
    } else {
        stop("Invalid axis. Use 1 for columns or 2 for rows.")
    }
    return(df)
}

#' Impute Missing Values with Mean, for CytoMethIC models ONLY
#' This function replaces missing values (NA) in a data frame or matrix, default is col means.
#'
#' @param df A dataframe or matrix
#' @return A data frame or matrix with missing values imputed.


impute_mean_cmi <- function(df) {
  if (length(colnames(df)) == 6636) {
    col_indices <- seq_len(length(colnames(df)))
    df[1, col_indices] <- mapply(function(x, y) {
      if(is.na(x)) y else x
    }, df[1, col_indices], BrainTumorClassifierMeanValues[1, col_indices])
  }
  else {
    col_indices <- seq_len(length(colnames(df)))
    df[1, col_indices] <- mapply(function(x, y) {
      if(is.na(x)) y else x
    }, df[1, col_indices], PanCancerClassifierMeanValues[1, col_indices])
  }
  df
}


#' The cmi_classify function takes in a model and a sample, and uses the model to classify it.
#' This function supports randomForest, e1071::svm, xgboost, and keras/tensorflow models. For xgboost and keras models,
#' the features used in classification as well as a label mapping must be provided for output.
#' 
#' @param betas DNA methylation beta
#' @param cmi_model Cytomethic model downloaded from ExperimentHub
#' @param source_platform source platform
#' If not given, will infer from probe ID.
#' @param lift_over whether to allow liftOver to convert probe IDs
#' @return predicted cancer type label
#' @examples
#' 
#' library(sesameData)
#' library(ExperimentHub)
#' basedir = "https://github.com/zhou-lab/CytoMethIC_models/raw/main/models/"
#' 
#' ## Cancer Type
#' model = ExperimentHub()[["EH8395"]]
#' cmi_classify(openSesame(sesameDataGet("EPICv2.8.SigDF")[[1]]), model)
#' cmi_classify(openSesame(sesameDataGet('EPIC.1.SigDF')), model)
#' cmi_classify(sesameDataGet("HM450.1.TCGA.PAAD")$betas, model)
#'
#' ## Sex
#' model = readRDS(url(sprintf("%s/Sex2_HM450.rds", basedir)))
#' cmi_classify(openSesame(sesameDataGet("EPICv2.8.SigDF")[[1]]), model)
#' cmi_classify(openSesame(sesameDataGet('EPIC.1.SigDF')), model)
#' cmi_classify(sesameDataGet("HM450.1.TCGA.PAAD")$betas, model)
#' 
#' \dontrun{
#' ## Ethnicity
#' model = readRDS(url(sprintf("%s/Race3_rfcTCGA_InfHum3.rds", basedir)))
#' cmi_classify(openSesame(sesameDataGet("EPICv2.8.SigDF")[[1]]), model)
#' cmi_classify(openSesame(sesameDataGet('EPIC.1.SigDF')), model)
#' cmi_classify(sesameDataGet("HM450.1.TCGA.PAAD")$betas, model)
#' 
#' }
#' @import stats
#' @import tools
#' @import sesameData
#' @import ExperimentHub
#' @importFrom tibble tibble
#' @importFrom sesame liftOver
#' @importFrom methods is
#' @export
cmi_classify <- function (betas, cmi_model, source_platform = NULL,
    lift_over = TRUE) { #Change model_list to cmi_model
    
    if(names(cmi_model)[[1]] == "model_serialized") {
        requireNamespace("keras")
        requireNamespace("tensorflow")
        cmi_model <- list(model = keras::unserialize_model(
            cmi_model[["model_serialized"]]),
            features = cmi_model[["features"]],
            label_levels = cmi_model[["label_levels"]])
    }

    if (!is.matrix(betas)) {
        betas <- cbind(betas)
    }

    features <- cmi_model$features
    idx <- match(features, rownames(betas))
    if (sum(is.na(idx)) > 0 || sum(is.na(betas[idx,])) > 0) {
        if (lift_over) {
            source_platform <- sesameData::sesameData_check_platform(
                source_platform, rownames(betas))
            target_platform <- cmi_model$feature_platform
            if (is.null(target_platform)) {
                target_platform <- "HM450"
            }
            betas <- liftOver(betas, source_platform = source_platform,
                target_platform = target_platform, impute=TRUE)
        } else {
            stop("Missing data. Consider turning on lift_over to do probe ID conversion and imputation.")
        }
    }
    betas <- betas[features,,drop=FALSE]
    if (cmi_model$model_name == "Threshold-based Sex Model") {
        vals <- mean(betas[cmi_model$model$hyperMALE,1], na.rm = TRUE) -
            betas[cmi_model$model$hypoMALE,1]
        dd <- density(na.omit(vals))
        if (dd$x[which.max(dd$y)] > 0.4) {
            res <- "MALE"
        } else {
            res <- "FEMALE"
        }
        tibble(response = res, prob = NA)
    } else if (is(cmi_model[["model"]], "randomForest")) {
        requireNamespace("randomForest")
        res <- sort(predict(cmi_model$model,
            newdata = t(betas), type = "prob")[1,], decreasing = TRUE)
        tibble(response = names(res)[1], prob = res[1])
    } else if (is(cmi_model[["model"]], "svm")) {
        betas <- t(as.data.frame(betas))
        if (!requireNamespace("e1071", quietly = TRUE)) stop("e1071 not installed")
        model <- cmi_model[["model"]]
        betas <- t(as.data.frame(betas[, attr(model$terms, "term.labels")]))
        betas <- impute_mean_cmi(betas)
        res <- as.character(predict(model, newdata = betas))
        probs <- attr(predict(model, newdata = betas, probability = TRUE), "probabilities")
        prob_max <- apply(probs, MARGIN = 1, FUN = max)[1]
        tibble(response = as.character(res), prob = prob_max)
    } else if (is(cmi_model[["model"]], "xgb")) {
        betas <- t(as.data.frame(betas))
        if (!requireNamespace("xgboost", quietly = TRUE)) stop("xgboost not installed")
        if (setequal(features, NULL)) stop("Must provide feature parameter with xgboost model")
        if (setequal(cmi_model$label_levels, NULL)) stop("Must provide label_levels parameter with xgboost model")
        betas <- impute_mean_cmi(t(as.data.frame(betas[, features])))
        model <- cmi_model[["model"]]
        betas <- xgboost::xgb.DMatrix(t(as.matrix(betas)))
        pred_probabilities <- predict(model, betas)
        num_classes <- length(pred_probabilities)
        pred_prob_matrix <- matrix(pred_probabilities, nrow = 1, ncol = num_classes, byrow = TRUE)
        max_probability <- apply(pred_prob_matrix, 1, max)
        pred_label <- cmi_model$label_levels[apply(pred_prob_matrix, 1, which.max)]
        tibble(response = pred_label, prob = max_probability)
    } else if (is(cmi_model[["model"]], "keras")) {
        betas <- t(as.data.frame(betas))
        if (setequal(features, NULL)) stop("Must provide feature parameter with Keras model")
        if (setequal(cmi_model$label_levels, NULL)) stop("Must provide label_levels parameter with Keras model")
        betas <- impute_mean_cmi(t(as.data.frame(betas[, features])))
        model <- cmi_model[["model"]]
        betas <- t(as.matrix(betas))
        pred_prob_matrix <- predict(model, betas)
        max_probability <- apply(pred_prob_matrix, 1, max)
        highest_prob_prediction <- apply(pred_prob_matrix, 1, function(x) which.max(x))
        pred_label <- cmi_model$label_levels[highest_prob_prediction]
        tibble(response = pred_label, prob = max_probability)
    } else {
        stop("Package not supported")
    }
}
