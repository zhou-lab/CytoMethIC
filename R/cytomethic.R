#' Master data frame for all model objects
#'
#' This is an internal object which will be updated on every new release
#'
#' @name cmi_models
#' @docType data
#' @return master sheet of CytoMethIC model objects
#' @examples print(cmi_models[,c("EHID","Title")])
#' @export
NULL

clean_features <- function(
    betas, cmi_model,
    source_platform = source_platform,
    lift_over = lift_over, verbose = verbose) {
    
    features <- cmi_model$features
    idx <- match(features, names(betas))
    if (sum(is.na(idx)) == 0 && sum(is.na(betas[idx])) == 0) {
        return(betas[idx])
    }
    
    if (verbose) {
        warning(sprintf("Missing %d/%d features.",
            sum(is.na(betas[idx])), length(idx)))
    }

    if (lift_over) {
        source_platform <- sesameData::sesameData_check_platform(
            source_platform, rownames(betas))
        target_platform <- cmi_model$feature_platform
        if (is.null(target_platform)) {
            target_platform <- "HM450"
        }
        betas <- mLiftOver(betas, source_platform = source_platform,
            target_platform = target_platform, impute=TRUE)
    }

    ## if still have missing values
    idx <- match(features, names(betas))
    if (sum(!is.na(idx)) == 0 && !lift_over) {
        stop("No overlapping probes. Consider lift_over=TRUE")
    }
    if (sum(!is.na(betas[idx])) == 0 ||  # all-NA
        (sum(is.na(betas[idx])) > 0 && ( # some NA and model requires all
            is.null(cmi_model$features_require_all) ||
            cmi_model$features_require_all))) {
        stop(sprintf("Missing %d/%d features. Consider lift_over=TRUE",
            sum(is.na(betas[idx])), sum(!is.na(idx))))
    }
    betas <- betas[features]
}

#' The cmi_predict function takes in a model and a sample, and uses the model
#' to predict it.  This function supports randomForest, e1071::svm, xgboost,
#' and keras/tensorflow models. For xgboost and keras models, the features used
#' in classification as well as a label mapping must be provided for output.
#' 
#' @param betas DNA methylation beta
#' @param cmi_model Cytomethic model downloaded from ExperimentHub
#' @param source_platform source platform
#' If not given, will infer from probe ID.
#' @param lift_over whether to allow mLiftOver to convert probe IDs
#' @param BPPARAM use MulticoreParam(n) for parallel processing
#' @param verbose be verbose with warning
#' @return predicted cancer type label
#' @examples
#' 
#' library(sesame)
#' library(ExperimentHub)
#' library(CytoMethIC)
#'
#' ## Cancer Type
#' model = ExperimentHub()[["EH8395"]]
#' cmi_predict(openSesame(sesameDataGet("EPICv2.8.SigDF")[[1]]), model, lift_over=TRUE)
#' cmi_predict(openSesame(sesameDataGet('EPIC.1.SigDF')), model, lift_over=TRUE)
#' cmi_predict(sesameDataGet("HM450.1.TCGA.PAAD")$betas, model, lift_over=TRUE)
#'
#' @import stats
#' @import tools
#' @import sesameData
#' @import ExperimentHub
#' @import BiocParallel
#' @importFrom tibble tibble
#' @importFrom sesame mLiftOver
#' @importFrom methods is
#' @export
cmi_predict <- function(betas, cmi_model,
    source_platform = NULL, lift_over = FALSE, verbose = FALSE,
    BPPARAM = SerialParam()) {
    
    if(names(cmi_model)[[1]] == "model_serialized") {
        requireNamespace("keras")
        requireNamespace("tensorflow")
        cmi_model <- list(
            model = keras::unserialize_model(cmi_model$model_serialized),
            features = cmi_model$features,
            label_levels = cmi_model$label_levels)
    }

    if (is.matrix(betas)) {
        return(do.call(rbind, bplapply(seq_len(ncol(betas)), function(i) {
            cmi_predict(betas[,i], cmi_model, source_platform = source_platform,
                lift_over = lift_over, verbose = verbose)
        }, BPPARAM = BPPARAM)))
    }
    stopifnot(is.numeric(betas))

    betas <- clean_features(betas, cmi_model,
        source_platform = source_platform,
        lift_over = lift_over, verbose = verbose)

    if (is(cmi_model$model, "function")) {
        cmi_model$model(betas)
    } else if (is(cmi_model$model, "randomForest")) {
        requireNamespace("randomForest")
        res <- sort(predict(cmi_model$model,
            newdata = t(betas), type = "prob")[1,], decreasing = TRUE)
        tibble(response = names(res)[1], prob = res[1])
    } else if (is(cmi_model$model, "svm")) {
        betas <- t(as.data.frame(betas))
        if (!requireNamespace("e1071", quietly = TRUE)) stop("e1071 not installed")
        model <- cmi_model$model
        betas <- t(as.data.frame(betas[, attr(model$terms, "term.labels")]))
        res <- as.character(predict(model, newdata = betas))
        probs <- attr(predict(model, newdata = betas, probability = TRUE), "probabilities")
        prob_max <- apply(probs, MARGIN = 1, FUN = max)[1]
        tibble(response = as.character(res), prob = prob_max)
    } else if (is(cmi_model$model, "xgb")) {
        betas <- t(as.data.frame(betas))
        if (!requireNamespace("xgboost", quietly = TRUE)) stop("xgboost not installed")
        if (is.null(cmi_model$features)) stop("Must provide feature parameter with xgboost model")
        if (is.null(cmi_model$label_levels)) stop("Must provide label_levels parameter with xgboost model")
        betas <- (t(as.data.frame(betas[, cmi_model$features])))
        model <- cmi_model$model
        betas <- xgboost::xgb.DMatrix(t(as.matrix(betas)))
        pred_probabilities <- predict(model, betas)
        num_classes <- length(pred_probabilities)
        pred_prob_matrix <- matrix(pred_probabilities, nrow = 1, ncol = num_classes, byrow = TRUE)
        max_probability <- apply(pred_prob_matrix, 1, max)
        pred_label <- cmi_model$label_levels[apply(pred_prob_matrix, 1, which.max)]
        tibble(response = pred_label, prob = max_probability)
    } else if (is(cmi_model$model, "keras")) {
        if (is.null(cmi_model$features)) stop("Must provide feature parameter with Keras model")
        if (is.null(cmi_model$label_levels)) stop("Must provide label_levels parameter with Keras model")
        betas <- t(as.data.frame(betas))
        betas <- t(as.data.frame(betas[, cmi_model$features]))
        betas <- t(as.matrix(betas))
        pred_prob_matrix <- predict(cmi_model$model, betas)
        max_probability <- apply(pred_prob_matrix, 1, max)
        highest_prob_prediction <- apply(pred_prob_matrix, 1, function(x) which.max(x))
        pred_label <- cmi_model$label_levels[highest_prob_prediction]
        tibble(response = pred_label, prob = max_probability)
    } else {
        stop("Package not supported")
    }
}



#' Check CytoMethIC versions
#'
#' print package verison of cytomethic and depended packages to help
#' troubleshoot installation issues.
#'
#' @return print the versions of cytomethic and dependencies
#' @importFrom utils packageVersion
#' @export
#' @examples
#' cmi_checkVersion()
cmi_checkVersion <- function() {
    rv <- R.Version()
    msg <- paste0(
        "CytoMethIC requires matched versions of ",
        "R, sesame, sesameData and ExperimentHub.\n",
        "Here is the current versions installed:\n",
        sprintf("R: %s.%s\n", rv$major, rv$minor),
        sprintf("Bioconductor: %s\n", BiocManager::version()),
        sprintf("CytoMethIC: %s\n", packageVersion("CytoMethIC")),
        sprintf("sesame: %s\n", packageVersion("sesame")),
        sprintf("sesameData: %s\n", packageVersion("sesameData")),
        sprintf("ExperimentHub: %s\n", packageVersion("ExperimentHub")))
    message(msg)
}
