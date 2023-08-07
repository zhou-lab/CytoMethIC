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
    for (i in seq_len(ncol(df))) {
      df[is.na(df[, i]), i] <- mean(df[, i], na.rm = TRUE)
    }
  } else if (axis == 2) {
    for (i in seq_len(nrow(df))){
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

#' Infer sex.
#'
#' @param betas DNA methylation beta
#' @return Inferred sex of sample
#' @examples
#' library(sesameData)
#' betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
#' inferSex(betas)
#' @export



inferSex <- function(betas) {
  hypoMALE <- c(
    "cg21983484","cg23696472","cg11673471","cg01742836","cg13574945",
    "cg08059778","cg24186901","cg26023405","cg15977272","cg13023833",
    "cg20766178","cg20455959","cg26584339","cg13130271","cg13244998",
    "cg05872808","cg21290550","cg05806018","cg07861180","cg20015269",
    "cg12576145","cg10991108","cg02333283","cg16357225","cg25206026",
    "cg20749341","cg03773146","cg04872051","cg16590821","cg09520212",
    "cg22221554","cg11152253","cg23429746","cg00813156","cg25132467",
    "cg16221895","cg09307104","cg15165114","cg18998000","cg00723973",
    "cg06041068","cg10860619","cg09514431","cg07912337","cg03334316",
    "cg17399684","cg05534333","cg23493872","cg12413138","cg05374090",
    "cg27501007","cg08855111","cg21159768","cg16488754","cg12075609",
    "cg07446674","cg01342901","cg02869694","cg12277627","cg19992190",
    "cg10717149","cg14191108","cg01869765","cg26505478","cg23685102",
    "cg02195366","cg06334238","cg02615131","cg15565409","cg15693668",
    "cg03505772","cg00845806","cg26439324","cg12935118","cg18932686",
    "cg24264679","cg08782677","cg13649400","cg06779802","cg23554546",
    "cg23951868","cg00337921","cg08479532","cg00114625","cg03391801",
    "cg22776211","cg07674503","cg22452543","cg18140045","cg15450782",
    "cg07674075","cg06510592","cg21137943","cg24479484","cg27501723",
    "cg20439892","cg18107314","cg08405463","cg09146364","cg16894263")

  hyperMALE <- c(
    "cg26359388","cg02540440","cg11049634","cg22874828","cg09182733",
    "cg01123965","cg15822015","cg05130312","cg17072671","cg22655232",
    "cg05695959","cg21010298","cg06143713","cg22759686","cg11143827",
    "cg04303560","cg11717280","cg14372935","cg05533223","cg16405492",
    "cg15765801","cg08156775","cg24183173","cg21797452","cg03161453",
    "cg10474871","cg11516614","cg18813691","cg08614574","cg08456555",
    "cg16440909","cg13326840","cg16822540","cg03801901","cg09039264",
    "cg01383599","cg14931238","cg04071644","cg22208280","cg05559023",
    "cg23317607","cg26327984","cg07801607","cg06870560","cg24156613",
    "cg04101819","cg07422795","cg14261068","cg12622895","cg09192294",
    "cg26695278","cg12653510","cg03554089","cg11166197","cg04032096",
    "cg25047306","cg07818713","cg21258987","cg07981033","cg14492530",
    "cg18157587","cg12030638","cg17498624","cg01816615","cg08723064",
    "cg05193067","cg27167763","cg15521097","cg25456959","cg16576300",
    "cg07318999","cg22417678","cg22671388","cg23644934","cg00267352",
    "cg22223709","cg23698976","cg06780606","cg13920260","cg15861835",
    "cg10039267","cg12454245","cg22067189","cg00150874","cg08401365",
    "cg13781721","cg02931660","cg01316390","cg14746118","cg21294096",
    "cg11871337","cg00408231","cg09641151","cg05226646","cg11291200",
    "cg01109660","cg23607813","cg04624564","cg07452499","cg18123612")
  vals <- mean(betas[hyperMALE], na.rm = TRUE) - betas[hypoMALE]
  dd <- density(na.omit(vals))
  if (dd$x[which.max(dd$y)] > 0.4) {
    "MALE"
  } else {
    "FEMALE"
  }
}