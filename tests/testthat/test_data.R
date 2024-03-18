context("data")

test_that("Impute mean functions properly", {
  df <- data.frame(a = c(NA, 2, 3, 4), b = c(1, 6, 5, NA), c = c(3, 6, 7, 8))
  df_imputed_cols <- CytoMethIC::impute_mean(df)
  df_imputed_rows <- CytoMethIC::impute_mean(df, axis = 2)
  expect_true(df_imputed_cols[1,1] == 3)
  expect_true(df_imputed_cols[4,2] == 4)
  expect_true(df_imputed_rows[1,1] == 2)
  expect_true((df_imputed_rows[4,2]) == 6)

})

test_that("test cmi_classify returns tibble for RFC", {
  library(tibble)
  library(sesameData)
  library(ExperimentHub)
  betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
  eh <- ExperimentHub()
  modelrfc <- eh[["EH8395"]]
  result <- CytoMethIC::cmi_classify(betas, modelrfc, lift_over=TRUE)
  expect_is(result, "tbl_df")
})


test_that("test cmi_classify returns tibble for SVM", {
  library(tibble)
  library(sesameData)
  library(ExperimentHub)
  betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
  eh <- ExperimentHub()
  modelrfc <- eh[["EH8396"]]
  result <- CytoMethIC::cmi_classify(betas, modelrfc, lift_over=TRUE)
  expect_is(result, "tbl_df")
})

