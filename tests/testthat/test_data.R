context("data")

test_that("test cmi_predict returns tibble for RFC", {
  library(tibble)
  library(sesameData)
  library(ExperimentHub)
  library(CytoMethIC)
  betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
  eh <- ExperimentHub()
  modelrfc <- eh[["EH8395"]]
  result <- cmi_predict(betas, modelrfc, lift_over=TRUE)
  expect_is(result, "tbl_df")
})


test_that("test cmi_predict returns tibble for SVM", {
  library(tibble)
  library(sesameData)
  library(ExperimentHub)
  library(CytoMethIC)
  betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
  eh <- ExperimentHub()
  modelrfc <- eh[["EH8396"]]
  result <- cmi_predict(betas, modelrfc, lift_over=TRUE)
  expect_is(result, "tbl_df")
})

