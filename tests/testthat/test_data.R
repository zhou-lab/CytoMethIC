context('data')

## test_that("test='HM27.address' gives correct data", {
##     dt <- sesameDataGet('HM27.address')
##     expect_is(dt, "list")
## })


library(sesameData)
setwd("/Users/jfanale/Documents/GitHub/askme")

betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
betas <- t(as.data.frame(readRDS("/Users/jfanale/HM450/GSM2403235.rds")))
colnames(betas) <- readr::read_tsv("/Users/jfanale/Desktop/Research/Data/LABEL/HM450.tsv")$HM450_Probe_ID


rfc_model <- readRDS("/Users/jfanale/Desktop/Final_Askme_Models/rfc_Capper_CNSTumor.rds")
askme_classify(betas, rfc_model)

svm_model <- readRDS("/Users/jfanale/Desktop/Final_Askme_Models/svm_Capper_CNSTumor.rds")
askme_classify(betas, svm_model)

load("/Users/jfanale/Desktop/Final_Askme_Models/xgb_Capper_CNSTumor66_pkg.rda")
askme_classify(betas, xgb_model, feature = xgb_feature_6636$Probe_ID,
                      label_levels = xgb_levels)


load("/Users/jfanale/Desktop/Final_Askme_Models/mlp_Capper_CNSTumor66_pkg.rda")
