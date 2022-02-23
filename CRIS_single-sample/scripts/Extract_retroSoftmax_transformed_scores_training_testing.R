CRIS_CLASSES

load('ENV_AFTER_RUNNING.RData')

pred_train <- matrix(0, nrow = 393, ncol = 6)
pred_train<-data.frame(pred_train)
colnames(pred_train)<-c(CRIS_CLASSES, 'class prediction')
for (i in 1:nrow(pred_train)){
  pred_train[i, CRIS_CLASSES] <- ohenery::inv_smax(as.matrix(
    thresholds[["svmLinear2"]][["pred"]][i,CRIS_CLASSES]))
  pred_train[i, 'class prediction']<-as.character(thresholds[["svmLinear2"]]
                                                  [["pred"]][["predict.label2"]][i])
}


min_cls<-apply(pred_train[,1:5], 2, min)
min_cls
max_cls<-apply(pred_train[,1:5], 2, max)
max_cls


pred_test_pdx <- matrix(0, nrow = 550, ncol = 6)
pred_test_pdx<-data.frame(pred_test_pdx)
colnames(pred_test_pdx)<-c(CRIS_CLASSES, 'class prediction')
for (i in 1:550){
  pred_test_pdx[i, CRIS_CLASSES] <- ohenery::inv_smax(as.matrix(
    ml_adapted_testing_pdx[["results"]][["svmLinear2"]][["pred"]][i,CRIS_CLASSES]))
  pred_test_pdx[i, 'class prediction']<-as.character(ml_adapted_testing_pdx[["results"]]
                                [["svmLinear2"]][["pred"]][["predict.label2"]][i])
}

min_cls_test_pdx<-apply(pred_test_pdx[,1:5], 2, min)
min_cls_test_pdx
max_cls_test_pdx<-apply(pred_test_pdx[,1:5], 2, max)
max_cls_test_pdx


pred_test_tcga <- matrix(0, nrow = 169, ncol = 6)
pred_test_tcga<-data.frame(pred_test_tcga)
colnames(pred_test_tcga)<-c(CRIS_CLASSES, 'class prediction')
for (i in 1:nrow(pred_test_tcga)){
  pred_test_tcga[i, CRIS_CLASSES] <- ohenery::inv_smax(as.matrix(
    ml_adapted_testing_tcga[["results"]][["svmLinear2"]][["pred"]][i,CRIS_CLASSES]))
  pred_test_tcga[i, 'class prediction']<-as.character(ml_adapted_testing_tcga[["results"]]
                                                     [["svmLinear2"]][["pred"]][["predict.label2"]][i])
}

min_cls_test_tcga<-apply(pred_test_tcga[,1:5], 2, min)
min_cls_test_tcga
max_cls_test_tcga<-apply(pred_test_tcga[,1:5], 2, max)
max_cls_test_tcga

pred_train_scaled<-pred_train
pred_test_tcga_scaled<-pred_test_tcga
pred_test_pdx_scaled<-pred_test_pdx

for (c in CRIS_CLASSES){
  min_not_inf <- min_cls[c]
  max_not_inf <- max_cls[c]
  pred_train_scaled[pred_train_scaled[,c] == Inf, c]  <- max_not_inf
  pred_train_scaled[pred_train_scaled[,c] == -Inf, c] <- min_not_inf
  pred_test_tcga_scaled[pred_test_tcga_scaled[,c] == Inf, c]  <- max_not_inf
  pred_test_tcga_scaled[pred_test_tcga_scaled[,c] == -Inf, c] <- min_not_inf
  pred_test_pdx_scaled[pred_test_pdx_scaled[,c] == Inf, c]  <- max_not_inf
  pred_test_pdx_scaled[pred_test_pdx_scaled[,c] == -Inf, c] <- min_not_inf
  if (max_not_inf == min_not_inf){
    pred_train_scaled[,c] <- 0
    pred_test_tcga_scaled[,c] <- 0
    pred_test_pdx_scaled[,c] <- 0
  }
  else{
    pred_train_scaled[,c] <- as.double(pred_train[,c] - min_not_inf)/as.double(max_not_inf - min_not_inf)
    pred_test_tcga_scaled[,c] <- as.double(pred_test_tcga[,c] - min_not_inf)/as.double(max_not_inf - min_not_inf)
    pred_test_pdx_scaled[,c] <- as.double(pred_test_pdx[,c] - min_not_inf)/as.double(max_not_inf - min_not_inf)
 }
}

write.xlsx(pred_train, 'pred_train_retroSoftmax.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
write.xlsx(pred_test_tcga, 'pred_test_tcga_retroSoftmax.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
write.xlsx(pred_test_pdx, 'pred_test_pdx_retroSoftmax.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)

write.xlsx(pred_train_scaled, 'pred_train_retroSoftmax_scaled.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
write.xlsx(pred_test_tcga_scaled, 'pred_test_tcga_retroSoftmax_scaled.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
write.xlsx(pred_test_pdx_scaled, 'pred_test_pdx_retroSoftmax_scaled.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
