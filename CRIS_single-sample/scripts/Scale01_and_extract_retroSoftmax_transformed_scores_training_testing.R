CRIS_CLASSES
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaleAlt<-function(x){(x-mean(as.numeric(x)))/sd(x)}
min_not_zero<-function(x){
  x[x == 0] <- 1
  min(x)}
max_not_one<-function(x){
  x[x == 1] <- 0
  max(x)}


differences<-function(x){
  diff<-rep(0,5)
  for (j in 1:length(x)-1)
  diff[j]<-as.numeric(x[j+1])-as.numeric(x[j])
  diff
  }


#load('ENV_AFTER_RUNNING.RData')

pred_train <- matrix(0, nrow = 393, ncol = 6)
pred_train<-data.frame(pred_train)
colnames(pred_train)<-c(CRIS_CLASSES, 'class prediction')
for (i in 1:nrow(pred_train)){
  pred_train[i, CRIS_CLASSES] <- ohenery::inv_smax(as.matrix(
    thresholds[["svmLinear2"]][["pred"]][i,CRIS_CLASSES]))
  pred_train[i, CRIS_CLASSES] <- scale01(pred_train[i, CRIS_CLASSES])
  pred_train[i, 'class prediction']<-as.character(thresholds[["svmLinear2"]]
                                                  [["pred"]][["predict.label2"]][i])
}



min_cls<-apply(pred_train[,1:5], 2, min_not_zero)
min_cls
max_cls<-apply(pred_train[,1:5], 2, max_not_one)
max_cls


pred_test_pdx <- matrix(0, nrow = 550, ncol = 6)
pred_test_pdx<-data.frame(pred_test_pdx)
colnames(pred_test_pdx)<-c(CRIS_CLASSES, 'class prediction')
for (i in 1:550){
  pred_test_pdx[i, CRIS_CLASSES] <- ohenery::inv_smax(as.matrix(
    ml_adapted_testing_pdx[["results"]][["svmLinear2"]][["pred"]][i,CRIS_CLASSES]))
  pred_test_pdx[i, CRIS_CLASSES] <- scale01(pred_test_pdx[i, CRIS_CLASSES])
  pred_test_pdx[i, 'class prediction']<-as.character(ml_adapted_testing_pdx[["results"]]
                                                     [["svmLinear2"]][["pred"]][["predict.label2"]][i])
}

min_cls_test_pdx<-apply(pred_test_pdx[,1:5], 2, min_not_zero)
min_cls_test_pdx
max_cls_test_pdx<-apply(pred_test_pdx[,1:5], 2, max_not_one)
max_cls_test_pdx


pred_test_tcga <- matrix(0, nrow = 169, ncol = 6)
pred_test_tcga<-data.frame(pred_test_tcga)
colnames(pred_test_tcga)<-c(CRIS_CLASSES, 'class prediction')
for (i in 1:nrow(pred_test_tcga)){
  pred_test_tcga[i, CRIS_CLASSES] <- ohenery::inv_smax(as.matrix(
    ml_adapted_testing_tcga[["results"]][["svmLinear2"]][["pred"]][i,CRIS_CLASSES]))
  pred_test_tcga[i, CRIS_CLASSES] <- scale01(pred_test_tcga[i, CRIS_CLASSES])
  pred_test_tcga[i, 'class prediction']<-as.character(ml_adapted_testing_tcga[["results"]]
                                                      [["svmLinear2"]][["pred"]][["predict.label2"]][i])
}



min_cls_test_tcga<-apply(pred_test_tcga[,1:5], 2, min_not_zero)
min_cls_test_tcga
max_cls_test_tcga<-apply(pred_test_tcga[,1:5], 2, max_not_one)
max_cls_test_tcga

pred_train_retroSoftmax_scaled01<-pred_train
pred_test_tcga_retroSoftmax_scaled01<-pred_test_tcga
pred_test_pdx_retroSoftmax_scaled01<-pred_test_pdx



write.xlsx(pred_train_retroSoftmax_scaled01, 'pred_train_retroSoftmax_scaled01.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
write.xlsx(pred_test_tcga_retroSoftmax_scaled01, 'pred_test_tcga_retroSoftmax_scaled01.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
write.xlsx(pred_test_pdx_retroSoftmax_scaled01, 'pred_test_pdx_retroSoftmax_scaled01.xlsx', col.names = TRUE,
           row.names = FALSE, append = FALSE)
