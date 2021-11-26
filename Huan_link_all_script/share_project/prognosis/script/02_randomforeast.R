# install.packages("randomForest")
library(randomForest)
library(pROC)
library(dplyr)
library(caret)

# data("iris")
# head(iris)
# iris.rf = randomForest(Species~.,data=iris,importance=TRUE,proximity=TRUE)

setwd("/home/huanhuan/prognosis/output/")

load("01_relate_fill.Rdata")
# org1 <-prog_roughfix%>%select(NO,pod.total,ECOG,Bsymptom,end.evaluation,R.course,progress)
org1 <-prog_roughfix%>%select(NO,pod.total,ECOG,Bsymptom,end.evaluation,progress)
org1$pod.total <-as.factor(org1$pod.total)
#---------------
load("time_and_status.Rdata")
times <-org_in%>%select(NO,status,times,FLIPI1,FLIPI2)
times$NO <-as.integer(times$NO)
org_times <-inner_join(times,org1,by="NO")
#--------
#--------------------------------------------- 
              
              #cv

#--------------------------------------------------
set.seed(1234)
folds <- createFolds(y=org1$pod.total,k=5)
rs <- data.frame()
for(i in 1:5){
  
  fold_test <- org1[folds[[i]],]   #取folds[[i]]作为测试集
  fold_train <- org1[-folds[[i]],]   # 剩下的数据作为训练集
  true_value1 =fold_test$pod.total
  num <-fold_test$NO
  fold_test<-fold_test%>%dplyr::select(-c(pod.total,NO))
  fold_train <- fold_train %>%dplyr::select(-NO)
  
  print("***组号***")
  # fold_pre <- glm(drug_repurposing ~.,family=binomial(link='logit'),data=fold_train)
  # fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
  # fold_pre <- svm(drug_repurposing ~.,data=fold_train)
#   fold_train_drug_re<-fold_train$drug_repurposing
#   fold_train<-fold_train%>%dplyr::select(-drug_repurposing)
  # fold_pre <- svm(x=fold_train,y=fold_train_drug_re,type = 'nu-regression', kernel = 'radial', gamma=0.01, cost =100)
#   fold_pre <- svm(x=fold_train,y=fold_train_drug_re,type = 'nu-regression', kernel = 'radial')
  fold_pre = randomForest(pod.total~.,data=fold_train,importance=TRUE,proximity=TRUE)
  fold_predict <- predict(fold_pre,fold_test,type="class")
  
  # fold_pre <- svm(drug_repurposing ~.,data=fold_train)
  #----------------------------------------------
  tmp<-data.frame(true_value1= true_value1,predict_value1=fold_predict,num=num)
  rs <- bind_rows(rs,tmp)
}

aaa <- multiclass.roc(rs$true_value1, rs$predict_value1)
t <-table(rs$true_value1,rs$predict_value1)
acc =sum(diag(t))/nrow(rs)*100
print(paste0("model acc:",round(acc,4),"%"))
#------------------------------------------------------

#                        predict 

#------------------------------------------------------

i=1
train_na <-org_times[which(is.na(org_times$times)),]
pre_test <-org_times[which(!is.na(org_times$times)),]
#-------------------------------
set.seed(122)
random_number <-sample(x=c(1:nrow(pre_test)), 250,replace = F)
fold_test_org<-pre_test[random_number,]
pre_train <-setdiff(pre_test,fold_test_org)
fold_train_org<-bind_rows(train_na)
# #-------------------------
# set.seed(12349)
# folds <- createFolds(y=org_times$pod.total,k=5)
# #--------------------------
# i=1
# fold_test_org <- org_times[folds[[i]],]   #取folds[[i]]作为测试集
# fold_train_org <- org_times[-folds[[i]],] 
#---------------------------------------
true_value1 =fold_test_org$pod.total
# num <-fold_test$NO
fold_test<-fold_test_org%>%dplyr::select(-c(pod.total,NO,status,times,FLIPI1,FLIPI2))
fold_train <- fold_train_org %>%dplyr::select(-c(NO,status,times,FLIPI1,FLIPI2))
fold_pre = randomForest(pod.total~.,data=fold_train,importance=TRUE,proximity=TRUE)
fold_predict <- predict(fold_pre,fold_test)
# fold_test <- org1[folds[[i]],]  
fold_test_org$predict_POD <-fold_predict

#-----------------------------------------------------------

#                             KM

#------------------------------------------------------------
library("survival")
library("survminer")

fit <- survfit(Surv(times, status) ~ predict_POD, data=fold_test_org)

# fit_km <- survfit(Surv(ECOG, Bsymptom,end.evaluation,R.course,progress) ~ predict_POD, data=fold_test)
# ,end_evaluation,R_course
# fit_km <- survfit(Surv(ECOG, Bsymptom,progress,end_evaluation) ~ predict_POD, data=fold_test)
pdf("random_foreast.pdf")
p1 <- ggsurvplot(fit)
p1 
dev.off()
#-------------------------------------
fit <- survfit(Surv(times, status) ~ FLIPI1, data=fold_test_org)
pdf("FLIPI1.pdf")
p1 <- ggsurvplot(fit)
p1 
dev.off()
fit <- survfit(Surv(times, status) ~ FLIPI2, data=fold_test_org)
pdf("FLIPI2.pdf")
p1 <- ggsurvplot(fit)
p1 
dev.off()
library(nnet)
fold_pre = multinom(pod.total~.,data=fold_train,probabilities = TRUE, model = TRUE)
fold_predict <- predict(fold_pre,fold_test)
# fold_test <- org1[folds[[i]],]  
fold_test_org$predict_POD <-fold_predict
#----------------------------------------
pdf("logistic.pdf")
p1 <- ggsurvplot(fit)
p1 
dev.off()


# plot.roc(aaa)
# dev.off()



# round(importance(iris.rf),2)



