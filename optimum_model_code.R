rm(list = ls())
library(survival)
library(caret)
library(plotrix)
library(survminer)
library(glmnet)
library(gbm)
library(randomForestSRC)

#data preprocessing

dat = read.csv("sarcoma.csv", stringsAsFactors = TRUE)
dat = subset(dat, select = -c(X,surtim,surind))
dat = na.omit(dat)
rownames(dat) <- 1:nrow(dat)
dat$sex =( as.numeric(as.factor(dat$sex)))
dat$grade = (as.numeric(as.factor(dat$grade)))
dat$tumor.subtype =( as.numeric(as.factor(dat$tumor.subtype)))

dim(dat)
data_subset = subset(dat,select=c(sum.var_GLCM,energy_GLCM,grade,max.grad_HIST,H0,pca.elongation,sex,Reg.grad.min,dfsind,dfstim))
repeats = 10
folds = 5
total = folds*repeats
threshold_limit = (60*total)/100
original_columns = ncol(dat)-2
subset_columns = ncol(data_subset)-2
subset_rows = nrow(data_subset)
subset_column_names = colnames(subset(data_subset, select = -c(dfstim, dfsind)))


# VAriable declaration

stp_pred_matrx = matrix(0,nrow=total,ncol=subset_rows)
full_step_model = vector("list",folds)
final_stp_train_con  = vector("list",total)
final_stp_test_con = vector("list",total)
final_stp_models = vector("list",total)
lasso_train_con = lasso_test_con = numeric(total)
lasso_pred_matrix  = matrix(0, nrow=total, ncol=nrow(data_subset))
lasso.feat.selected = matrix(0, nrow=total, ncol=subset_columns)
colnames(lasso.feat.selected) = subset_column_names
ridge_train_con = ridge_test_con =  numeric(total)
ridge_pred_matrix  = matrix(0, nrow=total, ncol=nrow(data_subset))
gbm_summary = vector("list",total) 
gbm.train.con = gbm.test.con = vector("double",total) 
gbm_pred_matrix = matrix(0,nrow=total,ncol=subset_rows)
rsf_model = rsf_pred = vector("list",total)
rsf_con_train = rsf_bs = vector("list",total)
rsf_con_test = vector("list",total)
rsf_con_test = vector("list",total)

## Kaplan Meier curve

split.kms <- function(zos,os,oevn,NQ=100,zqmin=.05,zqmax=.95){
  qmin = quantile(zos,zqmin,na.rm=TRUE)
  qmax = quantile(zos,zqmax,na.rm=TRUE)
  p0q = numeric(NQ)
  med0s = seq(qmin,qmax,length=NQ)
  izs = matrix(NA,ncol=NQ,nrow=length(zos))
  for(iq in 1:NQ){
    IZ = as.integer(zos<med0s[iq])
    p0q[iq] = 1-pchisq(survdiff(Surv(os,oevn)~factor(IZ))$chisq,1)
    izs[,iq] = IZ
  }
  best.med0 = med0s[which.min(p0q)]
  IZ = as.integer(zos<best.med0)
  return(list(iz=IZ,izs=izs,best.t=best.med0,
              ps=p0q,best.p=p0q[which.min(p0q)]))
}

flag = 1
set.seed(6)
for(a in 1:repeats){
  shuffled_dat = data_subset[sample(1:nrow(data_subset),nrow(data_subset)), ]
  
  cvIndex = createFolds(factor(shuffled_dat$dfsind), folds,returnTrain = T)
  for(b in 1:length(cvIndex)){
    train_data = data_subset[cvIndex[[b]],]
    test_data = data_subset[-cvIndex[[b]],]
    
    stp_cox_fit = coxph(Surv(dfstim,dfsind)~.,data=train_data)
    stp_cox_pred = predict(stp_cox_fit,newdata = test_data,type="lp")
    
    for(sw_pred_op in 1:length(stp_cox_pred)){
      stp_var = stp_cox_pred[sw_pred_op]
      stp_indx = as.numeric(names(stp_var))
      stp_pred_matrx[flag,stp_indx] = as.numeric(stp_var)[1]
    }
    
    final_stp_train_con[flag]<-concordance(stp_cox_fit, newdata = train_data)$concordance
    final_stp_test_con[flag] <- concordance(stp_cox_fit, newdata = test_data)$concordance
    
    x.train = subset(train_data, select = -c(dfstim, dfsind))
    y.train = subset(train_data, select = c(dfstim, dfsind))
    surv.train = Surv(y.train$dfstim,y.train$dfsind)
    x.train.m = model.matrix(surv.train~.+0,data=x.train)
    y.test = subset(test_data, select = c(dfstim, dfsind))
    surv.test = Surv(y.test$dfstim,y.test$dfsind)
    x.test = subset(test_data, select = -c(dfstim, dfsind))
    x.test.m = model.matrix(surv.test ~.+0,data=x.test)
    
    lasso.cvfit = cv.glmnet(x.train.m,  surv.train, family = "cox",type.measure ="C", folds=folds,maxit=100000)
    
    lasso.glmnet.fit = glmnet(x.train.m, surv.train,family ="cox", lambda = lasso.cvfit$lambda.min)
    lasso.predict.train = predict(lasso.glmnet.fit,newx=x.train.m,type="link")[,1]
    lasso.predict.test = predict(lasso.glmnet.fit,newx=x.test.m,type="link")[,1]
    for(lass_pred_op in 1:length(lasso.predict.test)){
      lasso_var = lasso.predict.test[lass_pred_op]
      lasso_indx = as.numeric(names(lasso_var))
      lasso_pred_matrix[flag,lasso_indx] = as.numeric(lasso_var)[1]
    }
    lasso_train_con[flag] =  Cindex(lasso.predict.train,y= surv.train)
    lasso_test_con[flag] =  Cindex(lasso.predict.test,y= surv.test)
    lasso.feat.selected[flag,] = as.numeric(coef(lasso.glmnet.fit)!=0)
    
    ridge.cvfit = cv.glmnet(x.train.m,  surv.train, family = "cox",alpha=0,type.measure ="C", folds=folds,maxit=100000)
    
    ridge.glmnet.fit = glmnet(x.train.m, surv.train,family ="cox", alpha=0,lambda = ridge.cvfit$lambda.min)
    ridge.predict.train = predict(ridge.glmnet.fit,newx=x.train.m,type="link")[,1]
    ridge.predict.test = predict(ridge.glmnet.fit,newx=x.test.m,type="link")[,1] 
    for(ridg_pred_op in 1:length(ridge.predict.test)){
      ridge_var = ridge.predict.test[ridg_pred_op]
      ridge_indx = as.numeric(names(ridge_var))
      ridge_pred_matrix[flag,ridge_indx] = as.numeric(ridge_var)[1]
    }
    ridge_train_con[flag] =  Cindex(ridge.predict.train,y= surv.train)
    ridge_test_con[flag] =  Cindex(ridge.predict.test,y= surv.test)
    
    gbm_model = gbm(surv.train~., data=x.train, distribution="coxph")
    gbm_summary[[flag]] = summary(gbm_model)
    gbm_pred_train = predict(gbm_model, newdata=train_data, type="link")
    gbm_pred_test = predict(gbm_model, newdata=test_data, type="link")
    gbm.train.con[flag] = Cindex(gbm_pred_train, y=surv.train)
    gbm.test.con[flag] = Cindex(gbm_pred_test, y=surv.test)
    gbm_indxs = rownames(test_data)
    for(indx_pos in 1:length(gbm_indxs)){
      gbm.indx.pos = as.numeric(gbm_indxs[indx_pos])
      gbm_pred_matrix[flag,gbm.indx.pos] = gbm_pred_train[indx_pos]
    }
    
    temp_rsf_model = rfsrc(Surv(dfstim, dfsind) ~ .,data = train_data, importance = TRUE)
    temp_rsf_pred = predict.rfsrc(temp_rsf_model,newdata = test_data,type="lp")
    rsf_con_train[flag] = Cindex(temp_rsf_model$predicted.oob,y=surv.train)
    rsf_con_test[flag] = Cindex(temp_rsf_pred$predicted,y=surv.test)
    flag = flag + 1
  }
}

set.seed(6)
rfe_train_con = vector("integer",subset_columns)
rfe_test_con = vector("integer",subset_columns)
variables_picked = list()

model_under_consideration = data_subset
rfe_flag = 1
while(length(model_under_consideration)>2){
  shuffled_orig_dat = model_under_consideration[sample(1:nrow(model_under_consideration),nrow(model_under_consideration)), ]
  cvIndex1 = createFolds(factor(shuffled_orig_dat$dfsind), folds,returnTrain = T)
  orig_train_data = shuffled_orig_dat[cvIndex1[[1]],]
  orig_test_data = shuffled_orig_dat[-cvIndex1[[1]],]
  orig.surv.train = Surv(orig_train_data$dfstim,orig_train_data$dfsind)
  orig.surv.test = Surv(orig_test_data$dfstim,orig_test_data$dfsind)
  model_rfsrc=rfsrc(Surv(dfstim,dfsind) ~ .,data = orig_train_data, importance = TRUE)
  pred = predict.rfsrc(model_rfsrc,newdata = orig_test_data,type="lp")
  rfe_train_con[rfe_flag] = Cindex(model_rfsrc$predicted.oob, y = orig.surv.train)
  rfe_test_con[rfe_flag] = Cindex(pred$predicted, y = orig.surv.test)
  response_subset=subset(model_under_consideration, select = c(dfstim,dfsind))
  var_imp=vimp.rfsrc(model_rfsrc)$importance
  sorted_vars= sort(var_imp,decreasing=TRUE)
  selected_vars = head(sorted_vars,-1)
  variables_picked[[rfe_flag]] = toString(names(selected_vars))
  new_subset=subset(model_under_consideration,select=c(names(selected_vars)))
  new_new_subset=cbind(new_subset,response_subset)
  model_under_consideration=new_new_subset
  rfe_flag = rfe_flag + 1
}

round(median(unlist(final_stp_train_con)),4)
round(std.error(unlist(final_stp_train_con)),4)
round(median(unlist(final_stp_test_con)),4)
round(std.error(unlist(final_stp_test_con)),4)
round(median(lasso_train_con),4)
round(std.error(lasso_train_con),4)
round(median(lasso_test_con),4)
round(std.error(lasso_test_con),4)

round(median(ridge_train_con),4)
round(std.error(ridge_train_con),4)
round(median(ridge_test_con),4)
round(std.error(ridge_test_con),4)

round(median(gbm.train.con),4)
round(std.error(gbm.train.con),4)
round(median(gbm.test.con),4)
round(std.error(gbm.test.con),4)

round(median(unlist(rsf_con_train)),4)
round(std.error(unlist(rsf_con_train)),4)
round(median(unlist(rsf_con_test)),4)
round(std.error(unlist(rsf_con_test)),4)




# boxplot

boxplot(unlist(final_stp_train_con),lasso_train_con,
        ridge_train_con,gbm.train.con,unlist(rsf_con_train),
        names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Train Concordances")

boxplot(unlist(final_stp_test_con),lasso_test_con,
        ridge_test_con,gbm.test.con,unlist(rsf_con_test),
        names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Test Concordances")

boxplot(unlist(final_stp_train_con),unlist(final_stp_test_con), lasso_train_con, lasso_test_con,
        ridge_train_con,ridge_test_con,gbm.train.con,gbm.test.con, unlist(rsf_con_train),
        unlist(rsf_con_test), col= c(0,5,0,5,0,5,0,5,0,5) ,
        names=c("Stepwise","","Lasso","","Ridge","", "GBM","","RSF",""), 
        main="Train - test concordances")
legend("topright",c("train concordance","test concordance"),fill=c(0,5),cex=0.5)



## KM plots


stp_mean_pred = as.array(apply(stp_pred_matrx,2,mean))
stp_kmo = split.kms(stp_mean_pred, os=data_subset$dfstim, oevn=data_subset$dfsind)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(data_subset$dfstim, data_subset$dfsind)~stp_temp_var, data=data_subset)
ggsurvplot(stp_km.split,title="Progression Free Survival for Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


lasso_mean_pred = as.array(apply(lasso_pred_matrix,2,mean))
lasso_kmo = split.kms(lasso_mean_pred, os=data_subset$dfstim, oevn=data_subset$dfsind)
lasso_temp_var = lasso_kmo$iz
lasso_km.split = survfit(Surv(data_subset$dfstim, data_subset$dfsind)~lasso_temp_var, data=data_subset)
ggsurvplot(lasso_km.split,title="Progression Free Suvival Penalised Cox model(lasso)",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


gbm_mean_pred = as.array(apply(gbm_pred_matrix,2,mean))
gbm_kmo = split.kms(gbm_mean_pred, os=data_subset$dfstim, oevn=data_subset$dfsind)
gbm_temp_var = gbm_kmo$iz
gbm_km.split = survfit(Surv(data_subset$dfstim, data_subset$dfsind)~gbm_temp_var, data=data_subset)
ggsurvplot(gbm_km.split,title="Progression Free Survival for Boosting Cox model",pval = TRUE,surv.median.line = NULL, conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))

ridge_mean_pred = as.array(apply(ridge_pred_matrix,2,mean))
ridge_kmo = split.kms(ridge_mean_pred, os=data_subset$dfstim, oevn=data_subset$dfsind)
ridge_temp_var = ridge_kmo$iz
ridge_km.split = survfit(Surv(data_subset$dfstim, data_subset$dfsind)~ridge_temp_var, data=data_subset)
ggsurvplot(ridge_km.split,title="Progression Free Survival for Penalised Cox model(ridge)",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))




