library(plotrix)
library(glmnet)
library(gbm)
library(randomForestSRC)
library(survival)
library(caret)
library(survminer)

## data pre-processing

dat = read.csv("sarcoma.csv", stringsAsFactors = TRUE)
dim(dat)
dat = subset(dat, select = -c(X,surtim,surind) )
dat = na.omit(dat)
rownames(dat) <- 1:nrow(dat)

table(dat$dfsind)
dat$sex =( as.numeric(as.factor(dat$sex)))
dat$grade = (as.numeric(as.factor(dat$grade)))
dat$tumor.subtype =( as.numeric(as.factor(dat$tumor.subtype)))

data_subset = dat

# Function for plotting Kaplan Meier Curve

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

#Threshold_limit Calculation

repeats = 10
folds = 5

total = folds*repeats
threshold_limit = (30*total)/100
original_cols = ncol(dat)-2
subset_columns = ncol(data_subset)-2
subset_rows = nrow(data_subset)
subset_col_names = colnames(subset(data_subset, select = -c(dfstim, dfsind)))

stp_feat_selected = matrix(0,nrow=total,ncol=subset_columns)
colnames(stp_feat_selected) = subset_col_names
stp_pred_matrx = matrix(0,nrow=total,ncol=subset_rows)
full_step_model = vector("list",folds)

final_stp_train_con  = vector("list",total)
final_stp_test_con = vector("list",total)
final_stp_models = vector("list",total)

lasso_train_con = lasso_test_con = numeric(total)
lasso.feat.selected = matrix(0, nrow=total, ncol=subset_columns)
colnames(lasso.feat.selected) = subset_col_names
lasso_pred_matrix  = matrix(0, nrow=total, ncol=nrow(data_subset))

ridge_cvIndex = createFolds(factor(data_subset$dfsind), folds,returnTrain = T)
ridge_train_con = ridge_test_con =  numeric(total)
ridge.coef.matrix = matrix(0, nrow=total, ncol=subset_columns)
colnames(ridge.coef.matrix) = subset_col_names
ridge_pred_matrix  = matrix(0, nrow=total, ncol=nrow(data_subset))

gbm_summary = vector("list",total) 
gbm.train.con = gbm.test.con = vector("double",total) 
gbm.coef.matrix = matrix(0, nrow=total, ncol= subset_columns)
colnames(gbm.coef.matrix) = subset_col_names
gbm_pred_matrix = matrix(0,nrow=total,ncol=subset_rows)

rsf_model = rsf_pred = vector("list",total)
rsf_con_train = rsf_bs = vector("list",total)
rsf_con_test = vector("list",total)
rsf_con_test = vector("list",total)
rsf_feat_sel = matrix(0,nrow=total,ncol=subset_columns)
colnames(rsf_feat_sel) = subset_col_names

## Fitting all the models Using Cross Validation

n=nrow(data_subset)
counter = 1
for(a in 1:repeats){
  shuffled_dat = data_subset[sample(1:nrow(data_subset),nrow(data_subset)), ]
  cvIndex = createFolds(factor(shuffled_dat$dfsind), folds,returnTrain = T)
  for(b in 1:length(cvIndex)){
    training_data = data_subset[cvIndex[[b]],]
    test_dat = data_subset[-cvIndex[[b]],]
    
    start_cox = coxph(Surv(dfstim,dfsind) ~ 1, data = training_data)
    full_cox = coxph(Surv(dfstim,dfsind) ~ ., data = training_data)
    fit_step = step(start_cox, direction = "both", scope = full_cox$formula)
    stp_names = names(fit_step$coefficients)
    for(c in 1:length(stp_names)){
      stp_feat_selected[counter,stp_names[c]] = 1
    }
    full_form = fit_step$formula[3]
    full_form_1 = formula(paste("Surv(dfstim,dfsind)~",full_form))
    
    stp_cox_fit = coxph(full_form_1,data=training_data)
    stp_cox_pred = predict(stp_cox_fit,newdata = test_dat,type="lp")
    
    for(d in 1:length(stp_cox_pred)){
      stp_var = stp_cox_pred[d]
      stp_indx = as.numeric(names(stp_var))
      stp_pred_matrx[counter,stp_indx] = as.numeric(stp_var)[1]
    }
    
    final_stp_models[counter] <- full_form
    final_stp_train_con[counter]<-concordance(stp_cox_fit, newdata = training_data)$concordance
    final_stp_test_con[counter] <- concordance(stp_cox_fit, newdata = test_dat)$concordance
    
    x.train = subset(training_data, select = -c(dfstim, dfsind))
    y.train = subset(training_data, select = c(dfstim, dfsind))
    
    surv.train = Surv(y.train$dfstim,y.train$dfsind)
    x.train.m = model.matrix(surv.train~.+0,data=x.train)
    
    y.test = subset(test_dat, select = c(dfstim, dfsind))
    surv.test = Surv(y.test$dfstim,y.test$dfsind)
    x.test = subset(test_dat, select = -c(dfstim, dfsind))
    x.test.m = model.matrix(surv.test ~.+0,data=x.test)
    
    lasso.cvfit = cv.glmnet(x.train.m,  surv.train, family = "cox",type.measure ="C", folds=folds,maxit=100000)
    
    lasso.glmnet.fit = glmnet(x.train.m, surv.train,family ="cox", lambda = lasso.cvfit$lambda.min)
    lasso.predict.train = predict(lasso.glmnet.fit,newx=x.train.m,type="link")[,1]
    lasso.predict.test = predict(lasso.glmnet.fit,newx=x.test.m,type="link")[,1] # the fitted relative-risk for "cox";
    for(g in 1:length(lasso.predict.test)){
      lasso_var = lasso.predict.test[g]
      lasso_indx = as.numeric(names(lasso_var))
      lasso_pred_matrix[counter,lasso_indx] = as.numeric(lasso_var)[1]
    }
    lasso_train_con[counter] =  Cindex(lasso.predict.train,y= surv.train)
    lasso_test_con[counter] =  Cindex(lasso.predict.test,y= surv.test)
    lasso.feat.selected[counter,] = as.numeric(coef(lasso.glmnet.fit)!=0)
    
    ridge.cvfit = cv.glmnet(x.train.m,  surv.train, family = "cox",alpha=0,type.measure ="C", folds=folds,maxit=100000)
    
    ridge.glmnet.fit = glmnet(x.train.m, surv.train,family ="cox", alpha=0,lambda = ridge.cvfit$lambda.min)
    ridge.predict.train = predict(ridge.glmnet.fit,newx=x.train.m,type="link")[,1] # the fitted relative-risk for "cox";
    ridge.predict.test = predict(ridge.glmnet.fit,newx=x.test.m,type="link")[,1] # the fitted relative-risk for "cox";
    for(j in 1:length(ridge.predict.test)){
      ridge_var = ridge.predict.test[j]
      ridge_indx = as.numeric(names(ridge_var))
      ridge_pred_matrix[counter,ridge_indx] = as.numeric(ridge_var)[1]
    }
    ridge_train_con[counter] =  Cindex(ridge.predict.train,y= surv.train)
    ridge_test_con[counter] =  Cindex(ridge.predict.test,y= surv.test)
    
    ridge.coefs = coef(ridge.glmnet.fit)[,1]
    for(indx in 1:length(ridge.coefs)){
      ridge.coef.matrix[counter,indx] = abs(as.numeric(ridge.coefs[indx]))
    }
    
    gbm_model = gbm(surv.train~., data=x.train, distribution="coxph")
    gbm_summary[[counter]] = summary(gbm_model)
    gbm_pred_train = predict(gbm_model, newdata=training_data, type="link")
    gbm_pred_test = predict(gbm_model, newdata=test_dat, type="link")
    gbm.train.con[counter] = Cindex(gbm_pred_train, y=surv.train)
    gbm.test.con[counter] = Cindex(gbm_pred_test, y=surv.test)
    gbm_indxs = rownames(test_dat)
    for(indx_pos in 1:length(gbm_indxs)){
      gbm.indx.pos = as.numeric(gbm_indxs[indx_pos])
      gbm_pred_matrix[counter,gbm.indx.pos] = gbm_pred_train[indx_pos]
    }
    
    for(m in 1:total ){
      for(n in 1:subset_columns){
        gbm_index_var=which(subset_col_names==gbm_summary[[m]][n,]$var)
        gbm.coef.matrix[m,gbm_index_var]=gbm_summary[[m]][n,]$rel.inf
      }
    }
    
    temp_rsf_model = rfsrc(Surv(dfstim, dfsind) ~ .,data = training_data, importance = TRUE)
    temp_rsf_pred = predict.rfsrc(temp_rsf_model,newdata = test_dat,type="lp")
    rsf_model[[counter]] = temp_rsf_model
    rsf_pred[[counter]] = temp_rsf_pred
    rsf_con_train[counter] = Cindex(temp_rsf_model$predicted.oob,y=surv.train)
    rsf_con_test[counter] = Cindex(temp_rsf_pred$predicted,y=surv.test)
    rsf_var_imp = subsample(temp_rsf_model)
    rsf_feats = head(subset_col_names[order(temp_rsf_model$importance, decreasing=TRUE)],10)
    for(q in 1:length(rsf_feats)){
      rsf_feat_sel[counter,rsf_feats[q]] = 1
    }
    
    counter = counter + 1
  }
}


# Random Forest survival(RFS) + RFE
rfe_train_con = vector("integer",)
rfe_test_con = vector("integer",subset_columns)
variables_picked = list()
subset_columns
model_under_consideration = data_subset
rfe_counter = 1
while(length(model_under_consideration)>2){
  shuffled_orig_dat = model_under_consideration[sample(1:nrow(model_under_consideration),nrow(model_under_consideration)), ]
  cvIndex1 = createFolds(factor(shuffled_orig_dat$dfsind), folds,returnTrain = T)
  orig_training_data = shuffled_orig_dat[cvIndex1[[b]],]
  orig_test_dat = shuffled_orig_dat[-cvIndex1[[b]],]
  orig.surv.train = Surv(orig_training_data$dfstim,orig_training_data$dfsind)
  orig.surv.test = Surv(orig_test_dat$dfstim,orig_test_dat$dfsind)
  model_rfsrc=rfsrc(Surv(dfstim,dfsind) ~ .,data = orig_training_data, importance = TRUE)
  pred = predict.rfsrc(model_rfsrc,newdata = orig_test_dat,type="lp")
  rfe_train_con[rfe_counter] = Cindex(model_rfsrc$predicted.oob, y = orig.surv.train)
  rfe_test_con[rfe_counter] = Cindex(pred$predicted, y = orig.surv.test)
  response_subset=subset(model_under_consideration, select = c(dfstim,dfsind))
  var_imp=vimp.rfsrc(model_rfsrc)$importance
  sorted_vars= sort(var_imp,decreasing=TRUE)
  selected_vars = head(sorted_vars,-1)
  variables_picked[[rfe_counter]] = toString(names(selected_vars))
  new_subset=subset(model_under_consideration,select=c(names(selected_vars)))
  new_new_subset=cbind(new_subset,response_subset)
  model_under_consideration=new_new_subset
  rfe_counter = rfe_counter + 1
}


## Finding the features which are important


variables_picked = unlist(variables_picked)
indx_with_max_train_con = which.max(rfe_train_con)
indx_with_max_test_con = which.max(rfe_test_con)
rfe_feat_sel = variables_picked[indx_with_max_test_con]

sw_feat_sel = which(colSums(stp_feat_selected)>threshold_limit) == TRUE

lasso_feat_sel = which(colSums(lasso.feat.selected)>threshold_limit) == TRUE
rfsrc_feat_sel = which(colSums(rsf_feat_sel)>threshold_limit) == TRUE
ridge_coef_sum = colSums(ridge.coef.matrix)
ridge_feat_sel = names(head(ridge_coef_sum[order(-ridge_coef_sum)],6))

gbm_coef_sum = colSums(gbm.coef.matrix)
gbm_feat_sel = names(head(gbm_coef_sum[order(-gbm_coef_sum)],6))


## finding the performance Metric


round(median(ridge_train_con),4)
round(std.error(ridge_train_con),4)
round(median(ridge_test_con),4)
round(std.error(ridge_test_con),4)

round(median(rfe_train_con),4)
round(std.error(rfe_train_con),4)
round(median(rfe_test_con),4)
round(std.error(rfe_test_con),4)

round(median(gbm.train.con),4)
round(std.error(gbm.train.con),4)
round(median(gbm.test.con),4)
round(std.error(gbm.test.con),4)

round(median(unlist(rsf_con_train)),4)
round(std.error(unlist(rsf_con_train)),4)
round(median(unlist(rsf_con_test)),4)
round(std.error(unlist(rsf_con_test)),4)



round(median(unlist(final_stp_train_con)),4)
round(std.error(unlist(final_stp_train_con)),4)
round(median(unlist(final_stp_test_con)),4)
round(std.error(unlist(final_stp_test_con)),4)
round(median(lasso_train_con),4)
round(std.error(lasso_train_con),4)
round(median(lasso_test_con),4)
round(std.error(lasso_test_con),4)


## Kaplan Meier plots

step_mean_predict = as.array(apply(stp_pred_matrx,2,mean))
step_kmo = split.kms(step_mean_predict, os=data_subset$dfstim, oevn=data_subset$dfsind)
step_tmp_var = step_kmo$iz
stp_km.split = survfit(Surv(data_subset$dfstim, data_subset$dfsind)~step_tmp_var, data=data_subset)
ggsurvplot(stp_km.split,title="Progression Free Survival for Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


gbm_mean_predict = as.array(apply(gbm_pred_matrix,2,mean))
gbm_kmo = split.kms(gbm_mean_predict, os=data_subset$dfstim, oevn=data_subset$dfsind)
gbm_tmp_var = gbm_kmo$iz
gbm_km.split = survfit(Surv(data_subset$dfstim, data_subset$dfsind)~gbm_tmp_var, data=data_subset)
ggsurvplot(gbm_km.split,title="Progression Free Survival for Boosting Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


lasso_mean_predict = as.array(apply(lasso_pred_matrix,2,mean))
lasso_kmo = split.kms(lasso_mean_predict, os=data_subset$dfstim, oevn=data_subset$dfsind)
lasso_tmp_var = lasso_kmo$iz
lasso_km.split = survfit(Surv(data_subset$dfstim, data_subset$dfsind)~lasso_tmp_var, data=data_subset)
ggsurvplot(lasso_km.split,title="Progression Free Survival for Penalised Cox model(lasso)",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))



##boxplots without pre-filtering

boxplot(unlist(final_stp_train_con),unlist(final_stp_test_con), lasso_train_con, lasso_test_con,
        ridge_train_con,ridge_test_con,gbm.train.con,gbm.test.con, unlist(rsf_con_train),
        unlist(rsf_con_test),rfe_train_con,rfe_test_con ,col= c(0,5,0,5,0,5,0,5,0,5,0,5) ,
        names=c("Stepwise","","Lasso","","Ridge","", "GBM","","RSF","","RFE",""), 
        main="Train - test concordances")
legend("topright",c("train concordance","test concordance"),fill=c(0,5),cex=0.4)

