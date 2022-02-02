rm(list=ls())
require(moonBook)
library(Boruta)
library(lme4)
library(pROC)
library(caret)

data <- read.csv("LMS_data_r.csv")


######## R2 ########

ResVar = "malignancy"
Model = as.formula(paste(ResVar, "~shape_R2+margin_R2+hemorrhage_R2+T2dark_R2+
                     necrosis_R2+ADC_R2+RCR_R2"))

## cross table analysis
table <- mytable(Model,data=data)
table

## Boruta
lm_fit = lm(Model, data = data)
Boruta_fit = Boruta(Model, data = data)
summary(Boruta_fit)
Boruta_signif <- names(Boruta_fit$finalDecision[Boruta_fit$finalDecision %in% c("Confirmed", "Tentative")])  # collect Confirmed and Tentative variables
print(Boruta_signif)
Boruta_fit$finalDecision
plot(Boruta_fit, cex.axis=.7, las=2, xlab="", main="Variable Importance")


## logistic regression model with OR calculation for qualitative features

ORtable=function(x,digits=4){
    suppressMessages(a<-confint(x))
    result=data.frame(exp(coef(x)),exp(a))
    result=round(result,digits)
    result=cbind(result,round(summary(x)$coefficient[,4],3))
    colnames(result)=c("OR","2.5%","97.5%","p")
    result
}

feature_qual = as.formula(paste(ResVar, "~shape_R2+margin_R2+T2dark_R2+
                     necrosis_R2"))
logistic_qual = glm(feature_qual, data = data, family = binomial)
summary(logistic_qual)
ORtable(logistic_qual)
ORplot(logistic_qual,type=2,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of Qualitative Features") 

## ROC curve for LR model

Log_odds=predict(logistic_qual, newdata= data)
Prob = predict(logistic_qual, newdata= data, type='response')

PREDICTED_C = ifelse(Prob > 0.5 , 1 , 0)
PREDICTED_C = factor(PREDICTED_C,levels = c(0,1))
data$malignancy = factor(data$malignancy, levels = c(0,1))
confusionMatrix(PREDICTED_C,data$malignancy)
ROC = roc(data$malignancy, Prob)
plot.roc(ROC, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)


## ROC curve for ADC only

ROC_adc = roc(data$malignancy, data$ADC_R2)
plot.roc(ROC_adc, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

## ROC curve for RCR only

ROC_rcr = roc(data$malignancy, data$RCR_R2)
plot.roc(ROC_rcr, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

## logistic regression model for qualitative + ADC

feature_adc = as.formula(paste(ResVar, "~shape_R2+margin_R2+T2dark_R2+
                     necrosis_R2+ADC_R2"))
logistic_adc = glm(feature_adc, data = data, family = binomial)
summary(logistic_adc)
ORtable(logistic_adc)
ORplot(logistic_adc,type=2,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of Qualitative Features + ADC") 

## ROC curve of LR model: qualitative + ADC

Log_odds=predict(logistic_adc, newdata= data)
Prob = predict(logistic_adc, newdata= data, type='response')

PREDICTED_C = ifelse(Prob > 0.5 , 1 , 0)
PREDICTED_C = factor(PREDICTED_C,levels = c(0,1))
confusionMatrix(PREDICTED_C,data$malignancy)
ROC_ADC = roc(data$malignancy, Prob)
plot.roc(ROC_ADC, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

## logistic regression model for qualitative + RCR

feature_rcr = as.formula(paste(ResVar, "~shape_R2+margin_R2+T2dark_R2+
                     necrosis_R2+RCR_R2"))
logistic_rcr = glm(feature_rcr, data = data, family = binomial)
summary(logistic_rcr)
ORtable(logistic_rcr)
ORplot(logistic_rcr,type=2,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of Qualitative Features + RCR") 

## ROC curve of LR model: qualitative + RCR

Log_odds=predict(logistic_rcr, newdata= data)
Prob = predict(logistic_rcr, newdata= data, type='response')

PREDICTED_C = ifelse(Prob > 0.5 , 1 , 0)
PREDICTED_C = factor(PREDICTED_C,levels = c(0,1))
confusionMatrix(PREDICTED_C,data$malignancy)
ROC_RCR = roc(data$malignancy, Prob)
plot.roc(ROC_RCR, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)





######## K ########

ResVar = "malignancy"
Model = as.formula(paste(ResVar, "~shape_K+margin_K+hemorrhage_K+T2dark_K+
                     necrosis_K+ADC_K+RCR_K"))

## cross table analysis
table <- mytable(Model,data=data)
table

## Boruta
lm_fit = lm(Model, data = data)
Boruta_fit = Boruta(Model, data = data)
summary(Boruta_fit)
Boruta_signif <- names(Boruta_fit$finalDecision[Boruta_fit$finalDecision %in% c("Confirmed", "Tentative")])  # collect Confirmed and Tentative variables
print(Boruta_signif)
Boruta_fit$finalDecision
plot(Boruta_fit, cex.axis=.7, las=2, xlab="", main="Variable Importance")


## logistic regression model with OR calculation for qualitative features

ORtable=function(x,digits=4){
    suppressMessages(a<-confint(x))
    result=data.frame(exp(coef(x)),exp(a))
    result=round(result,digits)
    result=cbind(result,round(summary(x)$coefficient[,4],3))
    colnames(result)=c("OR","2.5%","97.5%","p")
    result
}

feature_qual = as.formula(paste(ResVar, "~shape_K+margin_K+T2dark_K+
                     necrosis_K"))
logistic_qual = glm(feature_qual, data = data, family = binomial)
summary(logistic_qual)
ORtable(logistic_qual)
ORplot(logistic_qual,type=2,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of Qualitative Features") 

## ROC curve for LR model

Log_odds=predict(logistic_qual, newdata= data)
Prob = predict(logistic_qual, newdata= data, type='response')

PREDICTED_C = ifelse(Prob > 0.5 , 1 , 0)
PREDICTED_C = factor(PREDICTED_C,levels = c(0,1))
data$malignancy = factor(data$malignancy, levels = c(0,1))
confusionMatrix(PREDICTED_C,data$malignancy)
ROC = roc(data$malignancy, Prob)
plot.roc(ROC, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)


## ROC curve for ADC only

ROC_adc = roc(data$malignancy, data$ADC_K)
plot.roc(ROC_adc, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

## ROC curve for RCR only

ROC_rcr = roc(data$malignancy, data$RCR_K)
plot.roc(ROC_rcr, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

## logistic regression model for qualitative + ADC

feature_adc = as.formula(paste(ResVar, "~shape_K+margin_K+T2dark_K+
                     necrosis_K+ADC_K"))
logistic_adc = glm(feature_adc, data = data, family = binomial)
summary(logistic_adc)
ORtable(logistic_adc)
ORplot(logistic_adc,type=2,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of Qualitative Features + ADC") 

## ROC curve of LR model: qualitative + ADC

Log_odds=predict(logistic_adc, newdata= data)
Prob = predict(logistic_adc, newdata= data, type='response')

PREDICTED_C = ifelse(Prob > 0.5 , 1 , 0)
PREDICTED_C = factor(PREDICTED_C,levels = c(0,1))
confusionMatrix(PREDICTED_C,data$malignancy)
ROC_ADC = roc(data$malignancy, Prob)
plot.roc(ROC_ADC, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

## logistic regression model for qualitative + RCR

feature_rcr = as.formula(paste(ResVar, "~shape_K+margin_K+T2dark_K+
                     necrosis_K+RCR_K"))
logistic_rcr = glm(feature_rcr, data = data, family = binomial)
summary(logistic_rcr)
ORtable(logistic_rcr)
ORplot(logistic_rcr,type=2,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of Qualitative Features + RCR") 

## ROC curve of LR model: qualitative + RCR

Log_odds=predict(logistic_rcr, newdata= data)
Prob = predict(logistic_rcr, newdata= data, type='response')

PREDICTED_C = ifelse(Prob > 0.5 , 1 , 0)
PREDICTED_C = factor(PREDICTED_C,levels = c(0,1))
confusionMatrix(PREDICTED_C,data$malignancy)
ROC_RCR = roc(data$malignancy, Prob)
plot.roc(ROC_RCR, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

roc.test(ROC,ROC_RCR)
