rm(list=ls())
require(moonBook)
library(Boruta)
library(lme4)
library(pROC)
library(caret)

data <- read.csv("LMS_data_r.csv")
table <- mytable(malignancy~., data = data, show.all=TRUE)
table

ResVar = "malignancy"
Model_String = paste(ResVar, "~shape_R1+margin_R2+hemorrhage_R2+T2dark_R1+
                     cystic_R1+DR_R1+ADC_R1+RCR_R1")


Model = formula(Model_String)
lm_fit = lm(Model, data = data)
Boruta_fit = Boruta(Model, data = data)

summary(Boruta_fit)


Boruta_signif <- names(Boruta_fit$finalDecision[Boruta_fit$finalDecision %in% c("Confirmed", "Tentative")])  # collect Confirmed and Tentative variables
print(Boruta_signif)

Boruta_fit$finalDecision

plot(Boruta_fit, cex.axis=.7, las=2, xlab="", main="Variable Importance")

glm_fit = glm(Model, data = data, family = binomial)
summary(glm_fit)
ORtable(glm_fit)
reduced.model=step(glm_fit)
summary(reduced.model)
ORtable(reduced.model)
#par("mar")
#par(mar=c(1,1,1,1))

ORplot(glm_fit,type=1,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of All Variables") 

ORplot(reduced.model,type=2,show.OR=TRUE,show.CI=TRUE, 
       main="Odds Ratios of Selected Variables") 

ADC_roc <- roc(malignancy~ADC_R1,data=data)
RCR_roc <- roc(malignancy~RCR_R1,data=data)

plot.roc(ADC_roc, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)
plot.roc(RCR_roc, add=TRUE, col="blue", print.auc=TRUE, print.auc.adj=c(1.11,2.5),
         print.thres=TRUE,legacy.axes = TRUE)

legend("bottom", legend=c("ADC", "RCR"), col=c("red", "blue"), lwd=2) 

data2 = data[,c("malignancy","shape_R1","margin_R1","hemorrhage_R1","T2dark_R1",
                "cystic_R1","DR_R1","ADC_R1","RCR_R1")]

data2$feature_sum = rowSums(data2[,c(1,2,4,6)])
feature_roc <- roc(malignancy~feature_sum,data=data2)
plot.roc(feature_roc, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=c(1,2,3,4), print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

logistic_qual = glm(malignancy~feature_sum, data=data2, family = binomial())
summary(logistic_qual)


Log_odds=predict(logistic_qual, newdata= data2)
Prob = predict(logistic_qual, newdata= data2, type='response')


PREDICTED_C = ifelse(Prob > 0.5 , 1 , 0)
PREDICTED_C = factor(PREDICTED_C,levels = c(1,0))
data2$malignancy = factor(data2$malignancy, levels = c(1,0))
confusionMatrix(PREDICTED_C,data2$malignancy)
ROC = roc(data2$malignancy, Prob)
plot.roc(ROC, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)


logistic_qual2 = glm(malignancy~shape_R1+margin_R1+T2dark_R1+ADC_R1, data=data2, family = binomial())
summary(logistic_qual2)

Log_odds2=predict(logistic_qual2, newdata= data2)
Prob2 = predict(logistic_qual2, newdata= data2, type='response')

PREDICTED_2 = ifelse(Prob2 > 0.424 , 1 , 0)
PREDICTED_2 = factor(PREDICTED_2,levels = c(1,0))
data2$malignancy = factor(data2$malignancy, levels = c(1,0))
confusionMatrix(PREDICTED_2,data2$malignancy)

ROC2 = roc(data2$malignancy, Prob2)
plot.roc(ROC2, col="red", print.auc=TRUE, print.auc.adj=c(1.11,1.2),
         print.thres=TRUE, print.thres.adj=c(0.3,-1.0), legacy.axes = TRUE)

require(Epi)
require(ztable)
require(survival)

a1=ROC(form=malignancy~feature_sum,data=data2,plot="ROC")


glm_fit = glm(malignancy~shape_R1+T2dark_R1+DR_R1, data = data2, family = binomial)
summary(glm_fit)

ORtable=function(x,digits=4){
    suppressMessages(a<-confint(x))
    result=data.frame(exp(coef(x)),exp(a))
    result=round(result,digits)
    result=cbind(result,round(summary(x)$coefficient[,4],3))
    colnames(result)=c("OR","2.5%","97.5%","p")
    result
}

ORtable(glm_fit)

glm_fit = glm(malignancy~RCR_R1, data = data2, family = binomial)
summary(glm_fit)
ORtable(glm_fit)
