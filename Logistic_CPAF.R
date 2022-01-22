
rm(list=ls())
library("xlsx")
library(lme4)
library(Boruta)
library(moonBook)
library(car)


options(stringsAsFactors = FALSE)

CurDir = "C:/Users/hkkim/Desktop/"
FilePath = paste(CurDir, "CPAF_data_r.xlsx", sep="")

Data = read.xlsx(FilePath, 1, header=TRUE, encoding = "UTF-8")

ColTitle = colnames(Data)

DataMat = data.frame(Data)

ResVar = "Treatment"


Model_String = paste(ResVar, "~Origin_site+Age+Size+Aneurysm+DM+Dyslipidemia+Hypertension+CAD+Sex+Chest_pain+Smoking+Dyspnea+Palpitation")

#Model_String = paste(ResVar, "~Smoking+Sex+Size+Palpitation")

#Model_String = paste(ResVar, "~age")
Model = formula(Model_String)
# 
lm_fit = lm(Model, data = DataMat)

Boruta_fit = Boruta(Model, data = DataMat)

summary(Boruta_fit)


Boruta_signif <- names(Boruta_fit$finalDecision[Boruta_fit$finalDecision %in% c("Confirmed", "Tentative")])  # collect Confirmed and Tentative variables
print(Boruta_signif)

Boruta_fit$finalDecision

plot(Boruta_fit, cex.axis=.7, las=2, xlab="", main="Variable Importance")



############################# È¸±Í ¸ðµ¨ #############################################

glm_fit = glm(Model, data = DataMat, family = binomial)
#glm_fit = glm(Model, data = DataMat, family = poisson)
summary(glm_fit)


reduced.model=step(glm_fit)
summary(reduced.model)

ORplot(glm_fit,type=2,show.OR=TRUE,show.CI=TRUE,
       main="Odds Ratios of All Variables") 
ORplot(reduced.model,type=2,show.OR=TRUE,show.CI=TRUE, 
       main="Odds Ratios of Selected Variables") 


#ORtable=function(x,digits=2){
#  suppressMessages(a<-confint(x))
#  result=data.frame(exp(coef(x)),exp(a))
#  result=round(result,digits)
#  result=cbind(result,round(summary(x)$coefficient[,4],3))
#  colnames(result)=c("OR","2.5%","97.5%","p")
#  result
#}

#ORtable(reduced.model)
