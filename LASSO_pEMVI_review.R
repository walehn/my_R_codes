rm(list=ls())
library(dplyr)
library(heatmaply)
library(glmnet)
library(mRMRe)
library(caret)
library(tidyverse)
library(pROC)
library(ggplot2)
library(reshape2)
library(ggpubr)

## hyper parameter
f_count = 15
ratio = 0.6

setwd("~/DL/pyradiomics")
df <- read.csv("./results/results_v2.csv", header = T)
clinical <- read.csv("./results/radiomics_list_v4.csv", header = T)
df2 <-cbind(df, clinical$pEMVI_review)
df2 <- dplyr::rename(df2, pEMVI_review = 'clinical$pEMVI_review')
df3 <- na.omit(df2)

## remove environment variables
df3 <- df3[, -c(1:39)]

## sort by EMVI and numeric conversion
df3 <- df3[order(-df3$pEMVI_review),]
df3$pEMVI_review <- as.numeric(df3$pEMVI_review)
#df3$rEMVI <- as.numeric(df3$rEMVI)

## MRMR feature selection
radiomics_features <- heatmaply::normalize(df3[-1133])
radiomics_features <- cbind(radiomics_features, df3$pEMVI_review)
radiomics_features <- dplyr::rename(radiomics_features, pEMVI_review = 'df3$pEMVI_review')

dd <- mRMR.data(data = radiomics_features)
fs <- mRMR.ensemble(data=dd,target_indices=1133, feature_count = f_count, solution_count = 1)
sel_vals <- print(apply(solutions(fs)[[1]], 2, function(x, y) { return(y[x]) }, y=featureNames(dd)))

## drawing heatmap with selected variables

annotation <- radiomics_features["pEMVI_review"]
temp_list <- list()
temp_list <- sel_vals
temp <-  df3[,temp_list]
sel_fs <- temp
sel_fs <- heatmaply::normalize(sel_fs)
# 
# heatmaply(heatmaply::normalize(sel_fs),
#           # Rowv = row_dend,
#           # Colv = col_dend,
#           seriate = "OLO",
#           plot_method = "plotly",
#           dendrogram = "column",
#           row_side_colors = data.frame("pEMVI_review" = annotation, check.names = FALSE),
#           distfun = "spearman",
#           margins = c(60,130,10,10)) %>%
# 
#   colorbar(
#     tickfont = list(size = 10),
#     titlefont = list(size = 10,family = "Times"), which = 1) %>%
# 
#   colorbar(
#     tickfont = list(size = 10),
#     titlefont = list(size = 10), which = 2) %>%
# 
#   layout(
#     font = list(
#       family = "Times New Roman",
#       font ='italic')
#   )
# 
# ## correlation heatmap
# 
# r <- cor(heatmaply::normalize(sel_fs), method = "spearman")
# 
# cor.test.p <- function(x){
#   FUN <- function(x, y) cor.test(x, y)[["p.value"]]
#   z <- outer(
#     colnames(x),
#     colnames(x),
#     Vectorize(function(i,j) FUN(x[,i], x[,j]))
#   )
#   dimnames(z) <- list(colnames(x), colnames(x))
#   z
# }
# p <- cor.test.p(heatmaply::normalize(sel_fs))
# 
# heatmaply_cor(
#   r,
#   node_type = "scatter",
#   point_size_mat = -log10(p),
#   point_size_name = "-log10(p-value)",
#   label_names = c("x", "y", "Correlation")
# )
# 
# heatmaply_cor(
#   r,
#   xlab = "Features",
#   ylab = "Features",
#   k_col = 2,
#   k_row = 2
# )

## LASSO

# divide train/val set

sel_fs <- cbind(sel_fs, radiomics_features$pEMVI_review)
sel_fs <- dplyr::rename(sel_fs, pEMVI_review = 'radiomics_features$pEMVI_review')
set.seed(5)
inTrain <- createDataPartition(y=sel_fs$pEMVI_review, p=ratio, list=F)
sel_fs.train <- sel_fs[inTrain,]
sel_fs.test <- sel_fs[-inTrain,]


####### train / test data set difference estimation

t.test_p.value_df <- data.frame() # blank data.frame for saving

for (i in 1:(length(sel_fs)-1)){
  t.test_p.value <- t.test(sel_fs.train[i],sel_fs.test[i])$p.value
  t.test_p.value_df[i,1] <- names(sel_fs)[i]
  t.test_p.value_df[i,2] <- t.test_p.value
}
colnames(t.test_p.value_df) <- c("x_var_name", "p.value")
arrange(t.test_p.value_df, p.value)

aaa <- table(sel_fs.train$pEMVI_review)
bbb <- table(sel_fs.test$pEMVI_review)
table_1 <- rbind(aaa, bbb)
chisq.test(as.table(table_1))$p.value

# reconstruct data frame 
sel_fs.train_marking <- cbind(sel_fs.train,rep("train",nrow(sel_fs.train)))
sel_fs.train_marking <- rename(sel_fs.train_marking, dataset = 'rep("train", nrow(sel_fs.train))')
sel_fs.test_marking <- cbind(sel_fs.test,rep("test",nrow(sel_fs.test)))
sel_fs.test_marking <- rename(sel_fs.test_marking, dataset = 'rep("test", nrow(sel_fs.test))')
sel_fs_combine <- rbind(sel_fs.train_marking,sel_fs.test_marking)

sel_fs_melt <- melt(sel_fs_combine[-16])
sel_fs_melt_log <- transform(sel_fs_melt)
# visualization
ggplot(data=sel_fs_melt_log,aes(x=variable, y=value))+
  geom_boxplot(aes(fill=dataset))+coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(method = "t.test", ref.group = ".all.",label.x = 0, label.y = -0.2)


###### LASSO start

x <- model.matrix(~.,sel_fs.train[-(f_count+1)])[,-1]
y <- sel_fs.train$pEMVI_review

set.seed(1575)
# train = sample(1:nrow(x),nrow(x)/2)
# test = (-train)
# ytest = y[test]

# cv.lasso <- cv.glmnet(x[train,], y[train], alpha=1, family="binomial")
cv.lasso <- cv.glmnet(x, y, alpha=1, family="binomial", nfolds=10)
lasso.coef = predict(cv.lasso, type = "coefficients", s=cv.lasso$lambda.1se)


#plot
plot(cv.lasso)
plot(cv.lasso$glmnet.fit, xvar = "lambda", label = TRUE)
plot(cv.lasso$glmnet.fit, xvar = "norm", label = TRUE)

a = cv.lasso$lambda.min
b = cv.lasso$lambda.1se
c = coef(cv.lasso, s= cv.lasso$lambda.1se)
c

#coefficient plot
Name <- c@Dimnames[[1]]
Value <- c@x
coef_df <- data.frame(Name[c@i+1], Value)
coef_df <- coef_df[-1,]
coe <- coef_df
coef_df$Name <- as.factor(coef_df$Name)
g <- ggplot(coef_df,aes(reorder(Name,-Value),Value))
g + geom_col(fill="blue",position = position_stack(reverse = TRUE))+
  coord_flip()+theme_classic()+xlab("Feature") + ylab("Coefficient Value")


##factorize
sel_fs.train$pEMVI_review <-factor(sel_fs.train$pEMVI_review, levels = c(1,0), labels = c("Pos", "Neg"))
sel_fs.test$pEMVI_review <-factor(sel_fs.test$pEMVI_review, levels = c(1,0), labels = c("Pos", "Neg"))

####### Make prediction on training data


x.train <- model.matrix(~.,sel_fs.train[-(f_count+1)])[,-1]
lasso.prediction_train = predict(cv.lasso, s=cv.lasso$lambda.1se, newx=x.train, type="response")
predicted.classes_train <- ifelse(lasso.prediction_train > 0.406, 1, 0)
predicted.classes_train <-factor(predicted.classes_train, levels = c(1,0), labels = c("Pos", "Neg"))

cf_predict = confusionMatrix(data = predicted.classes_train, reference = sel_fs.train$pEMVI_review, positive = "Pos")

draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Class1', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Class2', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Class1', cex=1.2, srt=90)
  text(140, 335, 'Class2', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

draw_confusion_matrix(cf_predict)
#ROC curve
par(pty="m", cex=1)
result_train.roc <- roc(sel_fs.train$pEMVI_review,lasso.prediction_train)
plot(result_train.roc, col = "orangered", legacy.axes=T,
     print.auc=T, print.auc.adj=c(0.9, -13), #print.thres.adj=c(0.3,-1.0),
     max.auc.polygon=T, #print.thres="best", print.thres.pch=19, print.thres.col="black",
     #print.thres.best.method="closest.topleft",
     auc.polygon=T, auc.polygon.col="lightpink")


####### Make prediction on validation data

x.test <- model.matrix(~.,sel_fs.test[-(f_count+1)])[,-1]
lasso.prediction = predict(cv.lasso, s=cv.lasso$lambda.1se, newx=x.test, type="response")
predicted.classes <- ifelse(lasso.prediction > 0.314, 1, 0)
predicted.classes <-factor(predicted.classes, levels = c(1,0), labels = c("Pos", "Neg"))

cm_val <- confusionMatrix(data = predicted.classes, reference = sel_fs.test$pEMVI_review, positive = "Pos")
draw_confusion_matrix(cm_val)
# ROC curve
result.roc <- roc(sel_fs.test$pEMVI_review,lasso.prediction)
plot(add=TRUE, result.roc, #print.thres="best", print.thres.pch=19, print.thres.best.method="closest.topleft",
     print.thres.col="blue", auc.polygon=T, auc.polygon.col="slategray2",
     col="blue", print.auc=T, print.auc.adj=c(0.9, 3))

legend("bottomright", legend=c("training set", "validation set"), col=c("red", "blue"), lwd=2)  


##rad logistic regression
rad <- clinical[c(3,5)]
rad.train <- na.omit(rad[inTrain,])
rad.test <- na.omit(rad[-inTrain,])

fit.train <- glm(rad.train$pEMVI_review~ rEMVI,data = rad.train, family = "binomial")
summary(fit.train)
anova(fit.train,test="Chisq")
pred_rad <- predict(fit.train, newdata=rad.test,type="response")

rad.roc <- roc(rad.test$pEMVI_review,rad.test$rEMVI)
plot(rad.roc, add=TRUE, print.thres="best", print.thres.best.method="closest.topleft",
     col="blue", print.auc=T, print.auc.adj=c(1.11,1.2))
print(result.coords)#to get threshold and accuracy
rad.roc[["auc"]]

rad.roc2 <- roc(rad.train$pEMVI_review,rad.train$rEMVI)
plot(rad.roc2, add=TRUE, print.thres="best", print.thres.best.method="closest.topleft",
     col="black", print.auc=T, print.auc.adj=c(1.11,1.2))
print(result.coords)#to get threshold and accuracy
rad.roc2[["auc"]]

confusionMatrix(data = factor(rad.test$rEMVI), reference = factor(rad.test$pEMVI_review), positive = "1")
confusionMatrix(data = factor(rad.train$rEMVI), reference = factor(rad.train$pEMVI_review), positive = "1")
confusionMatrix(data = factor(rad$rEMVI), reference = factor(rad$pEMVI_review), positive = "1")

roc.test(result.roc, rad.roc)

# #### Logistic regression without LASSO
# full.model <- glm(rEMVI~., sel_fs.train, family = binomial)
# # Make predictions
# probabilities <- full.model %>% predict(sel_fs.test[-16], type = "response")
# predicted.classes_lr <- ifelse(probabilities > 0.024, 1, 0)
# predicted.classes_lr <-factor(predicted.classes_lr, levels = c(1,0), labels = c("Pos", "Neg"))
# 
# confusionMatrix(data = predicted.classes_lr, reference = sel_fs.test$rEMVI, positive = "Pos")
# 
# # ROC curve
# result_lr.roc <- roc(sel_fs.test$rEMVI, probabilities)
# plot(result_lr.roc, print.thres='best', print.thres.best.method="closest.topleft", 
#      col="red", print.auc=T)
# result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
# print(result.coords)#to get threshold and accuracy