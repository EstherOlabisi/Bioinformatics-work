#####1. Differences in expression between B and T immunophenotype----
#After studying Karlee's paper again and comparing it to this project, I realized that we never needed a KW plus post-hoc test. Karlee's paper was testing miRNA expression between three groups: healthy controls, Bcells, and T cells, hence the need for KW. It's also why they needed a post-hoc to identify which of the three groups differed per comparison. No wonder I remained stuck on this the whole time because I was trying to apply KW + post-hoc to this dataset, which was not applicable. Two-sample Kruskal Wallis test alone would work too because it’s basically the same as Mann-Whitney(MW). However, I chose the non-parametric MW test over KW since the former is designed for independent two sample testing. Following the testing steps for all 27 miRNAs, I corrected the p-values for multiple testing error using the Benjamin and Hochberg method (Yoav Benjamini & Hochberg, 1995). Then, I selected the significant miRNAs (n=17) with p < 0.05. 


#install_github("KKPMW/matrixTests")
library(tidyverse)
library(matrixTests)

#Read processed ct values into a dataframe. Make miRNA ids rownames. Select columns of ct values then transform the dataframe so that miRNA ids become column names and individuals become rows
ct_data <- read.csv("Processed_data.csv", check.names = F)
rownames(ct_data) <- ct_data$`Mature ID`
ct_data <- ct_data[c(4:39)]
ct_data <- as.data.frame(t(ct_data))

#Separate into B and T dataframes
B.df <- ct_data[1:22,]
T.df <- ct_data[23:36,]

#Run Mann Whitney U test by column (1 column: 1 miRNA )
results <- as.data.frame(
  col_wilcoxon_twosample(x = B.df, y = T.df,exact = TRUE))

#This step adjusts for multiple testing errors using the Benjamin and Hochberg correction method (Yoav Benjamini & Hochberg, 1995). Add adjusted p-values to the results df. 
results <- results %>%
  mutate(pvalue.adj = p.adjust(pvalue, method = "BH")) %>%
  relocate(pvalue.adj, .after = pvalue)

#Select significant results with p < 0.05. 
sig.res <- results[results$pvalue.adj < 0.05, 5:6]  #Final result. 

#Export sig res
write.csv(sig.res, file = "BT_updated_sigdiff.csv", row.names = T)



#References
#Yoav Benjamini, & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x 




####2. Survival Analysis Curves ----
library(survival)
library(survminer)
library(ggfortify)
library(ggpubr)

#Read pfs and os files as dataframes
pfs <- read.csv("Updated PFS.csv")
os <- read.csv("Updated OS.csv")

#Subset portions of dataframes (df) needed. Note: the PFS or OS column must remain as the first column of any df used.

#PFS
pfs.T <- pfs[1:14, ]
pfs.B <- pfs[16:37, ]

#OS
os.T <- os[1:14, ]
os.B <- os[16:37, ]

#Change survival input dataframe and mirna column range with the following variables.
#View for pfs.B mirnas for example
survdf <- pfs.B #change df here
range.col <- 6:27 #change mirna column range here

##Exploratory----
#You can run these lines to view all survival plots for all mirnas regardless of significance. Each plot will show on the viewer pane
surv.fits <- lapply(survdf[range.col], function(x) surv_fit(Surv(PFS.days, Censoring) ~ x, data = survdf))
surv.plots <- ggsurvplot_list(surv.fits, data = pfs.B, pval = T)
surv.plots
arrange_ggsurvplots(surv.plots[1:8], ncol = 4, nrow = 2)  #First 8 plots at once  


##Main analysis----
#---- function for plotting curves----
#Create a function for plotting km fits. Median values and pvalues are annotated onto the plot.

###ARGUMENTS
#surv.obj = fitted survival object
#df = full dataframe containing survival and censoring data
#legend.label = Preferred string for legend title e.g "B cell "or "T cell"
#xlabel and ylabel = labels for x and y axes
kmplot <- function(surv.obj, df, legend.label = "", xlabel = "Days", ylabel = "Percent survival") {
  
  xmax <- max(surv.obj$time)  #max survival time for x-axis range
  main.plot <- ggsurvplot(surv.obj,
                          data = df,
                          surv.scale = "percent",
                          size = 1.2,
                          pval = F,
                          conf.int = F,
                          xlab = xlabel,
                          xlim = c(0, xmax),
                          ylab = ylabel,
                          legend.title = legend.label,
                          ggtheme = theme_classic()) # Change ggplot2 theme
  
  pvalue <- surv_pvalue(surv.obj) #extract pvalues and round
  pvalue <- signif(pvalue$pval, digits = 3)
  median <- surv_median(surv.obj, combine = T) #extract median values for annotation
  high.m <- median$median[1] 
  low.m <- median$median[2]
  
  #paste0() concatenates each median value with the letter 'd' for days
  main.plot$plot+
    theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          axis.title = element_text(size = 8))+
    annotate("text",
             x = xmax-100,
             y = 1,
             label = paste0("p = ", pvalue))+
    annotate("text", 
             x = high.m,
             y = 0.4,
             label = paste0(high.m, "d"),
             size = 3)+
    annotate("text",
             x = low.m,
             y = 0.6,
             label = paste0(low.m, "d"),
             size = 3)
}



#---- function for fitting survival data ----
#This function extracts the pvalues of statistically significant miRNAs and plots those miRNA survival data only. The surv_fit() function is used instead of survfit() since the former supports surv lists. As.formula() was used with miRNA names to fit the final plots since we want those names to be displayed on the plots.

###ARGUMENTS
#df = full df including censoring data
#col.range = range of columns to be fitted (miRNAs only)
#Plot customization arguments applicable to km.plot()
plot.sig.fits <- function(df, col.range, legend.label = "", xlabel = "Days", ylabel = "Percent survival"){
  subset.df <- df[col.range] 
  time.col <- df[,1]  #PFS or OS time column
  
  #fit all miRNA survival data
  all.fits <- lapply(subset.df, function(x) surv_fit(Surv(time.col, Censoring) ~ x, data = df))
  pvalues <- surv_pvalue(all.fits, combine = T)
  pvalues$pval <- signif(pvalues$pval, digits = 3)

  #Obtain statistically significant miRNAs ids for fitting final plots.
  sig.mirna <- pvalues[pvalues$pval < 0.05,] 
  sig.ids <- sig.mirna$id
  num.sig <- length(sig.ids)
  
  if (num.sig >= 1) {
    print(paste("There are", num.sig, "significant miRNAs:", 
                paste(sig.ids, collapse = ", ")))
    
    sig.formulae <- lapply(sig.ids, function(x) as.formula(
      paste0("Surv(time.col, Censoring) ~ ", x) ))
    sig.fits <-  surv_fit(sig.formulae, data = df)
    
    #plot the fits for significant miRNAs 
    plot.list <- lapply(sig.fits,
                        function(x) kmplot(surv.obj = x,
                                           df = df,
                                           legend.label,
                                           xlabel,
                                           ylabel))
    
    #create options for the output
    out.list <- list(
      all.pvalues = pvalues,
      plot.list = plot.list)
    return(out.list)
    
  } else if (num.sig == 0) {
    cat("There are", num.sig, "significant miRNAs. These are all the p-values: \n\n")
    print(pvalues)
  }
  
}


#---- Results ----
#One can select pvalues or plots only for the output
#PFS for B cell
pfs.B.km <- plot.sig.fits(pfs.B, 4:25, legend.label = "B cell:  ") 
pfs.B.km
ggarrange(plotlist = pfs.B.km$plot.list, labels = "AUTO") #all PFS.B plots at once
pfs.B.km$all.pvalues

#PFS for T cell - no signif miRNAs
pfs.T.km <- plot.sig.fits(pfs.T, 3:25, legend.label = "T cell:  ") 

#OS for B cell
os.B.km <- plot.sig.fits(os.B, 3:25, legend.label = "B cell:  ") 
ggarrange(plotlist = os.B.km$plot.list)

#OS for T cell 
os.T.km <- plot.sig.fits(os.T, 3:25, legend.label = "T cell:  ") 
ggarrange(plotlist = os.T.km$plot.list) 

#Join OS B and T
test <- c(os.T.km$plot.list, os.B.km$plot.list)
ggarrange(plotlist = test[1:4], labels = "AUTO")
ggarrange(plotlist = test[5:7], labels = c("E", "F", "G"))

#Export all pvalues
write.csv(pfs.B.km$all.pvalues, file = "PFS_Bcell.csv", row.names = T)
write.csv(pfs.T.km, file = "PFS_Tcell.csv", row.names = T)
write.csv(os.B.km$all.pvalues, file = "OS_Bcell.csv", row.names = T)
write.csv(os.T.km$all.pvalues, file = "OS_Tcell.csv", row.names = T)



####3. Classification of B and T cell lymphoma ----
library(caret)
library(tidyverse)
library(pROC)

#Use transposed processed count dataframe from #1. Subset the required columns
df <- ct_data %>%
  select(c("dre-miR-18a", "hsa-miR-19a-3p", "cfa-miR-19b", "hsa-miR-21-5p", "hsa-miR-29b-3p", "hsa-miR-34a-5p", "cfa-miR-125a", "hsa-miR-130b-3p", "hsa-miR-146a-5p", "cfa-miR-181a", "cfa-miR-181b", "cfa-miR-181c", "hsa-miR-182-5p"))

#Add labels for B and T cells
df$Celltype <- as.factor(rep(c("Bcell", "Tcell"), c(22, 14)))

#Ensure that the values are "numeric" and celltype, "factor"
lapply(df, class)

#Split train and test data (75:25)
set.seed(2)
train.index <- createDataPartition(df$Celltype, p = .75, list = F)

train.data <- df[train.index, ] 
test.data <- df[-train.index, ] 
class.col <- 14
num.feats <- ncol(train.data) #number of features

trainX <- train.data[, -class.col]  #features
trainY <- train.data[,class.col]  #class



#Function to do RFE and training with hyperparameter tuning
rfe_train <- function(X, Y, algorithm, rfe.ctrl, tr.ctrl) {
  rfe(x = X, 
      y = Y, 
      sizes = seq(5, num.feats, 2),
      rfeControl = rfe.ctrl,
      
      #pass to train()
      method = algorithm,
      tuneLength = 6,
      trControl = tr.ctrl)
}



caretFuncs$summary <- prSummary #use precision-recall summary
rfectrl <- rfeControl(functions = caretFuncs,
                      method = "repeatedcv",
                      number = 10,
                      repeats = 5,
                      verbose = T)

trctrl <- trainControl(method = "LOOCV",
                        classProbs = T,
                       summaryFunction = prSummary)

#KNN 
set.seed(2)
rfe.knn <- rfe_train(trainX, trainY, algorithm = "knn", rfectrl, trctrl)
rfe.knn
rfe.knn$fit

#Linear SVM
set.seed(2)
rfe.svm <- rfe_train(trainX, trainY, "svmLinear", rfectrl, trctrl)
rfe.svm
rfe.svm$fit


#create a function to predict new data and return metrics
predict_new <- function (model, newdataX, newdataY) {
  pred.class <- predict(model, newdataX, 
                        type = ifelse(class(model) == "train", 
                                      yes = "raw", 
                                      no = "class")) #predict class ("raw" for svm)
  pred.prob <- predict(model, newdataX, type = "prob") #predict prob
  
  conf.mat <- confusionMatrix(pred.class, 
                              newdataY, 
                              mode = "prec_recall") #prediction metrics
  auroc <- roc(pred.prob[,1], response = newdataY) #auc
  ci.auroc <- ci.auc(auroc) #conf interval of auc

  return(
    list(conf.mat, auroc, ci.auroc) )
}

#Extract selected features for each model
knn.feats <- predictors(rfe.knn)
svm.feats <- predictors(rfe.svm)

#Test models with (already seen) training data 
set.seed(2)
predict_new(rfe.knn$fit$finalModel, trainX[knn.feats], newdataY = trainY) #knn
predict_new(rfe.svm$fit, trainX[svm.feats], newdataY = trainY) #svm


#Test with real test data 
testY <- test.data[, class.col]

set.seed(2)
predict_new(rfe.knn$fit$finalModel, test.data[knn.feats], newdataY = testY) #knn
predict_new(rfe.svm$fit, test.data[svm.feats], newdataY = testY) #knn


#KNN is chosen as the final model because metrics indicate less overfitting in KNN vs SVM 
final.trainX <- df[, -class.col] 
final.trainY <- df[,class.col]

set.seed(2)
final_model <- rfe_train(X = final.trainX, Y = final.trainY, algorithm = "knn", rfectrl, trctrl)
predictions <- final_model$fit$pred[final_model$fit$pred$k == final_model$fit$bestTune$k, ] #predictions at optimal k

confusionMatrix(predictions$pred, predictions$obs, mode = "prec_recall")
auroc <- roc(predictor = predictions$Bcell, response = predictions$obs)
plot(auroc, print.auc = T, 
     auc.polygon = T, 
     auc.polygon.col = "#B9D9EB", 
     max.auc.polygon = T, 
     max.auc.polygon.col = "beige")
ci.auc(auroc) 
        


####4. Use final model on new data- not part of the dataset ----
new_data <- readxl::read_excel("Test 1.1.xlsx")

#extract the ct values and transpose dataframe 
new_data <- new_data %>%
  slice(1:5) %>%
  column_to_rownames(var = "...1") %>%
  t(.)

#reorder to match the order in predictors(final_model) 
new_data <- as.data.frame(
  new_data[, c("hsa-miR-34a-5p", "hsa-miR-29b-3p", "cfa-miR-181b", "cfa-miR-181a", "hsa-miR-21-5p")] )
new.preds <- predict(final_model$fit$finalModel, newdata = new_data, type = "class")
new_data$Predictions <- new.preds

write.csv(new_data, file = "Test_data_with_predictions.csv")

####5. Differences in miRNA expression for patients receiving CHOP chemo (remission vs non-remission) ----
#B6, B10, and B17 did not reach remission in B 
#T2, T5, T6 in T

#B and T dataframes
B.df <- ct_data[1:22,]
T.df <- ct_data[23:36,]
 
B.remission <- B.df[-c(6, 10, 17), ]
B.nonremission <- B.df[c(6, 10, 17), ]

T.remission <- T.df[-c(2, 5, 6), ]
T.nonremission <- T.df[c(2, 5, 6), ]

#Mann Whitney test for remission vs remission in B and T with Benjamin and Hochberg corrected pvalues. df1 and df2 are the different dfs to be compared (e.g remission vs nonremission)
MW_test_adjusted <- function(df1, df2) { 
  results <- col_wilcoxon_twosample(df1, df2, exact = T)
  results %>% 
    mutate(pvalue.adj = p.adjust(pvalue, method = "BH")) %>%
    relocate(pvalue.adj, .after = pvalue)
}

rem.nonrem.B <- MW_test_adjusted(B.remission, B.nonremission)
rem.nonrem.T <- MW_test_adjusted(T.remission, T.nonremission)

#Export results (none are significant based on padjust)
write.csv(rem.nonrem.B[,5:6], file = "Rem_Nonrem_Bcell.csv")
write.csv(rem.nonrem.T[,5:6], file = "Rem_Nonrem_Tcell.csv")




####6. Differences in miRNA expression for patients receiving CHOP chemo (alive vs deceased at one year)

#Add OS days to B and T dfs
B.df$OSdays <- c(1052, 343, 1328,	1331,	744, 28,	261, 110,	308,	67,	288,	78,	489,	489,	328,	221,	123,	501,	594, 154, NA, 142)
#remove 21B since its OS days is NA
B.df <- na.omit(B.df)

T.df$OSdays <- c(152,	5, 160, 293,	22,	20,	243,	849,	120,	183,	103,	121,	128,	100)

#Separate dfs by alive vs not alive at one year
B.alive <- B.df[B.df$OSdays >= 365, ]
B.deceased <- B.df[B.df$OSdays < 365, ]

T.alive <- T.df[T.df$OSdays >= 365, ]
T.deceased <- T.df[T.df$OSdays < 365, ]

#Use the function from section #5 to compare miRNA expressn in alive vs deceased dataframes
OS.col <- ncol(B.df) #index of the OS column (last column of B or T df)

alive.vs.deceased.B <- MW_test_adjusted(B.alive[, -OS.col], B.deceased[, -OS.col])
alive.vs.deceased.T <- MW_test_adjusted(T.alive[, -OS.col], T.deceased[, -OS.col])

#Export results (none are significant based on padjust)
write.csv(alive.vs.deceased.B[, 5:6], "Bcell_alive_vs_deceased.csv")
write.csv(alive.vs.deceased.T[, 5:6], "Tcell_alive_vs_deceased.csv")
