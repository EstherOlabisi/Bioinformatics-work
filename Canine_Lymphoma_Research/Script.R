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

#This step adjusts for multiple testing errors using the Benjamin and Hochberg correction method (Yoav Benjamini & Hochberg, 1995). Add adjusted p-values to the results df. Place the new column right after the unadjusted pvalue column.
results <- results %>%
  mutate(pvalue.adj = p.adjust(pvalue, method = "BH")) %>%
  relocate(pvalue.adj, .after = pvalue)

#Select significant results with p < 0.05. 
sig.res <- results[results$pvalue.adj < 0.05, 5:6]  #Final result. 
unadj.sig.res <- results[results$pvalue < 0.05,] #For comparison sake. Only one additional miRNA here that's absent in the adjusted pvalue results

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

##Exploratory----
#To view all survival plots regardless of significance, change survival df and range.col as needed here.
survdf <- pfs.B #change df here
range.col <- 4:25 #change mirna column range here

#Then run these lines. Each plot will show on the viewer pane
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


