####HEATMAP ----
library(ggplot2)
library(gplots)


#Added comment.char argument to skip the comment line
df <- read.delim("../Sylvain_1(keywords).txt", header=T, stringsAsFactors=F, comment.char = "#") 

onedheatmap <- function(oned.annotation.enrichment.df, plot.title = "") {
  #some annotations have the same name but belong to a different database
  #add additional column combining those entries
  oned.annotation.enrichment.df$unique.annotation <- paste(oned.annotation.enrichment.df$Name, " (", oned.annotation.enrichment.df$Type , ")", sep="")
  
  #unique annotations and celltypes
  annotations <- sort(unique(oned.annotation.enrichment.df$unique.annotation))
  annotations <- annotations[!grepl("^\\+", annotations)]
  celltypes <- sort(unique(oned.annotation.enrichment.df$Column))
  
  #generate 1D score and pvalue(BH) matrix (takes long)
  M.score <- matrix(ncol=length(celltypes), nrow=length(annotations))
  M.pvalue <- matrix(ncol=length(celltypes), nrow=length(annotations))
  
  for(i in 1:length(celltypes)){
    for(j in 1:length(annotations)){
      score.value <- oned.annotation.enrichment.df$Score[oned.annotation.enrichment.df$Column==celltypes[i] & oned.annotation.enrichment.df$unique.annotation==annotations[j]]
      p.value <- oned.annotation.enrichment.df$Benj..Hoch..FDR[oned.annotation.enrichment.df$Column==celltypes[i] & oned.annotation.enrichment.df$unique.annotation==annotations[j]]
      if(length(score.value)==0){score.value <- NA}
      if(length(p.value)==0){p.value <- NA}
      M.score[j,i] <- score.value
      M.pvalue[j,i] <- p.value
    }
    print(paste(i, celltypes[i], sep="; "))
  }
  
  
  m <- as.matrix(M.score)
  rownames(m) <- annotations
  colnames(m) <- celltypes
  
  ##m <- m[,c(2,3,1)]
  #Sequence of samples ('cell types' as currently above) for columns of heat map
  m <- m[,c(1,2,3,4,5,6,7,8)]
  m[is.na(m)] <- 0
  
  return ( heatmap.2(m,
            col = bluered(256),
            # breaks=col_breaks,
            #ColSideColors=NULL,
            margins = c(10, 20),
            trace = "none", 
            main = plot.title,
            xlab = "T-test differences",
            ylab = "Annotations",
            # labRow = row.label,
            # labCol = col.label,
            #lmat=rbind(c(4,3,3), c(2,1,1), c(2,1,1)), 
            #lhei=c(1,2,3), 
            #lwid=c(1,2,3),
            #scale = c("none"),
            #symbreaks = min(18, na.rm=TRUE),
            na.color="grey",
            na.rm=F,
            #cexRow = .9, cexCol = .9,
            #main = header.heatmap, 
            dendrogram = "none", 
            density.info = "none",
            Colv = F,
            Rowv = T,
            cexRow = 1,
            cexCol = 1,
            # block sepration
            colsep = c(1:ncol(m)),
            rowsep = c(1:nrow(m)),
            sepcolor="grey",
            sepwidth=c(0.05,0.05))
            
           #cellnote=round(m.exclusive.proteins, 1),
           #notecol="black",
           #notecex=.7
  )
}


##Halos to pca plot
##change range of axes
#Versioning software 
##Use second table as input

####VOLCANO PLOT ----
#BiocManager::install(c("msmsTests", "msmsEDA"))
#library("msmsTests")
#library("msmsEDA")
require(tidyverse)
require(RColorBrewer)
require(ggnewscale)


##Features
#Can show levels of significance
#Can color code for gene ontology terms among sigs and nonsigs
#Working on size differences for protein intensities based on avgs
#Working on most sig labels
#Don't move app file out of #App-1 folder
#setwd to app location for easy access



#A function to process input data ---- 
#Read data containing t-test result. Remove all comment rows manually since the "comment.char" argument for read.delim() fails to do so. 
#Change column names to shorter generic names since the column names may slightly different among datasets. Columns of focus: -log10 p-values, Difference, and Significant proteins. Also, select the first GO terms using regex. 

process.df <- function(data = "") {
  
  dataset <- read.delim(data)
  comment_char <- grep("^#", dataset[,1])
  dataset <- dataset[-(comment_char),]
  lfq.cols <- grep("lfq", colnames(df), ignore.case = T) # lfq columns
   
  dataset <- dataset %>%
    setNames(sub(".+p.value.+", "minus.log10.pval", names(.))) %>% #-log10pval
    setNames(sub(".+Difference.+", "Difference", names(.))) %>% #Difference
    setNames(sub(".+Significant.+", "Significant", names(.))) %>% #signficant proteins
    mutate(Keywords =  sub(";.*$", "", Keywords)) %>% #select first keywords
    remove_rownames() %>%
    column_to_rownames(var = "Majority.protein.IDs") %>%
    mutate(Expression = case_when(Difference >= 0 & Significant == "+" ~ "High",
                                  Difference <= 0 & Significant == "+" ~ "Low",
                                  TRUE ~ "Not significant")) %>%
    mutate(Expression = factor(Expression,
                               levels = c("High", "Low", "Not significant"))) %>%
    mutate_at(c('Difference', 'minus.log10.pval'), as.numeric) 
  
  
  return(dataset)
}



#Creating a function for volcano plot: ----
#The function requires an input dataframe(df) while other arguments are optional.
#curves.df = dataframe for fdr curves
#go.terms = display GO terms for significant or non-significant proteins
##plot.title = custom plot title top left
##default.cols = colours for the default plot in the following order: high-expressing, low-expressing colours. 
##s0 = s0 value to show on the plot. Default = 0.1
##fdr = fdr value to show on the plot. Default = 0.05
##fdr.lines = 'yes' to show lines and 'no' to remove lines
##palette.col = based on hcl.color() palettes. "Viridis" is the default one
##left.st, left.end = index range for left group's LFQ columns - for calculating mean across the columns
##right.st, right.end = index range for right group's LFQ columns

volcano_plot <- function(df, curves.df, go.terms = paste(""), plot.title = "", s0 = 1, fdr = 0.05, fdr.lines = paste("yes"), palette.col = paste("Viridis"), default.cols = c("#FF6666", "#00B1B2")) {
  
  #Subsetting all significant and non-significant(ns) proteins into different dfs. Extract the number of unique GO terms in both sig and ns dfs. The numbers will be used for color-coding the GO terms downstream. Also, extract the x-axis label (in "Group1_Group2" format) from "significant" column
  sig.df <- df[df$Significant == "+",]
  ns.df <- df[df$Significant == "",] 
  num.sig.go <- length(unique(sig.df$Keywords)) 
  num.ns.go <- length(unique(ns.df$Keywords)) 
  label.col <- grep("significant", colnames(df)) #column index 
  xlab <- unique(sig.df[label.col]) 
  
  
  
  #Custom settings for ggplots
  a = 0.2 #alpha
  s = 21  #shape
  ns.sz = 2 #non-sig shape size  
  s.sz = 4  #sig shape size
  
  #Base volcano plot to plot points, fdr curves, axes boundaries, and other custom arguments.
  base_vplot <- ggplot(data = ns.df, 
                       aes(x = Difference,
                           y = minus.log10.pval))+
    geom_point(aes(fill = "grey"), size = ns.sz,
               alpha = a, shape = s)+
    labs(title = plot.title,
         caption = paste("s0 =", s0, "  ", "FDR =", fdr),
         x = paste0("Difference (", xlab, ")"), 
         y = expression(-log[10]("p-value")))+
    xlim(min(df$Difference)-1,
         max(df$Difference)+1)+
    ylim(min(df$minus.log10.pval),
         max(df$minus.log10.pval)+2)+
    theme_minimal()+
    theme(plot.title = element_text(face="bold", size = 20),
          plot.caption = element_text(
            colour = "darkviolet", size = 11))
  
  #Fdr curves layer 
  curve.plot <- geom_line(data = curves.df, aes(x, y), linetype=2)
  
  
  ###Adding layers to the base plot
  #Add GO terms for sig proteins
  sig.go <- base_vplot +
    scale_fill_identity()+  #hides grey legend from base plot
    new_scale_fill()+
    geom_point(data = sig.df,
               aes(fill = Keywords), size = s.sz,
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.sig.go, palette.col, rev = F))
  
  #Add GO terms for non-sig proteins. Grey out sig proteins. 
  nonsig.go <- base_vplot +
    scale_fill_identity()+
    new_scale_fill()+
    geom_point(aes(fill = Keywords),
               size = ns.sz, shape = s)+
    geom_point(data = sig.df, fill = "grey",
               size = s.sz, alpha = a,
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.ns.go,
      palette.col,
      rev = F))
  
  ###Add layer showing average protein abundance across samples 
  #Where right sample has the higher averages, the right averages are used for size variation. The same for where the left sample has the higher averages
  # rowwise() %>%
  #   mutate(left.grp = mean(c_across(left.st : left.end), na.rm = T)) %>%
  #   mutate(right.grp = mean(c_across(right.st : right.end), na.rm = T))
  # 
  # 
  #   vary.size.pl <- base_vplot + 
  #   scale_fill_identity()+
  #   new_scale_fill()+
  #   geom_point(data = sig.df[sig.df$Difference > 0, ],
  #              aes(size = right.grp),
  #              fill = alpha("#FF6666", 0.2),
  #              shape = 21)+
  #   scale_size_binned(range = c(1,8),
  #                     name = "Avg Intensity",
  #                     n.breaks = 4)+
  #   new_scale(new_aes = "size")+
  #   geom_point(data = sig.df[sig.df$Difference < 0, ],
  #              aes(size = left.grp), fill = alpha("#00B1B2", 0.2),
  #              shape = 21)+
  #   scale_size_binned(range = c(1,8),
  #                     name = "Avg Intensity",
  #                     n.breaks = 4)
  # 
  #Options: GO term plot for sig proteins, with fdr lines and without fdr lines
  if (go.terms == "sig" && fdr.lines == "yes") {
    return (sig.go + curve.plot)
    
  } else if (go.terms == "sig" && fdr.lines == "no") {
    return (sig.go)
  
  #Options: GO terms for non-sig proteins with fdr lines, without fdr lines, and default plot with color-coding of high and low expression proteins. 
  } else if (go.terms == "non-sig" && fdr.lines == "yes") {
    return (nonsig.go + curve.plot)
    
  } else if (go.terms == "non-sig" && fdr.lines == "no") {
    return (nonsig.go)
  
  #Sizes
  # } else if(vary.sizes == "yes" && fdr.lines == "yes") {
  #   return (vary.size.pl + curve.plot)
  # 
  # } else if (vary.sizes == "yes" && fdr.lines == "no") {
  #   return (vary.size.pl)
  #   
  } else if (fdr.lines == "yes") {
    return(
      base_vplot +
        scale_fill_identity(name = NULL,
                            labels = "Not significant",
                            guide = "legend")+
        new_scale_fill()+
        geom_point(data = sig.df, aes(fill = Expression),
                   size = s.sz, shape = s) +
        scale_fill_manual(values = default.cols)+
        guides(fill = guide_legend(order = 1))+
        theme(legend.spacing = unit(-0.5, "cm"))+
        curve.plot)
    
  } else if (fdr.lines == "no") {
    base_vplot +
      scale_fill_identity(name = NULL,
                          labels = "Not significant",
                          guide = "legend")+
      new_scale_fill()+
      geom_point(data = sig.df, aes(fill = Expression),
                 size = s.sz, shape = s) +
      scale_fill_manual(values = default.cols)+
      guides(fill = guide_legend(order = 1))+
      theme(legend.spacing = unit(-0.5, "cm"))
  }
  
}



#dt <- process.df(data = "ttest_table.txt")
#curvesdff <- read.table("curve_matrix.txt", header = T)
volcano_plot(dt, curvesdff)

#left.st, left.end, right.st, right.end
dt2 <- dt %>%
  mutate_at(c('Difference', 'minus.log10.pval'), as.numeric) %>%
  mutate_at(1:8, as.numeric) %>%
  rowwise() %>%
  mutate(left.grp = mean(c_across(5:8), na.rm = T)) %>%
  mutate(right.grp = mean(c_across(1:4), na.rm = T))
sig.df <- dt2[dt2$Significant == "+",]
ns.df <- dt2[dt2$Significant == "",]


#test
ggplot(data = ns.df,
       aes(x = Difference,
           y = minus.log10.pval))+
  geom_point(aes(fill = "grey"), size = 2,
             alpha = 0.2, shape = 21)+
  scale_fill_identity()+
  new_scale_fill()+
  geom_point(data = sig.df[sig.df$Difference > 0, ],
                          aes(size = right.grp),
                          fill = alpha("#FF6666", 0.2),
                          shape = 21)+
               scale_size_binned(range = c(1,8),
                                 name = "Avg Intensity",
                                 n.breaks = 4)+
               new_scale(new_aes = "size")+
               geom_point(data = sig.df[sig.df$Difference < 0, ],
                          aes(size = left.grp), fill = alpha("#00B1B2", 0.2),
                          shape = 21)+
               scale_size_binned(range = c(1,8),
                                 name = "Avg Intensity",
                                 n.breaks = 4)






##SHINY GUI ----
require(colourpicker)
require(shiny)
#flowchart 
#perseus vs my result
#these proteins associated with this go term
#GO terms on heatmap only
#make transparent circles


#Shiny UI ---- 
ui <- fluidPage(
  titlePanel("Visualization for Proteomics"),
  
  sidebarLayout(
    
    sidebarPanel(width = 4,
                 fileInput(inputId = "proteinfile", label = "Protein data input"),
                 
                 fileInput("curvesfile", label = "FDR curves data input"),
                 
                 fileInput("onedfile", label = "1D annotation data input"),
                 
                 br(),
                 h4("ggplot Options"),
                 
                 colourInput("high.col",
                             label = "High expression colour",
                             value = "#FF6666"),
                 
                 colourInput("low.col",
                             label = "Low expression colour",
                             value = "#00B1B2"),
                 
                 actionButton("reset", label = "Reset colours"),
                 br(),
                 
                 br(),
                 helpText("Note: Changing the s0 and FDR values is purely aesthetic for the
               captions at the base of the plot. The plotted points remain unaffected."),
                 
                 numericInput("s.knot",
                              label = "s0 value",
                              value = 1),
                 
                 numericInput("fdr.val",
                              label = "FDR value",
                              value = 0.05),
                 
                 radioButtons("fdr.lines",
                              label = "FDR lines",
                              choices = list("yes", "no"),
                              selected = "yes"),
                 
                 radioButtons("go.term",
                              label = "GO terms",
                              choices = list("sig", "non-sig", "reset"), selected = "reset"),
                 
                 selectInput("palette",
                             label = "GO terms color palette",
                             choices = list("Viridis", "Plasma", "Inferno", "Rocket", "Lajolla", "Turku", "Hawaii", "Batlow"),
                             selected = "Viridis"),
                 
    ),
    
    mainPanel(
      
      br(),
      h4("Volcano Plot",  align = "center"),
      plotOutput("volcanoplot"), #volcano plot
      textInput("vplot.title",
                label = "Volcano Plot Title (Optional)",
                value = NULL,
                placeholder = "Enter text..."),
      
      br(),
      br(),
      br(),
      h4("1D Annotation Heatmap", align = "center", ),
      plotOutput("onedheatmap"),
      textInput("hmplot.title",
                label = "Heat Map Title (Optional)",
                value = NULL,
                placeholder = "Enter text..."),
    )
  )
)


#Shiny Server ----
server <- function(input, output, session) {
  
  #process both protein and fdr curve input files
  df.prot <- eventReactive(input$proteinfile, {
    process.df(data = input$proteinfile$datapath)
  })
  df.curves <- eventReactive(input$curvesfile, {
    read.table(input$curvesfile$datapath, header = T)
  })
  df.1d <- eventReactive(input$onedfile, {
    read.delim(input$onedfile$datapath, header=T, 
               stringsAsFactors=F, comment.char = "#")
  })
  
  ##create volcano plot with UI input variables   
  output$volcanoplot <- renderPlot({
    
    volcano_plot(df = df.prot(),
                 curves.df = df.curves(),
                 go.terms =  input$go.term,
                 palette.col = input$palette,
                 plot.title = input$vplot.title,
                 s0 = input$s.knot,
                 fdr = input$fdr.val,
                 fdr.lines = input$fdr.lines,
                 default.cols = c(input$high.col, input$low.col)
    )
    
  })   
  #allow for the volcano plot to reset to the default figure
  observeEvent(req(input$go.term == "reset"), {
    volcano_plot(df = df.prot(),
                 curves.df = df.curves(),
                 default.cols = c("#FF6666", "#00B1B2"))
  })
  
  observeEvent(input$reset, {
    updateColourInput(session, inputId = "high.col", value = "#FF6666")
    updateColourInput(session, inputId = "low.col", value = "#00B1B2")
    #    output$volcanoplot <- renderPlot({
    #  volcano_plot(df = df.prot(),
    #            curves.df = df.curves(),
    #           default.cols = c("#FF6666", "#00B1B2"))
    
  })
  
  
  #Heatmap 
  output$onedheatmap <- renderPlot({
    onedheatmap(df.1d(), plot.title = input$hmplot.title)
  })
  
}

shinyApp(ui = ui, server = server)

#add image download as pdf/ png 