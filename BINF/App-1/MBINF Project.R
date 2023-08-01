####HEATMAP ----
# library(ggplot2)
# library(ggpubr)
# library(gplots)
# library(ggrepel)
# library(tidyverse)
# library(reshape2)
# require(colourpicker)
# require(shiny)
# require("shinyWidgets")


onedheatmap <- function(oned.df, plot.title = "") {
  
  oned.df <- oned.df %>%
    mutate_at(c("Score", "P.value", "Benj..Hoch..FDR"), as.numeric)
  
  #some annotations have the same name but belong to a different database
  #add additional column combining those entries
  oned.df$unique.annotation <- paste(oned.df$Name, " (", oned.df$Type , ")", sep="")
  
  #unique annotations and celltypes
  annotations <- sort(unique(oned.df$unique.annotation))
  annotations <- annotations[!grepl("^\\+", annotations)]
  celltypes <- sort(unique(oned.df$Column))
  
  #generate 1D score and pvalue(BH) matrix (takes long)
  M.score <- matrix(ncol=length(celltypes), nrow=length(annotations))
  M.pvalue <- matrix(ncol=length(celltypes), nrow=length(annotations))
  
  for(i in 1:length(celltypes)){
    for(j in 1:length(annotations)){
      score.value <- oned.df$Score[oned.df$Column==celltypes[i] & oned.df$unique.annotation==annotations[j]]
      p.value <- oned.df$Benj..Hoch..FDR[oned.df$Column==celltypes[i] & oned.df$unique.annotation==annotations[j]]
      if(length(score.value)==0){score.value <- NA}
      if(length(p.value)==0){p.value <- NA}
      M.score[j,i] <- score.value
      M.pvalue[j,i] <- p.value
    }
    #print(paste(i, celltypes[i], sep="; "))
  }
  
  m <- M.score
  rownames(m) <- annotations
  colnames(m) <- celltypes
  
  #Sequence of samples ('cell types' as currently above) for columns of heat map
  m <- m[,c(1:8)]
  m[is.na(m)] <- 0
  
  
  m <- melt(m)
  colnames(m) <- c("Annotations", "T-test differences", "value")
  
  return(
    ggplot(data = m, aes(x = `T-test differences`, 
                         y = reorder(Annotations, value), 
                         fill = value))+
      geom_tile(colour = "grey", linewidth = 1)+
      scale_fill_gradientn(colors = c(low = "blue", mid = "white", high = "red"),
                           na.value = "grey")+
      scale_y_discrete(position = "right")+
      guides(fill = guide_colourbar(title = "Colour Key"))+
      theme(axis.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 90),
            axis.text = element_text(colour = "black"),
            legend.position = "left")+
      labs(title = plot.title, 
           y = "Annotations")
  )

}


##Halos to pca plot
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
#Don't move app file out of #App-1 folder
#setwd to app location for easy access



#A function to process input data ---- 
#Read data containing t-test result. Remove all comment rows manually since the "comment.char" argument for read.delim() fails to do so. 
#Change column names to shorter generic names since the column names may slightly different among datasets. Columns of focus: -log10 p-values, Difference, and Significant proteins. Also, select the first GO terms using regex. 

process.df <- function(dataset, left.range, right.range) {
  comment_char <- grep("^#", dataset[,1])
  dataset <- dataset[-(comment_char),]
  
  dataset <- dataset %>%
    setNames(sub(".+p.value.+", "minus.log10.pval", names(.))) %>% #-log10pval
    setNames(sub(".+Difference.+", "Difference", names(.))) %>% #Difference
    setNames(sub(".+Significant.+", "Significant", names(.))) %>% #signficant proteins
    mutate(Keywords =  sub(";.*$", "", Keywords))  #select first keywords
  
  
  
  #Extract group labels from the "significant" column for use in legend
  sig.df <- dataset[dataset$Significant == "+",]  
  gr.name.col <- grep("significant", colnames(sig.df)) 
  group.labs <- unique(sig.df[, gr.name.col]) #extract label 
  group.labs <- str_split_1(group.labs, "_") #separate labels
  
  #Add the group labels based on positive or negative x axis then calculate mean across groups
  dataset <- dataset %>%
    mutate(Expression = case_when(Difference >= 0 & Significant == "+" ~ paste0("Higher abundance in ", group.labs[1]),
                                  Difference <= 0 & Significant == "+" ~ paste0("Higher abundance in ", group.labs[2]),
                                  TRUE ~ "Not significant")) %>%
    mutate(Expression = factor(Expression, levels = c(
      paste0("Higher abundance in ", group.labs[2]),
      paste0("Higher abundance in ", group.labs[1]),
      "Not significant"))) %>%  #reorder
    mutate_at(c('Difference', 'minus.log10.pval'), as.numeric) %>%
    mutate_at(vars(contains("lfq")), as.numeric) %>%
    rowwise() %>%
    mutate(left.group = mean(
      c_across(all_of(left.range)))) %>%
    mutate(right.group = mean(
      c_across(all_of(right.range)))) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Majority.protein.IDs")
  
  return(dataset)
}



#A function for the Volcano plot: ----
#The function requires an input dataframe(df) while other arguments are optional.
#curves.df = dataframe for fdr curves
#go.terms = display GO terms for significant or non-significant proteins
##plot.title = custom plot title top left
##group.cols = colours for the default plot in the following order: higher-left, higher-right colours. 
##s0 = s0 value to show on the plot. Default = 0.1
##fdr = fdr value to show on the plot. Default = 0.05
##fdr.lines = 'yes' to show lines and 'no' to remove lines
##palette.col = based on hcl.color() palettes. "Viridis" is the default one
##left.st, left.end = index range for left group's LFQ columns - for calculating mean across the columns
##right.st, right.end = index range for right group's LFQ columns

volcano_plot <- function(df, curves.df, 
                         go.terms = "", plot.title = "", 
                         s0 = 1, fdr = 0.05, fdr.lines = "yes", 
                         palette.col = "Viridis", group.cols = c("#00B1B2", "#FF6666"),
                         vary.sizes = "no",
                         select.pts = c()) {
  
  #Subset significant and non-significant(ns) proteins into their respective dfs. 
  #Extract the number of unique GO terms in both sig and ns dfs. The numbers will be used for color-coding the GO terms downstream. 
  #Subset the selected proteins so that we can use their protein ids(rownames) as point labels on the plot
  #Also, extract the x-axis label (in "Group1_Group2" format) from "significant" column
  sig.df <- df[df$Significant == "+",]
  ns.df <- df[df$Significant == "",] 
  num.sig.go <- length(unique(sig.df$Keywords)) 
  num.ns.go <- length(unique(ns.df$Keywords)) 
  select.pts <- df[select.pts, ]
  gr.name.col <- grep("significant", colnames(df)) #column index 
  xlab <- unique(sig.df[, gr.name.col]) 
  xlab <- sub("_", " - ", xlab)
  group.labs <- str_split_1(xlab, " - ")
  
  
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
    coord_cartesian(xlim = c(min(df$Difference),
                             max(df$Difference)),
                    ylim = c(min(df$minus.log10.pval),
                             max(df$minus.log10.pval)),
                    expand = T)+
    theme_minimal()+
    theme(plot.title = element_text(face="bold", size = 20),
          plot.caption = element_text(
            colour = "darkviolet", size = 11))
  
  
  #Fdr curves layer 
  curve.plot <- geom_line(data = curves.df, aes(x, y), linetype=2)
  
  #Layer for labelling points 
  label.pt.plot <- geom_label_repel(data = select.pts, 
                                    aes(label = rownames(select.pts)),
                                    nudge_y = 0.4,  fill = "grey")
  
  ###Adding layers to the base plot
  #Add GO terms for sig proteins
  sig.go <- base_vplot +
    scale_fill_identity()+  #hides grey legend that comes from the base plot
    new_scale_fill()+
    geom_point(data = sig.df,
               aes(fill = Keywords), 
               size = s.sz,
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.sig.go, palette.col, rev = F)) + 
    label.pt.plot
  
  #Add GO terms for non-sig proteins. Grey out sig proteins. 
  nonsig.go <- base_vplot +
    scale_fill_identity()+
    new_scale_fill()+
    geom_point(aes(fill = Keywords),
               size = ns.sz, 
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.ns.go,
      palette.col,
      rev = F)) + 
    geom_point(data = sig.df, fill = "grey",
               size = s.sz, alpha = a,
               shape = s) + 
    label.pt.plot
  
  ###Add layer showing average protein abundance across samples 
  #Where right sample has the higher averages, the right averages are used for size variation. The same for where the left sample has the higher averages
  vary.size.pl <- base_vplot +
    scale_fill_identity()+
    new_scale_fill()+
    geom_point(data = sig.df[sig.df$Difference > 0, ],
               aes(size = right.group),
               fill = alpha("#FF6666", a),
               shape = s)+
    scale_size_binned(range = c(1,8),
                      name = paste0("Average Intensity ", group.labs[1]),
                      n.breaks = 4)+
    new_scale(new_aes = "size")+
    geom_point(data = sig.df[sig.df$Difference < 0, ],
               aes(size = left.group), 
               fill = alpha("#00B1B2", a),
               shape = s)+
    scale_size_binned(range = c(1,8),
                      name = paste0("Average Intensity ", group.labs[2]),
                      n.breaks = 4) + 
    label.pt.plot
  
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
  } else if(vary.sizes == "yes" && fdr.lines == "yes") {
    return (vary.size.pl + curve.plot)
    
  } else if (vary.sizes == "yes" && fdr.lines == "no") {
    return (vary.size.pl)
    
  } else if (fdr.lines == "yes") {
    return(
      base_vplot +
        scale_fill_identity(name = NULL,
                            labels = "Not significant",
                            guide = "legend")+
        new_scale_fill()+
        geom_point(data = sig.df, aes(fill = Expression),
                   size = s.sz, shape = s) +
        scale_fill_manual(values = group.cols)+
        guides(fill = guide_legend(order = 1))+
        label.pt.plot + curve.plot + 
        theme(legend.spacing = unit(-0.5, "cm"))
    )
    
  } else if (fdr.lines == "no") {
    return(
      base_vplot +
        scale_fill_identity(name = NULL,
                            labels = "Not significant",
                            guide = "legend")+
        new_scale_fill()+
        geom_point(data = sig.df, aes(fill = Expression),
                   size = s.sz, shape = s) +
        scale_fill_manual(values = group.cols)+
        guides(fill = guide_legend(order = 1))+
        label.pt.plot +
        theme(legend.spacing = unit(-0.5, "cm"))
    )
  }
  
}


####test ----
# dt <- read.delim("../ttest_table.txt", na.strings = c("NA", "NaN"))
# dt <- process.df(dataset = dt, left.range = c("LFQ.intensity.WT.1", "LFQ.intensity.WT.2", "LFQ.intensity.WT.3", "LFQ.intensity.WT.4"), right.range = c("LFQ.intensity.401.1", "LFQ.intensity.401.2", "LFQ.intensity.401.3", "LFQ.intensity.401.4"))
#    curvesdff <- read.table("../curve_matrix.txt", header = T)
#volcano_plot(dt, curvesdff, go.terms = "non-sig", fdr.lines = "no", select.pts = c("A6T4E4"))
# 
#  df <- read.delim("../1d annot.txt", stringsAsFactors = F, header=T)
#  df<- df[-c(1),]
# unique.annot.types <- df %>%
#    group_by(Type) %>%
#    group_split()
#  plots <- lapply(unique.annot.types,
#                  function(x) onedheatmap(x,
#                                          plot.title = paste(unique(x$Type))
#                  ))
# plots
# ggarrange(plotlist = plots)
# onedheatmap(df)



##SHINY GUI ----
#flowchart 
#perseus vs my result
#these proteins associated with this go term
#GO terms on heatmap only
#make transparent circles
#ignore NA warning, comes form the ones left in dataframe

#Shiny UI ---- 
ui <- fluidPage(
  titlePanel("Visualization for Proteomics"),
  
  sidebarLayout(
    
    sidebarPanel(width = 4,
                 
                 fluidRow(
                   column(width = 4,
                          fileInput(inputId = "proteinfile", label = "Protein data input") ),
                   column(width = 4, 
                          fileInput("curvesfile", label = "Volcano curves data input") ),
                   column(width = 4,
                          fileInput("onedfile", label = "1D annotation data input") ),
                   uiOutput("enter.range"),
                   column(width = 6, selectInput("left.gr", 
                                                 label = "Column range for Left group",
                                                 choices = c(),
                                                 multiple = T) ),
                   column(width = 6, selectInput("right.gr", 
                                                 label = "Column range for Right group",
                                                 choices = c(),
                                                 multiple = T) ),
                   
                   h4("ggplot Options"),
                   column(width = 12, br()),
                   column(width = 6, colourInput("left.col",
                                                 label = "Left group colour",
                                                 value = "#00B1B2") ),
                   column(width = 6, colourInput("right.col",
                                                 label = "Right group colour",
                                                 value = "#FF6666") )
                   
                 ),
                 actionButton("reset", label = "Reset Volcano Plot"),
                 br(),
                 
                 br(),
                 fluidRow(
                   
                   column(width = 3, numericInput("xmin",
                                                  label = "x min", 
                                                  value = NULL) ),
                   column(width = 3, numericInput("xmax", 
                                                  label = "x max", 
                                                  value = NULL) ),
                   column(width = 3, numericInput("ymin", 
                                                  label = "y min", 
                                                  value = NULL) ),
                   column(width = 3, numericInput("ymax", 
                                                  label = "y max", 
                                                  value = NULL) ),
                   
                   column(width = 12, br()),
                   
                   column(width = 6, numericInput("s.knot",
                                                  label = "s0 value",
                                                  value = 1) ),
                   
                   column(width = 6, numericInput("fdr.val",
                                                  label = "FDR value",
                                                  value = 0.05) ),
                   
                   column(width = 12, 
                          helpText("Note: Changing the s0 and FDR values is purely aesthetic for the
               captions at the base of the plot. The plotted points remain unaffected.") ),
                   
                   column(width = 12, br()),
                   column(width = 6, radioButtons("fdr.lines",
                                                  label = "Display curves",
                                                  choices = list("yes", "no"),
                                                  selected = "yes", inline = T) ),
                   
                   column(width = 6, radioButtons("pt.sizes",
                                                  label = "Average protein intensities",
                                                  choices = list("yes", "no"),
                                                  selected = "no", inline = T) ),
                   
                   
                   column(width = 7, radioButtons("go.term",
                                                  label = "GO terms",
                                                  choices = list("sig", "non-sig", "reset"),
                                                  selected = "reset", inline = T) ),
                   
                   column(width = 5, selectInput("palette",
                                                 label = "GO terms color palette",
                                                 choices = list("Viridis", "Plasma", "Inferno", "Rocket", "Lajolla", "Turku", "Hawaii", "Batlow"),
                                                 selected = "Viridis") )
                   
                 )
    ),
    
    mainPanel(
      uiOutput("df.preview"),
      dataTableOutput("dataframe"),
      br(),
      h4("Volcano Plot",  align = "center"),
      
      ##volcano plot output with brush selection of points
      plotOutput("volcanoplot",
                 brush = "vplot.brush"),
      textInput("vplot.title",
                label = "Volcano Plot Title (Optional)",
                value = NULL,
                placeholder = "Enter text..."),
      fluidRow(
        ##initialize drop down list of protein labels
        column(width = 3, selectizeInput("protein.labs",
                                         label = "Manually select protein IDs:",
                                         choices = list(),
                                         multiple = T)),
        column(width = 6, verbatimTextOutput("brushed.ids"))
      ),
      
      br(),
      br(),
      br(),
      ##1D plots output
      h4("1D Annotation Heatmaps", align = "center", ),
      plotOutput("onedhm1.bio.proc"),
      br(),
      plotOutput("onedhm2.cell.comp"),
      br(),
      plotOutput("onedhm3.mol.func"),
      br(),
      plotOutput("onedhm4.keywords"),
    )
  )
)




#Shiny Server ----
server <- function(input, output, session) {
  
  #1D annotation file
  df.1d <- eventReactive(input$onedfile, {
    read.delim("../1d annot.txt", stringsAsFactors = F, header=T)
    df<- df[-c(1),]
  })
  
  #Curve points file
  df.curves <- eventReactive(input$curvesfile, {
    read.table(input$curvesfile$datapath, header = T)
  })
  
  #Protein input file
  df.prot <- eventReactive(input$proteinfile, {
    read.delim(input$proteinfile$datapath, na.strings = c("NA", "NaN"))
  })
  
  #Add column names for selection through the UI
  observe({
    req(df.prot())
    updateSelectInput(session,
                      inputId = "left.gr",
                      choices = colnames(df.prot()))
    updateSelectInput(session,
                      inputId = "right.gr",
                      choices = colnames(df.prot()))
    output$enter.range <- renderUI({
      helpText("To calculate average protein intensities across groups, select columns for each group:")
    })
  })
  
  #Process the protein dataframe and calculate means across group columns
  df.prot2 <- reactive({
    req(df.prot(), input$left.gr, input$right.gr)
    process.df(dataset = df.prot(),
               left.range = input$left.gr,
               right.range = input$right.gr )
  })
  
  ##Preview df
  observe({
    req(df.prot2())
    output$df.preview <- renderUI({
      helpText(h4("Dataframe preview"))
    })
    output$dataframe <- renderDataTable({
      df.prot2() %>% head(4) })
    #fill dropdown list with protein names (labels) from which to choose
    updateSelectizeInput(session,
                         inputId = "protein.labs",
                         choices = rownames(df.prot2()), 
                         server = T )
  })
  
  ##activate brush selection    
  brushed.pts <- reactive({
    brushedPoints(df.prot2(), input$vplot.brush)
  })
  ##display selected ids as text on the UI
  output$brushed.ids <- renderPrint({ 
    cat("Selected proteins: ", rownames(brushed.pts()) ) 
  }) 
  
  ##create volcano plot using UI input variables  
  vplot <- reactive({
    volcano_plot(df = df.prot2(),
                 curves.df = df.curves(),
                 go.terms =  input$go.term,
                 palette.col = input$palette,
                 plot.title = input$vplot.title,
                 s0 = input$s.knot,
                 fdr = input$fdr.val,
                 group.cols = c(input$left.col, input$right.col),
                 fdr.lines = input$fdr.lines,
                 vary.sizes = input$pt.sizes,
                 #show protein labels either by quick or fixed selection 
                 select.pts = c(rownames(brushed.pts()), 
                                input$protein.labs) 
    ) 
  })
  #send vplot to the UI
  output$volcanoplot <- renderPlot({
    vplot()
  })
  
  ##Axes limit adjustments
  observeEvent(req(input$xmin | input$xmax | input$ymin | input$ymax), {
    output$volcanoplot <- renderPlot({
      vplot()+
        coord_cartesian(xlim = c(input$xmin, input$xmax),
                        ylim = c(input$ymin, input$ymax),
                        expand = T)
    })  
  })
  
  ##Volcano plot reset 
  observeEvent(input$reset, {
    updateColourInput(session, inputId = "left.col", value = "#00B1B2")
    updateColourInput(session, inputId = "right.col", value = "#FF6666")
    output$volcanoplot <- renderPlot({ vplot() })
  })
  
  
  ##Heatmap plot list
  plotlist.1d <- reactive({
    unique.annot.types <- df.1d() %>%
      group_by(Type) %>%
      group_split()
    lapply(unique.annot.types, 
           function(x) onedheatmap(x, 
                                   plot.title = paste(unique( x$Type ))
           ))
  })
  
  #plot heatmaps individually
  output$onedhm1.bio.proc <- renderPlot({ plotlist.1d()[1] }) #render 1d plot1
  output$onedhm2.cell.comp <- renderPlot({ plotlist.1d()[2] }) #render 1d plot2
  output$onedhm3.mol.func <- renderPlot({ plotlist.1d()[3] }) #render 1d plot3
  output$onedhm4.keywords <- renderPlot({ plotlist.1d()[4] }) #render 1d plot4
  
  
  
}#do 1, 2, 3, 4

shinyApp(ui = ui, server = server)

#add image download as pdf/ png