#Esther Olabisi-Adeniyi
#Bioinformatics Master's Research Project - BINF6999
#August 8, 2023

#All packages used
# my.packages <- c("ggplot2", "ggpubr", "ggrepel", "tidyverse", "DT", "reshape2", "scales", "RColorBrewer", "colourpicker", "ggnewscale", "ggfortify", "ggrepel", "shiny", "shinyjs", "shinyWidgets", "eulerr", "webshot")
#Check if packages are installed
# installed.pkgs <- my.packages %in% rownames(installed.packages())
#Install packages that aren't already installed
# if(any(installed.pkgs == FALSE)) {
  # install.packages(my.packages[!installed.pkgs])
# }
#Load all packages
#invisible(lapply(my.packages, library, character.only = T))

#Option to run packages individually when running app on Shinyapps.io
library(ggplot2); library(ggpubr); library(ggrepel); library(tidyverse); library(DT); library(reshape2); library(scales); library(RColorBrewer); library(colourpicker); library(ggnewscale); library(ggfortify); library(ggrepel); library(shiny); library(shinyjs); library(shinyWidgets); library(ggvenn); library(shinyBS)


####1. PROCESS DATA MATRIX ----
#Read protein.txt matrix processed in Perseus for a volcano plot. Remove all comment rows. 
#Change column names to shorter generic names that are easier to work with. 
#prot.dataset: Imported from Perseus following quality control, two samples t-test, and GO enrichment
#left.range: Intensity columns belonging to the left group of the volcano plot / two samples test in Perseus
#right.range: The same as left.range but for the right group

process.df <- function(prot.dataset, left.range, right.range) {
  comment_char <- grep("^#", prot.dataset[,1])
  prot.dataset <- prot.dataset[-(comment_char),]
  col.names <- colnames(prot.dataset)
  intensity.cols <- grep("intensity.+", col.names, ignore.case = T, value = T)
  major.prot.col <- grep("majority", col.names, ignore.case = T, value = T)
  grp.name.col <- grep("significant", col.names, value = T) 
  
  #select relevant columns only
  prot.dataset <- prot.dataset %>% 
    setNames(sub(".*p.value.*", "Minus.log.pval", names(.))) %>% #-log10pval
    setNames(sub(".*Difference.*", "Difference", names(.))) %>% #Difference
    setNames(sub(".*Significant.*", "Significant", names(.))) %>%
    setNames(sub(".*biological.*", "Biological Process", ignore.case = T, names(.))) %>%
    setNames(sub(".*cellular.*", "Cellular Component", ignore.case = T, names(.))) %>%
    setNames(sub(".*molecular.*", "Molecular Function", ignore.case = T, names(.))) %>%
    setNames(sub(".*keywords.*", "Keywords", ignore.case = T, names(.))) %>%
    select(all_of(intensity.cols), 
           all_of(major.prot.col), 
           all_of(grp.name.col),
           Minus.log.pval, Difference, Significant, `Biological Process`, `Cellular Component`, `Molecular Function`, Keywords)
  
  
  #Extract group labels from the "significant" column for use in legend
  sig.df <- prot.dataset[prot.dataset$Significant == "+",]  
  grp.names <- unique(sig.df[, grp.name.col]) #extract label 
  grp.names <- str_split_1(grp.names, "_") #separate labels
  
  #Add the group labels based on positive or negative x axis then calculate mean across groups
  prot.dataset <- prot.dataset %>%
    mutate(Expression = case_when(Difference >= 0 & Significant == "+" ~ paste0("Higher abundance in ", grp.names[1]),
                                  Difference <= 0 & Significant == "+" ~ paste0("Higher abundance in ", grp.names[2]),
                                  TRUE ~ "Not significant")) %>%
    mutate(Expression = factor(Expression, levels = c(
      paste0("Higher abundance in ", grp.names[2]),
      paste0("Higher abundance in ", grp.names[1]),
      "Not significant"))) %>%  #reorder
    mutate_at(c('Difference', 'Minus.log.pval', intensity.cols), as.numeric) %>%
    rowwise() %>%
    mutate(`Left-Group` = mean(
      c_across(all_of(left.range)))) %>%
    mutate(`Right-Group` = mean(
      c_across(all_of(right.range)))) %>%
    remove_rownames() %>%
    column_to_rownames(var = major.prot.col)
  
  
  return(prot.dataset)
}





####2. 1D HEATMAPS ----
#This function accepts the dataframe produced after 1D annotation enrichment in Perseus. The title argument is optional
onedheatmap <- function(oned.df, plot.title = "") {
  
  oned.df <- oned.df %>%
    mutate_at(c("Score"), as.numeric)
  annotations <- oned.df$Name
  samples <- sort(unique(oned.df$Column))
  
  #generate 1D score matrix
  M.score <- matrix(ncol=length(samples), nrow=length(annotations))
  for(i in 1:length(samples)){
    for(j in 1:length(annotations)){
      score.value <- oned.df$Score[oned.df$Column==samples[i] & oned.df$Name==annotations[j]]
      if(length(score.value)==0){score.value <- NA}
      M.score[j,i] <- score.value
    }
  }
  
  m <- M.score
  rownames(m) <- annotations
  colnames(m) <- samples
  #Assign zero to NA values.
  m[is.na(m)] <- 0
  #collapse matrix for use in ggplot
  m <- melt(m)
  colnames(m) <- c("Annotations", "T-test differences", "value")
  n.columns = length(unique(m$`T-test differences`))
  hm.height = length(m$`Annotations`)
  
  return(
    ggplot(data = m, aes(x = `T-test differences`, 
                         y = reorder(Annotations, value), 
                         fill = value))+
      geom_tile(colour = "grey", linewidth = 1)+
      scale_fill_gradientn(colors = c(low = "blue", mid = "white", high = "red"),
                           na.value = "grey")+
      guides(fill = guide_colourbar(title = "Colour Key"))+
      theme(axis.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = ifelse(
              n.columns <= 4, 0, 90 
            )),
            axis.text = element_text(colour = "black"),
            text = element_text(size = 14),
            legend.position = "right")+
      coord_fixed(ratio = 0.5)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      labs(title = plot.title, 
           y = "Annotations")
  )
  
}





####3. VOLCANO Plot ----
#The function requires a dataframe that has been processed by process.df() (step 1). It also requires the curves.df exported from Perseus to create the lines on the volcano plot. The other arguments are optional.
#go.terms = display GO terms for significant or non-significant proteins
#plot.title = custom plot title top left
#group.cols = colours for the default plot in the following order: higher-left, higher-right colours. 
#s0 = s0 value to be shown on the plot. Default = 0.1
#fdr = fdr value to be shown on the plot. Default = 0.05
#fdr.lines = 'yes' to show lines and 'no' to remove lines
#palette.col = based on hcl.colors() palettes. "Viridis" is the default one

volcano_plot <- function(df, curves.df, 
                         go.terms = "", which.go = "Keywords", plot.title = "", 
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
  num.sig.go <- length(unique(sig.df[, which.go])) 
  num.ns.go <- length(unique(ns.df[, which.go])) 
  select.pts <- df[select.pts, ]
  grp.name.col <- grep("significant", colnames(sig.df)) #column index 
  xlab <- unique(sig.df[, grp.name.col]) 
  xlab <- sub("_", " - ", xlab) #for x axis label 
  grp.names <- str_split_1(xlab, " - ") #for legend label
  
  
  #Custom settings for ggplots
  a = 0.2 #alpha
  s = 21  #shape
  ns.sz = 2 #non-sig shape size  
  ns.col = "grey" #non-sig colour
  s.sz = 4  #sig shape size
  text.sz = 18
  
  #Base volcano plot to plot points, fdr curves, axes boundaries, and other custom arguments.
  base_vplot <- ggplot(data = ns.df, 
                       aes(x = Difference,
                           y = Minus.log.pval))+
    geom_point(aes(fill = ns.col), size = ns.sz,
               alpha = a, shape = s)+
    labs(title = plot.title,
         caption = paste("s0 =", s0, "  ", "FDR =", fdr),
         x = paste0("Difference (", xlab, ")"), 
         y = expression(-log[10]("p-value")))+
    coord_cartesian(xlim = c(min(df$Difference),
                             max(df$Difference)),
                    ylim = c(min(df$Minus.log.pval),
                             max(df$Minus.log.pval)),
                    expand = T)+
    theme_minimal()+
    theme(plot.title = element_text(face="bold", size = 20),
          plot.caption = element_text(
            colour = "darkviolet", size = 11),
          axis.line = element_line(colour = "black"),
          text = element_text(size = text.sz))
  
  
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
               aes(fill = !!sym(which.go)), 
               size = s.sz,
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.sig.go, 
      palette.col), 
      na.value = alpha(ns.col, a)) + 
    label.pt.plot 
  
  #Add GO terms for non-sig proteins. Grey out sig proteins. 
  nonsig.go <- base_vplot +
    scale_fill_identity()+
    new_scale_fill()+
    geom_point(aes(fill = !!sym(which.go)),
               size = ns.sz, 
               shape = s)+
    scale_fill_manual(values = hcl.colors(
      num.ns.go,
      palette.col),
      na.value = alpha(ns.col, a)) + 
    geom_point(data = sig.df, fill = "grey",
               size = s.sz, alpha = a,
               shape = s) + 
    label.pt.plot
  
  ###Add layer showing average protein abundance across samples 
  vary.size.pl <- base_vplot +
    scale_fill_identity()+
    new_scale_fill()+
    geom_point(data = sig.df[sig.df$Difference > 0, ],
               aes(size = `Right-Group`),
               fill = alpha("#FF6666", a),
               shape = s)+
    scale_size_binned(range = c(1,8),
                      name = paste0("Average Intensity ", grp.names[1]),
                      n.breaks = 4)+
    new_scale(new_aes = "size")+
    geom_point(data = sig.df[sig.df$Difference < 0, ],
               aes(size = `Left-Group`), 
               fill = alpha("#00B1B2", a),
               shape = s)+
    scale_size_binned(range = c(1,8),
                      name = paste0("Average Intensity ", grp.names[2]),
                      n.breaks = 4) + 
    label.pt.plot 
  
  #Options: GO term plot for sig proteins, with fdr lines and without fdr lines
  if (go.terms == "significant" && fdr.lines == "yes") {
    return (sig.go + curve.plot)
    
  } else if (go.terms == "significant" && fdr.lines == "no") {
    return (sig.go)
    
    #Options: GO terms for non-significant proteins with fdr lines, without fdr lines, and default plot with color-coding of high and low expression proteins. 
  } else if (go.terms == "non-significant" && fdr.lines == "yes") {
    return (nonsig.go + curve.plot)
    
  } else if (go.terms == "non-significant" && fdr.lines == "no") {
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





####4. PCA plot ---- 

#This function also accepts the processed dataframe from step 1 and creates a PCA of all intensity columns.
#Serial ID are removed from group names so that ellipses are drawn by group


pca_plot <- function(df, ellipse = "yes") {
  counts <- df[, grep("intensity.+", colnames(df), value = T)]
  counts <- na.omit(counts)
  counts <- as.data.frame(t(counts)) 
  pca <- prcomp(counts) 
  counts$Group.names <- sub("\\.\\d+$", "", rownames(counts)) #remove serial numbers to extract group names
  text.sz = 16
  point.sz = 3
  scl = 0
  
  pca_ellipse <- autoplot(pca, data = counts, size = point.sz, 
                          scale = scl, colour = "Group.names", frame = T, frame.type = "norm") + 
    guides(colour=guide_legend("Sample Type"), fill = "none")+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          text = element_text(size = text.sz),
          panel.border = element_blank())
  
  
  pca_plain <- autoplot(pca, data = counts, size = point.sz, 
                        scale = scl, colour = "Group.names") + guides(colour=guide_legend("Sample Type"), fill = "none")+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          text = element_text(size = text.sz),
          panel.border = element_blank())
  
  
  if (ellipse == "yes") {
    pca_ellipse
    
  } else {
    pca_plain
  }
}


####5. S-curve plot ----
#The function accepts the same dataframe produced after step 1 but only produces S-curves for the group pair (e.g WT vs MT) in the volcano plot
#The function create one s-curve plot for each group and puts them on the same image with ggarrange

scurve <- function(df) {
  ranked.df <- df %>%
    mutate(Left.rank = rank(desc(`Left-Group`)))%>%
    mutate(Right.rank = rank(desc(`Right-Group`)))
  
  sig.df <- df[df$Significant == "+", ]
  grp.names <- unique(
    sig.df[, grep("significant", colnames(sig.df))] ) 
  grp.names <- str_split_1(grp.names, "_")  
  text.sz = 16
  
  left.gr.pl <- ggplot(data = ranked.df,
                       aes(x = Left.rank, y = `Left-Group`))+
    geom_point(size = 2, colour = "royalblue")+
    labs(x = paste(grp.names[2], "proteins ranked by iBAQ Intensity"),
         y = bquote("iBAQ Intensity for " ~ .(grp.names[2]) * " " * (Log[10])))+
    theme_minimal()+
    theme(axis.line = element_line(colour = "black"),
          text = element_text(size = text.sz))
  
  
  
  right.gr.pl <- ggplot(data = ranked.df,
                        aes(x = Right.rank, y = `Right-Group`))+
    geom_point(size = 2, colour = "royalblue")+
    labs(x = paste(grp.names[1], "proteins ranked by iBAQ Intensity"),
         y = bquote("iBAQ Intensity for " ~ .(grp.names[1]) * " " * (Log[10])))+
    theme_minimal()+
    theme(axis.line = element_line(colour = "black"),
          text = element_text(size = text.sz))
  
  
  
  return(
    #plot with a tiny column in between to serve as a larger gap
    ggarrange(left.gr.pl, NULL, right.gr.pl, 
              nrow = 1, 
              widths = c(1, 0.05, 1))
  )
}  




####6. Venn Diagram ----
#This function takes the matrix from Perseus and replaces the column names with group names. The function returns a narrowed down dataframe with intensity and protein ID columns only
#add options for group names e,g experiment, type, etc**
groupnames_to_colnames <- function(venn.df) {
  intensity.cols <- grep("intensity.+", colnames(venn.df), ignore.case = T) #extract indexes of intensity columns 
  grp.names.row <- grep("group", venn.df[ ,1], ignore.case = T) #index of the row containing group-names
  grp.names <- venn.df[grp.names.row, intensity.cols] #obtain group names by row index
  colnames(venn.df)[intensity.cols] <- sub(".*}", "", grp.names, perl = T) #replace the intensity column headers with the group names
  comment_char <- grep("^#", venn.df[, 1]) #index of non-numeric rows
  venn.df <- venn.df[-(comment_char), ] #remove non-numeric rows
  venn.df <- type.convert(venn.df, as.is = T) #ensure that intensities are "numeric"
  id.col <- grep("majority", names(venn.df), ignore.case = T) #id column
  venn.df <- venn.df[, c(intensity.cols, id.col)]   #keep intensity + protein ID columns
}

#This functions takes the dataframe from the above step and evaluates the presence/absence of proteins based on minimum valid values
evaluate_valid_vals <- function(venn.df, group.names, min.percent = 75) {
  min.percent = min.percent/100
  df_list <- lapply(group.names, function(x) 
    select(venn.df, contains(x))
  ) #split the main df into the different groups
  binary.res <- lapply(df_list, \(x) ifelse(
    rowSums(!is.na(x))/ncol(x) >= min.percent,
    1, 0))  #binary results: If valid values > minimum percentage in each df/group, the protein is present and evaluates to 1. Else, it evaluates to 0 for absent proteins
  names(binary.res) <- paste0(unlist(group.names), " ")  #add the extracted labels to the binary results 
  #Insert binary results in venn.df
  merged.df <- data.frame(df_list, binary.res, check.names = F)
  
  #Return combined df or binary df only
  output <- list(binary.result = binary.res,
                 merged.df = merged.df)
  return(output)
}


####7. SHINY ----
#It accepts the protein matrix from Perseus and processes it using step 1 for the remaining program. It also accepts the 1D matrix from Perseus to create 1D heatmaps and to filter out less significant GO terms on the volcano plot. Thirdly, it accepts the curve matrix from Perseus for the lines on the volcano plot 
#A Shiny app's function has two sections as seen below: the UI and the Server

dl.unit = "in"
#Shiny UI ---- 
ui <- shinyUI(
  fluidPage(
    titlePanel("Visualization for Proteomics"),
    
    tabsetPanel(
      #Main Tab for Volcano plot, Scurve, and PCA
      tabPanel("Main", 
               sidebarLayout(
                 sidebarPanel(width = 4,
                              fluidRow(
                                column(width = 6,
                                       fileInput(inputId = "proteinfile", label = "Protein groups") ),
                                column(width = 6, 
                                       fileInput("curvesfile", label = "Threshold curve") ),
                                column(width = 12, selectInput("left.gr", 
                                                               label = "Column range for Left group",
                                                               choices = c(),
                                                               multiple = T)),
                                column(width = 12, selectInput("right.gr", 
                                                               label = "Column range for Right group",
                                                               choices = c(),
                                                               multiple = T)),
                                column(width = 12, br()),
                                column(width = 6, 
                                       actionButton("gen.pca", label = "Generate PCA plot")),
                                column(width = 6, radioButtons("ellipse", 
                                                               label = "PCA with ellipse",
                                                               choices = list("yes", "no"),
                                                               selected = "yes",
                                                               inline = T) ),
                                column(width = 6, 
                                       actionButton("gen.scurve", label = "Generate S-curve plot")),
                                column(width = 12, br())
                              ),
                              
                              useShinyjs(),
                              fluidRow(
                                div(id = "reset.grp",
                                    h4("Plot Options"), 
                                    column(width = 3, numericInput("xmin",
                                                                   label = "x min", 
                                                                   value = NA) ),
                                    column(width = 3, numericInput("xmax", 
                                                                   label = "x max", 
                                                                   value = NA) ),
                                    column(width = 3, numericInput("ymin", 
                                                                   label = "y min", 
                                                                   value = NA) ),
                                    column(width = 3, numericInput("ymax", 
                                                                   label = "y max", 
                                                                   value = NA) ),
                                    
                                    column(width = 12, actionButton("reset", label = "Reset Volcano Plot")),
                                    column(width = 12, br()),
                                    column(width = 6, colourpicker::colourInput("left.col",
                                                                                label = "Left group's colour",
                                                                                value = "#00B1B2") ),
                                    column(width = 6, colourpicker::colourInput("right.col",
                                                                                label = "Right group's colour",
                                                                                value = "#FF6666") ),
                                    column(width = 12, br()),
                                    column(width = 6, numericInput("s.knot",
                                                                   label = "s0 value",
                                                                   value = 0.1) ),
                                    column(width = 6, numericInput("fdr.val",
                                                                   label = "FDR value",
                                                                   value = 0.05) ),
                                    column(width = 12, 
                                           helpText("Note: Changing the s0 and FDR values does not change the plot. The caption only indicates the parameters used for creating the plot.") ),
                                    column(width = 12, br()),
                                    column(width = 6, radioButtons("fdr.lines",
                                                                   label = "Display curves",
                                                                   choices = list("yes", "no"),
                                                                   selected = "yes", inline = T) ),
                                    column(width = 6, radioButtons("pt.sizes",
                                                                   label = "Average protein intensities",
                                                                   choices = list("yes", "no"),
                                                                   selected = "no", inline = T) ),
                                    column(width = 12, radioButtons("go.term",
                                                                    label = "Gene Ontology terms",
                                                                    choices = list("significant", "non-significant", "reset"),
                                                                    selected = "reset", inline = T) ),
                                    column(width = 7, selectInput("go.types",
                                                                  label = "Gene Ontology annotation category for volcano plot",
                                                                  choices = list("Keywords", "Cellular Component", "Molecular Function", "Biological Process"),
                                                                  selected = "Keywords", 
                                                                  width = '100%') ),
                                    column(width = 5, selectInput("palette",
                                                                  label = "GO terms color palette",
                                                                  choices = list("Viridis", "Plasma", "Inferno", "Rocket", "Lajolla", "Turku", "Hawaii", "Batlow", "Spectral", "Blue-Red", "Green-Orange", "RdYlBu", "Zissou 1", "Roma"),
                                                                  selected = "Viridis") )
                                )
                              )
                 ),
                 mainPanel(
                   img(src = "proteovis.png"),
                   br(),
                   
                   h4("Volcano Plot",  align = "center"),
                   ##Volcano plot output with brush selection of points
                   plotOutput("volcanoplot",
                              brush = "vplot.brush"),
                   textInput("vplot.title",
                             label = "Volcano Plot Title (Optional)",
                             value = NULL,
                             placeholder = "Enter text..."),
                   fluidRow(
                     ##initialize drop down list of protein labels
                     column(width = 6, selectizeInput("protein.labs",
                                                      label = "Manually select protein IDs:",
                                                      choices = list(),
                                                      multiple = T)),
                     column(width = 6, verbatimTextOutput("brushed.ids")),
                     column(width = 6, downloadBttn("dlvplot", "Download Volcano plot", size = "xs"))
                   ),
                   br(),
                   br(),
                   
                   ##S-curve plot output
                   h4("S-curve Plots",  align = "center"),
                   plotOutput("scurveplot"),
                   downloadBttn("dlscurve", "S-curves", size = "xs"),
                   br(),
                   br(),
                   
                   #PCA plot output
                   h4("PCA Plot", align = "center", ),
                   plotOutput("pcaplot"),
                   downloadBttn("dlpca", "PCA", size = "xs")
                 )
               )
      ),
      
      #Panel for 1D Heatmap
      tabPanel("1D Annotation Heatmap", 
               sidebarLayout(
                 sidebarPanel(width = 4,
                              fluidRow(
                                column(width = 12,
                                       fileInput("onedfile", label = "1D annotation enrichment file")),
                                p("Dimensions for download (inches):"),
                                column(width = 4, numericInput("width.val",
                                                               label = "Width ", 
                                                               value = 9, min = 5, max = 25)),
                                column(width = 4, numericInput("height.val",
                                                               label = "Height", 
                                                               value = 9, min = 5, max = 25)),
                                # column(width = 6, radioButtons("hmtype",
                                #                                label = "Sort by Gene Ontology groups",
                                #                                choices = list("yes", "no"),
                                #                                selected = "yes", inline = T) )
                              )
                 ),
                 mainPanel(
                   br(),
                   br(),
                   ##1D plots output
                   h4("1D Annotation Heatmaps", align = "center", ),
                   
                   plotOutput("onedhm1.bio.proc"),
                   downloadBttn("dlhm1", "Heatmap 1", size = "xs"),
                   br(),
                   br(),
                   plotOutput("onedhm2.cell.comp"),
                   downloadBttn("dlhm2", "Heatmap 2", size = "xs"),
                   br(),
                   br(),
                   plotOutput("onedhm3.mol.func"),
                   downloadBttn("dlhm3", "Heatmap 3", size = "xs"),
                   br(),
                   br(),
                   plotOutput("onedhm4.keywords"),
                   downloadBttn("dlhm4", "Heatmap 4", size = "xs"), 
                 )
               )
      ),
      
      #Panel for Venn diagram
      tabPanel("Venn Diagram", 
               sidebarLayout(
                 sidebarPanel(width = 4,
                              fluidRow(
                                column(width = 12,
                                       fileInput("vennfile", label = "Venn Diagram Matrix")),
                                column(width = 6, selectInput("venngroups", 
                                                            label = "Select groups:", 
                                                            choices = c(), multiple = T)),
                                bsTooltip(id = "venngroups", title = "Select groups to be compared."),
                                column(width = 6, numericInput("min.valid.vals",
                                                               label = "Minimum % valid values",
                                                               min = 0, max = 100, value = 75)),
                                column(width = 6, actionButton("eval.valid.vals", label = "Evaluate valid values")),
                                column(width = 6, actionButton("gen.venn", label = "Generate Venn diagram")),
                                column(width = 12, br()),
                                column(width = 12, uiOutput("venn.colour.picker")),
                                column(width = 3, selectInput("linetype", 
                                                             label = "Stroke (type)",
                                                             choices = c("solid", "dashed", "dotted", "longdash"),
                                                             selected = "solid")),
                                column(width = 3, numericInput("linewidth", 
                                                              label = "Stroke (width)",
                                                              min = 0, max = 10, value = 1)),
                                column(width = 3, numericInput("labelsize", 
                                                               label = "Size (labels)",
                                                               min = 4, max = 25, value = 12)),
                                column(width = 3, numericInput("numsize", 
                                                               label = "Size (values)",
                                                               min = 3, max = 20, value = 10)),
                                p("Dimensions for download (inches):"),
                                column(width = 4, numericInput("width.val",
                                                               label = "Width ", 
                                                               value = 9, min = 5, max = 25)),
                                column(width = 4, numericInput("height.val",
                                                               label = "Height", 
                                                               value = 9, min = 5, max = 25))
                              )
                              
                 ),
                 mainPanel(
                   tags$head(
                     tags$style(
                       HTML("#venndf\\.preview 
                                     {font-size: 15px; color: red; font-family: 'Courier New';}")
                     )),
                   htmlOutput("venndf.preview", align="center"),
                   br(),
                   DTOutput("venntable"),
                   br(), br(), br(), br(), br(),
                   plotOutput("vennplot", height = '300px'),
                   downloadBttn("dlvennplot", "Venn plot", size = "xs")
                 )
               )
      )
    )
  )
)




#Shiny Server ----
server <- function(input, output, session) {
  
  #Curves file
  df.curves <- eventReactive(input$curvesfile, {
    read.table(input$curvesfile$datapath, header = T)
  })
  
  #Protein input file
  df.prot <- eventReactive(input$proteinfile, {
    read.delim(input$proteinfile$datapath, na.strings = c("NA", "NaN"))
  })
  
  #1D annotation file
  df.1d <- eventReactive(input$onedfile, {
    df <- read.delim(input$onedfile$datapath, stringsAsFactors = F, header=T)
    df <- df[-c(1),]
    df <- df %>% mutate(Column = sub(".*Difference ", "", Column))
  })
  
  #Update column selection on the UI
  observe({
    req(df.prot())
    updateSelectInput(session,
                      inputId = "left.gr",
                      choices = colnames(df.prot()))
    updateSelectInput(session,
                      inputId = "right.gr",
                      choices = colnames(df.prot()))
    output$enter.range <- renderUI({
      helpText("Select columns for each group to calculate average protein intensities across groups")
    })
  })
  
  #Process the protein dataframe and calculate means across group columns
  df.prot2 <- reactive({
    req(df.prot(), input$left.gr, input$right.gr)
    df <- process.df(prot.dataset = df.prot(),
                     left.range = input$left.gr,
                     right.range = input$right.gr )
    
    if (is.null(input$onedfile)) {
      df
      
    } else {
      # Retain only those GO terms found in the 1D Annotation file
      df$Keywords[-which(df$Keywords %in% df.1d()$Name)] <- NA
      df$`Cellular Component`[-which(df$`Cellular Component` %in% df.1d()$Name)] <- NA
      df$`Molecular Function`[-which(df$`Molecular Function` %in% df.1d()$Name)] <- NA
      df$`Biological Process`[-which(df$`Biological Process` %in% df.1d()$Name)] <- NA
      
      df
    }    
  })
  
  observe({
    req(df.prot2())
    #fill dropdown list with protein names (labels) from which to choose
    updateSelectizeInput(session,
                         inputId = "protein.labs",
                         choices = rownames(df.prot2()), 
                         server = T )
  })
  
  #activate brush selection    
  brushed.pts <- reactive({
    brushedPoints(df.prot2(), input$vplot.brush)
  })
  #display selected ids as text on the UI
  output$brushed.ids <- renderPrint({ 
    cat("Selected proteins: ", rownames(brushed.pts()) ) 
  }) 
  
  
  
  ####create volcano plot using UI input variables  
  vplot <- reactive({
    vplot <- volcano_plot(df = df.prot2(),
                          curves.df = df.curves(),
                          go.terms =  input$go.term,
                          palette.col = input$palette,
                          plot.title = input$vplot.title,
                          which.go = input$go.types,
                          s0 = input$s.knot,
                          fdr = input$fdr.val,
                          group.cols = c(input$left.col, input$right.col),
                          fdr.lines = input$fdr.lines,
                          vary.sizes = input$pt.sizes,
                          #show protein labels either by quick or fixed selection 
                          select.pts = c(rownames(brushed.pts()), 
                                         input$protein.labs)
    )
    
    if (is.na(input$xmin) && is.na(input$xmax) && is.na(input$ymin) && is.na(input$ymax)) {
      vplot
      
    } else {
      #Axes limits adjustment
      vplot+
        coord_cartesian(xlim = c(input$xmin, input$xmax),
                        ylim = c(input$ymin, input$ymax),
                        expand = T)
    }
  })
  #send vplot to the UI
  output$volcanoplot <- renderPlot({
    vplot()
  })
  #Reset to default plot
  observeEvent(input$reset, {
    shinyjs::reset("reset.grp")
  })
  output$dlvplot <- downloadHandler( 
    filename = "Volcano_plot.png",
    content = function(file) {
      ggsave(file, plot = vplot(),
             width = input$width.val, height = input$height.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  
  
  
  ####PCA plot
  pca.pl <- reactive({req(input$gen.pca) 
    pca_plot( df.prot2(), ellipse = input$ellipse )
  })
  output$pcaplot <- renderPlot({
    pca.pl()
  })
  output$dlpca <- downloadHandler( 
    filename = "PCA.png",
    content = function(file) {
      ggsave(file, plot = pca.pl(),
             width = input$width.val, height = input$height.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  
  
  
  
  ####S-curve plot 
  scurve.pl <- reactive({
    req(input$gen.scurve)
    scurve(df.prot2()) 
  })
  output$scurveplot <- renderPlot({
    scurve.pl() 
  })
  #download
  output$dlscurve <- downloadHandler( 
    filename = function() { paste("scurve.png") },
    content = function(file) {
      ggsave(file, plot = scurve.pl(),
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  
  
  
  
  ####Heatmap plots
  #create heatmap list
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
  bio.proc <- reactive(plotlist.1d()[[1]])
  cell.comp <- reactive(plotlist.1d()[[2]])
  mol.func <- reactive(plotlist.1d()[[3]])
  keywords <- reactive(plotlist.1d()[[4]])
  output$onedhm1.bio.proc <- renderPlot({ 
    bio.proc()
  }) 
  output$onedhm2.cell.comp <- renderPlot({ 
    cell.comp()
  }) 
  output$onedhm3.mol.func <- renderPlot({ 
    mol.func()
  }) 
  output$onedhm4.keywords <- renderPlot({ 
    keywords()
  })
  
  #downloads
  output$dlhm1 <- downloadHandler( 
    filename = "Bio_process.png",
    content = function(file) {
      ggsave(file, plot = bio.proc(), 
             width = input$width.val, height = input$height.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  output$dlhm2 <- downloadHandler( 
    filename = "Cell_component.png",
    content = function(file) {
      ggsave(file, plot = cell.comp(),
             width = input$width.val, height = input$height.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  output$dlhm3 <- downloadHandler( 
    filename = "Mol_function.png",
    content = function(file) {
      ggsave(file, plot = mol.func(),
             width = input$width.val, height = input$height.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  output$dlhm4 <- downloadHandler( 
    filename = "Keywords.png",
    content = function(file) {
      ggsave(file, plot = keywords(), 
             width = input$width.val, height = input$height.val,
             device = "png", bg = "white", units = dl.unit)
    }
  )
  
  
  
  ####Venn Diagram
  #Read venn matrix 
  df.venn <- eventReactive(input$vennfile, {
    df <- read.delim(input$vennfile$datapath, check.names = F)
    df <- groupnames_to_colnames(df)
  })

  #Evaluate valid values using the column indexes provided as input i.e "venngroups"
  observe({req(df.venn())
    venn.groups <- select(df.venn(), !contains("majority"))
    venn.groups <- sub("\\..*", "", names(venn.groups))
    venn.groups <- unique(venn.groups)
    updateSelectInput(session,
                      inputId = "venngroups",
                      choices = venn.groups,
                      selected = venn.groups)
  })

  #Evaluate valid values
  venn.input <- reactive({req(df.venn(), input$eval.valid.vals, input$venngroups)
    evaluate_valid_vals(venn.df = df.venn(), 
                        group.names = input$venngroups,
                        min.percent = input$min.valid.vals)
  })

  #Extract venn group names and generate dynamic colour pickers for each venn group
  
  
  observe({req(input$venngroups, input$gen.venn)
    output$venn.colour.picker <- renderUI({
      lapply(input$venngroups, function(x)
        column(width = 4, colourpicker::colourInput(inputId = x,
                                                    label = paste0(x), value = "grey", 
                                                    allowTransparent = T)
        )
      )
    })
  })
  
  #Input ids for colour pickers
  colour.input.ids <- reactive({req(input$venngroups) 
    sapply(input$venngroups, function(x) { 
      input[[x]] 
    }) 
  })
  
  #Create venn diagram with binary values and colourinput variables
  venn.pl <- reactive({req(input$gen.venn, venn.input())
    venn.fit <- lapply(venn.input()[["binary.result"]],
                     \(x) which(x==1))
    venn.fit <- eulerr::venn(venn.fit)
    plot(venn.fit, 
         quantities = list(type = c("counts", "percent"),
                           fontsize = input$numsize),
         fills= list(fill = colour.input.ids()),
         edges = list(lwd = input$linewidth, lty = input$linetype),
         labels = list(fontsize = input$labelsize))
  })
  #Generate Venn diagram and add download handler
  output$vennplot <- renderPlot({ venn.pl() })
  
  output$dlvennplot <- downloadHandler(
    filename = "Venn.png",
    content = function(file) {
      ggsave(file, plot = venn.pl(),
             width = input$width.val, height = input$height.val,
             device = "png", bg = 'white', units = dl.unit)
    }
  )
  

  # if (!interactive()) {
  #   session$onSessionEnded(function() {
  #     stopApp()
  #     q("no")
  #   })
  # }
    
}

####RUN PROTEOMICS APP -----
shinyApp(ui, server)
