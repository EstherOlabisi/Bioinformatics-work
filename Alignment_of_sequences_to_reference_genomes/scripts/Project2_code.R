####Project 2
####Alignment of burbot sequences to reference genomes
####February 20, 2023


####RCODE####
library(tidyverse)
library(reshape2)
library(ggplot2)

#Create a function that reads the txt file as a table and adds a column containing % of raw data aligned
create.table <- function(txt_file) {
  txt_table <- read.table(txt_file, header = T)
  txt_table$percent_aligned <- (txt_table$assembled/txt_table$raw) *100
  txt_table$percent_aligned <- round(txt_table$percent_aligned, digits = 2)
  return(txt_table)
}

BWA_to_Bur <- create.table("b_bwa_assembled_per_ind.txt")#bwa alignment to burbot genome
BWA_to_Cod <- create.table("c_bwa_assembled_per_ind.txt") #bwa alignment to cod genome
BT2_to_Bur <- create.table("b_bt2_assembled_per_ind.txt") #bowtie2 alignment to burbot genome
BT2_to_Cod <- create.table("c_bt2_assembled_per_ind.txt") ##bowtie2 alignment to cod genome
  

#I created a table starting with BWA_to_Bur, selected the ID and percent-aligned columns to be kept, and added % aligned columns from the other tables - each with their respective column titles.

Merged_table <- BWA_to_Bur %>%
  select(ind, BWA_to_Burbot=percent_aligned) %>%
  add_column(BWA_to_Cod= BWA_to_Cod$percent_aligned)  %>%
  add_column(BT2_to_Burbot = BT2_to_Bur$percent_aligned) %>%
  add_column(BT2_to_Cod = BT2_to_Cod$percent_aligned) 

#view
Merged_table

#remove first copies
rm(BWA_to_Bur, BWA_to_Cod, BT2_to_Bur, BT2_to_Cod)

#A function to create bar plots
create_bar_plot <- function(column1, column2, column_for_x_axis, colour_vector, plot_title){
  joined_data <- t(cbind(column1, column2))
  y <- as.matrix(joined_data)
  x <- barplot(joined_data, ylim=c(0,100), beside=T,
               names.arg=column_for_x_axis, cex.names=0.5, col=colour_vector,
               xlab="Burbot IDs", ylab="Raw data aligned (%)",
               main=plot_title)
  legend(x="top", legend = c("Burbot ref genome", "Cod ref genome"),
         fill = colour_vector, cex = 0.7,
         horiz=T, bty="n", , inset=-0.3)
  text(x, y+2, labels=y, cex=0.7)
}


#Mapped reads: Alignment of raw burbot data to burbot and cod reference genomes with Burrows-Wheeler Aligner (BWA)
create_bar_plot(Merged_table$BWA_to_Bur, Merged_table$BWA_to_Cod, Merged_table$ind, c("brown2", "deeppink3"), plot_title = "Burrows-Wheeler Aligner (BWA) alignment rates by burbot individuals")

create_bar_plot(Merged_table$BT2_to_Bur, Merged_table$BT2_to_Cod, Merged_table$ind, c("cyan4", "darkgoldenrod2"), plot_title = "Bowtie2 alignment rates by burbot individuals")



# Function to make a modified dataframe and box plot
Box_plot <- function(df,column_names_vector){
  modified <- melt(df, id.vars=NULL, measure.vars = column_names_vector)
  colnames(modified) <- c("Key" ,"Raw data aligned (%)")
  
  # creating a plot
  ggplot(modified, aes(x=NULL, y=`Raw data aligned (%)`, fill=Key)) +
    
    geom_boxplot(linewidth=1, outlier.size = 3,
                 outlier.colour="red") +
    
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom", legend.title=element_blank(),
          legend.text = element_text(size=15))
}

Box_plot(Merged_table, c('BWA_to_Burbot','BT2_to_Burbot'))
Box_plot(Merged_table, c('BWA_to_Cod','BT2_to_Cod'))


  
  