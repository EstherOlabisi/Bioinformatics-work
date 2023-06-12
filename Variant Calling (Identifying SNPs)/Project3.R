#Project3: March 28, 2023

############Compare the number of SNPs and the overlap -----#############
library("VennDiagram")
library("gridExtra")
library("ggplot2")
library("tidyverse")

#Import table containing number of SNPs. Add row names according to the descriptions in README.txt.  
no_of_SNPs <- read.table("compared_SNPs.txt", header = F)
row.names(no_of_SNPs) <- c("BCFtools", "Freebayes", "Shared1", "Shared2")

#Manually typing in the values is more efficient in this case than selecting table cells.
BCFtools <-  3998 + 8997
Freebayes <- 37725 + 8997
Shared <- 8997

#Draw venn diagram with number of shared and unique number of SNPs. The option to represent the plot in percentages is not working as expected (see commented line of code) so I went with raw values, unfortunately. 
grid.newpage()
venn <- draw.pairwise.venn(area1=BCFtools, area2=Freebayes, cross.area=Shared, category=c("BCFtools","Freebayes"),fill=c("blue","green"), col='transparent')
#draw.pairwise.venn(area1=BCFtools, area2=Freebayes, cross.area=Shared, category=c("BCFtools","Freebayes"),fill=c("blue","green"), col='transparent', print.mode = 'percent')

#Add title
grid.arrange(gTree(children = venn), top="Venn diagram: Number of SNPs detected by Freebayes and BCFtools mpileup")


############Plot the minor allele frequency ---- ############

#import allele frequency tables and add column names
bcf_allele_freq<- read.table("bcf_freq_bi.txt", header = F)
colnames(bcf_allele_freq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_ALLELE_FREQ", "ALT_ALLELE_FREQ")

freebayes_allele_freq <- read.table("freebayes_freq_bi.txt", header = F)
colnames(freebayes_allele_freq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_ALLELE_FREQ", "ALT_ALLELE_FREQ")

#The first frequency column will be plotted for each VCF. Edit the first frequency column to remove allele alphabets and colons, leaving numbers only. Convert values from 'character' class to 'numeric.'
bcf_allele_freq$REF_ALLELE_FREQ <- as.numeric(gsub("[AGCT]+:", "", bcf_allele_freq$REF_ALLELE_FREQ))
freebayes_allele_freq$REF_ALLELE_FREQ <- as.numeric(gsub("[AGCT]+:", "", freebayes_allele_freq$REF_ALLELE_FREQ))

#histogram for bcftools allele frequency
ggplot(bcf_allele_freq, aes(x=REF_ALLELE_FREQ))+
  geom_histogram(aes(y=..count..), color="black", fill="darkblue", alpha=0.6)

#histogram for freebayes allele frequency
ggplot(freebayes_allele_freq, aes(x=REF_ALLELE_FREQ))+
  geom_histogram(aes(y=..count..), color="black", fill="darkgreen", alpha=0.7)



#############Plot the depth of all SNPs (summed across individuals)----#####
#import tables containing depths. Add a new column to each table indicating the source of the results. This distinction will be important for creating a histogram in the next step. Finally, merge both tables into one dataframe. 
bcf_snp_dp <- read.table("all_BCF_SNP_dp.txt", header = T)
bcf_snp_dp$variant_caller <- "BCFtools"

freebayes_snp_dp <- read.table("all_freebayes_SNP_dp.txt", header = T)
freebayes_snp_dp$variant_caller <- "Freebayes"

merged_snp_dp <- merge.data.frame(bcf_snp_dp, freebayes_snp_dp, all=T)

#Plot a histogram of SNP depth results for freebayes and bcftools mpileup. Red represents BCFtools while yellow represents freebayes.
ggplot(merged_snp_dp, aes(x=SUM_DEPTH))+
  geom_histogram(data=subset(merged_snp_dp, variant_caller=="BCFtools"),
                 fill="red", position="identity",
                 color="black")+
  geom_histogram(data=subset(merged_snp_dp, variant_caller=="Freebayes"),
                 fill="yellow", alpha=0.4,
                 position="identity", color="black")+
  xlab("Depth of SNPs per locus")+
  ggtitle("Depth of SNPs per locus for Freebayes and BCFtools mpileup results")+
  theme(plot.title = element_text(hjust = 0.5)) 


###########Mean depth for each VCF (reads per locus per individual)----
#Import files, add variant caller names, and merge tables as in the previous steps. I also used gsub to shorten the burbot IDs so that the IDs fit in the axis as labels.
mean_dp_bcf <- read.table("mean_dp_ind_BCF.txt", header = T)
mean_dp_bcf$variant_caller <- "BCFtools"

mean_dp_freebayes <- read.table("mean_dp_ind_freebayes.txt", header = T)
mean_dp_freebayes$variant_caller <- "Freebayes"

combined_mean_dp <- merge.data.frame(mean_dp_bcf, mean_dp_freebayes, all=T)
combined_mean_dp$INDV <- gsub(".sorted.bam", "", combined_mean_dp$INDV)


#Plot bar graphs for mean depth of coverage.
ggplot(combined_mean_dp, aes(x=INDV, y=MEAN_DEPTH, fill=variant_caller))+
  geom_bar(stat="identity", position="dodge", color="black")+
  xlab("Burbot individuals")+ ylab("Mean depth")+
  ggtitle("Mean depth of coverage per individual: Freebayes and BCFtools mpileup results")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 7))

