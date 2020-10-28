####Info####
"This file describes the work process for obtaining the SNP-distribution 
 (Chapter 6.1 \"SNP-distribution analysis along chromosome 25 of two cattle populations\"), 
 which includes:
 
1. Running VCF-tools on the orgininal VCF-file.
2. Loading the data
3. Transforming
4. Visualizing in R
"
####Setting working directory and loading library#####
setwd("~/SNP-Distribution")
library(ggplot2)
library(gridExtra)
library(evobiR)
library(devtools)

####Running VCFTools####
#VCFTools:
"To obtain the density of SNPs, VCFTools has to be run in the terminal via:"
#vcftools --vcf braunvieh.vcf --SNPdensity 4000 --out SNP_freqs_25_4k
#vcftools --vcf fleckvieh.vcf --SNPdensity 4000 --out SNP_freqs_25_4k

####Functions####
#Function to plot distribution
snp_frequency_fnc <- function(snpdens_file, what_species, what_subset){ 
  snpdens <- read.table(snpdens_file, header = T)
  barfreqs <- as.data.frame(cbind(as.numeric(names(table(snpdens$SNP_COUNT)[1:120])),table(snpdens$SNP_COUNT)[1:120]))
  colnames(barfreqs) <- c("Segment", "SNP_Count")
  p_snp_frequency <- ggplot2::ggplot(barfreqs, aes(x=Segment, y=SNP_Count)) +
    ggtitle(paste0(what_species," - ", what_subset, "\nSNP distribution of chromosome 25 per 4kb segment"))+
    geom_bar(stat = "identity")+
    xlab("SNP count per 4kb segment") +
    ylab("Frequency")+
    scale_y_continuous(breaks = seq(0,300,by=25))+
    scale_x_continuous(breaks = seq(0, 120, by=20))+
    coord_cartesian(ylim=c(0, 303))+
    annotate("label", x = 105, y = 300, label = paste0("Total SNP count:\n", sum(snpdens$SNP_COUNT)))+
    theme(axis.text = element_text(size = 8),
          plot.title=element_text(hjust=0.5))
  return(p_snp_frequency)
}

#Braunvieh: 338,122 SNPs
#Fleckvieh: 428,439 SNPs

####Run function, plot results####
bv_all_freq <- snp_frequency_fnc(snpdens_file = "bv_all_SNP_freqs_4k.snpden", what_species = "Braunvieh" , what_subset = "All")
fv_all_freq <- snp_frequency_fnc(snpdens_file = "fv_all_SNP_freqs_4k.snpden", what_species = "Fleckvieh" , what_subset = "All")

prow <- cowplot::plot_grid(bv_all_freq,
                           fv_all_freq,
                           labels = c("A", "B"),
                           hjust = -1,
                           ncol = 2)

####Save plot####
ggsave(
  filename="SNP-Frequency-BV-FV.png",
  plot = prow,
  device = "png",
  path = NULL,
  scale = 1,
  width = 32 ,
  height = 16,
  units = c("cm"),
  dpi = 300,
  limitsize = TRUE,
)










