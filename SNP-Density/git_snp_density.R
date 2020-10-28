####Info####
"This file describes the work process for obtaining the SNP-density 
 (Chapter 6.2 \"Identification of the highest and lowest SNP-density regions of chromosome
25 for two cattle populations \")
 
1. Running VCF-tools on the orgininal VCF-file.
2. Loading the data
3. Transforming
4. Visualizing in R
"

####Setting working directory and loading library#####
setwd("~/SNP-Density")
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
snp_density_fnc <- function(snpdens_file, what_species, what_subset){
  snpdens <- read.table(snpdens_file, header = T)
  size_window <- 500
  window_size2mb <- evobiR::SlidingWindow(FUN= sum, data = snpdens$SNP_COUNT, window= size_window, step = 1)
  ourmax2mb <- which(window_size2mb==max(window_size2mb)) 
  ourmin2mb <- which(window_size2mb==min(window_size2mb)) 
  
  print(paste0("Average SNP-density: ", mean(snpdens$SNP_COUNT)))
  
  print(paste0("Maximum: ", ourmax2mb))
  print(paste0("Start: ", snpdens$BIN_START[ourmax2mb]))
  print(paste0("End: ", snpdens$BIN_START[ourmax2mb+size_window-1]+4000))
  
  
  print(paste0("Minimum: ", ourmin2mb))
  print(paste0("Start: ", snpdens$BIN_START[ourmin2mb]))
  print(paste0("End: ", snpdens$BIN_START[ourmin2mb+size_window-1]+4000))
  
  p_snp_density <- ggplot2::ggplot(snpdens[ourmin2mb:(ourmax2mb+499),], aes(x=BIN_START, y=SNP_COUNT, fill = SNP_COUNT)) +
    ggtitle(paste0(what_species," - ", what_subset, "\nHighest and lowest SNP-density regions"))+
    geom_vline(xintercept = c(snpdens$BIN_START[ourmax2mb], snpdens$BIN_START[ourmax2mb+499]), linetype="longdash")+
    geom_vline(xintercept = c(snpdens$BIN_START[ourmin2mb],snpdens$BIN_START[ourmin2mb+499]), linetype="longdash")+
    coord_cartesian(ylim=c(0, 175))+
    scale_x_continuous(breaks = seq(0, 14e6, by=1e6))+
    scale_y_continuous(breaks = seq(0,175,by=25))+
    geom_bar(stat = "identity")+
    scale_fill_gradient(low="blue", high="red", limits = c(0,160))+
    #  scale_color_gradientn(colours = rainbow(2), limits=c(0,7000000000), labels = comma)
    xlab("Chromosome 25") +
    ylab("SNP count per 4kb segment") +
    labs(fill = "SNP Count")+
    theme(axis.text = element_text(size = 8),
          plot.title=element_text(hjust=0.5)) +
    annotate("label", x=(snpdens$BIN_START[ourmin2mb] + 
                           (snpdens$BIN_START[ourmin2mb+499] - snpdens$BIN_START[ourmin2mb])/2),
             y=160, label="LDR", size = 8) +
    annotate("label", x=(snpdens$BIN_START[ourmax2mb] +
                           (snpdens$BIN_START[ourmax2mb+499] - snpdens$BIN_START[ourmax2mb])/2),
             y=160, label="HDR", size = 8) +
    geom_hline(yintercept = mean(snpdens$SNP_COUNT), linetype="solid")
  
  return(p_snp_density)
}

#Braunvieh: 338,122 SNPs
#Fleckvieh: 428,439 SNPs

####Run function, plot results####
bv_all_dens <- snp_density_fnc(snpdens_file = "bv_all_SNP_freqs_4k.snpden", what_species = "Braunvieh" , what_subset = "All")
fv_all_dens <- snp_density_fnc(snpdens_file = "fv_all_SNP_freqs_4k.snpden", what_species = "Fleckvieh" , what_subset = "All")

prow <- cowplot::plot_grid(bv_all_dens,
                           fv_all_dens,
                           labels = c("A", "B"),
                           hjust = -1,
                           ncol = 1)

####Save plot####
ggsave(
  filename="SNP-Density-BV-FV.png",
  plot = prow,
  device = "png",
  path = NULL,
  scale = 1,
  width = 25 ,
  height = 20,
  units = c("cm"),
  limitsize = TRUE,
)



