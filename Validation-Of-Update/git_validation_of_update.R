####Info####
"This file describes the work process for the validation of LDJump (Chapter 5.1.3. \"Validation of the update\"), which includes
1. Loading the data
2. Transforming
3. Visualizing in R
"
####Running LDJump####

#Loading the packages:
library(parallel)
library(seqinr)
library(vcfR)
library(ape)
library(LDJump)

#Setting the working directory:
"Please consider setting your working directory in the folder in which you have all the files that you will be using (VCF-file, Reference-File)"
setwd("~/Validation-Of-Update")

#The files we are using:
#Group 1:
ref_seq = "Reference_CH21_41187000_41290679.fasta"
vcf_file = "TSI_21_41187000_41290679.vcf"
fasta_seq = "TSI_CH21_41187000_41290679.fa"

startofseq = 41187000
endofseq = 41290679

#Group 2:
ref_seq = "41m_41m10k.fa"
vcf_file = "41m_41m10k.vcf"
fasta_seq = "my_fasta.fasta"

startofseq = 41000000
endofseq = 41010000

#Running the code:

##Using our VCF-File:
start_time_vcf <- Sys.time()
results_testing = LDJump(vcf_file, chr = 21 , segLength = 1000, cores = 6, pathPhi = "/home/roots/anaconda3/bin/Phi", format = "vcf", refName = ref_seq, lengthofseq = endofseq-startofseq, startofseq = startofseq, endofseq = endofseq)
end_time_vcf <- Sys.time()
time_taken_vcf <- end_time_vcf - start_time_vcf

##Using a FASTA-File of the same DNA sequence for comparison:
start_time_fasta <- Sys.time()
fasta_only_result = LDJump(fasta_seq , segLength = 1000, cores = 3, pathPhi = "/home/roots/anaconda3/bin/Phi", format = "fasta", refName = ref_seq)
end_time_fasta <- Sys.time()
time_taken_fasta <- end_time_fasta - start_time_fasta

#Comparison of the results:
postscript("ResultsVCF.pdf", horiz = F)
plot(results[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump using VCF")
dev.off()

postscript("ResultsFASTA.pdf", horiz = F)
plot(fasta_only_result[[1]], xlab = "Segments", ylab = "Estimated Recombination Rate", main = "Estimated recombination map with LDJump using FASTA")
dev.off()



####Plotting####
#Load dataset
load("results_comparison_vcf_fasta.RData")

#Transform VCF result
demoF <- results
mappingvalsF <- data.frame(steps = sort(c(demoF[[1]]$leftEnd, demoF[[1]]$rightEnd)), 
                           values = as.vector(matrix(c(demoF[[1]]$value, demoF[[1]]$value), nrow=2, byrow=TRUE)),
                           Demography = FALSE)
#Plot
rmpvcf <- ggplot2::ggplot() +
  geom_step(data = mappingvalsF, 
            mapping= aes(x = steps, y=values),
            direction = "vh", linetype = 1) +
  scale_y_continuous(breaks = seq(0,0.15,by=0.05))+
  coord_cartesian(ylim=c(0, 0.15))+
  ggtitle(expression(atop("Estimated recombination map using" , paste(italic("LDJump"), " with VCF format")))) +
  xlab("Jumps between 1kb segments") +
  ylab("Recombination rate per 1kb segment") +
  theme(axis.text = element_text(size = 8),
        plot.title=element_text(hjust=0.5),
        legend.key.size = unit(0.5, "cm")) 

#Transform FASTA result
demoF <- fasta_only_result
mappingvalsF <- data.frame(steps = sort(c(demoF[[1]]$leftEnd, demoF[[1]]$rightEnd)), 
                           values = as.vector(matrix(c(demoF[[1]]$value, demoF[[1]]$value), nrow=2, byrow=TRUE)),
                           Demography = FALSE)
#Plot
rmpfasta <- ggplot2::ggplot() +
  geom_step(data = mappingvalsF, 
            mapping= aes(x = steps, y=values),
            direction = "vh", linetype = 1) +
  scale_y_continuous(breaks = seq(0,0.15,by=0.05))+
  coord_cartesian(ylim=c(0, 0.15))+
  ggtitle(expression(atop("Estimated recombination map using" , paste(italic("LDJump"), " with FASTA format")))) +
  xlab("Jumps between 1kb segments") +
  ylab("Recombination rate per 1kb segment") +
  theme(axis.text = element_text(size = 8),
        plot.title=element_text(hjust=0.5),
        legend.key.size = unit(0.5, "cm"),
        axis.title.y = element_blank()) 

#Combine plots
combined <- cowplot::plot_grid(
  rmpvcf,
  rmpfasta,
  ncol = 2,
  align = "h",
  labels = c("A", "B"))

#Save plot
ggsave(
  filename="Comparison-VCF-FASTA.png",
  plot = combined,
  device = "png",
  path = NULL,
  scale = 1,
  width = 20,
  height = 10,
  units = c("cm"),
  dpi = 300,
  limitsize = TRUE,
)