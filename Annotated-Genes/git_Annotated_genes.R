#####Info####
"This file describes the work process for 
 (Chapter 6.9 \"Annotated genes in the HDR and LDR\").

 Here, we present the annotated genes in the HDR and the LDR.

 The genes are downloaded from NCBI (search for genes in the genomic regions of the HDR and LDR and save the file)"

####Libraries#####
setwd("~/Annotated-Genes")

####Function####
compare_breeds_rmp <- function(snpdens_file, what_subset, what_dr, my_label=my_label, braunvieh, fleckvieh){
  
  snpdens <- read.table(snpdens_file, header = T)
  
  print(paste0("Average SNP-density: ", mean(snpdens$SNP_COUNT)))
  
  if(what_dr == "High"){
    regionstart = 10692000
    regionend = 12728000
  } else if (what_dr == "Low") {
    regionstart = 36000
    regionend = 2036000
  } else stop("Wrong input for density region.")
  
  dr_s <- which(snpdens$BIN_START == regionstart)
  dr_e <- which(snpdens$BIN_START == regionend)
  
  
  braunvieh = readRDS(braunvieh)
  fleckvieh = readRDS(fleckvieh)
  
  
  mappingvals_bv <- data.frame(steps = sort(c(braunvieh[[1]]$leftEnd, braunvieh[[1]]$rightEnd)) + 9, 
                               values = as.vector(matrix(c(braunvieh[[1]]$value, braunvieh[[1]]$value), nrow=2, byrow=TRUE)), Group = "Braunvieh")
  mappingvals_fv <- data.frame(steps = sort(c(fleckvieh[[1]]$leftEnd, fleckvieh[[1]]$rightEnd)), 
                               values = as.vector(matrix(c(fleckvieh[[1]]$value, fleckvieh[[1]]$value), nrow=2, byrow=TRUE)), Group = "Fleckvieh")
  
  
  bv_median <- cbind(mappingvals_bv, Freqs=c(diff(mappingvals_bv$steps),1)); 
  bv_median <- data.frame(Recombination = rep(bv_median$values, bv_median$Freqs)); 
  bv_median <- median(bv_median$Recombination)
  fv_median <- cbind(mappingvals_fv, Freqs=c(diff(mappingvals_fv$steps),1)); 
  fv_median <- data.frame(Recombination = rep(fv_median$values, fv_median$Freqs));
  fv_median <- median(fv_median$Recombination)
  
  
  df_all <- rbind(mappingvals_bv, mappingvals_fv)
  
  
  if(what_dr == "High"){dr = "HDR"} else if (what_dr == "Low"){dr = "LDR"} else stop("Wrong input for density region.")
  
  rmp <- ggplot2::ggplot(data = df_all, aes(x=(steps*4000)+snpdens[dr_s,]$BIN_START, y=values, color = Group)) + geom_step(direction = "vh", linetype = 1, size = 1) +
    scale_y_continuous(breaks = seq(0,0.20,by=0.05))+
    coord_cartesian(ylim=c(0, 0.165), xlim=c(snpdens[dr_s,]$BIN_START-4000, snpdens[dr_e,]$BIN_START))+
    scale_x_continuous(breaks = seq(snpdens[dr_s,]$BIN_START, snpdens[dr_e,]$BIN_START, 2e5), labels=function(x) format(x, scientific = T))+
    ggtitle(paste0("Recombination Map - ", dr, " - ", what_subset))+
    xlab("Chromosome 25") +
    ylab("Recombination rate per 4kb segment") +
    #linetype 2 == dashed
    #linetype 3 == dotted
    geom_hline(yintercept = bv_median*3, linetype = 2, color = "black", size = 1) +
    geom_hline(yintercept = fv_median*3, linetype = 3, color = "black", size = 1) +
    theme(axis.text = element_text(size = 15),
          plot.title=element_text(hjust=0.5, size=19),
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size = 19),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 19),
          axis.title.x = element_text(size = 19)) + scale_color_manual(labels = my_label, values = c("red", "deepskyblue"))
  
  
  return(rmp)
}


####Reading into genes files####
#Either one by one:
all_genes <- rbind(
cbind(read.table("/home/roots/Desktop/projects/LD/LDJump_Thesis/runningthesis/writing_thesis/gene_annotation/HDR_noncoding_genes", header = T, stringsAsFactors = F), Region = "HDR", Function = "Non-Coding"),
cbind(read.table("/home/roots/Desktop/projects/LD/LDJump_Thesis/runningthesis/writing_thesis/gene_annotation/HDR_protein_coding_genes", header = T, stringsAsFactors = F), Region = "HDR", Function = "Protein-Coding"),
cbind(read.table("/home/roots/Desktop/projects/LD/LDJump_Thesis/runningthesis/writing_thesis/gene_annotation/LDR_noncoding_genes", header = T, stringsAsFactors = F), Region = "LDR", Function = "Non-Coding"),
cbind(read.table("/home/roots/Desktop/projects/LD/LDJump_Thesis/runningthesis/writing_thesis/gene_annotation/LDR_protein_coding_genes", header = T, stringsAsFactors = F), Region = "LDR", Function = "Protein-Coding"),
cbind(read.table("/home/roots/Desktop/projects/LD/LDJump_Thesis/runningthesis/writing_thesis/gene_annotation/LDR_pseudogenes_genes", header = T, stringsAsFactors = F), Region = "LDR", Function = "Pseudo-Genes")
)
#Or at once:
all_genes <- read.table("all_genes.txt", header = T, stringsAsFactors = F)


####Plotting Recombination Maps####
#Braunvieh & Fleckvieh, highest SNP-density region:
snpdens_file = "bv_fv_combined_all_SNP_freqs_4k.snpden" 
what_subset = "No Cut-off"
what_dr = "High"
braunvieh = "bv_all_high_demoT.rds"
fleckvieh = "fv_all_high_demoT.rds"
my_label = c("Braunvieh\n(91 cows)\n", "\nFleckvieh\n(161 cows)")
comallh <- compare_breeds_rmp(snpdens_file = snpdens_file, what_subset = what_subset, what_dr = what_dr, my_label=my_label, braunvieh = braunvieh, fleckvieh = fleckvieh)

comallhh <- comallh + ggtitle("Annotated Genes - HDR - No Cut-off")+
  geom_segment(data = all_genes[which(all_genes$Region=="HDR"), ], 
               aes(x = start_position_on_the_genomic_accession, y = c(0.14, 0.16, 0.15, 0.15)),  
               xend = all_genes[which(all_genes$Region=="HDR"), ]$end_position_on_the_genomic_accession, 
               yend=c(0.14, 0.16, 0.15, 0.15), inherit.aes = FALSE, size = 4)+
  geom_point(data = all_genes[which(all_genes$Region=="HDR" & all_genes$Function == "Non-Coding"), ], 
             aes(x = all_genes[which(all_genes$Region=="HDR" & all_genes$Function == "Non-Coding"), ]$start_position_on_the_genomic_accession,
                 y = c(0.16)), inherit.aes = FALSE, size = 2, shape = 15, color = "red") 
#comallhh


#Braunvieh & Fleckvieh, lowest SNP-density region:
snpdens_file = "bv_fv_combined_all_SNP_freqs_4k.snpden" 
what_subset = "No Cut-off"
what_dr = "Low"
braunvieh = "bv_all_low_demoT.rds"
fleckvieh = "fv_all_low_demoT.rds"
my_label = c("Braunvieh\n(91 cows)\n", "\nFleckvieh\n(161 cows)")
comalll <- compare_breeds_rmp(snpdens_file = snpdens_file, what_subset = what_subset, what_dr = what_dr, my_label=my_label, braunvieh = braunvieh, fleckvieh = fleckvieh)


####Plotting gene segments####
#Plotting lines
poscod <- c(rep(c(0.150, 0.13, 0.14, 0.13, 0.13), 23), c(0.150, 0.14))
posncod <- c(rep(c(0.15, 0.13, 0.14, 0.13, 0.13), 3), c(0.150, 0.13, 0.14, 0.13))

comalllh <- comalll + ggtitle("Annotated Genes - LDR - No Cut-off") +
  geom_segment(data = all_genes[which(all_genes$Region=="LDR"& all_genes$Function == "Protein-Coding"), ], 
               aes(x = start_position_on_the_genomic_accession, y =poscod), 
               xend = all_genes[which(all_genes$Region=="LDR"& all_genes$Function == "Protein-Coding"), ]$end_position_on_the_genomic_accession, 
               yend=poscod, inherit.aes = FALSE, size = 4) +
  geom_segment(data = all_genes[which(all_genes$Region=="LDR"& all_genes$Function == "Non-Coding"), ], 
               aes(x = start_position_on_the_genomic_accession, y =posncod), color = "red", 
               xend = all_genes[which(all_genes$Region=="LDR"& all_genes$Function == "Non-Coding"), ]$end_position_on_the_genomic_accession, 
               yend=posncod, inherit.aes = FALSE, size = 4) 

#comalllh


####Combining plots####
#Combining the plots
legend <- cowplot::get_legend(comallhh + theme(legend.box.margin = margin(0, 0, 0, 0),             
                                               legend.key.size = unit(0.8, "cm",),
                                               legend.text = element_text(size = 19),
                                               legend.title = element_text(size = 19)))
prow <- cowplot::plot_grid(comallhh + theme(legend.position = "none"), comalllh + theme(legend.position = "none"), align = "vh", ncol = 1,
                           labels = c("A", "B"))

com_genes <- cowplot::plot_grid(prow, 
                                legend, 
                                rel_widths = c(3,.4))
#com_genes



####Saving plot####
ggsave(
  filename="Portrait-Annotated-RMPs.png",
  plot = com_genes,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 26,
  units = c("cm"),
  limitsize = TRUE,
)



