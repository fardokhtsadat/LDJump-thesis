#####Info####
"This file describes the work process for 
 (Chapter 6.6 \"Comparison of recombination patterns with varying levels of inbreeding
and their correlation to SNP density\").

 Here, we plot the recombination maps of all categories (all, 0.125, 0.0625) for each cattle-population,
 in addition to the SNP density of the highest density region.

 Also, the same type of plot is shown for a random region of chromosome 25."

####Libraries#####
library(ggplot2)
library(cowplot)
setwd("RMPs-and-SNP-density")

####Functions for plot####
combined_rm_snpdens <- function(snpdens_file, what_species, what_dr, my_label, all, sub125, sub625, dr_s, dr_e){
  snpdens <- read.table(snpdens_file, header = T)
  
  print(paste0("Average SNP-density: ", mean(snpdens$SNP_COUNT)))
  
  dr_s <- which(snpdens$BIN_START == dr_s)
  dr_e <- which(snpdens$BIN_START == dr_e)
  
  p_snp_density <- ggplot2::ggplot(snpdens[dr_s:dr_e,], aes(x=BIN_START, y=SNP_COUNT, fill = SNP_COUNT)) +
    ggtitle(paste0(what_dr, " SNP-density region"))+
    coord_cartesian(ylim=c(0, 175), xlim=c(snpdens[dr_s,]$BIN_START-4000, snpdens[dr_e,]$BIN_START +4000))+
    scale_x_continuous(breaks = seq(snpdens[dr_s,]$BIN_START, snpdens[dr_e,]$BIN_START, 2e5), labels=function(x) format(x, scientific = T))+
    scale_y_continuous(breaks = seq(0,175,by=25))+
    geom_bar(stat = "identity")+
    scale_fill_gradient(low="blue", high="red", limits = c(0,200))+
    xlab("Chromosome 25") +
    ylab("SNP count per 4kb segment") +
    labs(fill = "SNP Count")+
    theme(axis.text = element_text(size = 15),
          plot.title=element_text(hjust=0.5, size=19),
          legend.text = element_text(size = 19),
          axis.title.x = element_text(size = 19),
          axis.title.y = element_text(size = 19),
          legend.title = element_text(size = 19),
          legend.key.size = unit(1, "cm"))+
    annotate("label", x=(snpdens$BIN_START[dr_s+457]),
             y = 168, label = paste0("Total SNP count:\n", sum(snpdens$SNP_COUNT[dr_s:dr_e])), size = 6)        
  

  all = readRDS(all)
  sub125 = readRDS(sub125)
  sub625 = readRDS(sub625)
  
  mappingvals_all <- data.frame(steps = sort(c(all[[1]]$leftEnd, all[[1]]$rightEnd)), 
                                values = as.vector(matrix(c(all[[1]]$value, all[[1]]$value), nrow=2, byrow=TRUE)), Group = "No-Cutoff")
  mappingvals_sub125 <- data.frame(steps = sort(c(sub125[[1]]$leftEnd, sub125[[1]]$rightEnd)), 
                                   values = as.vector(matrix(c(sub125[[1]]$value, sub125[[1]]$value), nrow=2, byrow=TRUE)), Group = "0.125")
  mappingvals_sub625 <- data.frame(steps = sort(c(sub625[[1]]$leftEnd, sub625[[1]]$rightEnd)), 
                                   values = as.vector(matrix(c(sub625[[1]]$value, sub625[[1]]$value), nrow=2, byrow=TRUE)), Group = "0.0625")
  
  df_all <- rbind(mappingvals_all, mappingvals_sub125, mappingvals_sub625)
  
  
  if(what_dr == "Highest"){dr = "HDR"} else if (what_dr == "Lowest"){dr = "LDR"} else if (what_dr == "Random"){dr = "Random"} else stop("Wrong input for density region.")
  
  
  rmp <- ggplot2::ggplot(data = df_all, aes(x=(steps*4000)+snpdens[dr_s,]$BIN_START, y=values, color = Group)) + geom_step(direction = "vh", linetype = 1, size = 1) +
    scale_y_continuous(breaks = seq(0,0.20,by=0.05))+
    coord_cartesian(ylim=c(0, 0.165), xlim=c(snpdens[dr_s,]$BIN_START-4000, snpdens[dr_e,]$BIN_START))+
    scale_x_continuous(breaks = seq(snpdens[dr_s,]$BIN_START, snpdens[dr_e,]$BIN_START, 2e5), labels=function(x) format(x, scientific = T))+
    ggtitle(paste0(what_species, "\nRecombination Map - ", dr))+
    xlab("Chromosome 25") +
    ylab("Recombination rate per 4kb segment") +
    theme(axis.text = element_text(size = 15),
          plot.title=element_text(hjust=0.5, size=19),
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size = 19),
          axis.title.y = element_text(size = 19),
          axis.title.x = element_text(size = 19),
          legend.text = element_text(size = 19),
          legend.text.align = 0) + scale_color_manual(name = "Group", labels = my_label, values = c("#0000CC", "#00CC00", "#CC6600"))
  
  return(list(rmp, p_snp_density))
}


####Plots####

#Within group BRAUNVIEH
snpdens_file = "bv_all_SNP_freqs_4k.snpden" 
what_species = "Braunvieh - Comparison" 
what_dr = "Highest"
all = "bv_all_high_demoT.rds"
sub125 = "bv_125_high_demoT.rds"
sub625 =  "bv_625_high_demoT.rds"
dr_s = 10692000 
dr_e = 12692000
my_label = c("No Cut-off\n(91 cows)", "\n0.125\n(57 cows)", "\n0.0625\n(42 cows)")
comall <- combined_rm_snpdens(snpdens_file = snpdens_file, what_species = what_species, what_dr = what_dr, my_label = my_label, all = all, sub125 = sub125, sub625 = sub625, dr_s = dr_s, dr_e = dr_e)
comall <- cowplot::plot_grid(comall[[1]],
                             comall[[2]],
                             labels = c("A", "B"),
                             nrow = 2,
                             align = "vh")
#comall

ggsave(
  filename="BV-RMP-WithinPop.png",
  plot = comall,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 32,
  units = c("cm"),
  limitsize = TRUE,
)

#Within group FLECKVIEH
snpdens_file = "fv_all_SNP_freqs_4k.snpden" 
what_species = "Fleckvieh - Comparison" 
what_dr = "Highest"
all = "fv_all_high_demoT.rds"
sub125 = "fv_125_high_demoT.rds"
sub625 =  "fv_625_high_demoT.rds"
dr_s = 10728000 
dr_e = 12728000
my_label = c("No Cut-off\n(161 cows)", "\n0.125\n(77 cows)", "\n0.0625\n(53 cows)")
comall <- combined_rm_snpdens(snpdens_file = snpdens_file, what_species = what_species, what_dr = what_dr, my_label=my_label, all = all, sub125 = sub125, sub625 = sub625, dr_s = dr_s, dr_e = dr_e)
comall <- cowplot::plot_grid(comall[[1]],
                             comall[[2]],
                             labels = c("A", "B"),
                             nrow = 2,
                             align = "vh")
#comall

ggsave(
  filename="FV-RMP-WithinPop.png",
  plot = comall,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 32,
  units = c("cm"),
  limitsize = TRUE,
)


####Random Region####

#Braunvieh Random:
snpdens_file = "bv_all_SNP_freqs_4k.snpden" 
what_species = "Braunvieh - Comparison" 
what_dr = "Random"
all = "random_bv_all_demoT.rds"
sub125 = "random_bv_125_demoT.rds"
sub625 =  "random_bv_625_demoT.rds"
#20074456
#22074456
dr_s = 20076000 
dr_e = 22076000
my_label = c("No Cut-off\n(91 cows)", "\n0.125\n(57 cows)", "\n0.0625\n(42 cows)")
comall <- combined_rm_snpdens(snpdens_file = snpdens_file, what_species = what_species, what_dr = what_dr, my_label = my_label, all = all, sub125 = sub125, sub625 = sub625, dr_s = dr_s, dr_e = dr_e)
comall <- cowplot::plot_grid(comall[[1]],
                             comall[[2]],
                             labels = c("A", "B"),
                             nrow = 2,
                             align = "vh")
#comall

ggsave(
  filename="BV-RMP-Random-WithinPop.png",
  plot = comall,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 32,
  units = c("cm"),
  limitsize = TRUE,
)



#Fleckvieh Random:
snpdens_file = "fv_all_SNP_freqs_4k.snpden" 
what_species = "Fleckvieh - Comparison" 
what_dr = "Random"
all = "random_fv_all_demoT.rds"
sub125 = "random_fv_125_demoT.rds"
sub625 =  "random_fv_625_demoT.rds"
dr_s = 20076000 
dr_e = 22076000
my_label = c("No Cut-off\n(161 cows)", "\n0.125\n(77 cows)", "\n0.0625\n(53 cows)")
comall <- combined_rm_snpdens(snpdens_file = snpdens_file, what_species = what_species, what_dr = what_dr, my_label=my_label, all = all, sub125 = sub125, sub625 = sub625, dr_s = dr_s, dr_e = dr_e)
comall <- cowplot::plot_grid(comall[[1]],
                             comall[[2]],
                             labels = c("A", "B"),
                             nrow = 2,
                             align = "vh")
#comall

ggsave(
  filename="FV-RMP-Random-WithinPop.png",
  plot = comall,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 32,
  units = c("cm"),
  limitsize = TRUE,
)



