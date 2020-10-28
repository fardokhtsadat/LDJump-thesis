#####Info####
"This file describes the work process for 
 (Chapter 6.5 \"Shared recombination patterns between breeds\").

 Here, we plot shared recombination patterns by using the combined dataset. 
 The combined recombination map is plotted along with the individual recombination maps."

####Libraries#####
library(ggplot2)
library(cowplot)
setwd("~/RMPs-Shared")

####Functions for plot####
combined_rm_snpdens <- function(snpdens_file, what_species, what_subset, what_dr, my_label, braunvieh, fleckvieh, combined){
  
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
  combined = readRDS(combined)
  
  mappingvals_bv <- data.frame(steps = sort(c(braunvieh[[1]]$leftEnd, braunvieh[[1]]$rightEnd)) + 9, 
                               values = as.vector(matrix(c(braunvieh[[1]]$value, braunvieh[[1]]$value), nrow=2, byrow=TRUE)), Group = "Braunvieh")
  mappingvals_fv <- data.frame(steps = sort(c(fleckvieh[[1]]$leftEnd, fleckvieh[[1]]$rightEnd)), 
                               values = as.vector(matrix(c(fleckvieh[[1]]$value, fleckvieh[[1]]$value), nrow=2, byrow=TRUE)), Group = "Fleckvieh")
  mappingvals_com <- data.frame(steps = sort(c(combined[[1]]$leftEnd, combined[[1]]$rightEnd)), 
                                values = as.vector(matrix(c(combined[[1]]$value, combined[[1]]$value), nrow=2, byrow=TRUE)), Group = "Combined")
  
  com_median <- cbind(mappingvals_com, Freqs=c(diff(mappingvals_com$steps),1)); 
  com_median <- data.frame(Recombination = rep(com_median$values, com_median$Freqs)); 
  com_median <- median(com_median$Recombination)
  
  df_all <- rbind(mappingvals_bv, mappingvals_fv, mappingvals_com)
  
  
  if(what_dr == "High"){dr = "HDR"} else if (what_dr == "Low"){dr = "LDR"} else stop("Wrong input for density region.")
  
  rmp <- ggplot2::ggplot(data = df_all, aes(x=(steps*4000)+snpdens[dr_s,]$BIN_START, y=values, color = Group)) + geom_step(direction = "vh", linetype = 1, size = 1) +
    scale_y_continuous(breaks = seq(0,0.20,by=0.05))+
    coord_cartesian(ylim=c(0, 0.165), xlim=c(snpdens[dr_s,]$BIN_START-4000, snpdens[dr_e,]$BIN_START))+
    scale_x_continuous(breaks = seq(snpdens[dr_s,]$BIN_START, snpdens[dr_e,]$BIN_START, 2e5), labels=function(x) format(x, scientific = T))+
    ggtitle(paste0("Recombination Map - ", dr, " - ", what_subset))+
    xlab("Chromosome 25") +
    ylab("Recombination rate per 4kb segment") +
    geom_hline(yintercept = com_median*3, linetype = 3, color = "black", size = 1) +
    theme(axis.text = element_text(size = 15),
          plot.title=element_text(hjust=0.5, size=19),
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size = 19),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 19),
          axis.title.x = element_text(size = 19)) + scale_color_manual(labels = my_label, values = c("red", "deepskyblue", "black"))
  
  
  return(rmp)
}


####Braunvieh and Fleckvieh dataset combined####
#125 High combined
snpdens_file = "bv_fv_combined_125_SNP_freqs_4k.snpden" 
what_species = "Combined Data" 
what_subset = "0.125"
what_dr = "High"
braunvieh = "bv_125_high_demoT.rds"
fleckvieh = "fv_125_high_demoT.rds"
combined =  "combined_125_high_demoT.rds"
my_label = c("Braunvieh\n(57 cows)\n", "Fleckvieh\n(77 cows)", "\nCombined\n(134 cows)")
com125 <- combined_rm_snpdens(snpdens_file = snpdens_file, what_species = what_species, what_subset = what_subset, what_dr = what_dr, my_label = my_label, braunvieh = braunvieh, fleckvieh = fleckvieh, combined = combined)


#625 High combined
snpdens_file = "bv_fv_combined_625_SNP_freqs_4k.snpden" 
what_species = "Combined Data" 
what_subset = "0.0625"
what_dr = "High"
braunvieh = "bv_625_high_demoT.rds"
fleckvieh = "fv_625_high_demoT.rds"
combined =  "combined_625_high_demoT.rds"
my_label = c("Braunvieh\n(42 cows)\n", "Fleckvieh\n(53 cows)", "\nCombined\n(95 cows)")
com625 <- combined_rm_snpdens(snpdens_file = snpdens_file, what_species = what_species, what_subset = what_subset, what_dr = what_dr, my_label = my_label, braunvieh = braunvieh, fleckvieh = fleckvieh, combined = combined)


####Full plot

com_625_125 <- cowplot::plot_grid(com125, 
                                  com625, 
                                  labels = c("A", "B"),
                                  align = "vh",
                                  ncol=1)
#com_625_125


ggsave(
  filename="Portrait-Combined-RMPs-625125.png",
  plot = com_625_125,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 26,
  units = c("cm"),
  limitsize = TRUE,
)

