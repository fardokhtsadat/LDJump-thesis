#####Info####
"This file describes the work process for 
 (Chapter 6.4 \"Comparison of recombination patterns between two cattle breeds\").

 Here, we compare the recombination maps of both populations with each other for all three subsets."

####Libraries#####
library(ggplot2)
library(cowplot)
setwd("~/RMPs-Comparison")

####Functions for plot####
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


####Braunvieh and Fleckvieh dataset combined####

#All High compared
snpdens_file = "bv_fv_combined_all_SNP_freqs_4k.snpden" 
what_subset = "No Cut-off"
what_dr = "High"
braunvieh = "bv_all_high_demoT.rds"
fleckvieh = "fv_all_high_demoT.rds"
my_label = c("Braunvieh\n(91 cows)\n", "\nFleckvieh\n(161 cows)")
comall <- compare_breeds_rmp(snpdens_file = snpdens_file, what_subset = what_subset, what_dr = what_dr, my_label=my_label, braunvieh = braunvieh, fleckvieh = fleckvieh)
#comall


#125 High compared
snpdens_file = "bv_fv_combined_125_SNP_freqs_4k.snpden" 
what_subset = "0.125"
what_dr = "High"
braunvieh = "bv_125_high_demoT.rds"
fleckvieh = "fv_125_high_demoT.rds"
my_label = c("Braunvieh\n(57 cows)\n", "\nFleckvieh\n(77 cows)")
com125 <- compare_breeds_rmp(snpdens_file = snpdens_file, what_subset = what_subset, what_dr = what_dr, my_label=my_label, braunvieh = braunvieh, fleckvieh = fleckvieh)
#com125


#625 High compared
snpdens_file = "bv_fv_combined_625_SNP_freqs_4k.snpden" 
what_subset = "0.0625"
what_dr = "High"
braunvieh = "bv_625_high_demoT.rds"
fleckvieh = "fv_625_high_demoT.rds"
my_label = c("Braunvieh\n(42 cows)\n", "\nFleckvieh\n(53 cows)")
com625 <- compare_breeds_rmp(snpdens_file = snpdens_file, what_subset = what_subset, what_dr = what_dr, my_label=my_label, braunvieh = braunvieh, fleckvieh = fleckvieh)
#com625

####Full plot
com_allthree <- cowplot::plot_grid(comall,
                                   com125, 
                                   com625,
                                   labels = c("A", "B", "C"),
                                   align = "vh",
                                   ncol=1)
#com_allthree

#portrait
ggsave(
  filename="Portrait-Compared-RMPs.png",
  plot = com_allthree,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 38,
  units = c("cm"),
  limitsize = TRUE,
)
