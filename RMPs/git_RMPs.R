#####Info####
"This file describes the work process for:
 (Chapter 6.3 \"Estimation of recombination rates under neutrality and demography using genotyped cattle data\").

 We will plot the recombination maps for all three groups for both cow populations 
 on top of each other, where we want to compare demography settings equal to TRUE and FALSE."

####Libraries and working directory#####
library(ggplot2)
library(cowplot)
setwd("~/RMPs")

####Functions for plot####
collapsed_rmps <- function(snpdens_file, what_species, what_subset, what_dr, demoTresult, demoFresult, dr_s, dr_e){
  snpdens <- read.table(snpdens_file, header = T)
  
  print(paste0("Average SNP-density: ", mean(snpdens$SNP_COUNT)))
  
  dr_s <- which(snpdens$BIN_START == dr_s)
  dr_e <- which(snpdens$BIN_START == dr_e)
  demoT = readRDS(demoTresult)
  demoF = readRDS(demoFresult)
  
  mappingvalsT <- data.frame(steps = sort(c(demoT[[1]]$leftEnd, demoT[[1]]$rightEnd)), 
                             values = as.vector(matrix(c(demoT[[1]]$value, demoT[[1]]$value), nrow=2, byrow=TRUE)),
                             Demography = TRUE)
  
  
  my_median <- cbind(mappingvalsT, Freqs=c(diff(mappingvalsT$steps),1))
  my_median <- data.frame(Recombination = rep(my_median$values, my_median$Freqs))
  
  my_median <- median(my_median$Recombination)
  
  mappingvalsF <- data.frame(steps = sort(c(demoF[[1]]$leftEnd, demoF[[1]]$rightEnd)), 
                             values = as.vector(matrix(c(demoF[[1]]$value, demoF[[1]]$value), nrow=2, byrow=TRUE)),
                             Demography = FALSE)
  
  df_all <- rbind(mappingvalsT, mappingvalsF)
  
  if(what_dr == "High"){dr = "HDR"} else if (what_dr == "Low"){dr = "LDR"} else stop("Wrong input for density region.")
  
  rmp <- ggplot2::ggplot(data = df_all, aes(x=(steps*4000)+snpdens[dr_s,]$BIN_START, y=values, color = Demography)) + geom_step(direction = "vh", linetype = 1, size = 1) +
    scale_y_continuous(breaks = seq(0,0.20,by=0.05))+
    coord_cartesian(ylim=c(0, 0.165), xlim=c(snpdens[dr_s,]$BIN_START-4000, snpdens[dr_e,]$BIN_START))+
    scale_x_continuous(breaks = seq(snpdens[dr_s,]$BIN_START, snpdens[dr_e,]$BIN_START, 2e5), labels=function(x) format(x, scientific = T))+
    ggtitle(paste0("Recombination Map - ", what_species, "\n", what_subset, " Cut-off - ", dr))+
    xlab("Chromosome 25") +
    ylab("Recombination rate per 4kb segment") +
    theme(axis.text = element_text(size = 11),
          plot.title=element_text(hjust=0.5, size=17),
          legend.key.size = unit(1.2, "cm"),
          legend.title = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          legend.text = element_text(size = 13),
          axis.title.x = element_text(size = 17)) 
  return(rmp)
}

####Plotting correlated recombination maps with SNP-density####

####High-Density####

bv_all_high <- collapsed_rmps(snpdens_file = "bv_all_SNP_freqs_4k.snpden", 
                              what_species = "Braunvieh" , 
                              what_subset = "No",
                              what_dr = "High", 
                              demoT = "bv_all_high_demoT.rds",
                              demoF = "bv_all_high_demoF.rds",
                              dr_s = 10728000,
                              dr_e = 12728000)

bv_125_high <- collapsed_rmps(snpdens_file = "bv_125_SNP_freqs_4k.snpden", 
                              what_species = "Braunvieh" , 
                              what_subset = "0.125",
                              what_dr = "High", 
                              demoT = "bv_125_high_demoT.rds",
                              demoF = "bv_125_high_demoF.rds",
                              dr_s = 10728000,
                              dr_e = 12728000)

bv_625_high <- collapsed_rmps(snpdens_file = "bv_625_SNP_freqs_4k.snpden", 
                              what_species = "Braunvieh" , 
                              what_subset = "0.0625",
                              what_dr = "High", 
                              demoT = "bv_625_high_demoT.rds",
                              demoF = "bv_625_high_demoF.rds",
                              dr_s = 10728000,
                              dr_e = 12728000)

fv_all_high <- collapsed_rmps(snpdens_file = "fv_all_SNP_freqs_4k.snpden", 
                              what_species = "Fleckvieh" , 
                              what_subset = "No",
                              what_dr = "High", 
                              demoT = "fv_all_high_demoT.rds",
                              demoF = "fv_all_high_demoF.rds",
                              dr_s = 10692000,
                              dr_e = 12692000)

fv_125_high <- collapsed_rmps(snpdens_file = "fv_125_SNP_freqs_4k.snpden", 
                              what_species = "Fleckvieh" , 
                              what_subset = "0.125",
                              what_dr = "High", 
                              demoT = "fv_125_high_demoT.rds",
                              demoF = "fv_125_high_demoF.rds",
                              dr_s = 10692000,
                              dr_e = 12692000)

fv_625_high <- collapsed_rmps(snpdens_file = "fv_625_SNP_freqs_4k.snpden", 
                              what_species = "Fleckvieh" , 
                              what_subset = "0.0625",
                              what_dr = "High", 
                              demoT = "fv_625_high_demoT.rds",
                              demoF = "fv_625_high_demoF.rds",
                              dr_s = 10692000,
                              dr_e = 12692000)


my_leg <- cowplot::get_legend(bv_all_high)

plot_noleg <- cowplot::plot_grid(bv_all_high + theme(legend.position = "none", axis.title.x = element_blank()),
                                 fv_all_high + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()),
                                 bv_125_high + theme(legend.position = "none", axis.title.x = element_blank()),
                                 fv_125_high + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()),
                                 bv_625_high + theme(legend.position = "none"),
                                 fv_625_high + theme(legend.position = "none", axis.title.y = element_blank()),
                                 align = 'vh',
                                 labels = c("A", "B", "C", "D", "E", "F"),
                                 hjust = -1,
                                 nrow = 3,
                                 ncol = 2)
full_plot <- cowplot::plot_grid(plot_noleg, my_leg, rel_widths = c(3,.2))


ggsave(
  filename="collapsed-rmp-high.png",
  plot = full_plot,
  device = "png",
  path = NULL,
  scale = 1,
  width = 56,
  height = 36,
  units = c("cm"),
  limitsize = TRUE,
)

####Low-Density####

bv_all_low <- collapsed_rmps(snpdens_file = "bv_all_SNP_freqs_4k.snpden", 
                             what_species = "Braunvieh" , 
                             what_subset = "No",
                             what_dr = "Low", 
                             demoT = "bv_all_low_demoT.rds",
                             demoF = "bv_all_low_demoF.rds",
                             dr_s = 36000,
                             dr_e = 2036000)

bv_125_low <- collapsed_rmps(snpdens_file = "bv_125_SNP_freqs_4k.snpden", 
                             what_species = "Braunvieh" , 
                             what_subset = "0.125",
                             what_dr = "Low", 
                             demoT = "bv_125_low_demoT.rds",
                             demoF = "bv_125_low_demoF.rds",
                             dr_s = 36000,
                             dr_e = 2036000)

bv_625_low <- collapsed_rmps(snpdens_file = "bv_625_SNP_freqs_4k.snpden", 
                             what_species = "Braunvieh" , 
                             what_subset = "0.0625",
                             what_dr = "Low", 
                             demoT = "bv_625_low_demoT.rds",
                             demoF = "bv_625_low_demoF.rds",
                             dr_s = 36000,
                             dr_e = 2036000)
bv_625_low <- bv_625_low +
  annotate("label", x = 1.838e+06, y = 0.15, label = paste0("Recombination-rate value:\n", 0.270907))

fv_all_low <- collapsed_rmps(snpdens_file = "fv_all_SNP_freqs_4k.snpden", 
                             what_species = "Fleckvieh" , 
                             what_subset = "No",
                             what_dr = "Low", 
                             demoT = "fv_all_low_demoT.rds",
                             demoF = "fv_all_low_demoF.rds",
                             dr_s = 36000,
                             dr_e = 2036000)

fv_125_low <- collapsed_rmps(snpdens_file = "fv_125_SNP_freqs_4k.snpden", 
                             what_species = "Fleckvieh" , 
                             what_subset = "0.125",
                             what_dr = "Low", 
                             demoT = "fv_125_low_demoT.rds",
                             demoF = "fv_125_low_demoF.rds",
                             dr_s = 36000,
                             dr_e = 2036000)

fv_625_low <- collapsed_rmps(snpdens_file = "fv_625_SNP_freqs_4k.snpden", 
                             what_species = "Fleckvieh" , 
                             what_subset = "0.0625",
                             what_dr = "Low", 
                             demoT = "fv_625_low_demoT.rds",
                             demoF = "fv_625_low_demoF.rds",
                             dr_s = 36000,
                             dr_e = 2036000)


my_leg <- cowplot::get_legend(bv_all_low)

plot_noleg <- cowplot::plot_grid(bv_all_low + theme(legend.position = "none", axis.title.x = element_blank()),
                                 fv_all_low + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()),
                                 bv_125_low + theme(legend.position = "none", axis.title.x = element_blank()),
                                 fv_125_low + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()),
                                 bv_625_low + theme(legend.position = "none"),
                                 fv_625_low + theme(legend.position = "none", axis.title.y = element_blank()),
                                 align = 'vh',
                                 labels = c("A", "B", "C", "D", "E", "F"),
                                 hjust = -1,
                                 nrow = 3,
                                 ncol = 2)
full_plot <- cowplot::plot_grid(plot_noleg, my_leg, rel_widths = c(3,.2))



ggsave(
  filename="collapsed-rmp-low.png",
  plot = full_plot,
  device = "png",
  path = NULL,
  scale = 1,
  width = 56,
  height = 36,
  units = c("cm"),
  limitsize = TRUE,
)



