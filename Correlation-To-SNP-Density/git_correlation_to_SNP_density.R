#####Info####
"This file describes the work process for 
 (Chapter 6.7 \"Correlation between recombination rate and SNP density\").

 Here, we plot and test the correlation between SNP density and recombination rates. 
 
 Also, the same type of plot is shown for a random region of chromosome 25."

####Libraries#####
library(ggplot2)
library(cowplot)
setwd("Correlation-To-SNP-Density")

####Functions#####
cor_data <- function(all, snpdens, dr_s, dr_e, sortit = T){
  
  all <- readRDS(all)
  snpdens <- read.table(snpdens, header=T)
  
  dr_s <- which(snpdens$BIN_START == dr_s)
  dr_e <- which(snpdens$BIN_START == dr_e)
  
  mappingvals_all <- data.frame(steps = sort(c(all[[1]]$leftEnd, all[[1]]$rightEnd)), 
                                values = as.vector(matrix(c(all[[1]]$value, all[[1]]$value), nrow=2, byrow=TRUE)), Group = "No-Cutoff")
  
  mappingvals_all <- cbind(mappingvals_all, Freqs=c(diff(mappingvals_all$steps),1))
  
  mappingvals_all <- data.frame(Recombination = rep(mappingvals_all$values, mappingvals_all$Freqs), SNP = snpdens[dr_s:(dr_e-1),]$SNP_COUNT) 
  
  if(sortit){
    mappingvals_all <- mappingvals_all[order(mappingvals_all$SNP),]
  }
  return(mappingvals_all)
}

cor_p <- function(mappingvals_all, my_title){
  cortest <- cor.test(mappingvals_all$Recombination, mappingvals_all$SNP, method="pearson")
  pval <- round(cortest$p.value, 12)
  pmethod <- cortest$method
  pcor <- round(cortest$estimate, 3)

  cor_p <- ggplot(data = mappingvals_all, aes(x = SNP, y = Recombination)) + 
    geom_point() + geom_smooth(method = 'loess', formula = y ~ x) +
    ggtitle(my_title)+
    xlab("SNP count") +
    ylab("Recombination rate") +
    scale_y_continuous(breaks = seq(0,0.20,by=0.05))+
    scale_x_continuous(breaks = seq(0,160,20))+
    coord_cartesian(ylim=c(0, 0.160), xlim=c(0, 161))+
    annotate("label", x= 140, y=0.15, label=paste0("Correlation: ", pcor, "\nP-value: ", format(pval, scientific = TRUE)), size = 5) +
    theme(axis.text = element_text(size = 12),
          plot.title=element_text(hjust=0.5, size=15),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15)) 
  return(cor_p)
}

correlation_plot <- function(all, snpdens, dr_s, dr_e, my_title){
  my_data <- cor_data(all = all, snpdens = snpdens, dr_s = dr_s, dr_e = dr_e)
  cor_p(mappingvals_all = my_data, my_title = my_title)
}

####Correlation individually, recombination rate and SNP-density####
bv_all <- correlation_plot(all = "bv_all_high_demoT.rds", snpdens = "bv_all_SNP_freqs_4k.snpden", dr_s = 10692000, dr_e = 12692000, my_title = "Braunvieh - No Cut-off")
bv_125 <- correlation_plot(all = "bv_125_high_demoT.rds", snpdens = "bv_125_SNP_freqs_4k.snpden", dr_s = 10692000, dr_e = 12692000, my_title = "Braunvieh - 0.125")
bv_625 <- correlation_plot(all = "bv_625_high_demoT.rds", snpdens = "bv_625_SNP_freqs_4k.snpden", dr_s = 10692000, dr_e = 12692000, my_title = "Braunvieh - 0.0625")

fv_all <- correlation_plot(all = "fv_all_high_demoT.rds", snpdens = "fv_all_SNP_freqs_4k.snpden", dr_s = 10728000, dr_e = 12728000, my_title = "Fleckvieh - No Cut-off")
fv_125 <- correlation_plot(all = "fv_125_high_demoT.rds", snpdens = "fv_125_SNP_freqs_4k.snpden", dr_s = 10728000, dr_e = 12728000, my_title = "Fleckvieh - 0.125")
fv_625 <- correlation_plot(all = "fv_625_high_demoT.rds", snpdens = "fv_625_SNP_freqs_4k.snpden", dr_s = 10728000, dr_e = 12728000, my_title = "Fleckvieh - 0.0625")

cor_individ <- cowplot::plot_grid(bv_all, fv_all,
                                  bv_125, fv_125,
                                  bv_625, fv_625,
                                  nrow = 3,
                                  labels = c("A", "B", "C", "D", "E", "F"),
                                  align = "vh")
ggsave(
  filename="Correlation-Individually.png",
  plot = cor_individ,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 38,
  units = c("cm"),
  limitsize = TRUE,
)

####RANDOM REGION: Correlation individually, recombination rate and SNP-density####
bv_all <- correlation_plot(all = "random_bv_all_demoT.rds", snpdens = "bv_all_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, my_title = "Braunvieh - No Cut-off - Random")
bv_125 <- correlation_plot(all = "random_bv_125_demoT.rds", snpdens = "bv_125_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, my_title = "Braunvieh - 0.125 - Random")
bv_625 <- correlation_plot(all = "random_bv_625_demoT.rds", snpdens = "bv_625_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, my_title = "Braunvieh - 0.0625 - Random")

fv_all <- correlation_plot(all = "random_fv_all_demoT.rds", snpdens = "fv_all_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, my_title = "Fleckvieh - No Cut-off - Random")
fv_125 <- correlation_plot(all = "random_fv_125_demoT.rds", snpdens = "fv_125_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, my_title = "Fleckvieh - 0.125 - Random")
fv_625 <- correlation_plot(all = "random_fv_625_demoT.rds", snpdens = "fv_625_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, my_title = "Fleckvieh - 0.0625 - Random")

cor_individ <- cowplot::plot_grid(bv_all, fv_all,
                                  bv_125, fv_125,
                                  bv_625, fv_625,
                                  nrow = 3,
                                  labels = c("A", "B", "C", "D", "E", "F"),
                                  align = "vh")
ggsave(
  filename="Correlation-RANDOM-Individually.png",
  plot = cor_individ,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 38,
  units = c("cm"),
  limitsize = TRUE,
)

