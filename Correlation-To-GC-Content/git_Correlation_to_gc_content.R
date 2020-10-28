#####Info####
"This file describes the work process for 
 (Chapter 6.8 \"Correlation between recombination rate and GC content\").

 Here, we plot and test the correlation between GC content and recombination rates. 
 
 Also, the same type of plot is shown for a random region of chromosome 25."

####Libraries#####
library(ggplot2)
library(cowplot)
setwd("~/Correlation-To-GC-Content")

####Functions#####
gc_cor <- function(data, snpdens, gcfile, dr_s, dr_e, sortit, my_title){
  my_data <- cor_data(all = data, snpdens = snpdens, dr_s = dr_s, dr_e = dr_e, sortit = sortit)
  gc_bv <- read.table(file="gc_content_reference_genome.txt", header = T)
  mappingvals_all <- data.frame(Recombination = my_data$Recombination, GC = gc_bv$GC., SNP = my_data$SNP)
  mappingvals_all <- mappingvals_all[order(mappingvals_all$GC),]
  cortest <- cor.test(mappingvals_all$Recombination, mappingvals_all$GC, method="pearson")
  pval <- round(cortest$p.value, 2)
  pmethod <- cortest$method
  pcor <- round(cortest$estimate, 2)
  
  
  cor_p <- ggplot(data = mappingvals_all, aes(y = GC, x = Recombination)) + 
    geom_point() + geom_smooth(method = 'loess', formula = y ~ x) +
    ggtitle(my_title)+
    ylab("GC content") +
    xlab("Recombination rate") +
    scale_x_continuous(breaks = seq(0,0.20,by=0.05))+
    scale_y_continuous(breaks = seq(0.2,0.6,0.2))+
    coord_cartesian(xlim=c(0, 0.160), ylim=c(0, 0.6))+
    annotate("label", y= 0.07, x=0.142,label=paste0("Correlation: ", pcor, "\nP-value: ", format(pval, scientific = TRUE)), size = 5) +
    theme(axis.text = element_text(size = 12),
          plot.title=element_text(hjust=0.5, size=15),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15)) 
  
  return(cor_p)
}

####Running Functions####

bv_all <- gc_cor(data = "bv_all_high_demoT.rds", snpdens = "bv_all_SNP_freqs_4k.snpden", dr_s = 10692000, dr_e = 12692000, sortit = F, my_title = "Braunvieh - No Cut-off")
bv_125 <- gc_cor(data = "bv_125_high_demoT.rds", snpdens = "bv_125_SNP_freqs_4k.snpden", dr_s = 10692000, dr_e = 12692000, sortit = F, my_title = "Braunvieh - 0.125")
bv_625 <- gc_cor(data = "bv_625_high_demoT.rds", snpdens = "bv_625_SNP_freqs_4k.snpden", dr_s = 10692000, dr_e = 12692000, sortit = F, my_title = "Braunvieh - 0.0625")

fv_all <- gc_cor(data = "fv_all_high_demoT.rds", snpdens = "fv_all_SNP_freqs_4k.snpden", dr_s = 10728000, dr_e = 12728000, sortit = F, my_title = "Fleckvieh - No Cut-off")
fv_125 <- gc_cor(data = "fv_125_high_demoT.rds", snpdens = "fv_125_SNP_freqs_4k.snpden", dr_s = 10728000, dr_e = 12728000, sortit = F, my_title = "Fleckvieh - 0.125")
fv_625 <- gc_cor(data = "fv_625_high_demoT.rds", snpdens = "fv_625_SNP_freqs_4k.snpden", dr_s = 10728000, dr_e = 12728000, sortit = F, my_title = "Fleckvieh - 0.0625")

cor_individ <- cowplot::plot_grid(bv_all, fv_all,
                                  bv_125, fv_125,
                                  bv_625, fv_625,
                                  nrow = 3,
                                  labels = c("A", "B", "C", "D", "E", "F"),
                                  align = "vh")


#cor_individ

ggsave(
  filename="Correlation-GCContent-Individually.png",
  plot = cor_individ,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 38,
  units = c("cm"),
  limitsize = TRUE,
)


#GC Content correlation random: 
bv_all <- gc_cor(data = "random_bv_all_demoT.rds", snpdens = "bv_all_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, sortit = F, my_title = "Braunvieh - No Cut-off - Random")
bv_125 <- gc_cor(data = "random_bv_125_demoT.rds", snpdens = "bv_125_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, sortit = F, my_title = "Braunvieh - 0.125 - Random")
bv_625 <- gc_cor(data = "random_bv_625_demoT.rds", snpdens = "bv_625_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, sortit = F, my_title = "Braunvieh - 0.0625 - Random")

fv_all <- gc_cor(data = "random_fv_all_demoT.rds", snpdens = "fv_all_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, sortit = F, my_title = "Fleckvieh - No Cut-off - Random")
fv_125 <- gc_cor(data = "random_fv_125_demoT.rds", snpdens = "fv_125_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, sortit = F, my_title = "Fleckvieh - 0.125 - Random")
fv_625 <- gc_cor(data = "random_fv_625_demoT.rds", snpdens = "fv_625_SNP_freqs_4k.snpden", dr_s = 20076000, dr_e = 22076000, sortit = F, my_title = "Fleckvieh - 0.0625 - Random")

cor_individ <- cowplot::plot_grid(bv_all, fv_all,
                                  bv_125, fv_125,
                                  bv_625, fv_625,
                                  nrow = 3,
                                  labels = c("A", "B", "C", "D", "E", "F"),
                                  align = "vh")


cor_individ

ggsave(
  filename="Correlation-GCContent-Random-Individually.png",
  plot = cor_individ,
  device = "png",
  path = NULL,
  scale = 1,
  width = 40,
  height = 38,
  units = c("cm"),
  limitsize = TRUE,
)






