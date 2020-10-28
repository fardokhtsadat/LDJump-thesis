####Info####
"This file describes the work process for subsetting the cattle-populations 
 (Chapter 6.3 \"Estimation of recombination rates under neutrality and demography using genotyped cattle data\"), 
 which includes:
 
1. Running PLINK on the orgininal VCF-file.
2. Loading the data
3. Subsetting groups

The workflow is shown only for one cow dataset. 
The other dataset is commented out, the procedure stays the same. 
"
####Running PLINK####
"In order to obtain a relationship matrix, PLINK has to be applied onto the VCF-files
 of both data sets (braunvieh.vcf.gz, fleckvieh.vcf.gz).
 
 PLINK was applied using the following options:
 
For Braunvieh:
  --maf 0.01
  --make-rel square
  --vcf braunvieh.vcf.gz
For Fleckvieh:
  --maf 0.01
  --make-rel square
  --vcf fleckvieh.vcf.gz
"
####Setting Working Directory####
setwd("~\Cattle-Subsets")

####Retrieving data#####
#In case PLINK was not applied as described above, the datasets are provided:
cows_rel <- as.data.frame(readxl::read_xlsx("braunvieh_plink.xlsx"))
#cows_rel <- read.csv2("fleckvieh_plink.csv")

####Preparing Data frame####
rownames(cows_rel) <- cows_rel[,1]
cows_rel <- cows_rel[,-1]
for(i in 1:ncol(cows_rel)){
  cows_rel[i,i] <- 0
}

#####Excluding-Relatives#####
testcowsrel <- cows_rel

#Depending on your subset, you comment one of them out, in our example we use the 0.125 threshold.
#df <- apply(testcowsrel, 1, function(x) x > 0.0625)
df <- apply(testcowsrel, 1, function(x) x > 0.125)

#We set a seed
set.seed(1234567)

#We start excluding
while(max(colSums(df))!=0){
  themost <- which(colSums(df) == max(colSums(df)))
  print(themost)
  if(length(themost)>1){themost <- sample(themost, 1)}
  testcowsrel <- testcowsrel[, -themost]
  testcowsrel <- testcowsrel[-themost,]
  #df <- apply(testcowsrel, 1, function(x) x > 0.0625)
  df <- apply(testcowsrel, 1, function(x) x > 0.125)
}


nrow(testcowsrel)
ncol(testcowsrel)
ncol(cows_rel) - ncol(testcowsrel)




#####Write out matrix, retrieve list of individuals#####

#Braunvieh - Write out files:

#We write out the file:
#write.csv2(testcowsrel, "0.0625_42cows.csv")
#write.csv2(testcowsrel, "0.125_57cows.csv")

#Fleckvieh - Write out files:
#write.csv2(testcowsrel, "0.0625_53cows_apr.csv")
#write.csv2(testcowsrel, "0.125_77cows_apr.csv")


#Retrieve list of individuals:
#Braunvieh
cows125 <- read.csv2("bv_0.125_57.csv")
rownames(cows125) <- cows125[,1]
cows125 <- cows125[,-1]

cows625 <- read.csv2("bv_0.0625_42.csv")
rownames(cows625) <- cows625[,1]
cows625 <- cows625[,-1]


fileConn<-file("bv_cows125.txt")
writeLines(rownames(cows125), fileConn)
close(fileConn)

fileConn<-file("bv_cows625.txt")
writeLines(rownames(cows625), fileConn)
close(fileConn)


#Fleckvieh
cows125 <- read.csv2("fv_0.125_77.csv")
rownames(cows125) <- cows125[,1]
cows125 <- cows125[,-1]

cows625 <- read.csv2("fv_0.0625_53.csv")
rownames(cows625) <- cows625[,1]
cows625 <- cows625[,-1]


fileConn<-file("fv_125_ids.txt")
writeLines(rownames(cows125), fileConn)
close(fileConn)

fileConn<-file("fv_625_ids.txt")
writeLines(rownames(cows625), fileConn)
close(fileConn)

