##########
#OLA data#
##########
Ola <- read.csv("~/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/FileS1.csv")
Accessions <- Ola[1:2]


#Convert VCF to plink
#plink --vcf NGR_607_nsc_maf_003_GWAS.vcf --make-bed --out NGR

##################
#Process SNP data#
##################
library(BGLR)
map.file <- read.table("/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/NGR.bim", 
                       header = F, sep = "\t")
fam.file <- read.table("/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/NGR.fam", 
                       header = F, sep = " ")
PED <- read_bed(bed_file = "/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/NGR.bed",
                fam_file = "/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/NGR.fam",
                bim_file = "/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/NGR.bim")
m <- PED$p
n <- PED$n
PED <- PED$x

##SNPs in PED are coded as 0, 1, 2, 3. 2 is missing data. 1 are heterozygous, 0 and 3 are homozygous for 1/1 and 2/2 for major allele and minor allele respectively
PED[PED == 2] <- NA 
PED[PED == 0] <- 0
PED[PED == 1] <- 1
PED[PED == 3] <- 2

W <- t(matrix(PED, nrow=m, ncol=n, byrow = T)) #n x m
colnames(W) <- map.file$V2
rownames(W) <- fam.file$V2

saveRDS(W, "/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/NGR_SNPmat.rds")



######################
#Check IDs in P and G#
######################

SNPs <- readRDS("/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/NGR_SNPmat.rds")
Lines <- rownames(SNPs)

#For each of the 648 accessions with feild data count the number of matches in the SNP data. This is to determine if they are names according to PI or IS numbers.
NameIndx <- NULL
for(i in 1:nrow(Accessions)){
  NoPIMatches <- sum(Accessions[i,1] %in% Lines)
  NoISMatches <- sum(Accessions[i,2] %in% Lines)
  NameIndx <- rbind(NameIndx, c(NoPIMatches, NoISMatches))
}
Accessions <- cbind(Accessions, NameIndx) # all names according to Taxa

#Order Ola to match SNP data
#Partition Ola into climate, phenotype data
Ola <- Ola[match(row.names(SNPs), Ola$Taxa) ,]
Clim <- Ola[c(12,21:29)]
Pheno <- Ola[c(14,16,18)]

saveRDS(list(Pheno = Pheno,
             Clim = Clim,
             SNPmat = SNPs), "/Users/malachycampbell/Documents/Dropbox/Work/Sorghum_vQTL/OlaData/OlaData.rds")