AccClim <- read.csv("~/Downloads/SOO/Lasky2012/AccessionClimateData.csv")
AccClim$ecotype_id <- paste0("X", AccClim$ecotype_id)

CommClim <- read.csv("~/Downloads/SOO/Lasky2012/CommonGard_Clim.csv")
Pheno <- read.csv("~/Downloads/SOO/Fournier-Level.csv")
Pheno$ecotype_id <- paste0("X", Pheno$ecotype_id)
SNPs <- read.csv("~/Downloads/SOO/call_method_75/call_method_75_TAIR9.csv", header = T, skip = 1)

AccClim <- AccClim[AccClim$ecotype_id %in% Pheno$ecotype_id ,] #114
Pheno <- Pheno[Pheno$ecotype_id %in% AccClim$ecotype_id ,] #137 accessions

SeqIndx <- which(colnames(SNPs) %in% Pheno$ecotype_id)
SNPs <- SNPs[c(1,2, SeqIndx)]

#Order the Phenotypes, Climate data and SNPs by ecotype ID
Pheno <- Pheno[match(colnames(SNPs[3:ncol(SNPs)]), Pheno$ecotype_id) ,]
AccClim <- AccClim[match(colnames(SNPs[3:ncol(SNPs)]), AccClim$ecotype_id) ,]

dim(Pheno); dim(AccClim); dim(SNPs)


#Recode SNPs
SNPinfo <- SNPs[1:2]
SNPs <- t(SNPs[3:ncol(SNPs)])

#Recode ATCG with 0,1,2
#minor allele gets coded as 0, major coded as 2
SNPmat <- matrix(0, ncol = ncol(SNPs), nrow = nrow(SNPs))
for (i in 1:ncol(SNPs)){
  tmpCnts <- table(SNPs[,i])
  MajAllele <- names(tmpCnts[which(tmpCnts == max(tmpCnts))])
  MinAllele <- names(tmpCnts[which(tmpCnts == min(tmpCnts))])
  SNPmat[which(SNPs[,i] %in% MajAllele),i] <- 2
  SNPmat[which(SNPs[,i] %in% MinAllele),i] <- 0
}
colnames(SNPmat) <- colnames(SNPs)
row.names(SNPmat) <- row.names(SNPs)

#filter low MAF (0.1) SNPs
freq <- colMeans(SNPmat) / 2
maf <- ifelse(freq > 0.5, 1-freq, freq)
maf.index <- which(maf < 0.1)
length(maf.index)

SNPmat <- SNPmat[, -maf.index] #137 x 200682
SNPinfo <- SNPinfo[-maf.index ,] #200682 x 2

dim(SNPmat); dim(SNPinfo) #114 168407

FullData <- list(SNPmat = SNPmat,
                 MAP = SNPinfo,
                 AccessionClimData = AccClim,
                 CommonGardenClimData = CommClim,
                 Phenotypes = Pheno)

saveRDS(FullData, "~/Downloads/SOO/FullData.rds")
