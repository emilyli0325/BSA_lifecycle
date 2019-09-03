### Figure 5 G value, P value, deltaSNP ####
BSA <- read.table("BSA.SNP.filter.2.table", sep = '\t',header = TRUE)

########################################################################################################################
### Oocyst #####
dfOo <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_Oo",lowBulk = "OocystD10")
df_filtOo <-filterSNPs(SNPset = dfOo,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filtOo <- runQTLseqAnalysis(SNPset = df_filtOo,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filtOo <- runGprimeAnalysis(df_filtOo,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filtOo, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)


#plotQTLStats(SNPset = df_filtOo, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#pOo <- plotQTLStats(SNPset = df_filtOo, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()

########################################################################################################################
#### spz ######
dfSpz <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_Spz",lowBulk = "SGSpz")
df_filtSpz <-filterSNPs(SNPset = dfSpz,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filtSpz <- runQTLseqAnalysis(SNPset = df_filtSpz,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filtSpz <- runGprimeAnalysis(df_filtSpz,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filtSpz, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)

#plotQTLStats(SNPset = df_filtSpz, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#pSpz <- plotQTLStats(SNPset = df_filtSpz, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()

########################################################################################################################
#### Liver ######
dfLiver <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_liver",lowBulk = "Liver")
df_filtLiver <-filterSNPs(SNPset = dfLiver,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 30,minGQ = 99)
df_filtLiver <- runQTLseqAnalysis(SNPset = df_filtLiver,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filtLiver <- runGprimeAnalysis(df_filtLiver,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filtLiver, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)


#plotQTLStats(SNPset = df_filtLiver, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#pLiver <- plotQTLStats(SNPset = df_filtLiver, var = "deltaSNP", line=FALSE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()

########################################################################################################################
#### in vivo ######
dfvivo <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_vivo",lowBulk = "BSA_020")
df_filtvivo <-filterSNPs(SNPset = dfvivo,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filtvivo <- runQTLseqAnalysis(SNPset = df_filtvivo,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filtvivo <- runGprimeAnalysis(df_filtvivo,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filtvivo, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)


#plotQTLStats(SNPset = df_filtvivo, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#pVivo <- plotQTLStats(SNPset = df_filtvivo, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()

########################################################################################################################
#### day 23 ######
df21 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_23",lowBulk = "BSA_021")
df_filt21 <-filterSNPs(SNPset = df21,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt21 <- runQTLseqAnalysis(SNPset = df_filt21,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt21 <- runGprimeAnalysis(df_filt21,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt21, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)


#plotQTLStats(SNPset = df_filt21, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#pd23 <- plotQTLStats(SNPset = df_filt21, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()

########################################################################################################################
#### day 27 ######
df24 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_27",lowBulk = "BSA_024")
df_filt24 <-filterSNPs(SNPset = df24,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt24 <- runQTLseqAnalysis(SNPset = df_filt24,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt24 <- runGprimeAnalysis(df_filt24,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt24, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)


#plotQTLStats(SNPset = df_filt24, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#plotQTLStats(SNPset = df_filt24, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()

df25 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_27",lowBulk = "BSA_025")
df_filt25 <-filterSNPs(SNPset = df25,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt25 <- runQTLseqAnalysis(SNPset = df_filt25,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt25 <- runGprimeAnalysis(df_filt25,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt25, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)

#plotQTLStats(SNPset = df_filt25, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#pd27 <- plotQTLStats(SNPset = df_filt25, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()+geom_line(data = df_filt24, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)


########################################################################################################################
#### day 32 ######
df28 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_32",lowBulk = "BSA_028")
df_filt28 <-filterSNPs(SNPset = df28,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt28 <- runQTLseqAnalysis(SNPset = df_filt28,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt28 <- runGprimeAnalysis(df_filt28,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt28, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)


#plotQTLStats(SNPset = df_filt28, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#plotQTLStats(SNPset = df_filt28, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()

df29 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_32",lowBulk = "BSA_029")
df_filt29 <-filterSNPs(SNPset = df29,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt29 <- runQTLseqAnalysis(SNPset = df_filt29,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt29 <- runGprimeAnalysis(df_filt29,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt29, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)### p value


#plotQTLStats(SNPset = df_filt29, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
#pd32 <- plotQTLStats(SNPset = df_filt29, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()+geom_line(data = df_filt28, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)

########################################################################################################################
#### day 36 ######
df30 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_36",lowBulk = "BSA_030")
df_filt30 <-filterSNPs(SNPset = df30,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt30 <- runQTLseqAnalysis(SNPset = df_filt30,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt30 <- runGprimeAnalysis(df_filt30,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt30, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 

df31 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_36",lowBulk = "BSA_031")
df_filt31 <-filterSNPs(SNPset = df31,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt31 <- runQTLseqAnalysis(SNPset = df_filt31,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt31 <- runGprimeAnalysis(df_filt31,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt31, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 

########################################################################################################################
#### day 42 ######
df34 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_42",lowBulk = "BSA_034")
df_filt34 <-filterSNPs(SNPset = df34,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt34 <- runQTLseqAnalysis(SNPset = df_filt34,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt34 <- runGprimeAnalysis(df_filt34,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt34, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 

df35 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_42",lowBulk = "BSA_035")
df_filt35 <-filterSNPs(SNPset = df35,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt35 <- runQTLseqAnalysis(SNPset = df_filt35,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt35 <- runGprimeAnalysis(df_filt35,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt35, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 

########################################################################################################################
#### day 46 ######
df38 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_42",lowBulk = "BSA_038")
df_filt38 <-filterSNPs(SNPset = df38,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt38 <- runQTLseqAnalysis(SNPset = df_filt38,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt38 <- runGprimeAnalysis(df_filt38,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt38, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 

df39 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_42",lowBulk = "BSA_039")
df_filt39 <-filterSNPs(SNPset = df39,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt39 <- runQTLseqAnalysis(SNPset = df_filt39,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(99))
df_filt39 <- runGprimeAnalysis(df_filt39,windowSize = 1e5,outlierFilter = "Hampel",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filt39, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 

########################################################################################################################
#### day 50 ######
df42 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_50",lowBulk = "BSA_042")
df_filt42 <-filterSNPs(SNPset = df42,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt42 <- runQTLseqAnalysis(SNPset = df_filt42,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt42 <- runGprimeAnalysis(df_filt42,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)### Threshold good
plotQTLStats(SNPset = df_filt42, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 
write.csv(df_filt42, "df_filt42.csv")


df44 <-importFromGATK(file = "BSA.SNP.filter.2.table",highBulk = "Mock_50",lowBulk = "BSA_044")
df_filt44 <-filterSNPs(SNPset = df44,minTotalDepth = 120,maxTotalDepth = 1000,minSampleDepth = 20)
df_filt44 <- runQTLseqAnalysis(SNPset = df_filt44,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filt44 <- runGprimeAnalysis(df_filt44,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)### Threshold too low
plotQTLStats(SNPset = df_filt44, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+ylim(0, 50) 
write.csv(df_filt44, "df_filt44.csv")

write.csv(df_filt28, "df_filt28.csv")
write.csv(df_filt29, "df_filt29.csv")
write.csv(df_filt30, "df_filt30.csv")
write.csv(df_filt31, "df_filt31.csv")
write.csv(df_filt34, "df_filt34.csv")
write.csv(df_filt35, "df_filt35.csv")
write.csv(df_filt38, "df_filt38.csv")
write.csv(df_filt39, "df_filt39.csv")

########################################################################################################################
### deltaSNP compare to mean
#library(ggpubr)

pOo <- plotQTLStats(SNPset = df_filtOo, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())
pSpz <- plotQTLStats(SNPset = df_filtSpz, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())
pLiver <- plotQTLStats(SNPset = df_filtLiver, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())
pVivo <- plotQTLStats(SNPset = df_filtvivo, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd23 <- plotQTLStats(SNPset = df_filt21, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd27 <- plotQTLStats(SNPset = df_filt25, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+geom_line(data = df_filt24, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd32 <- plotQTLStats(SNPset = df_filt29, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+geom_line(data = df_filt28, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd36 <- plotQTLStats(SNPset = df_filt30, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+geom_line(data = df_filt31, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd42 <- plotQTLStats(SNPset = df_filt34, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+geom_line(data = df_filt35, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd46 <- plotQTLStats(SNPset = df_filt38, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+geom_line(data = df_filt39, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd50 <- plotQTLStats(SNPset = df_filt42, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = FALSE)+theme_classic()+geom_line(data = df_filt44, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())

pdf('deltaSNP-2.pdf', width=15, height=20)
ggarrange(pOo, pSpz, pLiver, pVivo, pd23, pd27, pd32, pd36, pd42, pd46, pd50,
          labels = c("Oocyst", "Spz", "Liver","in vivo", "d23","d27","d32","d36","d42","d46","d50"),
          ncol = 1, nrow = 11)
dev.off()

pdf('deltaSNP-3.pdf', width=15, height=40)
ggarrange(pOo, pSpz, pLiver, pVivo, pd23, pd27, pd32, pd36, pd42, pd46, pd50,
          labels = c("Oocyst", "Spz", "Liver","in vivo", "d23","d27","d32","d36","d42","d46","d50"),
          ncol = 1, nrow = 11)
dev.off()

pd50 <- plotQTLStats(SNPset = df_filt42, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()+geom_line(data = df_filt44, aes_string(x = "POS", y = "tricubeDeltaSNP"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank(),legend.position="none")

pdf('deltaSNP-day50.pdf', width=8, height=3)
ggarrange(pd50,
          ncol = 1, nrow = 1)
dev.off()


########################################################################################################################
### Gprime compare to mean
pOo <- plotQTLStats(SNPset = df_filtOo, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pSpz <- plotQTLStats(SNPset = df_filtSpz, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pLiver <- plotQTLStats(SNPset = df_filtLiver, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pVivo <- plotQTLStats(SNPset = df_filtvivo, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd23 <- plotQTLStats(SNPset = df_filt21, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd27 <- plotQTLStats(SNPset = df_filt24, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+ geom_line(data = df_filt25, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd32 <- plotQTLStats(SNPset = df_filt28, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+ geom_line(data = df_filt29, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd36 <- plotQTLStats(SNPset = df_filt30, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+ geom_line(data = df_filt31, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd42 <- plotQTLStats(SNPset = df_filt34, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+ geom_line(data = df_filt35, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd46 <- plotQTLStats(SNPset = df_filt38, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+ geom_line(data = df_filt39, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())
pd50 <- plotQTLStats(SNPset = df_filt42, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0, 50)+ geom_line(data = df_filt44, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.5)+theme(strip.background = element_blank(),strip.text.x = element_blank())


pdf('G.pdf', width=8, height=24)
ggarrange(pOo, pSpz, pLiver, pVivo, pd23, pd27, pd32, pd36, pd42, pd46, pd50,
          labels = c("Oocyst", "Spz", "Liver","in vivo", "d23","d27","d32","d36","d42","d46","d50"),
          ncol = 1, nrow = 11)
dev.off()

pdf('G_3.pdf', width=8, height=24)
ggarrange(pOo, pSpz, pLiver, pVivo, pd23, pd27, pd32, pd36, pd42, pd46, pd50,
          labels = c("Oocyst", "Spz", "Liver","in vivo", "d23","d27","d32","d36","d42","d46","d50"),
          ncol = 1, nrow = 11)
dev.off()

########################################################################################################################
### selection coefficient 

setwd("C:/Users/xli/Documents/Emily/P01/BSA_cross_lifeCycle/1.Test_Run_Marina/NGS_analysis/refNHP1337")
#### load refFre.AD
BSA <- read.table("BSA.SNP.filter.table", sep = '\t',header = TRUE)
> dim(BSA)
[1] 7115  120

AD <- BSA[,seq(5,120,4)]
DP <- BSA[,seq(6,120,4)]
ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}
dim(AD)
[1] 7115   29
refFre.AD <- matrix(ncol=29,nrow=7115)
for(j in 1:29){
	refFre.AD[,j]<-sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
	}
refFre.AD[DP<30]<-NA
colnames(refFre.AD)<-colnames(AD)	
refFre.AD <- cbind(BSA[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

### load SC
Ln.ratio <- read.csv("Ln.ratio.csv", header = TRUE, row.names=1)
Ln.ratio.0 <- Ln.ratio[3:24,]
setA <- Ln.ratio.0[Ln.ratio.0$RepeatID == "A", ]
setB <- Ln.ratio.0[Ln.ratio.0$RepeatID == "B", ]

sc <- matrix(ncol=2,nrow=7115)
for(j in 1:7115){
    ifelse(
	       all(is.na(setA[,j+2])), 
		   sc[j,1] <- NA, 
		   sc[j,1] <- -coef(lm( as.numeric(setA[,j+2]) ~ setA$ALS))[2]
	      )
	ifelse(
	       all(is.na(setB[,j+2])), 
		   sc[j,2] <- NA, 
		   sc[j,2] <- -coef(lm( as.numeric(setB[,j+2]) ~ setB$ALS))[2]
	      )
	}	
sc <- cbind(refFre.AD[,1:4],sc)
colnames(sc)[5:6] <- c("setA","setB")

outliersMAD <- function(data, MADCutOff = 2.5, replace = NA, values = FALSE, bConstant = 1.4826, digits = 2) {
    absMADAway <- abs(   (data - median(data, na.rm = T))  /  mad(data, constant = bConstant, na.rm = T)  )
   
    data[absMADAway > MADCutOff] <- replace
    
    if (values == TRUE) { 
        return(round(absMADAway, digits)) #if values == TRUE, return number of mads for each value
    } else {
        return(round(data, digits)) #otherwise, return values with outliers replaced
    }
}

outlierByMAD <- function (x, k){
	n <- length(x)
     y <- x
     
 for (i in (k + 1):(n - k)) {
 	
    data <- x[(i - k):(i + k)]
    y[i] <- outliersMAD(data)[k+1]}
    
   return(y)
                               }
							   
SetA.filter <- outlierByMAD(sc$setA,50)
SetB.filter <- outlierByMAD(sc$setB,50)					   
sc.filter <- cbind(sc, SetA.filter, SetB.filter)
	
tricubeStat <- function(POS, Stat, windowSize = 1e5, ...)
{
    if (windowSize <= 0)
        stop("A positive smoothing window is required")
    stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0), ...), POS)
}

SAnalysis <-
    function(SNPset,
        windowSize = 1e5,
        ...)
    {
 
        SNPset <- SNPset %>%
            dplyr::group_by(CHROM) %>%
            dplyr::mutate(tricubeSetA = tricubeStat(POS = POS, Stat = SetA.filter, windowSize, ...))
		SNPset <- SNPset %>%
			dplyr::mutate(tricubeSetB = tricubeStat(POS = POS, Stat = SetB.filter, windowSize, ...))
        
        return(as.data.frame(SNPset))
    }

	
sc.filter <- SAnalysis(sc.filter)


> mean(c(sc.filter$tricubeSetA,sc.filter$tricubeSetB))
[1] 0.196668


SetA.format <- sc.filter$tricubeSetA - 0.197
SetB.format <- sc.filter$tricubeSetB - 0.197
sc.filter.2 <- cbind(sc.filter, SetA.format, SetB.format)

write.csv(sc.filter.2, "sc.filter.2.csv")

pSC.1 <- ggplot(data = sc.filter.2,aes(x = POS, y = SetA.format)) +scale_x_continuous(breaks = seq(from = 0,to = max(sc.filter.2$POS), by = 10^(floor(log10(max(sc.filter.2$POS))))), name = "Genomic Position (Mb)") +ylab("selection coefficient")+ facet_grid(~ CHROM, scales = "free_x", space = "free_x")+theme(plot.margin = margin(b = 10,l = 20,r = 20,unit = "pt"))+
geom_line(aes_string(x = "POS", y = "SetA.format"),color = "orange",size=1) +
geom_line(aes_string(x = "POS", y = "SetB.format"),color = "black",size=0.5)+geom_hline(yintercept = 0,color = "red",size = 1,alpha = 0.4)+theme_classic()
 
pSC.2 <- ggplot(data = sc.filter.2,aes(x = POS, y = tricubeSetA)) +scale_x_continuous(breaks = seq(from = 0,to = max(sc.filter.2$POS), by = 10^(floor(log10(max(sc.filter.2$POS))))), name = "Genomic Position (Mb)") +ylab("selection coefficient")+ facet_grid(~ CHROM, scales = "free_x", space = "free_x")+theme(plot.margin = margin(b = 10,l = 20,r = 20,unit = "pt"))+
geom_line(aes_string(x = "POS", y = "tricubeSetA"),color = "orange",size=1) +
geom_line(aes_string(x = "POS", y = "tricubeSetB"),color = "black",size=0.5)+geom_hline(yintercept = 0.197,color = "red",size = 1,alpha = 0.4)+theme_classic()

pdf('SC.pdf', width=8, height=6)
ggarrange(pSC.1,pSC.2,
          ncol = 1, nrow = 2)
dev.off()
