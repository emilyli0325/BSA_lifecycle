### Figure 5 G value, P value, deltaSNP ####
BSA <- read.table("BSA.SNP.filter.2.table", sep = '\t',header = TRUE)

########################################################################################################################
### deltaSNP compare to mean

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
