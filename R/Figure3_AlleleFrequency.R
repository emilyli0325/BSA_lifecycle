### Figure 3 Allele frequency ##### 

### remove outliers ###

### adopted from QTLseqr (Mansfeld BN, Grumet R. QTLseqr: An R package for bulk segregant analysis with next-generation sequencing. The Plant Genome. 2018)
format_genomic <- function(...) {
      function(x) {
            limits <- c(1e0,   1e3, 1e6)
            #prefix <- c("","Kb","Mb")
            # Vector with array indices according to position in intervals
            i <- findInterval(abs(x), limits)
            # Set prefix to " " for very small values < 1e-24
            i <- ifelse(i==0, which(limits == 1e0), i)
            paste(format(round(x/limits[i], 1),
                         trim=TRUE, scientific=FALSE, ...)
                #  ,prefix[i]
            )
      }
}

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

####### tricube smooth #######		
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
            dplyr::mutate(tricubeSetA = tricubeStat(POS = POS, Stat = setA, windowSize, ...))
		SNPset <- SNPset %>%
			dplyr::mutate(tricubeSetB = tricubeStat(POS = POS, Stat = setB, windowSize, ...))
        
        return(as.data.frame(SNPset))
    }

	
SAnalysis.1 <-
    function(SNPset,
        windowSize = 1e5,
        ...)
    {
 
        SNPset <- SNPset %>%
            dplyr::group_by(CHROM) %>%
            dplyr::mutate(tricubeBSA_039 = tricubeStat(POS = POS, Stat = BSA_039, windowSize, ...))
		return(as.data.frame(SNPset))	
    }
refFre.AD <- SAnalysis.1(refFre.AD)	

write.csv(refFre.AD, "Figure 3 Allele frequency.csv")
	
p0 <- ggplot(data = refFre.AD) +scale_x_continuous(breaks = seq(from = 0,to = max(refFre.AD$POS), by = 10^(floor(log10(max(refFre.AD$POS))))), labels=format_genomic())+ylim(0, 1.3)+ facet_grid(~ CHROM, scales = "free_x", space = "free_x")+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())+theme(axis.title.x=element_blank(),axis.line.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.line.y=element_blank())

p1 <- p0 + geom_point(aes_string(x = "POS", y = "OocystD10"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeOocystD10"),color = "red",size=1)+geom_hline(yintercept = 0.8135593,color = "black",linetype="dashed",alpha = 1,size=1)

p2 <- p0 + geom_point(aes_string(x = "POS", y = "SGSpz"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeSGSpz"),color = "red",size=1)+geom_hline(yintercept = 0.8,color = "black",linetype="dashed",alpha = 1,size=1)

p3 <- p0 + geom_point(aes_string(x = "POS", y = "Liver"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeLiver"),color = "red",size=1)+geom_hline(yintercept = 0.9021277,color = "black",linetype="dashed",alpha = 1,size=1)

pA0 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_020"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_020"),color = "red",size=1)+geom_hline(yintercept = 0.8484848,color = "black",linetype="dashed",alpha = 1,size=1)

pA1 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_021"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_021"),color = "red",size=1)+geom_hline(yintercept = 0.8366337,color = "black",linetype="dashed",alpha = 1,size=1)

pA3 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_024"),color = "black",size=0.5) + geom_point(aes_string(x = "POS", y = "BSA_025"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_024"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "tricubeBSA_025"),color = "red",size=1)+geom_hline(yintercept = 0.8717201,color = "black",linetype="dashed",alpha = 1,size=1)


pA5.5 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_028"),color = "black",size=0.5) + geom_point(aes_string(x = "POS", y = "BSA_029"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_028"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "tricubeBSA_029"),color = "red",size=1)+geom_hline(yintercept = 0.8142857,color = "black",linetype="dashed",alpha = 1,size=1)

pA7.5 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_030"),color = "black",size=0.5) + geom_point(aes_string(x = "POS", y = "BSA_031"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_030"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "tricubeBSA_031"),color = "red",size=1)+geom_hline(yintercept = 0.6590909,color = "black",linetype="dashed",alpha = 1,size=1)

pA10.5 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_034"),color = "black",size=0.5) + geom_point(aes_string(x = "POS", y = "BSA_035"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_034"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "tricubeBSA_035"),color = "red",size=1)+geom_hline(yintercept = 0.5454545,color = "black",linetype="dashed",alpha = 1,size=1)

pA12.5 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_038"),color = "black",size=0.5) + geom_point(aes_string(x = "POS", y = "BSA_039"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_038"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "tricubeBSA_039"),color = "red",size=1)+geom_hline(yintercept = 0.4761905,color = "black",linetype="dashed",alpha = 1,size=1)

pA14.5 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_042"),color = "black",size=0.5) + geom_point(aes_string(x = "POS", y = "BSA_044"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_042"),color = "blue",size=1)+geom_line(aes_string(x = "POS", y = "tricubeBSA_044"),color = "red",size=1)+geom_hline(yintercept = 0.4742268,color = "black",linetype="dashed",alpha = 1,size=1)

pdf('Figure3.pdf', width=10, height=15)
ggarrange(p1, p2, p3, pA0, pA1, pA3, pA5.5, pA7.5, pA10.5, pA12.5, pA14.5, ncol = 1, nrow = 11)
dev.off()
