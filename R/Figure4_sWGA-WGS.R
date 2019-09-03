##  sWGA vs WGS  #######

Mixture <- read.table("Mixture_sWGA_WGS.SNP.filter.table", sep = '\t',header = TRUE)

AD <- Mixture[,seq(5,36,4)]
DP <- Mixture[,seq(6,36,4)]
ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}

dim(AD)
[1] 9045    8

refFre.AD <- matrix(ncol=8,nrow=9045)

for(j in 1:8){
	refFre.AD[,j]<-sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
	}

refFre.AD[DP<30]<-NA
	
colnames(refFre.AD)<-colnames(AD)	
refFre.AD <- cbind(Mixture[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

### remove outliers ###

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
			   

BSA_044_sWGA.filter <- outlierByMAD(refFre.AD$BSA_044_sWGA,50)
BSA_044.filter <- outlierByMAD(refFre.AD$BSA_044,50)

refFre.AD.2 <- cbind(refFre.AD[,1:4],BSA_044_sWGA.filter,BSA_044.filter)

SAnalysis.1 <-
    function(SNPset,
        windowSize = 1e5,
        ...)
    {
        SNPset <- SNPset %>%
            dplyr::group_by(CHROM) %>%
            dplyr::mutate(tricubeBSA_044 = tricubeStat(POS = POS, Stat = BSA_044.filter, windowSize, ...))
		SNPset <- SNPset %>%
            dplyr::group_by(CHROM) %>%
            dplyr::mutate(tricubeBSA_044_sWGA = tricubeStat(POS = POS, Stat = BSA_044_sWGA.filter, windowSize, ...))
		
		return(as.data.frame(SNPset))	
    }
refFre.AD.2 <- SAnalysis.1(refFre.AD.2)	

p0 <- ggplot(data = refFre.AD.2) +scale_x_continuous(breaks = seq(from = 0,to = max(refFre.AD.2$POS), by = 10^(floor(log10(max(refFre.AD$POS))))), labels=format_genomic())+ylim(0, 1.3)+ facet_grid(~ CHROM, scales = "free_x", space = "free_x")+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())

p1 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_044.filter"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_044"),color = "black",size=0.5)+ylab("Before sWGA")

p2 <- p0 + geom_point(aes_string(x = "POS", y = "BSA_044_sWGA.filter"),color = "orange",size=0.5)+geom_line(aes_string(x = "POS", y = "tricubeBSA_044_sWGA"),color = "black",size=0.5)+ylab("After sWGA")

p3 <- ggplot(refFre.AD, aes(x = BSA_044.filter, y = BSA_044_sWGA.filter))+ geom_point(size = 2, shape = 20, color="grey") + labs(x="Before sWGA",  y = "After sWGA") + theme_bw(base_size = 14) + theme(legend.text = element_text(face = "italic"),axis.text.x = element_text(size=14),legend.position = c(0.8, 0.2),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ geom_smooth(method = "lm", se = FALSE, color="red", formula=y~x-1, size = 0.6)


library(gridExtra)

pdf('Figure4_sWGA-WGS-2.pdf', width=8, height=7)
grid.arrange(p1,p2,p3,
widths = c(1.5, 1, 1.5),heights = c(1,1,1.5), layout_matrix = rbind(c(1, 1,1),c(2, 2,2),c(3,NA, NA))
)
dev.off()
