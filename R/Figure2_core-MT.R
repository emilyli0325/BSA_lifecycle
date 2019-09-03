#### Figure 2 whole genome and MT shifting ####### 

# data processing
data <- read.table("BSA.SNP.filter.table", sep = '\t',header = TRUE)
> dim(data)
[1] 7115  120

AD <- data[,seq(5,120,4)]
DP <- data[,seq(6,120,4)]

ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}

dim(AD)
[1] 7115   29
refFre.AD <- matrix(ncol=29,nrow=7115)

for(j in 1:29){
	refFre.AD[,j]<-sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
	}

refFre.AD[DP<30]<-NA
	
colnames(refFre.AD)<-colnames(AD)	
refFre.AD <- cbind(data[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))
write.csv(refFre.AD, "refFre.AD.csv")

###### Figure 2B ######
Figure2B <- ggplot(refFre.AD, aes(x = BSA_042, y = BSA_043))+ geom_point(size = 2, shape = 20, color="grey") + labs(x="Repeat A",  y = "Repeat B") + theme_bw(base_size = 14) + theme(legend.text = element_text(face = "italic"),axis.text.x = element_text(size=14),legend.position = c(0.8, 0.2),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ geom_smooth(method = "lm", se = FALSE, color="red", formula=y~x-1)


###### Figure2A ###### 

refFre.AD <- refFre.AD[,c(1:4,31,33,30,5:29)]
library(reshape2)
refFre.BSAs <- setNames(melt(refFre.AD[,5:32]), c('BSAs', 'Allele frequency of NHP1337'))
Repeats <- c(rep(rep("A",7115),5), rep(c(rep("A",7115),rep("B",7115)),11), c(rep("C",7115)))
LifeCycle <- c(rep("A.Ost.D10",7115), rep("B.Spz",7115),rep("C.Liver",7115),rep("D.ALS_0",7115),rep("E.ALS_1",7115),rep("F.ALS_2",14230),rep("G.ALS_3",14230),rep("H.ALS_4.5",14230),rep("I.ALS_5.5",14230),rep("J.ALS_7.5",14230),rep("K.ALS_9.5",14230),rep("L.ALS_10.5",14230),rep("M.ALS_11.5",14230),rep("N.ALS_12.5",14230),rep("O.ALS_13.5",14230),rep("P.ALS_14.5",21345))
refFre.BSAs.2 <- cbind(Repeats, LifeCycle, refFre.BSAs)

###### Figure2C ###### 

# compare core genome and MT
refFre.BSAs.4<-refFre.BSAs.2[complete.cases(refFre.BSAs.2), ]
colnames(refFre.BSAs.4)[4]<- c("AlleleFrequency")
library(doBy)
sum <- summaryBy(AlleleFrequency ~ BSAs, data = refFre.BSAs.4, FUN = list(mean, median))

comp <- read.csv("Genome_MT.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)

# seperated by Location + 08052019 ### compare to Sanger sequencing result
myColors <- c("Orange","Black","red","blue","green")
names(myColors) <- levels(comp$Location)
colScale <- scale_colour_manual(name = "Location",values = myColors)
shapeScale <- scale_shape_manual(name = "Location",values = c(16, 3,19,15,17))

p2 <- ggplot(comp, aes(x = AsexualCycle, y = Percent.NHP1337, colour = Location, shape = Location)) + geom_point(size = 4) + labs(x="Parasite life cycle",  y = "Allele frequency of NHP1337") + theme_bw(base_size = 16) + theme(legend.text = element_text(face = "italic"),axis.text.x = element_text(size=14),legend.position = c(0.1, 0.3))+ colScale + shapeScale + scale_y_continuous(limits = c(0, 1))

Figure2C <- p2 + geom_smooth(se=FALSE)+scale_y_continuous(limits = c(0, 1))


###### Figure2D ###### 

sco <- read.csv("Genome_MT_s.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)
sco <- cbind(sco, log(sco$Percent.NHP1337/sco$Percent.MKK2835))
names(sco)[9] <- "Ln"

myColors2 <- c("Orange","Black")
names(myColors2) <- levels(sco$Location)
colScale <- scale_colour_manual(name = "Location",values = myColors2)
shapeScale <- scale_shape_manual(name = "Location",values = c(21, 21))

myColors <- c("Orange","Black")
names(myColors) <- levels(sco$Location)
colScale <- scale_colour_manual(name = "Location",values = myColors)
shapeScale <- scale_shape_manual(name = "Location",values = c(16, 3))

Figure2D <- ggplot(sco, aes(x = AsexualCycle, y = Ln, colour = Location))+ geom_point(size = 2) + labs(x="48h asexual cycles",  y = "Ln(genotype ratio)") + theme_bw(base_size = 14) + theme(legend.text = element_text(face = "italic"),axis.text.x = element_text(size=14),legend.position = c(0.8, 0.2),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ colScale + shapeScale + geom_smooth(method = "lm", se = FALSE) + geom_smooth(method = "lm", se = FALSE)


# Figure2A <- ggplot(refFre.BSAs.3, aes(x = `Allele frequency of NHP1337`, y = LifeCycle, fill= Repeats)) +
geom_density_ridges(scale = 2, alpha=0.5)+ theme_ridges() + theme(axis.text.x = element_blank(),legend.position = c(0.8, 0.2))+scale_y_discrete(expand = c(0.01, 0)) +  scale_x_continuous(expand = c(0, 0)) + coord_flip()

# Figure2B <- ggplot(refFre.AD, aes(x = BSA_042, y = BSA_043))+ geom_point(size = 2, shape = 20, color="grey") + labs(x="Repeat A",  y = "Repeat B") + theme_bw(base_size = 14) + theme(legend.text = element_text(face = "italic"),axis.text.x = element_text(size=14),legend.position = c(0.8, 0.2),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ geom_smooth(method = "lm", se = FALSE, color="red", formula=y~x-1, size = 0.6)

# Figure2C <- ggplot(comp, aes(x = AsexualCycle, y = Percent.NHP1337, colour = Location, shape = Location)) + geom_point(size = 3) + labs(x="Parasite life cycle",  y = "Allele frequency of NHP1337") + theme_bw(base_size = 16) + theme(legend.text = element_text(face = "italic"),axis.text.x = element_text(size=14),legend.position = c(0.8, 0.2))+ colScale + shapeScale + geom_smooth(se=FALSE, size = 0.6)+scale_y_continuous(limits = c(0, 1))

# Figure2D <- ggplot(sco, aes(x = AsexualCycle, y = Ln, colour = Location, shape = Location))+ geom_point(size = 3) + labs(x="48h asexual cycles",  y = "Ln(genotype ratio)") + theme_bw(base_size = 14) + theme(legend.text = element_text(face = "italic"),axis.text.x = element_text(size=14),legend.position = c(0.8, 0.2),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ colScale + shapeScale + geom_smooth(method = "lm", se = FALSE, size = 0.6)


library(gridExtra)

pdf('Figure2.pdf', width=12, height=8)
grid.arrange(Figure2A,Figure2B,Figure2C,Figure2D,
widths = c(1, 0.5, 1),layout_matrix = rbind(c(1, 1, 2),c(3, 3, 4))
)
dev.off()



