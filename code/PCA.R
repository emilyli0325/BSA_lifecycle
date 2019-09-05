### Sanger data clean-up, vcf download from: ftp://ngs.sanger.ac.uk/production/pf3k/release_5/

java -jar GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L ./Core_Genome.intervals \
-V T2.SangerR5.SNP.rm.MOI.and.lowCov.0.8GTLoci.MAF0.05.vcf.recode.vcf \
--sample_file SangerSEasia.list -nt 10 \
-o SangerR5.SEasia.vcf ### 678 samples

java -jar GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R ./PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
   -L ./Core_Genome.intervals \
   -V ./Parent.08272019.11Samples.SNP.hardfilter.di.vcf \
   --concordance SangerR5.SEasia.vcf \
   -o BSAParents.vcf 

# merge genotypes
bgzip -c SangerR5.SEasia.vcf > SangerR5.SEasia.vcf.gz
bgzip -c BSAParents.vcf  > BSAParents.vcf.gz
for f in *.vcf.gz; do tabix -p vcf $f; done
bcftools merge --merge none  -R ./Known_sites/Core.bed --output-type z --output Sanger.SEAsia_BSAparents.vcf.gz SangerR5.SEasia.vcf.gz BSAParents.vcf.gz
tabix -p vcf Sanger.SEAsia_BSAparents.vcf.gz

# HetFilter
java -jar GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ./PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L ./Core_Genome.intervals \
--variant Sanger.SEAsia_BSAparents.vcf.gz \
--genotypeFilterExpression "isHet == 1" \
--genotypeFilterName "HetFilter" \
--setFilteredGtToNocall \
--out Sanger.SEAsia_BSAparents.isHet.vcf.gz

gunzip Sanger.SEAsia_BSAparents.isHet.vcf.gz

# remove loci with low GT rate
java -jar GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L ./Core_Genome.intervals \
-V Sanger.SEAsia_BSAparents.isHet.vcf \
--excludeFiltered \
--excludeNonVariants \
--exclude_sample_name AB-BSA-220-NF54/ --exclude_sample_name LA476-1/ --exclude_sample_name LA485-1/ --exclude_sample_name MALawi31/ --exclude_sample_name AB-BSA-222-NHP4026/ \
--maxNOCALLfraction 0.2 \
--removeUnusedAlternates \
-nt 10 \
-o Sanger.SEAsia_BSAparents.0.8GT.isHet.vcf


# MAF 0.05 T4
vcftools --vcf Sanger.SEAsia_BSAparents.0.8GT.isHet.vcf --max-missing 0.8 --maf 0.05 --max-maf 0.95 --recode --out Sanger.SEAsia_BSAparents.clean.vcf

### dataset summary
VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf Sanger.SEAsia_BSAparents.0.8GT.isHet.vcf
        --maf 0.05
        --max-maf 0.95
        --max-missing 0.8
        --out Sanger.SEAsia_BSAparents.clean.vcf
        --recode

After filtering, kept 684 out of 684 Individuals
Outputting VCF file...
After filtering, kept 13329 out of a possible 17316 Sites

### PCA ###
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(ggplot2)
library(ggridges)
library(gridExtra)

setwd("")

seqVCF2GDS("Sanger.SEAsia_BSAparents.clean.vcf.recode.vcf", "Sanger.SEAsia_BSAparents.gds")
genofile <- seqOpen("Sanger.SEAsia_BSAparents.gds")
pca <- snpgdsPCA(genofile, autosome.only=FALSE, num.thread=4)
> pc.percent <- pca$varprop*100
> head(round(pc.percent, 2))
[1] 9.27 4.38 3.68 2.38 1.88 1.49

tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    
    EV2 = pca$eigenvect[,2],
    EV3 = pca$eigenvect[,3],
    EV4 = pca$eigenvect[,4],	
	EV5 = pca$eigenvect[,5],
    EV6 = pca$eigenvect[,6],
    stringsAsFactors = FALSE)

write.csv(tab, "PCA1-6.csv")
### add country information 

tab.C <- read.csv("PCA1-6.csv", header = T)

# seperated by Location
myColors <- c("yellow","grey","lightgreen", "pink", "purple","#FFDB6D", "lightblue")
names(myColors) <- levels(tab.C$Country)
colScale <- scale_colour_manual(name = "Country",values = myColors)
shapeScale <- scale_shape_manual(name = "Country",values = c(20,20,21,22,24,25,23))
filScale <- scale_fill_manual(name = "Country",values = myColors)

MKK2835.2G: geom_point(data=tab.C[2, ], aes(x=EV1,y=EV2), shape = 24, color = "black", fill ="red", size=4)
NHP1337.12C: geom_point(data=tab.C[3, ], aes(x=EV1,y=EV2), shape = 24,  color = "black", fill ="blue", size=4)

P.PCA1.2 <- ggplot(tab.C, aes(x=EV1,y=EV2, fill = Country, shape = Country)) + geom_point(size = 2) + labs(x="PC1 (9.27%)",  y = "PC2 (4.38%)") + theme_bw(base_size = 16) + shapeScale + filScale + geom_point(data=tab.C[2, ], aes(x=EV1,y=EV2), shape = 24, color = "black", fill ="red", size=4) + geom_point(data=tab.C[3, ], aes(x=EV1,y=EV2), shape = 24,  color = "black", fill ="blue", size=4)

P.PCA2.3 <- ggplot(tab.C, aes(x=EV3,y=EV2, fill = Country, shape = Country)) + geom_point(size = 2) + labs(x="PC3 (3.68%)",  y = "PC2 (4.38%)") + theme_bw(base_size = 16) + shapeScale + filScale + geom_point(data=tab.C[2, ], aes(x=EV1,y=EV2), shape = 24, color = "black", fill ="red", size=4) + geom_point(data=tab.C[3, ], aes(x=EV1,y=EV2), shape = 24,  color = "black", fill ="blue", size=4)

pdf('PCA.pdf', width=16, height=6)
grid.arrange(P.PCA1.2,P.PCA2.3,
widths = c(1, 1),layout_matrix = cbind(1,2)
)
dev.off()










