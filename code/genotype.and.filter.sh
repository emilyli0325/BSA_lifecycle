## GT and filter
java -jar GenomeAnalysisTK.jar \
-T GenotypeGVCFs -R ./PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V BSA.g.vcf \
--useNewAFCalculator \
--sample_ploidy 1 \
-o cross.052018.GATK.vcf

## select good loci according to parent's genotype 
# Seperate SNP and Indel
java -jar GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V cross.052018.GATK.vcf \
-selectType SNP -o cross.052018.SNP.vcf
	
java -jar GenomeAnalysisTK.jar \
-T SelectVariants \
-R ./Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V cross.052018.GATK.vcf  \
-selectType INDEL -o cross.052018.INDEL.vcf

# find SNPs called both parents and BSA-bulks    
java -jar GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R ./Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
   -V cross.052018.SNP.vcf \
   --concordance MKK2835.NHP1337.core.SNP.vcf \
   -o cross.052018.filter.vcf

java -jar GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R ./Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
   -V MKK2835.NHP1337.core.SNP.vcf \
   --concordance cross.052018.SNP.vcf \
   -o MKK2835.NHP1337.filter.vcf   

   
# vcf to table 
   
java -jar GenomeAnalysisTK.jar \
-T VariantsToTable \
-R ./Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V cross.052018.filter.vcf \
-F CHROM -F POS -F REF -F ALT \
-GF AD -GF DP -GF GQ -GF PL \
-o cross.052018.filter.table

