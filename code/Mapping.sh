#!/bin/bash
#mkdir SAM
#mkdir sorted.bam
#mkdir dedup.sorted.bam 
#mkdir BQSR
#mkdir recal.bam
#mkdir metrics
#mkdir g.vcf

#files=(...)


for file in ${files[*]}
do

echo "mapping $file..."

bwa mem /master/xli/Index/Pfal32_bwa_index/Pfal32 fq/$file.R1.fastq.gz fq/$file.R2.fastq.gz -t 16 -M -R "@RG\tID:$file\tLB:$file\tPL:ILLUMINA\tPM:HISEQ\tSM:$file/" > SAM/$file.sam

java -jar picard.jar SortSam \
     INPUT=SAM/$file.sam \
     OUTPUT=sorted.bam/$file.sorted.bam \
     SORT_ORDER=coordinate

java -jar picard.jar MarkDuplicates \
     INPUT=sorted.bam/$file.sorted.bam \
     OUTPUT=dedup.sorted.bam/$file.dedup.sorted.bam \
     METRICS_FILE=metrics/$file.metrics.txt

cd dedup.sorted.bam

java -jar picard.jar BuildBamIndex \
     INPUT=$file.dedup.sorted.bam

cd ..

echo "BQSR $file..."	 
	java -jar GenomeAnalysisTK.jar\
		 -T BaseRecalibrator\
		 -R ~/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta\
		 -I dedup.sorted.bam/$file.dedup.sorted.bam\
		 -knownSites /master/xli/Index/Known_sites/3d7_hb3.combined.final.karo.sort.vcf\
		 -knownSites /master/xli/Index/Known_sites/7g8_gb4.combined.final.karo.sort.vcf\
		 -knownSites /master/xli/Index/Known_sites/hb3_dd2.combined.final.karo.sort.vcf\
		 -o BQSR/$file.recal.table
	java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar\
		 -T PrintReads\
		 -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta\
		 -I dedup.sorted.bam/$file.dedup.sorted.bam\
		 -BQSR BQSR/$file.recal.table\
		 -o recal.bam/$file.recal.bam 
		 
echo "Variant calling $file..."
	java -jar GenomeAnalysisTK.jar\
         -T HaplotypeCaller\
         -R ~/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta\
         -I recal.bam/$file.recal.bam\
         --emitRefConfidence BP_RESOLUTION\
         -o g.vcf/$file.g.vcf\
         --useNewAFCalculator \
         --sample_ploidy 1 \
         --dontUseSoftClippedBases \
         --num_cpu_threads_per_data_thread 16
        echo "$file Variant calling done!"
		
done
