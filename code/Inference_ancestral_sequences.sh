### Orthologous group determination and alignment
### blast index
export PATH=$PATH:/master/xli/software/blast/ncbi-blast-2.9.0+/bin
makeblastdb -in PlasmoDB-32_Pfalciparum3D7_AnnotatedCDSs.fasta -out databaseBLAST.3D7.R32 -dbtype nucl

### blast
blastn -db /master/xli/Index/AncestralAlleles/3D7.r32/databaseBLAST.3D7.R32 -query /master/xli/Index/AncestralAlleles/P.reichenowi.G01/PlasmoDB-45_PreichenowiG01_AnnotatedCDSs.fasta -out PreichenowiG01_3D7.out -num_threads 10 -outfmt 6 -max_target_seqs 1

blastn -db /master/xli/Index/AncestralAlleles/3D7.r32/databaseBLAST.3D7.R32 -query /master/xli/Index/AncestralAlleles/P.adleri.G01/PlasmoDB-45_PadleriG01_AnnotatedCDSs.fasta -out P.adleri.G01_3D7.out -num_threads 10 -outfmt 6 -max_target_seqs 1

blastn -db /master/xli/Index/AncestralAlleles/3D7.r32/databaseBLAST.3D7.R32 -query /master/xli/Index/AncestralAlleles/P.billcollinsi.G01/PlasmoDB-45_PbillcollinsiG01_AnnotatedCDSs.fasta -out P.billcollinsi.G01_3D7.out -num_threads 10 -outfmt 6 -max_target_seqs 1

blastn -db /master/xli/Index/AncestralAlleles/3D7.r32/databaseBLAST.3D7.R32 -query /master/xli/Index/AncestralAlleles/P.blacklocki.G01/PlasmoDB-45_PblacklockiG01_AnnotatedCDSs.fasta -out P.blacklocki.G01_3D7.out -num_threads 10 -outfmt 6 -max_target_seqs 1

blastn -db /master/xli/Index/AncestralAlleles/3D7.r32/databaseBLAST.3D7.R32 -query /master/xli/Index/AncestralAlleles/P.gaboni.G01/PlasmoDB-45_PgaboniG01_AnnotatedCDSs.fasta -out P.gaboni.G01_3D7.out -num_threads 10 -outfmt 6 -max_target_seqs 1

### extract target sequences
xargs samtools faidx 3D7.r32/PlasmoDB-32_Pfalciparum3D7_AnnotatedCDSs.fasta < 3d7.list >> BSA1_QTLregion.Merge.fasta
xargs samtools faidx P.adleri.G01/PlasmoDB-45_PadleriG01_AnnotatedCDSs.fasta < PadleriG01.list >> BSA1_QTLregion.Merge.fasta
xargs samtools faidx P.billcollinsi.G01/PlasmoDB-45_PbillcollinsiG01_AnnotatedCDSs.fasta < PbillcollinsiG01.list >> BSA1_QTLregion.Merge.fasta
xargs samtools faidx P.blacklocki.G01/PlasmoDB-45_PblacklockiG01_AnnotatedCDSs.fasta < PblacklockiG01.list >> BSA1_QTLregion.Merge.fasta
xargs samtools faidx P.gaboni.G01/PlasmoDB-45_PgaboniG01_AnnotatedCDSs.fasta < PgaboniG01.list >> BSA1_QTLregion.Merge.fasta
xargs samtools faidx P.reichenowi.G01/PlasmoDB-45_PreichenowiG01_AnnotatedCDSs.fasta < PreichenowiG01.list >> BSA1_QTLregion.Merge.fasta

samtools faidx BSA1_QTLregion.Merge.fasta PF3D7_1227300.1 PRG01_1230600-t36_1 PADL01_1227900-t36_1 PBILCG01_1227300-t36_1 PBLACG01_1226700-t36_1 PGABG01_1226100-t36_1 > PF3D7_1227300.fasta

### change fasta header
conda create -n BSA1
source activate BSA1
conda install python==3.6
conda install -c conda-forge biopython
conda install -c bioconda seqkit

seqkit replace -p "\_" -r ' ' PF3D7_1462400.fasta | seqkit replace -p "\s.+" > PF3D7_1462400_new.fasta

## sequence alignment : http://wasabiapp.org/software/prank/#Methods
/master/xli/software/PRANK/development/src/prank -d=PF3D7_1227300_new.fasta -o=PF3D7_1227300_prank -t=Tree.dnd +F -codon -once -f=nexus -showall
/master/xli/software/PRANK/development/src/prank -d=PF3D7_1227400_new.fasta -o=PF3D7_1227400_prank -t=Tree.1227400.dnd +F -codon -once -f=nexus -showall

#### vcf to allele frequency

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
-V /data/infectious/malaria_XUE/sWGA/mergeThai.Myan.Sangerrelease5/filter/Sanger.release5.core.SNP.vcf.gz \
--sample_file SangerSEasia.list -nt 10 \
-o Sanger.release5.SEasia.vcf ### 678 samples

vcftools --vcf Sanger.release5.SEasia.vcf --freq --out Sanger.release5.SEasia.alleleFrequency

