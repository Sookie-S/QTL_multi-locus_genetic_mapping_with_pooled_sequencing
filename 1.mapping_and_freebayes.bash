#!/bin/bash -l
#$ -cwd

# Mapping and variant calling method from https://www.pnas.org/doi/full/10.1073/pnas.2101242118

#Requires:

#trimmomatic 0.36
#bwa 0.7.15 gcc-4.8.5
#samtools 1.4
#picard 2.1.0
#bcftools 1.11
#gatk 3.8.0
#freebayes 1.1.0

# Make a directory 'fastqs' containing the R1 and R2 fastq files per sample
# Each sample was processed individually such that there was one resultant VCF file per sample
 
for R1 in fastqs/*_R1_*
  do
    R1=$(echo ${R1} | sed 's,^[^/]*/,,')
    R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz} 
    R1paired=${R1//.fastq.gz/_paired.fastq}
    R1unpaired=${R1//.fastq.gz/_unpaired.fastq} 
    R2paired=${R2//.fastq.gz/_paired.fastq}
    R2unpaired=${R2//.fastq.gz/_unpaired.fastq}
    R1bam=${R1//R1_001.fastq.gz/R1.bam}
    
   
    ## Trimming
       
    echo ${R1//R1_001.fastq.gz}
    trimmomatic PE -phred33 -threads 4 \
    fastqs/$R1 \
    fastqs/$R2 \
    $R1paired $R1unpaired \
    $R2paired $R2unpaired \
    ILLUMINACLIP:/allprimers.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:35
    
	
    ref='/sacCer_YPS128_plus_sacKud_IFO1802T_US/fasta/saccerYPS128_SkudIFO1802TUS'  #path to reference

    ## align with bwa mem
    
    bwa mem -t 4 -M /mnt/data-sets/bcf/genomeIndexes/sacCer_YPS128_plus_sacKud_IFO1802T_US/bwa/saccerYPS128_SkudIFO1802TUS $R1paired $R2paired  | samtools view -@ 4 -f2 -F256 -bS - > $R1bam 
     
    ## BAM flagstat with samtools
    samtools flagstat $R1bam > $R1bam.raw.flagstat
    
   ## Fixmates with samtools
    samtools fixmate -@4 $R1bam $R1bam.fm.bam
 
 
    ## Add read groups with picard
    picard AddOrReplaceReadGroups I=$R1bam.fm.bam O=$R1bam.fm.RG.bam RGID=${R1//R1_001.fastq.gz} RGLB=${R1//R1_001.fastq.gz} RGPL=Illumina RGPU=unit1 RGSM=${R1//R1_001.fastq.gz}
    
    
    ## Sort by coordinates with samtools
    samtools sort -@ 4 $R1bam.fm.RG.bam > $R1bam.fm.RG.sorted.bam
    samtools index $R1bam.fm.RG.sorted.bam

     
    ##RealignerTargetCreator
    gatk -I $R1bam.fm.RG.sorted.bam -T RealignerTargetCreator -R $ref.fa -o $R1bam.fm.RG.sorted.bam.targets.intervals

    ##IndelRealigner:re-aligns target regions
    gatk -I $R1bam.fm.RG.sorted.bam -R $ref.fa -T IndelRealigner -targetIntervals $R1bam.fm.RG.sorted.bam.targets.intervals -o $R1bam.fm.RG.sorted.lr.bam

     
    ##Duplicate removal with picard
    picard MarkDuplicates I=$R1bam.fm.RG.sorted.lr.bam O=$R1bam.fm.RG.sorted.lr.deduped.bam M=${R1//R1_001.fastq.gz}_marked_dup_metrics.txt

    ##sort and index
    samtools sort $R1bam.fm.RG.sorted.lr.deduped.bam -o $R1bam.fm.RG.sorted.lr.deduped.sorted.bam
    samtools index $R1bam.fm.RG.sorted.lr.deduped.sorted.bam


    ##samtools flagstat
    samtools flagstat $R1bam.fm.RG.sorted.lr.deduped.sorted.bam > $R1bam.fm.RG.sorted.lr.deduped.sorted.flagstat.txt

    ##coverage
    gatk -T DepthOfCoverage -R $ref.fa -I $R1bam.fm.RG.sorted.lr.deduped.sorted.bam > $R1bam.fm.RG.sorted.lr.deduped.sorted.bam.depth.txt


   ## VARIANT CALLING   
   freebayes -k --ploidy 1 --min-mapping-quality 30 --min-base-quality 20 --min-alternate-count 2 --min-coverage 2 --min-alternate-fraction 0.1 -X --no-mnps -f $ref.fa $R1bam.fm.RG.sorted.lr.deduped.sorted.bam > freebayes_ploidy1.vcf
  
done

   