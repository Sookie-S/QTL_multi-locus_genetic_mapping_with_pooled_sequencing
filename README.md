# QTL_multi-locus_genetic_mapping_with_pooled_sequencing

Naseeb S, Visinoni F, Hu Y, Hinks Roberts AJ, Maslowska A, Walsh T, Smart KA, Louis EJ, Delneri D. Restoring fertility in yeast hybrids: Breeding and quantitative genetics of beneficial traits. Proc Natl Acad Sci U S A. 2021 Sep 21;118(38):e2101242118. doi: 10.1073/pnas.2101242118. PMID: 34518218; PMCID: PMC8463882.

The 'MultipoolQTLanalysis' directory contains various wrapper scripts from Yue Hu that are required to run the analysis (as well as Multipool itself). Navigate into the 'MultipoolQTLanalysis' directory.


## 1. Mapping and variant calling with freebayes

A bash script was run over each set of fastq files (R1 and R2). The samples were processed individually such that there is one resultant VCF file per sample. For ech VCF file the following columns were extracted and saved as csv files: Chrom, Pos, Ref, Alt, AO and RO (NB it is important the columns are in this order). If the sample was a founder (only one species) and mapped against a hybrid genome, remove irrelevent genome positions.


## 2. Process alleles in R using R scripts

The parental allele frequencies for bi-allelic SNVs were obtained for each of the segregant pools

Change the directory to 'Example_LowTemp'

R command:
```
setwd('Example_LowTemp')
```

The scripts require making the following sub directories:

Example_LowTemp/Input/Founder (Put the four founder csv files in here)

Example_LowTemp/Input/Pool    (Put the two pool csv files in here)

make empty directories to store results:

Example_LowTemp/Output_Scer/SN1_S68      (This will contain the output for s_cer alleles high fitness pool)

Example_LowTemp/Output_Scer/SN2_S69      (This will contain the output for s_cer alleles low fitness pool)

Example_LowTemp/Output_Skud/SN1_S68      (This will contain the output for s_kud alleles high fitness pool)

Example_LowTemp/Output_Skud/SN2_S69      (This will contain the output for s_kud alleles high fitness pool)

## 3. Run MULTIPOOL Analysis

This requires python 2.7 environment 

The MultipoolQTLanalysis directory contains the scripts required to run the analysis (as well as MULTIPOOL itself).
To run, navigate into the 'MultipoolQTLanalysis' directory. 

Bash command lines:
```
mkdir /LowTemp
cd LowTemp
mkdir MPresults_scer
mkdir MPresults_skud
```
You will also need to copy the multipool input files, that were made using the R scripts earlier (the R output), into this directory

Bash command lines:
```
cp -r Example_LowTemp/Output_Scer/ . 
cp -r Example_LowTemp/Output_Skud/ .
```

From inside the 'LowTemp' directory, proceed to run the analysis for each condition on Bash command line as follows:

Bash command lines:
```
for num in chrI    chrII   chrIII  chrIV   chrV    chrVI   chrVII  chrVIII chrIX   chrX    chrXI   chrXII  chrXIII chrXIV  chrXV  chrXVI;
  do
    .././mp_inference.py -n 150  Output_scer/SN1_S68/$num.txt Output_scer/SN2_S69/$num.txt -c 3300 -r 100 -m contrast -o MPresults_scer/scerLT$num.out --plotFile MPresults_scer/scerLT_$num |& tee MPresults_scer/scerLT_$num.log;
  done

for num in Skud_10  Skud_12  Skud_14  Skud_16  Skud_2  Skud_4  Skud_6  Skud_8 Skud_11  Skud_13  Skud_15  Skud_1   Skud_3  Skud_5  Skud_7  Skud_9;
  do
    .././mp_inference.py -n 150  Output_skud/SN1_S68/$num.txt Output_skud/SN2_S69/$num.txt -c 3300 -r 100 -m contrast -o MPresults_skud/skudLT$num.out --plotFile MPresults_skud/skudLT_$num |& tee MPresults_skud/skudLT_$num.log
  done
```


