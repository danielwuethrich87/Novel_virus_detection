Novel virus detection
=======================

This pipeline is searching for novel viruses based on protein alignment.<br />

scripts:<br />
step 0:get_dbs.sh: installation of databases<br />
step 1.1:mapping_DNA.sh: removal of host reads from DNA samples<br />
step 1.2:mapping_RNA.sh: removal of host reads from RNA samples<br />
step 2:metagenomic_assembly.sh: assembly of reads from metagnomics sample<br />
step 3:novel_virus_identification.sh: searching for Virus domains in scaffolds<br />
		-searching for Virus protein<br />
		-comparing to nt database, to excluded false positive<br />

#Requirements:

-Linux 64 bit system<br />

-python (version 2.7)<br />
-samtools (version 1.3)<br />
-bowtie2 (version 2.3.0)<br />
-STAR (version 2.3.0)<br />
-SPAdes (version 3.10.1)<br />

#Installation:

wget https://github.com/danielwuethrich87/Reference_based_virus_detection/archive/master.zip<br />
unzip master.zip<br />
cd Reference_based_virus_detection-master<br />
sh get_dbs.sh<br />

#Usage:
<br />
sh mapping_DNA.sh <Sample_ID> <Reads_R1> <Reads_R2> <host_bowtie2_index> <Number_of_cores> <br />
<br />
<Sample_ID>               Unique identifier for the sample<br />
<Reads_R1>                Foreward read file<br />
<Reads_R2>                Reversed read file<br />
<host_bowtie2_index>      Bowtie2 index of host genome<br />
<Number_of_cores>         number of parallel threads to run (int)<br />
<br />
sh mapping_RNA.sh <Sample_ID> <Reads_R1> <Reads_R2> <host_STAR_index> <Number_of_cores> <br />
<br />
<Sample_ID>               Unique identifier for the sample<br />
<Reads_R1>                Foreward read file<br />
<Reads_R2>                Reversed read file<br />
<host_STAR_index>         Bowtie2 index folder of host genome<br />
<Number_of_cores>         number of parallel threads to run (int)<br />
<br />
sh metagenomic_assembly.sh <Sample_ID> <Reads_R1> <Reads_R2> <Number_of_cores> <br />
<br />
<Sample_ID>               Unique identifier for the sample<br />
<Reads_R1>                Foreward read file<br />
<Reads_R2>                Reversed read file<br />
<Number_of_cores>         number of parallel threads to run (int)<br />
<br />
sh novel_virus_identification.sh <Sample_ID> <metagenomics_contigs> <Number_of_cores> <br />
<br />
<Sample_ID>               Unique identifier for the sample<br />
<metagenomics_contigs>    Contigs of metagenomics assembly<br />
<Number_of_cores>         number of parallel threads to run (int)<br />
<br />
#example mapping:

#!/bin/sh<br />
#$ -q all.q<br />
#$ -e $JOB_ID.cov.err<br />
#$ -o $JOB_ID.cov.out<br />
#$ -cwd <br />
#$ -pe smp 16<br />

<br />
module add UHTS/Aligner/bowtie2/2.3.0;<br />

for i in 23871a<br />

do<br />

# Tipp: create index: bowtie2-build reference.fa reference.fa<br />

sh /home/dwuethrich/Application/novel_virus_detection/mapping_DNA.sh "$i" ../reads_20150103/Project_Neurocenter_TS/Sample_"$i"/"$i"_*_R1_*.fastq.gz ../reads_20150103/Project_Neurocenter_TS/Sample_"$i"/"$i"_*_R2_*.fastq.gz ../mapping/reference_genome/Bos_taurus.UMD3.1.dna_sm.toplevel "$NSLOTS"<br />

done<br />



module add UHTS/Aligner/STAR/2.3.0;<br />

for i in 23871a<br />

do<br />

# Tipp: create index: STAR --runMode genomeGenerate --genomeDir "$working_dir"/star_index --genomeFastaFiles "$working_dir"/reference.fa --runThreadN 8<br />

sh /home/dwuethrich/Application/novel_virus_detection/mapping_RNA.sh "$i" ../reads_20150103/Project_Neurocenter_TS/Sample_"$i"/"$i"_*_R1_*.fastq.gz ../reads_20150103/Project_Neurocenter_TS/Sample_"$i"/"$i"_*_R2_*.fastq.gz ../mapping_RNA/alignment/genome/start_index/ "$NSLOTS"<br />

done<br />

#example assembly:

#!/bin/sh<br />
#$ -q all.q<br />
#$ -e $JOB_ID.assembly.err<br />
#$ -o $JOB_ID.assembly.out<br />
#$ -cwd<br />
#$ -pe smp 16<br />

module add UHTS/Assembler/SPAdes/3.10.1;<br />
module add UHTS/Analysis/samtools/1.3;<br />
module add UHTS/Aligner/bowtie2/2.3.0;<br />

for i in 23871a<br />

do<br />

sh /home/dwuethrich/Application/novel_virus_detection/metagenomic_assembly.sh "$i" Test_reads_R1.gz Test_reads_R2.gz "$NSLOTS"<br />

done<br />

#example virus search:

#!/bin/sh<br />
#$ -q all.q<br />
#$ -e $JOB_ID.assembly.err<br />
#$ -o $JOB_ID.assembly.out<br />
#$ -cwd #executes from the current directory and safes the ouputfiles there<br />
#$ -pe smp 16<br />


module add Blast/ncbi-blast/2.6.0+;<br />

for i in 23871a<br />

do<br />

sh /home/dwuethrich/Application/novel_virus_detection/novel_virus_identification.sh 23871a /data/projects/p187_cattle_and_sheep_encephalitis/cattle/28_LT_29_LT_high_depth/assembly/assemblies/cov_selection/28_LT/28_LT.fasta Bos_taurus "$NSLOTS"<br />

done<br />











