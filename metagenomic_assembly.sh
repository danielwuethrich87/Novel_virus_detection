#!/bin/bash

export working_dir=$PWD
export cores=$4
export reads_R2=$3
export reads_R1=$2
export sample_id=$1
export software_location=$(dirname $0)


echo
echo "Input:"
echo

echo number_of_cores:$cores
echo sample_id:$sample_id
echo read_file_R1:$reads_R1
echo read_file_R2:$reads_R2
echo genus:$genus
echo species:$species


echo
echo "Checking software ..."
echo

is_command_installed () {
if which $1 &>/dev/null; then
    echo "$1 is installed in:" $(which $1)
else
    echo
    echo "ERROR: $1 not found."
    echo
    exit
fi
}



is_command_installed python
is_command_installed spades.py
is_command_installed samtools
is_command_installed bowtie2
echo

if [ -r "$reads_R1" ] && [ -r "$reads_R2" ] && [ -n "$sample_id" ] && [ "$cores" -eq "$cores" ]

then

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#trimmomatic-------------------------------------------------------------------------------------------------------------------------


mkdir -p "$working_dir"/results/"$sample_id"/assembly/0_read_trimming

java -jar "$software_location"/software/trimmomatic-0.36.jar PE -threads "$cores" -phred33 "$reads_R1" "$reads_R2" "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r1.fastq.gz "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r1.not-paired.fastq.gz "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r2.fastq.gz "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r2.not-paired.fastq.gz SLIDINGWINDOW:4:15 MINLEN:127 2> "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/"$sample_id".read_trimm_info


#SPades-------------------------------------------------------------------------------------------------------------------------

mkdir -p "$working_dir"/results/"$sample_id"/assembly/1_spades_assembly

spades.py --meta -t "$cores" -k 21,33,55,77,99,127 -1 "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r1.fastq.gz -2 "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r2.fastq.gz -o "$working_dir"/results/"$sample_id"/assembly/1_spades_assembly

#coverage_calculation-------------------------------------------------------------------------------------------------------------------------


mkdir -p "$working_dir"/results/"$sample_id"/assembly/2_cov_selection

cp "$working_dir"/results/"$sample_id"/assembly/1_spades_assembly/scaffolds.fasta "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.fa

bowtie2-build "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.fa "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds
bowtie2 -p "$cores" --un-conc-gz "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/not_aligned_reads.fastq.gz -x "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds -1 "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r1.fastq.gz -2 "$working_dir"/results/"$sample_id"/assembly/0_read_trimming/r2.fastq.gz -S "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sam 2> "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/mapping_Info."$sample_id"

samtools sort -@ "$cores" -T "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/temp_sort -o "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sorted.bam "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sam
samtools rmdup "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sorted.bam "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sorted.removed_duplicates.bam
samtools index "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sorted.removed_duplicates.bam
samtools faidx "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.fa

samtools idxstats "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sorted.removed_duplicates.bam > "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.idxstats

python "$software_location"/software/filter_contigs_by_samtools_idxstats.py  "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.fa "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.idxstats 0.0 "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/Short_scaffolds_"$sample_id".fasta > "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/"$sample_id".fasta
#the 0.1 indicates the min coverage a scaffold must have, compared to large scaffolds

#clean-up-------------------------------------------------------------------------------------------------------------------------

rm "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sam "$working_dir"/results/"$sample_id"/assembly/2_cov_selection/scaffolds.sorted.bam

echo "$sample_id" finished `date`

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

else

echo " "
echo "ERROR: Incorrect input!"
echo "Metagenomic assembly pipe line version 0.1 by Daniel Wüthrich (danielwue@hotmail.com)"
echo " "
echo "Usage: "
echo "  sh metagenomic_assembly.sh <Sample_ID> <Reads_R1> <Reads_R2> <Number_of_cores>"
echo " "
echo "  <Sample_ID>               Unique identifier for the sample"
echo "  <Reads_R1>                Foreward read file"
echo "  <Reads_R2>                Reversed read file"
echo "  <Number_of_cores>         number of parallel threads to run (int)"
echo " "

if ! [ -n "$sample_id" ];then
echo Incorrect input: "$sample_id"
fi
if ! [ -r "$reads_R1" ];then
echo File not found: "$reads_R1"
fi
if ! [ -r "$reads_R2" ];then
echo File not found: "$reads_R2"
fi
if ! [ "$cores" -eq "$cores" ] ;then
echo Incorrect input: "$cores"
fi


fi

