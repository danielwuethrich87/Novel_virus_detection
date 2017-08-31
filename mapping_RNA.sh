#!/bin/bash

export working_dir=$PWD
export cores=$5
export reference=$4
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
echo reference STAR index folder:$reference


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


is_command_installed STAR
echo

if [ -r "$reads_R1" ] && [ -r "$reads_R2" ] && [ -n "$sample_id" ] && [ "$cores" -eq "$cores" ]

then

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Mapping reads-------------------------------------------------------------------------------------------------------------------------

mkdir -p "$working_dir"/results/"$sample_id"/mapping_RNA

STAR --genomeDir "$reference" --readFilesCommand zcat --readFilesIn "$reads_R1" "$reads_R2" --runThreadN "$cores" --outReadsUnmapped Fastx --outFileNamePrefix "$working_dir"/results/"$sample_id"/mapping_RNA/

mv "$working_dir"/results/"$sample_id"/mapping_RNA/Unmapped.out.mate1 "$working_dir"/results/"$sample_id"/mapping_RNA/"$sample_id"_RNA_not_aligned_R1.fastq
mv "$working_dir"/results/"$sample_id"/mapping_RNA/Unmapped.out.mate2 "$working_dir"/results/"$sample_id"/mapping_RNA/"$sample_id"_RNA_not_aligned_R2.fastq

gzip "$working_dir"/results/"$sample_id"/mapping_RNA/"$sample_id"_RNA_not_aligned_R1.fastq
gzip "$working_dir"/results/"$sample_id"/mapping_RNA/"$sample_id"_RNA_not_aligned_R2.fastq

rm "$working_dir"/results/"$sample_id"/mapping_RNA/Aligned.out.sam 

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

else

echo " "
echo "ERROR: Incorrect input!"
echo "RNA mapping pipeline version 0.1 by Daniel Wüthrich (danielwue@hotmail.com)"
echo " "
echo "Usage: "
echo "  sh mapping_RNA.sh <Sample_ID> <Reads_R1> <Reads_R2> <host_STAR_index> <Number_of_cores>"
echo " "
echo "  <Sample_ID>               Unique identifier for the sample"
echo "  <Reads_R1>                Foreward read file"
echo "  <Reads_R2>                Reversed read file"
echo "  <host_STAR_index>         Bowtie2 index folder of host genome"
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

