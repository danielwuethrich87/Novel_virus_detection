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
echo reference bowtie2 index:$reference


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


is_command_installed bowtie2
echo

if [ -r "$reads_R1" ] && [ -r "$reads_R2" ] && [ -n "$sample_id" ] && [ "$cores" -eq "$cores" ]

then

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Mapping reads-------------------------------------------------------------------------------------------------------------------------

mkdir -p "$working_dir"/results/"$sample_id"/mapping_DNA

bowtie2 -p "$cores" --un-conc-gz "$working_dir"/results/"$sample_id"/mapping_DNA/not_aligned_reads.fastq.gz -x "$reference" -1 "$reads_R1" -2 "$reads_R2" -S "$working_dir"/results/"$sample_id"/mapping_DNA/"$sample_id"_alignment.sam 2> "$working_dir"/results/"$sample_id"/mapping_DNA/"$sample_id"_mapping_Info.txt

rm "$working_dir"/results/"$sample_id"/mapping_DNA/"$sample_id"_alignment.sam

mv "$working_dir"/results/"$sample_id"/mapping_DNA/not_aligned_reads.fastq.1.gz "$working_dir"/results/"$sample_id"/mapping_DNA/"$sample_id"_DNA_not_aligned_R1.fastq.gz
mv "$working_dir"/results/"$sample_id"/mapping_DNA/not_aligned_reads.fastq.2.gz "$working_dir"/results/"$sample_id"/mapping_DNA/"$sample_id"_DNA_not_aligned_R2.fastq.gz

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

else

echo " "
echo "ERROR: Incorrect input!"
echo "DNA mapping pipeline version 0.1 by Daniel Wüthrich (danielwue@hotmail.com)"
echo " "
echo "Usage: "
echo "  sh mapping_DNA.sh <Sample_ID> <Reads_R1> <Reads_R2> <host_bowtie2_index> <Number_of_cores>"
echo " "
echo "  <Sample_ID>               Unique identifier for the sample"
echo "  <Reads_R1>                Foreward read file"
echo "  <Reads_R2>                Reversed read file"
echo "  <host_bowtie2_index>      Bowtie2 index of host genome"
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

