#!/bin/bash

export working_dir=$PWD
export cores=$4
export host=$3
export contigs=$2
export sample_id=$1
export software_location=$(dirname $0)


echo
echo "Input:"
echo

echo number_of_cores:$cores
echo sample_id:$sample_id
echo Assembly contigs:$contigs


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

#is_command_installed diamond
is_command_installed blastn

echo

if [ -r "$contigs" ] && [ -n "$sample_id" ] && [ "$cores" -eq "$cores" ]

then

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Search proteins reads-------------------------------------------------------------------------------------------------------------------------

mkdir -p "$working_dir"/results/"$sample_id"/virus_search/blastx

"$software_location"/software/./diamond blastx --sensitive --max-target-seqs 250 -p "$NSLOTS" -q "$contigs" -o "$working_dir"/results/"$sample_id"/virus_search/blastx/"$sample_id".results.tab -t "$working_dir"/results/"$sample_id"/virus_search/blastx/ -d "$software_location"/databases/uniprot_viruses.fasta.dmnd

python "$software_location"/software/search_virus_domains.py "$working_dir"/results/"$sample_id"/virus_search/blastx/"$sample_id".results.tab "$contigs" "$sample_id" > "$working_dir"/results/"$sample_id"/virus_search/blastx/"$sample_id"_virus_candidates.fasta

#Search compare to nt ref database-------------------------------------------------------------------------------------------------------------------------

mkdir -p "$working_dir"/results/"$sample_id"/virus_search/blastn/

blastn -db "$software_location"/databases/nt -max_hsps 1 -max_target_seqs 10 -evalue 1e-01 -num_threads "$cores" -query "$working_dir"/results/"$sample_id"/virus_search/blastx/"$sample_id"_virus_candidates.fasta -out "$working_dir"/results/"$sample_id"/virus_search/blastn/"$sample_id".results.csv -outfmt "6 qseqid sseqid stitle qlen slen length pident nident mismatch gaps evalue bitscore"

echo "Sample	Scaffold ID	Scaffold length	read per bp (read_depth)	fraction with homology to virus protein	Homologous virus proteins	best BlastN hit	BlastN homology to best hit	Ref genome was homolog	sequence" > "$working_dir"/results/"$sample_id"/virus_search/blastn/"$sample_id".final_candidates.tab

python "$software_location"/software/exclude_host.py "$working_dir"/results/"$sample_id"/virus_search/blastn/"$sample_id".results.csv "$working_dir"/results/"$sample_id"/virus_search/blastx/"$sample_id"_virus_candidates.fasta "$host" >> "$working_dir"/results/"$sample_id"/virus_search/blastn/"$sample_id".final_candidates.tab

#actual analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

else

echo " "
echo "ERROR: Incorrect input!"
echo "Novel virus detection pipeline version 0.1 by Daniel Wüthrich (danielwue@hotmail.com)"
echo " "
echo "Usage: "
echo "  sh novel_virus_identification.sh <Sample_ID> <metagenomics_contigs> <Number_of_cores>"
echo " "
echo "  <Sample_ID>               Unique identifier for the sample"
echo "  <metagenomics_contigs>    Contigs of metagenomics assembly"
echo "  <Number_of_cores>         number of parallel threads to run (int)"
echo " "

if ! [ -n "$sample_id" ];then
echo Incorrect input: "$sample_id"
fi
if ! [ -r "$contigs" ];then
echo File not found: "$contigs"
fi

if ! [ "$cores" -eq "$cores" ] ;then
echo Incorrect input: "$cores"
fi


fi

