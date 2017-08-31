mkdir databases
cd databases

echo
echo Dowloading databases ...
echo

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.dat.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.dat.gz

echo
echo Unpacking databases ...
echo

for i in nt*.tar.gz
do
tar xvzf "$i"
done


echo
echo Preparing databases ...
echo

python ../software/convert_uniprot_to_fasta.py uniprot_trembl_viruses.dat.gz > uniprot_viruses.fasta
python ../software/convert_uniprot_to_fasta.py uniprot_sprot_viruses.dat.gz >> uniprot_viruses.fasta

../software/./diamond makedb --in uniprot_viruses.fasta --db uniprot_viruses.fasta
