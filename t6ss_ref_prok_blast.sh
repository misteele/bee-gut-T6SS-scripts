t6ss_ref_prok_blast.sh - search genomes in the NCBI ref_prok_rep_genomes database for sequences to use as an outgroup when constructing T6SS phylogenies
####################
# Make BLAST databases for all proteins encoded by each genome of interest - one database per genome
for i in $SEQ_DIR/*.faa; do
        id=$(echo $i | sed 's/.*\///g' | sed 's/.faa//g')
        echo "Name for new blast directory is ${id}_prot_db"
        makeblastdb -in ${i} -out ${id}_prot_db -dbtype prot -parse_seqids
done


# Remove linebreaks in the reference FASTA files
more $SEQS | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' > $SEQS2

# Double check that files turned out as expected before looping through them
echo "Protein IDs"
head $LIST
echo " "
echo "Protein Seqs"
head $SEQS2

# Loop through proteins from reference T6SS locus and use tblastn to search for homologous sequences in the NCBI ref_prok_rep_genomes database
while read -r line; do
        echo "Reading line $line"

        # Extract the sequence ID and sequence for one of the protein sequences in the FASTA file (create a FASTA file with only one sequence)
        grep -A1 "$line" $SEQS2 > $DIR/${line}.fa
#		echo "$DIR/${line}.fa"
#		head $DIR/${line}.fa

        # Blast extracted protein sequence against the NCBI ref_prok_rep_genomes database (downloaded from NCBI FTP site)
         echo "Blast analysis for $line started at `date`"

        tblastn -db $REFDIR/ref_prok_rep_genomes -query $DIR/${line}.fa -outfmt "6 sseqid qcovs pident bitscore stitle" \
        -out ${query}_${line}_ref_prok_tblastn.txt -num_threads 16 -max_target_seqs 500 -evalue 0.0001;
done<$LIST

# Loop through the tblastn output files and filter out hits with a minimum 40% percent identity over 70% of the query length and take the 25 best hits that meet those criteria
for x in *_ref_prok_tblastn.txt; do
        id=$(echo $x | sed 's/_ref_prok_tblastn.txt//g' | sed 's/.*\///g')
        echo "id is $id"
        more $x | awk '$2>=70&&$3>=40{print $0}' | sort -k4,4nr | sort -u -k1,1 | sort -k4,4nr | awk 'FNR <= 25 {print $1}' > \
			${id}_ref_prok_tblastn.70cov40pident.top25.txt
done

# Identify genomes that have at least three decent hits to reference T6SS proteins - setting the bar low
cat *_ref_prok_tblastn.70cov40pident.top25.txt | sed 's/.*.ref.//g' | sed 's/\.1./.1/g' | sort | uniq -c | awk '$1>=3{print $2}' > \
${query}_possible_relatives.2.txt

mkdir $OUTDIR
mv *_ref_prok_tblastn.txt $OUTDIR
mv *_ref_prok_tblastn.70cov40pident.top25.txt $OUTDIR
