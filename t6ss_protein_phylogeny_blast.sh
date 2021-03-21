t6ss_protein_phylogeny_blast.sh - find homologs to T6SS proteins and prepare alignments to be concatenated to make phylogenies
########################
OUTGROUP_DIR=/Protein_Files/possible_relatives_ncbi-genomes-2020-01-28
SEQ_DIR=/Protein_Files/Bee_Gut_Prot_Seqs
DB=all_bee_out_prot_db											## Local BLAST protein database containing all proteins from genomes of interest and outgroup genomes
T6SS_PROT_DIR=/Protein_Files/T6SS_Ref_Prot_Sequences
T6SS_ALIGN_DIR=/Protein_Files/T6SS_Prot_Alignment_Files4

mkdir $T6SS_ALIGN_DIR

### Make a BLAST database that contains all of the proteins encoded by Gram-negative bee gut bacteria
### IMPORTANT: Download protein sequences from NCBI genome assembly database, NOT the protein database
### IMPORTANT: Append genome assembly accession number (or other strain identifier) to each sequence. Makeblastdb won't make a database that includes sequences with the same identifier.

### Make the BLAST database
makeblastdb -in all_bee_out_assembly_proteins.0213.faa -out $DB -dbtype prot -parse_seqids


for query in Bimp_T6SS_1 PEB0191_T6SS_1 PEB0191_T6SS_2 wkB180_T6SS_1; do	### Loop through a set of reference T6SS loci
        echo "Making T6SS protein alignments for T6SS gene cluster $query"
        LIST=Representative_T6SS_Protein_IDs/${query}_prot_IDs.2.txt		### A list of all the IDs for proteins encoded within a reference T6SS locus
        SEQS=Representative_T6SS_Protein_Sequences/${query}_prot_seqs.2.fa	### A fasta file with sequences of proteins from the reference T6SS

        mkdir $T6SS_ALIGN_DIR/${query}_Alignments

		### This loops through each of the proteins encoded within the reference T6SS locus
        while read line; do		
             echo "Reading line $line"

                # Use grep to extract FASTA sequence for protein in wkB2
                grep -A1 "$line" $SEQS > $T6SS_PROT_DIR/${line}.fa
                echo "Query file: $T6SS_PROT_DIR/${line}.fa Sequence:" >&1
                head $T6SS_PROT_DIR/${line}.fa >&1

                # Blast protein from reference T6SS locus against BLAST protein database
                echo "Blast analysis for $line started at `date`" >&1

                blastp -db $DB -query $T6SS_PROT_DIR/${line}.fa \
                -outfmt "6 sseqid qseqid qcovs pident bitscore" -out $T6SS_ALIGN_DIR/${line}_${query}_blastp.txt \
                -num_threads 16 -max_target_seqs 500
                echo "blastp output:" >&1
                head $T6SS_ALIGN_DIR/${line}_${query}_blastp.txt >&1

                ### Filter out any BLAST hits with less than 20% identity over 20% of length. I had to make this very low to identify related T6SSs in non-bee gut bacteria.
                ### Put a comma and tab between the genome assembly accession number and protein accession number to make it possible to filter out the best matching protein from each genome (in case bacteria have more than one T6SS)

                echo "Filtering blastp output:" >&1
                echo "Split genome accession and protein accession:" >&1 
                awk '$3>=20&&$4>=20{print $1","$5}' $T6SS_ALIGN_DIR/${line}_${query}_blastp.txt | sed 's/^\(GCF_[0-9]*\.[0-9]\)_/\1,\t/g' | \
                sed 's/OCG/,\tOCG/g' | sort -t, -k1,1 -k3,3nr | head >&1 # print to output to make sure sorting worked as expected
                echo " " >&1
                
                ### Next, filter BLAST hits and split qseqids as above, then filter out the best hit per genome (based on bitscore)
                ### The protein accession number is then rejoined to the genome assembly accession number
                ### The sseqid identifiers have to be returned to their original state to use blastdbcmd in the next step
                
                echo "Fuse genome accession and protein accession after filtering best match for genome:" >&1
                awk '$3>=20&&$4>=20{print $1","$5}' $T6SS_ALIGN_DIR/${line}_${query}_blastp.txt | sed 's/^\(GCF_[0-9]*\.[0-9]\)_/\1,\t/g' | \
                sed 's/OCG/,\tOCG/g' | sort -t, -k1,1 -k3,3nr | sort -t, -u -k1,1 | sed 's/,\tOCG/OCG/g' | sed 's/,\t/_/g' | \
                awk -F, '{print $1}' > $T6SS_ALIGN_DIR/${line}_${query}_blastp.hits.txt

                echo "Filtered blastp hits:" >&1
                head $T6SS_ALIGN_DIR/${line}_${query}_blastp.hits.txt >&1 # checking to make sure that the filtered hits file looks okay

				### Use blastdbcmd to extract the protein sequences for the top blast hits
				blastdbcmd -db $DB -entry_batch $T6SS_ALIGN_DIR/${line}_${query}_blastp.hits.txt \
                        -out $T6SS_ALIGN_DIR/${query}_${line}_seqs.fa
                        echo "blastdbcmd output:" >&1
                        head $T6SS_ALIGN_DIR/${query}_${line}_seqs.fa >&1
 				
 				### Use muscle to align proteins, then remove the protein accession numbers from the alignment (keeping only the genome accession number, which will be used to concatenate the alignments of all proteins encoded by genes within the reference T6SS locus (i.e., a single T6SS)
                /work/02883/ms63587/stampede2/OrthoApps/muscle3.8.31_i86linux64 \
                -in $T6SS_ALIGN_DIR/${query}_${line}_seqs.fa \
                -out $T6SS_ALIGN_DIR/${query}_Alignments/${query}_${line}_seqs.aligned.fa
                sed 's/>\(GCF_[0-9]*\.[0-9]\).*$/>\1/g' $T6SS_ALIGN_DIR/${query}_Alignments/${query}_${line}_seqs.aligned.fa > \
                $T6SS_ALIGN_DIR/${query}_Alignments/${query}_${line}_seqs.aligned.renamed.fa
                
                ### Remove unneeded files
                rm $T6SS_ALIGN_DIR/${line}_${query}_blastp.txt $T6SS_ALIGN_DIR/${line}_${query}_blastp.hits.txt \
                $T6SS_ALIGN_DIR/${query}_${line}_seqs.fa
                rm $T6SS_ALIGN_DIR/${query}_Alignments/${query}_${line}_seqs.aligned.fa
        done<$LIST
done