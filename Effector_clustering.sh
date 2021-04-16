rhs_clustering.NEW.sh.save - most recent version of the script that I used to do most of my work with Rhs sequences
#####################

### Cluster Rhs toxin domain, Rhs core domain, and VgrG protein sequences
CDHITDIR= # path to directory where cd-hit is installed
vers= # version ID, usually the date of the analysis
KEY= # CSV file linking genome assembly accession numbers to genome metadata
core= # FASTA file containing Rhs core domain amino acid sequences for clustering
tox= # FASTA file containing Rhs C-terminal toxin domain amino acid sequences for clustering
vgrg= # FASTA file containing VgrG amino acid sequences for clustering

# Use CD-HIT to cluster protein sequences with a 40% identity cutoff
$CDHITDIR/cd-hit -i $core -o rhs_core_${vers}_40perc.fa -c 0.40 -n 2 -d 0 -T 16
$CDHITDIR/cd-hit -i $tox -o rhs_tox_${vers}_40perc.fa -c 0.40 -n 2 -d 0 -T 16
$CDHITDIR/cd-hit -i $vgrg -o vgrg_${vers}_40perc.fa -c 0.40 -n 2 -d 0 -T 16

# Assign numbers to representative sequences from each cluster and rename sequences
more rhs_core_${vers}_40perc.fa | awk '/^>/{print ">lcl|Rhs_core_" ++i; next}{print}' | tail -n +4 > rhs_core_${vers}_40perc_clusters.fa
more rhs_tox_${vers}_40perc.fa | awk '/^>/{print ">lcl|Rhs_tox_" ++i; next}{print}' | tail -n +4 > rhs_tox_${vers}_40perc_clusters.fa
more vgrg_${vers}_40perc.fa | awk '/^>/{print ">lcl|VgrG_" ++i; next}{print}' | tail -n +4 > vgrg_${vers}_40perc_clusters.fa

### Remove Rhs C-terminus domains that share homology with core domains
DB=rhs_core_prot_db
REF=rhs_core_${vers}_40perc_clusters.fa
QUERY=rhs_tox_${vers}_40perc_clusters.fa
search=rhs_tox_core_overlap_second_check_${vers}

# Make a BLAST database containing representative core domains
makeblastdb -in $REF -out $DB -dbtype prot -parse_seqids
# Use protein BLAST to search for homology between representative toxin domains and representative core domains (using the same length and %identity cutoffs that will be used to assign proteins to clusters)
blastp -db $DB4 -query $QUERY -outfmt "6 sseqid qseqid qcovs pident length qseq" -out ${search4}_blastp.txt -num_threads 16
echo "${search4}_blastp.txt >40 %identity >90 aa alignment:"
more ${search4}_blastp.txt | awk '$4>=40&&$5>=90{print $0}' | sort -u

# Manually remove detected sequences 

#######################
name=all_bee_gut
DB=${name}_prot_${vers}_db # BLAST protein database containing all proteins from all genome assemblies to be analyzed
query=rhs_tox # or "rhs_core" or "vgrg"
QUERY=${query}_${vers}_40perc_clusters.fa
var=40pi90len # to change the percent identity and alignment length cutoffs, change the awk command below (awk '$4>=40&&$5>=90{print $0}'). $4 is the percent identity, $5 is the alignment length.
MUSCLEDIR= # Path to directory where MUSCLE is installed

# Use protein BLAST to identify proteins from genome assemblies that are homologous to representative Rhs toxin domains
blastp -db $DB -query $QUERY -outfmt "6 qseqid sseqid qcovs pident length bitscore" \
        -out ${name}_${query}_blastp_${vers}.txt -num_threads 16 -max_target_seqs 500

# Filter BLAST output and insert a tab between assembly and protein accession numbers in the sequence IDs 
more ${name}_${query}_blastp_${vers}.txt | sort -k2,2 -k6,6nr | sort -u -k2,2 | awk '$4>=40&&$5>=90{print $0}' | \
        sed 's/_WP_/\tWP_/g' | sed 's/OCG/\tOCG/g' | sed 's/lcl.//g' > ${name}_${query}_blastp_${vers}.${var}.txt

# Double check that filtered output looks the way it should
echo "Lines in ${name}_${query}_blastp_${vers}.${var}.txt:"
wc -l ${name}_${query}_blastp_${vers}.${var}.txt
head ${name}_${query}_blastp_${vers}.${var}.txt

# Make CSV file from filtered output
more ${name}_${query}_blastp_${vers}.${var}.txt | awk '{print $2"_"$3","$1","$5}' | \
        sort -t, -k1,1 | sed '/.txt/d' | sed '/:::/d' > ${name}_${query}_blastp_${vers}.${var}.csv

echo "Lines in ${name}_${query}_blastp_${vers}.${var}.csv:"
wc -l ${name}_${query}_blastp_${vers}.${var}.csv
head ${name}_${query}_blastp_${vers}.${var}.csv

# Print sequence IDs for BLAST hits - this is necessary to use blastdbcmd to extract the amino acid sequences
more ${name}_${query}_blastp_${vers}.txt | sort -k2,2 -k6,6nr | sort -u -k2,2 | \
        awk '$4>=40&&$5>=90{print $2}' > ${name}_${query}_blastp_${vers}.${var}.seqids.txt

echo "Lines in ${name}_${query}_blastp_${vers}.${var}.seqids.txt:"
wc -l ${name}_${query}_blastp_${vers}.${var}.seqids.txt

# Use blastdbcmd to extract protein sequences with homology to representative toxin domains from the BLAST database
blastdbcmd -db $DB5 -entry_batch ${name}_${query}_blastp_${vers}.${var}.seqids.txt \
        -out ${name}_${query}_blastp_${vers}.${var}.seqs.fa

echo "Sequences in ${name}_${query}_blastp_${vers}.${var}.seqs.fa:"
grep -c '>' ${name}_${query}_blastp_${vers}.${var}.seqs.fa

# Calculate lengths of the extracted protein sequences
more ${name}_${query}_blastp_${vers}.${var}.seqs.fa | \
        awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | tail -n +2 | \
        sed 's/ .*\t/,/g' | sed 's/OCG/_OCG/g' | sort -t, -k1,1 | sed '/.txt/d' | sed '/:::/d' \
> ${name}_${query}_blastp_${vers}.${var}.protein_lengths.csv

echo "Lines in ${name}_${query}_blastp_${vers}.${var}.protein_lengths.csv:"
wc -l ${name}_${query}_blastp_${vers}.${var}.protein_lengths.csv
head ${name}_${query}_blastp_${vers}.${var}.protein_lengths.csv

# Use muscle to align extracted sequences for BLAST hits
$MUSCLEDIR/muscle3.8.31_i86linux64 -in ${name}_${query}_blastp_${vers}.${var}.seqs.fa -out ${name}_${query}_blastp_${vers}.${var}.seqs.aligned.fa

# Count the number of proteins belonging to each Rhs toxin domain cluster identified in each genome
more ${name}_${query}_blastp_${vers}.${var}.txt | awk '{print $1"__"$2}' | sort | uniq -c | sed 's/__/\t/g' | \
        awk '{print $3,$2,$1}' | sort -t, -k1,1 | sed '/:::/d' | sed '/.txt/d' > ${name}_${query}_${vers}_counts_per_genome.txt

more ${name}_${query}_${vers}_counts_per_genome.txt | awk '{print $1"_"$2","$3}' | sed 's/\.[0-9],/,/g' | \
        sort -t, -k1,1 | sed '/:::/d' | sed '/.txt/d' > ${name}_${query}_${vers}_counts_per_genome.sorted.csv
echo "Lines in ${name}_${query}_${vers}_counts_per_genome.sorted.csv:"
wc -l ${name}_${query}_${vers}_counts_per_genome.sorted.csv
head ${name}_${query}_${vers}_counts_per_genome.sorted.csv

# Count the total number of Rhs toxin domains identified in each genome
more ${name}_${query}_${vers}_counts_per_genome.txt | awk '{a[$2] += $3} END {for (i in a) print i","a[i]}' |
        sort -t, -k1,1 > ${name}_${query}_${vers}.${var}.total_counts.csv
echo "Lines in ${name}_${query}_${vers}.${var}.total_counts.csv:"
wc -l ${name}_${query}_${vers}.${var}.total_counts.csv
head ${name}_${query}_${vers}.${var}.total_counts.csv

###########################################
common=${name}_${query}_blastp_${vers}.${var} # Just to keep all this information in the file names without having to type it out over and over

# Join the filtered BLAST output CSV file, which contains the % identity and alignment length for each blast hit to the protein length CSV file, which contains the lengths of the proteins extracted from the BLAST database
# First, verify that the two files have the same number of lines
echo "Lines in ${common}.csv:"
wc -l ${common}.csv
head ${common}.csv
echo " "
echo "Lines in ${common}.protein_lengths.csv:"
wc -l ${common}.protein_lengths.csv
head ${common}.protein_lengths.csv

echo "Join percent identity and protein lengths tables:"
join -t, -o 1.1,1.2,1.3,2.2 ${common}.csv ${common}.protein_lengths.csv > ${common}.pident.length.csv
echo "Lines in ${common}.pident.length.csv:"
wc -l ${common}.pident.length.csv
head ${common}.pident.length.csv
echo " "
echo " "

# Insert a comma between genome assembly and protein accession numbers, drop protein accession number and link genome accession numbers to Rhs cluster names
more ${common}.pident.length.csv | sed 's/^\(GCF_[0-9]*\.[0-9]\)_/\1,/g' | awk -F, '{print $1"_"$3","$4";"$5}' | sort -t, -k1,1 > ${common}.pident.length.2.csv
wc -l ${common}.pident.length.2.csv
head ${common}.pident.length.2.csv
echo " "
echo " "

# Add Rhs toxin domain counts per genome to CSV file
echo "Lines in ${name}_${query}_${vers}_counts_per_genome.sorted.csv:"
wc -l ${name}_${query}_${vers}_counts_per_genome.sorted.csv
head ${name}_${query}_${vers}_counts_per_genome.sorted.csv

echo "Join counts-per-genome table:"
join -t, -o 1.1,1.2,2.2 ${common}.pident.length.2.csv ${name}_${query}_${vers}_counts_per_genome.sorted.csv | \
        sed 's/_Rhs_tox/,Rhs_tox/g' | sort -t, -k2,2 >  ${common}.pident.length.counts.csv
echo "Lines in ${common}.pident.length.counts.csv:"
wc -l ${common}.pident.length.counts.csv
head ${common}.pident.length.counts.csv
echo " "
echo " "

# Add total Rhs toxin domain counts per genome to CSV file
echo "Lines in ${name}_${query}_${vers}.${var}.total_counts.csv:"
wc -l ${name}_${query}_${vers}.${var}.total_counts.csv
head ${name}_${query}_${vers}.${var}.total_counts.csv

echo "Join total counts table:"
join -t, -1 2 -2 1 -o 1.1,1.2,1.3,1.4,2.2 ${common}.pident.length.counts.csv ${name}_${query}_${vers}.${var}.total_counts.csv | \
        sed 's/\(GCF_[0-9]*\)\.[0-9],/\1,/g' | sort -t, -k1,1 > ${common}.pident.length.counts.total_counts.csv
        echo "Lines in  ${common}.pident.length.counts.total_counts.csv:"
wc -l ${common}.pident.length.counts.total_counts.csv
head ${common}.pident.length.counts.total_counts.csv
echo " "
echo " "

# Add genome assembly metadata to CSV file
echo "Check metadata table:"
head $KEY
echo " "
echo "Join metadata table:"
join -t, -a1 -a2 -o 1.1,1.2,1.3,1.4,1.5,2.2 ${common}.pident.length.counts.total_counts.csv $KEY | sed 's/;/,/g' | \
        sed 's/:/,/g' > ${common}.pident.length.counts.total_counts.metadata.csv
echo "Lines in ${common}.pident.length.counts.total_counts.metadata.csv:"
wc -l ${common}.pident.length.counts.total_counts.metadata.csv
head ${common}.pident.length.counts.total_counts.metadata.csv
echo " "
echo " "


##########################
# Assign C-term clusters to N-term CORE
tox_clust=rhs_tox_${vers}_40perc_clusters.fa
DB_tox=rhs_tox_${vers}_prot_db
core_seqs=${name}_rhs_core_blastp_${vers}.${var}.seqs.fa
tox_core=assign_tox_domain_rhs_CORE_proteins

makeblastdb -in $tox_clust -out $DB_tox -dbtype prot -parse_seqids

blastp -db $DB_tox -query $core_seqs -outfmt "6 sseqid qseqid pident length bitscore" -out ${tox_core}_blastp_${vers}.txt -num_threads 16
more ${tox_core}_blastp_${vers}.txt | awk '$3>=40&&$4>=90{print $0}' | sort -k2,2 -k5,5nr | sort -u -k2,2 | \
        sed 's/lcl.//g' | awk '{print $2","$1":"$3":"$4}' | sort -t, -k1,1 > ${tox_core}_blastp_${vers}.${var}.csv

join -t, -a1 -e "no_core_assigned:NA:NA" -o 1.1,1.2,1.3,1.4,2.2 ${common}.pident.length.csv <( sort -t, -k1,1 ${tox_core}_blastp_${vers}.${var}.csv) | \
	sed 's/^\(GCF_[0-9]*\.[0-9]\)_/\1,/g' | awk -F, '{print $1"_"$3","$4";"$5";"$6}' | sort -t, -k1,1 | sed '/.txt/d' | sed  '/:::/d' > ${common}.pident.length.tox.csv
echo "Lines in ${common}.pident.length.tox.csv:"
wc -l ${common}.pident.length.tox.csv
head ${common}.pident.length.tox.csv

###########################
