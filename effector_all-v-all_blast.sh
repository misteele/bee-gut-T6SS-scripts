##############################################
effector_all-v-all_blast.sh
######################
DIR= # Main work directory
FAA_DIR= # Path to directory containing protein FASTA files
PROT_DB_DIR= # Path to directory containing BLAST databases

vers= # Usually the date
name=rhs_tox # or vgrG
KEY= # CSV file containing genome assembly metadata
SEQS_full= # FASTA file containing Rhs toxins
SEQS_cterm= # FASTA file containing C-terminal domains from Rhs toxins
DB1=${name}_${vers}_full_length_prot_db
CLUST1=${name}_${vers}_40perc_clusters.fa
LIST1=${name}_${vers}_40perc_clusters.list.txt

#######################
name1= #name for analysis
var1=40pi90len90cov # %identity, length, and coverage cutoffs can be changed by altering awk '$3>=40&&$4>=90&&$5>=90{print $0}'
effectorkey=${name}_${query}_blastp_${vers}.${var}.pident.length.csv # Output of effector clustering script

# Make blast database from all identified Rhs proteins
makeblastdb -in $SEQS_full -out $DB1 -dbtype prot -parse_seqids

# Use protein BLAST to identify homologous Rhs toxin and Rhs toxin domain sequences (all v all BLAST)
blastp -db $DB1 -query $SEQS_cterm -outfmt "6 sseqid qseqid pident length qcovs bitscore" -out ${name1}_blastp_${vers}.txt -num_threads 16
more ${name1}_blastp_${vers}.txt | awk '$3>=40&&$4>=90&&$5>=90{print $0}' | sed 's/\t/,/g' | sed 's/1OCG/1_OCG/g' | \
        sort -t, -k1,1 > ${name1}_blastp_${vers}.${var1}.csv

# Join BLAST output to CSV file containing toxin domain names and protein lengths
echo "Lines in ${name1}_blastp_${vers}.${var1}.csv:"
wc -l ${name1}_blastp_${vers}.${var1}.csv
head ${name1}_blastp_${vers}.${var1}.csv
echo " "
echo "Lines in $effectorkey:"
wc -l $effectorkey
head $effectorkey
echo " "
join -t, -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.4,2.1 ${name1}_blastp_${vers}.${var1}.csv $effectorkey | \
        sed 's/\(Rhs_.*GCF_[0-9]*\)\.[0-9].*$/\1/g' | sort -t, -k9,9 > ${name1}_blastp_${vers}.${var1}.clust.length.csv
echo "Lines in ${name1}_blastp_${vers}.${var1}.clust.length.csv:"
wc -l ${name1}_blastp_${vers}.${var1}.clust.length.csv
head ${name1}_blastp_${vers}.${var1}.clust.length.csv
#echo " "

# Add genome assembly metadata to CSV file
join -t, -1 9 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2 ${name1}_blastp_${vers}.${var1}.clust.length.csv $KEY | \
        sed 's/:/,/g' > ${name1}_blastp_${vers}.${var1}.clust.length.metadata.csv

echo "Lines in ${name1}_blastp_${vers}.${var1}.clust.length.metadata.csv:"
wc -l ${name1}_blastp_${vers}.${var1}.clust.length.metadata.csv
head ${name1}_blastp_${vers}.${var1}.clust.length.metadata.csv

# Add column names to CSV file
echo "Subject,Query,Perc_identity,Align_length,Align_cov,Bitscore,Cluster,Length,Strain,Species,Taxa,Host,Host_taxa,Accession,Location,Location2,T6SS" > header.temp
cat header.temp ${name1}_blastp_${vers}.${var1}.clust.length.metadata.csv > ${name1}_blastp_${vers}.${var1}.clust.length.metadata.header.csv


###########################
# Similar approach for VgrG all v all BLAST

SEQS2=all_bee_gut_vgrg_blastp_0213.40pi300len.seqs.fa
DB2=vgrG_0213_prot_db
name2=vgrG_0213
var2=40pi300len90cov

blastp -db $DB2 -query $SEQS2 -outfmt "6 qseqid sseqid pident length qcovs bitscore" -out ${name2}_all-v-all_blastp_${vers}.txt -num_threads 16
more ${name2}_all-v-all_blastp_${vers}.txt | awk '$3>=40&&$4>=300&&$5>=90{print $0}' | sed 's/\t/,/g' | sed 's/1OCG/1_OCG/g' | \
        sort -t, -k1,1 > ${name2}_all-v-all_blastp_${vers}.${var2}.csv
echo "Lines in ${name2}_all-v-all_blastp_${vers}.${var2}.csv:"
wc -l ${name2}_all-v-all_blastp_${vers}.${var2}.csv
head ${name2}_all-v-all_blastp_${vers}.${var2}.csv
echo " "
echo "Lines in all_bee_gut_vgrg_blastp_0213.40pi300len.pident.length.csv:"
wc -l all_bee_gut_vgrg_blastp_0213.40pi300len.pident.length.csv
head all_bee_gut_vgrg_blastp_0213.40pi300len.pident.length.csv
echo " "
join -t, -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.4,2.1 ${name2}_all-v-all_blastp_${vers}.${var2}.csv all_bee_gut_vgrg_blastp_0213.40pi300len.pident.length.csv | \
        sed 's/\(VgrG_.*GCF_[0-9]*\)\.[0-9].*$/\1/g' | sort -t, -k9,9 > ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.csv
echo "Lines in ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.csv:"
wc -l ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.csv
head ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.csv
echo " "

echo "bee_gut_accession_key_with_T6SS.7.csv:"
head bee_gut_accession_key_with_T6SS.7.csv
echo " "
join -t, -1 9 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2 ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.csv bee_gut_accession_key_with_T6SS.csv | \
        sed 's/:/,/g' > ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.metadata.csv

echo "Lines in ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.metadata.csv:"
wc -l ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.metadata.csv
head ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.metadata.csv

echo "Query,Subject,Perc_identity,Align_length,Align_cov,Bitscore,Cluster,Length,Strain,Species,Taxa,Host,Host_taxa,Accession,Location,Location2,T6SS" > header.temp
cat header.temp ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.metadata.csv > ${name2}_all-v-all_blastp_${vers}.${var2}.clust.length.metadata.header.csv

