# bee-gut-T6SS-scripts

# t6ss_protein_phylogeny_blast.sh
# This script performs the first part of the process to make a phylogeny from concatenated alignments of multiple proteins that belong to a complex (in this case, a Type VI secretion system)
# This script requires a FASTA file containing sequences of each protein in the reference T6SS locus, a list of the sequence identifiers in that FASTA file, and a local BLAST database containing protein sequences for every genome of interest, as well as any outgroup genomes
# It also requires BLAST+ and MUSCLE 


# t6ss_ref_prok_blast.sh
# This script identifies possible outgroup sequences for making T6SS phylogenies by searching the NCBI Reference Prokaryote Representative Genomes database for any bacterial genomes that have homology to at least three proteins from the reference T6SS locus
