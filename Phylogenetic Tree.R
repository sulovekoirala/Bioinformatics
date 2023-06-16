# Loading Packages
library(Biostrings)
library(ape)
library(msa)

# Download the 16S rRNA sequences of some streptococcus bacteria from GenBank
seq1 <-  readDNAStringSet("zoo16s.fasta")
seq2 <-  readDNAStringSet("equi16s.fasta")
seq3 <-  readDNAStringSet("oralis16s.fasta")
seq4 <-  readDNAStringSet("salivarius16s.fasta")
# Create a DNAStringSet object from the sequences
seqs <- DNAStringSet(c(seq1, seq2, seq3, seq4))

# Align the sequences using msa() function from seqinr
aln <- msa(seqs)

# Convert the alignment to a DNAbin object
aln <- as.DNAbin(aln)

# Calculate the pairwise distances using dist.dna() function from ape
dist <- dist.dna(aln)

# Construct a phylogenetic tree using nj() function from ape
tree <- nj(dist)

# Plot the tree using plot.phylo() function from ape
plot.phylo(tree)
