# The following script will reproduce the principal result of the study, i.e. the clustering of flanking genes and calculation of Cas association scores
# For running the code, it is required that FASTA v36 and MCL be installed on the system

# In addition a file containing protein sequences for genes flanking CRISPR-Cas systems is needed
# This file is assumed to be called flank5.faa, because in our study we used the five genes immediately flanking Type III systems up and downstream
# Sequence headers must consist only of the sequence ID, and no spaces in particular
# The file should ideally not contain any core Cas genes, as these will be mistakenly classified as accessory


# Also, a file containing protein sequences from the whole genomes of each of the organisms used to generate the above file is needed
# I.e. this file should contain all protein sequences from all genomes used above, including the above sequences, in addition to Cas proteins, as well as everything else on each genome like house-keeping genes etc.
# The file is used to estimate the specificity and degree of coevolution of the found accessory genes by estimating their background distribution
# we call this file genomes.faa. Sequence headers must be in the same format as flank5.faa
# We used locus IDs as headers in both files, so that works. E.g. >SSO1493 or >CLC_2119
# In addition to the two files above, corresposing tsv files called flank5.lengths.tab and genomes.lengths.tab are needed which contain lengths for each protein sequence in those files identified by the sequence ID in the following format:
# CLC_2119<tab>569<newline>
# and so on for all sequences in each file

