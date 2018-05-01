#!/bin/bash
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

# Once all these five requirements are met the code can be run:
fasta36 flank5.faa flank5.faa -m 8 > flank5.fasta
cat flank5.fasta | join -1 1 -2 1 <(sort) <(sort flank5.lengths.tab) | join -1 2 -2 1 <(sort -k2,2) <(sort flank5.lengths.tab) | sed 's/ /\t/g' | awk '{print $2 "\t" $1 "\t" $11 "\t" $13/$14 "\t" ($8-$7)/(2*$13)+($10-$9)/(2*$14) "\t" ($7+$8-$9-$10)/($13+$14)}' | awk '{if ($3 <= 1) print}' | awk '{if ($5 >= 0.4) print}' | awk '{if (sqrt(($4-1)^2) - (sqrt(sqrt($5))-.8) + sqrt($6^2) <= 0.1) print $1 "\t" $2}' | mcl - -o - --abc -I 1.2 | perl -lane '$j++; $i = 0; while ($i <=$#F and $#F > 1) {$hash{$F[$i]} = $j; $i++} END {foreach $key (keys %hash) {print "$key\t$hash{$key}"}}' | sort -k2,2n -k1,1 > flank5.clusters.tab
# the clustering result is in the above output file

# In order to proceed through the next steps it will be necessary to create individual subsets of flank5.faa corresponding to each cluster in the above file and place them in a folder called "clusters/"
# I.e. for each cluster 1 to n in flank5.clusters.tab, make a file clusters/cluster1.faa to clusters/clustern.faa containing sequences for each protein cluster. The sequences can be extracted directly from flank5.faa.
# Once this is done you're ready to proceed to the next steps which will calculate the self similarity and coevolution/specificity scores:
seq 1 $(tail -1 flank5.clusters.tab | cut -f2) | while read line; do echo $line && fasta36 clusters/cluster$line.faa clusters/cluster$line.faa -m 8 > clusters/cluster$line.fasta; done
seq 1 $(tail -1 flank5.clusters.tab | cut -f2) | while read line; do echo -ne "$line\t" && cat clusters/cluster$line.fasta | join -1 1 -2 1 <(sort) <(sort <(cat clusters/cluster$line.fasta | awk '{if ($1 == $2) print $1 "\t" $12}')) | awk '{print $1 "\t" $2 "\t" int(100*$12/$13+.5)}' | awk '{if ($1 != $2) print}' | cut -f3 | awk '{i ++; sum += $1} END {print int(100*sum/i+.5)/100}'; done > flank5.clusters.selfsimilarity.tab
# the self similarities are contained in the above output file
# so that's half of the association score

# as for the other half:
seq 1 $(tail -1 flank5.clusters.tab | cut -f2) | while read line; do echo $line && fasta36 clusters/cluster$line.faa genomes.faa -m 8 > clusters/cluster$line.genomes.fasta; done
seq 1 $(tail -1 flank5.clusters.tab | cut -f2) | while read line; do cat clusters/cluster$line.genomes.fasta; done | join -1 2 -2 1 <(cat -n | sort -k2,2) <(sort genomes.lengths.tab) | join -1 3 -2 1 <(sort -k3,3) <(sort genomes.lengths.tab) | awk '{print $1, $2, $3, $12, $14/$15, ($9-$8)/(2*$14)+($11-$10)/(2*$15), ($8+$9-$10-$11)/($14+$15)}' | awk '{if ($4 <= 1) print}' | awk '{if ($6 >= 0.4) print}' | awk '{if (sqrt(($5-1)^2) - (sqrt(sqrt($6))-.8) + sqrt($7^2) <= 0.1) print $1, $2, $3}' | join -1 1 -2 1 <(cat) <(cat genomes.faa flank5.faa | grep '^>' | sed 's/^>//' | sort | uniq -c | awk '{print $2, $1}') | join -1 2 -2 1 <(sort -k2,2) <(sort flank5.clusters.tab) | join -1 5 -2 1 <(sort -k5,5n) <(cut -f2 flank5.clusters.tab | uniq -c | awk '{print $2, $1}') | sort -k4,4n | awk '{print $2, $3, $5-1, $1, $6}' | perl -lane '$i = 0 if $f ne $F[0]; $i++; print "$F[2]\t$F[3]" if $i <= $F[4]; $f = $F[0]' | perl -lane 'unless ($f eq $F[1]) {print "$s\t$i"; $i = 0; $s = 0} $i++; $s += $F[0]; $f = $F[1]; END {print "$s\t$i"}' | tail -n +2 | awk '{i++; print i "\t" int(10000*$1/$2+.5)/100}' > flank5.clusters.cospecificity.tab
# the combined Cas specificity and coevolution scores for each cluster are in the above output file

# To calculate the final Cas association score:
join -1 1 -2 1 flank5.clusters.cospecificity.tab flank5.clusters.selfsimilarity.tab | awk '{print $1, $2-$3}'

