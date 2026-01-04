# insect_gypsy
This repository includes code for processing and analysing small RNA sequencing data from ovaries and embryos for each species.

It additionaly contains the following files:

gypsy_tblastn-baits_${ORF}.fasta: bait sequences used for tBLASTn to identify gypsy insertions in the genome.

insect_PIWI_aa.with_AGO1.mafft.fasta: multiple sequence alignment of the insect Aubergine/PIWI protein sequences (N-PAZ-MID-PIWI domains) used for the phylogenetic analysis presented in Figure 5

insect_PIWI_aa.with_AGO1.mafft.fasta.treefile: the treefile used to build the phylogenetic tree presented in Figure 5

jellyfish_linkage.R: R code for measuring the in trans ping-pong linkage presented in Figure 8

measure_viral_sRNAs.R: R code for making a heatmap of the abundance of piRNAs (>22nt) mapping to the Liao ning virus genome, segments 1 to 12 for Figure 2.

phasing_linkage_minus.R, phasing_linkage_plus.R: R code for measuring piRNA phasing linkage for Figure 6

pingpong_cluster.R: R code for measuring the ping-pong linkage presented in Figure 7 and S10.

plotting_piRNAs_genome.unique.R: R code for extracting 5' and 3' end coverage of piRNAs mapping to a specified genomic region.

reference-mosquito-viral-genomes.fasta: nucleotide sequences of 680 mosquito viruses used to map piRNAs in Figure 2 and S3.

tile_analysis_plots.R: and example R code for making scatter plots, comparing the 0.5kb tile coverage of ovarian and embryonic piRNAs.
