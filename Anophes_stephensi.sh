### the code was used to analyse the sequencing data for Anopheles stephensi.

### software used
### fastx_toolkit
http://hannonlab.cshl.edu/fastx_toolkit/
### bowtie-1.2.3-linux-x86_64
https://github.com/BenLangmead/bowtie/releases/tag/v1.2.3
### samtools/1.10
https://github.com/samtools/samtools/releases/tag/1.10
### bedtools/2.28.0
https://github.com/arq5x/bedtools2/releases/tag/v2.28.0
### kent-ucsc tools
https://github.com/ucscGenomeBrowser/kent
### R/4.0.0
https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz
### weblogo 3.7.8
https://github.com/WebLogo/weblogo

### define the variables
references="<directory where the genome sequence file and the gtf file are kept>"
${references}/tblastn-baits="<directory where the bait sequences of gypsy GAG POL ENV are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_sRNA="<name of the small RNA seq library>"
mkdir -p ${analysis}/${lib_sRNA}


### processing and analysing the small RNA sequencing libraries --- START ---

### below are the small RNA libraries processed in this part of the script:
Anopheles_stephensi_ovaries_sRNA_rep1
Anopheles_stephensi_ovaries_sRNA_rep2
Anopheles_stephensi_eggs_sRNA_rep1
Anopheles_stephensi_eggs_sRNA_rep1

### genome assembly used for the analysis
ASSEMBLY="GCF_013141755.1_UCI_ANSTEP_V1.0"
SPECIES="Aste"

### preparing for the analyses, do this before processing individual libraries --- START ---
### step 0-1: generating a bowtie index for miscellaneous RNAs
### download the gtf file and the genome
mkdir -p ${references}/${SPECIES}/indices
wget --directory-prefix="${references}/${SPECIES}" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/141/755/GCF_013141755.1_UCI_ANSTEP_V1.0/GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.gtf.gz
wget --directory-prefix="${references}/${SPECIES}" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/141/755/GCF_013141755.1_UCI_ANSTEP_V1.0/GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna.gz
gunzip ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.gz

### making fasta including tRNA, miRNA, rRNA, snoRNA and snRNA sequences
zcat ${references}/${SPECIES}/${ASSEMBLY}_genomic.gtf.gz |\
grep -e 'transcript_biotype "tRNA"' -e 'transcript_biotype "miRNA"' -e 'transcript_biotype "rRNA"' -e 'transcript_biotype "snoRNA"' -e 'transcript_biotype "snRNA"' |\
awk -F";" '{print $1,$(NF-2)}' | grep transcript_biotype |\
tr -d '"' | awk '{print $1,$4-1,$5,"Sgre_"$NF":"$10,".",$7}' | tr ' ' '\t' |\
bedtools getfasta -name -s -fi ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna -bed - |\
fasta_formatter - -t | awk '{print $1,toupper($2)}' | awk '!seen[$2]++' | awk '!seen[$1]++' |\
awk '{print ">"$1"\n"$2}' > ${references}/${SPECIES}/${ASSEMBLY}_genomic.miscRNA.fasta

### making a bowtie index of miscellaneous infrastructural RNAs
bowtie-build ${references}/${SPECIES}/${ASSEMBLY}_genomic.miscRNA.fasta ${references}/${SPECIES}/indices/${ASSEMBLY}_genomic.miscRNA


### step 0-2: generating a bowtie index for the genome
bowtie-build ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna ${references}/${SPECIES}/indices/UCI_ANSTEP_V1.0

### generate the size file
fasta_formatter -i ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna -t |\
awk '{print $1,length($NF)}' | tr ' ' '\t' > ${references}/${SPECIES}/${ASSEMBLY}_genomic.sizes


### step 0-3: generating a bed tile for genome unique 0.5kb tiles
### obtain every 25mer from the genome, only sense
### convert lowercase to uppercase first
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna |\
fasta_formatter -  -t | awk '{for(i=0;i<=length($NF)-25;i++) print ">"$1"@"i"_+:"substr($NF,i+1,25)"\n"substr($NF,i+1,25)
}' > ${references}/${SPECIES}/${ASSEMBLY}_25mers.plus.fasta

### obtain genome unique mappers
index_genome="${references}/${SPECIES}/indices/UCI_ANSTEP_V1.0"
bowtie -p 24 -f -v 0 -m 1 -S ${index_genome} ${references}/${SPECIES}/${ASSEMBLY}_25mers.plus.fasta |\
samtools view -@ 24 -bS - | bamToBed -i - > ${references}/${SPECIES}/${ASSEMBLY}_25mers.plus.unique-mappers.bed

### count the unique mappers in 0.5kb bins
### take if the centre of 25mer sits between the coordinate of 1 to 500 as the first 0.5kb bin
### take as the second 0.5kb bin if it sits between 501 to 1000
awk '{split($4,a,"@"); split(a[2],b,"_"); TILE[$1"_"int((b[1]+12)/500)]++
} END {for(var in TILE) print var,TILE[var]}' ${references}/${SPECIES}/${ASSEMBLY}_25mers.plus.unique-mappers.bed > ${references}/${SPECIES}/${ASSEMBLY}_0.5kbtiles.25mers.counts

### make a bed file that contains all kb tiles that are gt 85% mappability
awk '{split($1,a,"_"); if($2>425) print a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_+ . +""\n"a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_- . -"
}' ${references}/${SPECIES}/${ASSEMBLY}_0.5kbtiles.25mers.counts |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/${SPECIES}/${ASSEMBLY}_85percent_0.5kbtiles.bed


### step 0-4: extract transposon insertions --- START ---
### run RepeatMasker 4.1.0 by specifying "Arthropoda"
fastafile="${references}/${SPECIES}/${ASSEMBLY}_genomic.fna"
RepeatMasker -species "Arthropoda" -parallel 10 -xsmall ${fastafile} -dir ${references}/${SPECIES}

### extract sense and antisense gypsy insertions, these bed files are used in the IGV browser tracks
cat ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="+" ) print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.gypsy.plus.bed
cat ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.gypsy.minus.bed

### extract all transposon insertions
### remove "Simple_repeat", "Low_complexity", "rRNA" and "tRNA"
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print
}' ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9}' | tr ' ' '\t' > ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.allTE.plus.bed
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print
}' ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out |\
awk '{if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"}' | tr ' ' '\t' > ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.allTE.minus.bed


### step 0-5: running tblastn using baits of gypsy GAG ENV POL
### make blast database
mkdir -p ${references}/${SPECIES}/blastdb
fastafile="${references}/${SPECIES}/${ASSEMBLY}_genomic.fna"
makeblastdb -in ${fastafile} -parse_seqids -dbtype nucl -title "${ASSEMBLY}_NCBI_assembly" \
-out ${references}/${SPECIES}/blastdb/${ASSEMBLY}_NCBI_blast

### run tblastn
mkdir -p ${references}/${SPECIES}/tblastn_results
for ORF in GAG ENV POL; do
tblastn -db ${references}/${SPECIES}/blastdb/${ASSEMBLY}_NCBI_blast -outfmt 10 \
-query ${references}/tblastn-baits/gypsy_tblastn-baits_${ORF}.fasta \
-out ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out
done

### summarising tBlastn results
for ORF in GAG ENV POL; do
### extracting high-score insertions and merge entries, discarding entries with fewer than five hits.
### longer than 300nt for ENV and GAG, and longer than 1500nt for POL
if [[ ${ORF} == "POL" ]]; then
  ### extract ORFs tblastn hits
  cat ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out | tr ',' ' ' |\
  awk '{if ($9<$10 && $12>50) print $2,$9-1,$10,$1":"$7"-"$8,$5":"$6":"$11":"$12,"+";
else if ($9>$10 && $12>50) print $2,$10-1,$9,$1":"$7"-"$8,$5":"$6":"$11":"$12,"-"}' | tr ' ' '\t' | sort -k1,1 -k2,2n > ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out.bed

  bedtools merge -d 300 -s -c 5,6 -o count,distinct -i ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out.bed |\
  awk '{if ($3-$2>1500 && $4>4) print $1,$2,$3,$4,".",$5}' | tr ' ' '\t' > ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out.high-score.bed
else
  ### extract ORFs tblastn hits
  cat ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out | tr ',' ' ' |\
  awk '{if ($9<$10 && $12>30) print $2,$9-1,$10,$1":"$7"-"$8,$5":"$6":"$11":"$12,"+";
else if ($9>$10 && $12>30) print $2,$10-1,$9,$1":"$7"-"$8,$5":"$6":"$11":"$12,"-"}' | tr ' ' '\t' | sort -k1,1 -k2,2n > ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out.bed

  bedtools merge -d 300 -s -c 5,6 -o count,distinct -i ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out.bed |\
  awk '{if ($3-$2>300 && $4>4) print $1,$2,$3,$4,".",$5}' | tr ' ' '\t' > ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_${ORF}_tblastn_results.out.high-score.bed
fi
done

### preparing for the analyses, do this before processing individual libraries --- END ---


### processing small RNA sequencing libraries --- common part per library START ---

### step 1: trim the adapter sequence and collapse reads by sequence
mkdir -p ${analysis}/${lib_sRNA}
### we take R1 reads if they are from the paired end sequencing library.
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
# collapse reads of the same sequence to one entry while retaining the number of reads
fastx_clipper -c -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 18 -i - | fasta_formatter - -t |\
awk '{if(length($NF)>25 && length($NF)<49) READS[substr($NF,5,length($NF)-8)]++} END {
for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa


### step 2: run bowtie to map reads to the miscRNA and take unmapped reads
misc_index="${references}/${SPECIES}/indices/${ASSEMBLY}_genomic.miscRNA"
bowtie -f -v 1 -a -S ${misc_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.misc-mapped.bed


### step 3: obtain genome unique mappers
### run bowtie to map reads to the genome, allowing up to 1MM, unique mappers
index_genome="${references}/${SPECIES}/indices/UCI_ANSTEP_V1.0"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed


### step 4: generate the bigwig files for the genome unique mappers
### make bedgraph files
### only count reads that are greater than 22nt
### normalise read counts to one million genome unique mappers
TOTAL=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed)`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/${SPECIES}/${ASSEMBLY}_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.plus.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/${SPECIES}/${ASSEMBLY}_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,-$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.minus.bg

### make bigwig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.${strand}.bg \
${references}/${SPECIES}/${ASSEMBLY}_genomic.sizes ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.${strand}.bw
done


### step 5: measure coverage per 0.5kb tile for the genome unique mappers
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed \
-b ${references}/${SPECIES}/${ASSEMBLY}_85percent_0.5kbtiles.bed |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA}  '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/TOTAL*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.gt22nt.0.5kbtiles.counts


### step 6: measuring nucleotide frequencies around 3' ends of piRNAs mapping to ovarian somatic clusters in Anophles ovary sRNAseq reads
### details of ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed are found in the next block.
bedtools intersect -f 0.5 -s -wa -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-unique-mappers.UCI_ANSTEP_V1.0.bed -b ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); for(i=1;i<=b[1];i++) {if($3-$2>22 && $6=="+") print $1,$3-5,$3+6,$4,$5,$6;
else if($3-$2>22 && $6=="-") print $1,$2-6,$2+5,$4,$5,$6}}' | tr ' ' '\t' |\
bedtools getfasta -s -fi ${fastafile} -tab -bed - | awk '{print ">"$1"\n"toupper($2)}' | tr 'T' 'U' > ${analysis}/${lib_sRNA}/${lib_sRNA}_UCI_ANSTEP_V1.0_piRNA-clusters.3end_11nt_window.fasta

### running weblogo
TYPE="UCI_ANSTEP_V1.0_piRNA-clusters.3end_11nt_window"
weblogo -U probability -A rna -f ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.fasta -F pdf -n 50 -c classic -o ${analysis}/${lib_sRNA}/${lib_sRNA}_${TYPE}.gt22_logo_prob.pdf
weblogo -U probability -A rna -f ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.fasta -F logodata -n 50 -c classic -s large -t ${lib_sRNA}"_VanRij" -o ${analysis}/${lib_sRNA}/${lib_sRNA}_${TYPE}.gt22_logo_prob.txt

### processing small RNA sequencing libraries --- common part per library END ---


### post-processing of the tile coverage data for small RNA libraries --- START ---

### step 7-1: print out the piRNA cluster coordinate
printf "NC_050202.1 37699046 37730483 NC_050202_37 . -""\n""NC_050202.1 40888767 40903473 NC_050202_40 . -""\n""NC_050203.1 46572904 46613151 NC_050203 . +""\n" |\
tr ' ' '\t' > ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed

### step 7-2: add annotations of piRNA cluster tiles to the tiles
bedtools intersect -wo -s -a ${references}/${SPECIES}/${ASSEMBLY}_85percent_0.5kbtiles.bed \
-b ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed |\
awk '{print $4,$NF/500,$10}' > ${references}/${SPECIES}/Aste_0.5kb_tiles.counts

### step 7-3: add counts from small RNA libraries
for lib_sRNA in Anopheles_stephensi_ovaries_sRNA_rep1 Anopheles_stephensi_ovaries_sRNA_rep2 Anopheles_stephensi_eggs_sRNA_rep1 Anopheles_stephensi_eggs_sRNA_rep2; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.gt22nt.0.5kbtiles.counts >> ${references}/${SPECIES}/Aste_0.5kb_tiles.counts
done
### step 7-4: tile_analysis_plots.R was used to make scatter plots

### post-processing of the tile coverage data for small RNA libraries --- END ---


### for revision --- START ---
### measure proportions of sense and antisense gypsy piRNAs and ovarian somatic cluster piRNAs out of total genome mappers --- START ---
### step 8: run bowtie to map reads to the genome, allowing up to 1MM, all mappers with --all --best --strata option
### genome all mappers
index_genome="${references}/${SPECIES}/indices/UCI_ANSTEP_V1.0"
bowtie -p 12 -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort --parallel=12 -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.bed

### making a fasta file using genome all mappers and mapping information
awk '{READS[$4]++} END {split($4,a,"@"); {for(var in READS) print var,READS[var]}
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.bed |\
awk '{split($1,a,"@"); if(length(a[1])>22) print ">"$1"@"$2"\n"a[1]}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mappings.fa

### genome all mappers using the mapping information
index_genome="${references}/${SPECIES}/indices/UCI_ANSTEP_V1.0"
bowtie -p 12 -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mappings.fa |\
samtools view -@ 12 -bS - | bamToBed -i - | sort --parallel 12 -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed

### step 8-1: measure piRNA mappers intersecting with respective genomic regions. 
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed | awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count}' )`

TOTAL_Gypsy_S=`(bedtools intersect -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed \
-b <(cat ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.gypsy.*.bed ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_*_tblastn_results.out.high-score.bed) |\
awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]/a[3]} END {print count}')`

TOTAL_Gypsy_AS=`(bedtools intersect -S -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed \
-b <(cat ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.gypsy.*.bed ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_*_tblastn_results.out.high-score.bed) |\
awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]/a[3]} END {print count}')`

Cluster=`(bedtools intersect -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed -b ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]/a[3]} END {print count}' )`

Cluster_Gypsy_S=`(bedtools intersect -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed -b ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed |\
bedtools intersect -s -a - -b <(cat ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.gypsy.*.bed ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_*_tblastn_results.out.high-score.bed) |\
awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]/a[3]} END {print count}')`

Cluster_Gypsy_AS=`(bedtools intersect -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed -b ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed |\
bedtools intersect -S -a - -b <(cat ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.out.gypsy.*.bed ${references}/${SPECIES}/tblastn_results/${ASSEMBLY}_${SPECIES}_gypsy_*_tblastn_results.out.high-score.bed) |\
awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]/a[3]} END {print count}')`

printf "total_genome_mappers "$TOTAL" "${lib_sRNA}"\n""total_Gypsy_S "$TOTAL_Gypsy_S" "${lib_sRNA}"\n""total_Gypsy_AS "$TOTAL_Gypsy_AS" "${lib_sRNA}"\n""cluster_mappers "$Cluster" "${lib_sRNA}"\n""cluster_Gypsy_S "$Cluster_Gypsy_S" "${lib_sRNA}"\n""cluster_Gypsy_AS "$Cluster_Gypsy_AS" "${lib_sRNA}"\n" > ${analysis}/${lib_sRNA}/${lib_sRNA}_piRNAs_stats.txt
### measure proportions of sense and antisense gypsy piRNAs and ovarian somatic cluster piRNAs out of total genome mappers --- END ---


### phasing linkage analysis for Figure 6C and C' --- START ---
### step 9-1: generate pseudo piRNAs per every coordinate of VanRij_Aaeg_piRNA-clusters
cat ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed | awk '{print $2,$3-$2,$4,$6}' | while read START_position length TE strand; do
for ((i=0; i<length; i++)); do
j=$((i+1))
printf "%s\t%d\t%d\tAAAAAAAAAAAAAAAAAAAAAAA@1\t.\t+\n" "$TE" "$i" "$j"
printf "%s\t%d\t%d\tAAAAAAAAAAAAAAAAAAAAAAA@1\t.\t-\n" "$TE" "$i" "$j"
done > "${references}/${SPECIES}/${TE}_pseudo.piRNAs.txt"
done

### step 9-2: count 5end and 3end of plus and minus mapping reads
### total number of piRNA genome mappers
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed | awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count}' )`

### looping per clusters
cat ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed | awk '{print $2,$3-$2,$4,$6}' | while read START_position length TE strand; do

### extract reads mapping to the cluster
bedtools intersect -wa -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.bed \
-b <(awk -v TE=${TE} '{if($4==TE) print }' ${references}/${SPECIES}/UCI_ANSTEP_V1.0_piRNA-clusters.bed | tr ' ' '\t' ) |\
awk -v TE=${TE} -v START=${START_position} '{split($4,a,"@"); print TE,$2-START,$3-START,a[1]"@"a[2]/a[3],$5,$6}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.${TE}.bed

## count 5end and 3end of plus and minus mapping reads (>22nt)
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-all-mappers.UCI_ANSTEP_V1.0.mappings.${TE}.bed ${references}/${SPECIES}/${TE}_pseudo.piRNAs.txt |\
awk -v TOTAL=${TOTAL} '{split($4,a,"@"); if($6=="+" && length(a[1])>22) {PLUS5[$2]+=a[2]; PLUS3[$3-1]+=a[2]
} else if($6=="-" && length(a[1])>22) {MINUS5[$3-1]+=a[2]; MINUS3[$2]+=a[2]
}} END {for(var in PLUS5) print var,(PLUS5[var]-1)/TOTAL*1000000,"plus_5end""\n"var,(PLUS3[var]-1)/TOTAL*1000000,"plus_3end""\n"var,(MINUS5[var]-1)/TOTAL*1000000,"minus_5end""\n"var,(MINUS3[var]-1)/TOTAL*1000000,"minus_3end"
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_${TE}_3MM_mappers.counts

DIRECTORY_lib_sRNA="${analysis}/${lib_sRNA}"
if [[ ${strand} == "-" ]]; then
### if the cluster is on minus strand
Rscript ${scripts}/phasing_linkage_minus.R DIRECTORY=${DIRECTORY_lib_sRNA} LIB=${lib_sRNA} TE=${TE}
else
### if the cluster is on plus strand
Rscript ${scripts}/phasing_linkage_plus.R DIRECTORY=${DIRECTORY_lib_sRNA} LIB=${lib_sRNA} TE=${TE}
fi
done
### phasing linkage analysis for Figure 6C and C' --- END ---


### measure in-trans ping-pong linkage for Figure 8C--- START ---
### Soma enriched piRNAs (>7 fold enrichment) and piRNAs derived from ovarian somatic piRNA clusters were separately analysed.

### focusing on soma enriched piRNAs (>7 fold enrichment)
### step 10-1: extract g1g9, g2g10_revComp and last9 sequences from soma enriched piRNA
### ${analysis}/mosquitoes/Anopheles_stephensi_soma-enriched-tiles.bed includes tiles whose coverage in the ovarian sRNAseq library was more than 7 times greater than that of the embryonic sRNAseq library in both replicates.
mkdir -p ${analysis}/${lib_sRNA}/jellyfish/
bedtools intersect -wa -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed \
-b ${analysis}/mosquitoes/Anopheles_stephensi_soma-enriched-tiles.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)
}' > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_soma-enriched_g1g9.fasta

bedtools intersect -wa -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed \
-b ${analysis}/mosquitoes/Anopheles_stephensi_soma-enriched-tiles.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)
}' | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_soma-enriched_g2g10_revComp.fasta

bedtools intersect -wa -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed \
-b ${analysis}/mosquitoes/Anopheles_stephensi_soma-enriched-tiles.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)
}' > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_soma-enriched_last9.fasta

### step 10-2: count the occurrences of 9mers by jellyfish
CLASS="soma-enriched"
if [[ -f ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab ]]; then
rm ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab
fi
for TYPE in g2g10_revComp g1g9 last9; do
### jellyfish
jellyfish count -m 9 -s 100M -t 10 -o ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9 ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}.fasta
jellyfish dump ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9 > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa
### collect all counts
TOTAL=`(fasta_formatter -i ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa -t | awk '{count+=$1} END {print count}')`
### get rid of simple repeats
fasta_formatter -i ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa -t |\
awk -v TYPE=${TYPE} -v TOTAL=${TOTAL} '{if($2!~"AAAAAA" && $2!~"CCCCCC" && $2!~"GGGGGG" && $2!~"TTTTTT") print $2,$1/TOTAL*1000,TYPE}' >> ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab
done

### step 10-3: measure linkage using R
CLASS="soma-enriched"
Rscript ${scripts}/jellyfish_linkage.R DIRECTORY="${analysis}" LIB=${lib_sRNA} CLASS=${CLASS}


### focusing on somatic piRNA clusters
### step 11-1: extract g1g9, g2g10_revComp and last9 sequences from ovarian somatic cluster piRNAs
mkdir -p ${analysis}/${lib_sRNA}/jellyfish/
bedtools intersect -wa -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed \
-b <(cat ${DIRECTORY}/Diptera/Aste/UCI_ANSTEP_V1.0_piRNA-clusters.bed | grep -v "cluster_NC_050201") |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)
}' > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_soma_clusters_g1g9.fasta

bedtools intersect -wa -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed \
-b <(cat ${DIRECTORY}/Diptera/Aste/UCI_ANSTEP_V1.0_piRNA-clusters.bed | grep -v "cluster_NC_050201") |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)
}' | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_soma_clusters_g2g10_revComp.fasta

bedtools intersect -wa -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCI_ANSTEP_V1.0.bed \
-b <(cat ${DIRECTORY}/Diptera/Aste/UCI_ANSTEP_V1.0_piRNA-clusters.bed | grep -v "cluster_NC_050201") |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)
}' > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_soma_clusters_last9.fasta

### step 11-2: count the occurrences of 9mers by jellyfish
CLASS="soma_clusters"
if [[ -f ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab ]]; then
rm ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab
fi
for TYPE in g2g10_revComp g1g9 last9; do
### jellyfish
jellyfish count -m 9 -s 100M -t 10 -o ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9 ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}.fasta
jellyfish dump ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9 > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa
### collect all counts
TOTAL=`(fasta_formatter -i ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa -t | awk '{count+=$1} END {print count}')`
### get rid of simple repeats
fasta_formatter -i ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa -t |\
awk -v TYPE=${TYPE} -v TOTAL=${TOTAL} '{if($2!~"AAAAAA" && $2!~"CCCCCC" && $2!~"GGGGGG" && $2!~"TTTTTT") print $2,$1/TOTAL*1000,TYPE}' >> ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab
done

### step 11-3: measure linkage using R
CLASS="soma_clusters"
Rscript ${scripts}/jellyfish_linkage.R DIRECTORY="${analysis}" LIB=${lib_sRNA} CLASS=${CLASS}
### measure in-trans ping-pong linkage for Figure 8C--- END ---

### for revision --- END ---
