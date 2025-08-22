### the code was used to analyse the sequencing data for Tetragonula carbonaria.

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

### define the variables
references="<directory where the genome sequence file and the gtf file are kept>"
${references}/tblastn-baits="<directory where the bait sequences of gypsy GAG POL ENV are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_sRNA="<name of the small RNA seq library>"
mkdir -p ${analysis}/${lib_sRNA}


### processing and analysing the small RNA sequencing libraries --- START ---

### below are the small RNA libraries processed in this part of the script:
Tetragonula_carbonaria_ovaries_sRNA_rep1
Tetragonula_carbonaria_ovaries_sRNA_rep2
Tetragonula_carbonaria_eggs_sRNA_rep1
Tetragonula_carbonaria_eggs_sRNA_rep1

### genome assembly used for the analysis
ASSEMBLY="GCA_032399595.1_TetCarb_2.0"
SPECIES="Tcar"

### preparing for the analyses, do this before processing individual libraries --- START ---
### step 0-1: generating a bowtie index for miscellaneous RNAs
### download the gtf file and the genome asseembly of Apis mellifera
mkdir -p ${references}/${SPECIES}/indices
wget --directory-prefix="${references}/${SPECIES}" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gtf.gz
wget --directory-prefix="${references}/${SPECIES}" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
gunzip ${references}/${SPECIES}/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz

### making fasta including tRNA, miRNA, rRNA, snoRNA and snRNA sequences
zcat ${references}/${SPECIES}/GCF_003254395.2_Amel_HAv3.1_genomic.gtf.gz |\
grep -e 'transcript_biotype "tRNA"' -e 'transcript_biotype "miRNA"' -e 'transcript_biotype "rRNA"' -e 'transcript_biotype "snoRNA"' -e 'transcript_biotype "snRNA"' |\
awk -F";" '{print $1,$(NF-2)}' | grep transcript_biotype |\
tr -d '"' | awk '{print $1,$4-1,$5,"Sgre_"$NF":"$10,".",$7}' | tr ' ' '\t' |\
bedtools getfasta -name -s -fi ${references}/${SPECIES}/GCF_003254395.2_Amel_HAv3.1_genomic.fna -bed - |\
fasta_formatter - -t | awk '{print $1,toupper($2)}' | awk '!seen[$2]++' | awk '!seen[$1]++' |\
awk '{print ">"$1"\n"$2}' > ${references}/${SPECIES}/GCF_003254395.2_Amel_HAv3.1_genomic.miscRNA.fasta

### making a bowtie index of miscellaneous infrastructural RNAs
bowtie-build ${references}/${SPECIES}/GCF_003254395.2_Amel_HAv3.1_genomic.miscRNA.fasta ${references}/${SPECIES}/indices/GCF_003254395.2_Amel_HAv3.1_genomic.miscRNA

### step 0-2: generating a bowtie index for the genome
### download the genome asseembly of Tetragonula carbonaria
wget --directory-prefix="${references}/${SPECIES}" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/399/595/GCA_032399595.1_TetCarb_2.0/GCA_032399595.1_TetCarb_2.0_genomic.fna.gz
gunzip ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna.gz
### making a bowtie index of the Tetragonula carbonaria genome
bowtie-build ${references}/${SPECIES}/${ASSEMBLY}_genomic.fna ${references}/${SPECIES}/indices/Tcarb.CSA7p2

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
index_genome="${references}/${SPECIES}/indices/Tcarb.CSA7p2"
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
### use sequences from Apis mellifera
misc_index="${references}/${SPECIES}/indices/GCF_003254395.2_Amel_HAv3.1_genomic.miscRNA"
bowtie -f -v 1 -a -S ${misc_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.misc-mapped.bed


### step 3: obtain genome unique mappers
### run bowtie to map reads to the genome, allowing up to 1MM, unique mappers
index_genome="${references}/${SPECIES}/indices/Tcarb.CSA7p2"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.bed


### step 4: generate the bigwig files for the genome unique mappers
### make bedgraph files
### only count reads that are greater than 22nt
### normalise read counts to one million genome unique mappers
TOTAL=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.bed)`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/${SPECIES}/${ASSEMBLY}_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.plus.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/${SPECIES}/${ASSEMBLY}_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,-$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.minus.bg

### make bigwig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.${strand}.bg \
${references}/${SPECIES}/${ASSEMBLY}_genomic.sizes ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.${strand}.bw
done


### step 5: measure coverage per 0.5kb tile for the genome unique mappers
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.bed \
-b ${references}/${SPECIES}/${ASSEMBLY}_85percent_0.5kbtiles.bed |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA}  '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/TOTAL*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.gt22nt.0.5kbtiles.counts

### processing small RNA sequencing libraries --- common part per library END ---


### post-processing of the tile coverage data for small RNA libraries --- START ---

### step 6-1: print out the piRNA cluster coordinate
printf "R2_1 792355 1038850 R2_1_+ . +""\n""R2_1 792355 1038850 R2_1_- . -""\n" |\
tr ' ' '\t' > ${references}/${SPECIES}/Tcarb.CSA7p2_piRNA-clusters.bed

### step 6-2: add annotations of piRNA cluster tiles to the tiles
bedtools intersect -wo -s -a ${references}/${SPECIES}/${ASSEMBLY}_85percent_0.5kbtiles.bed \
-b ${references}/${SPECIES}/Tcar_idEriTena2_piRNA-clusters.bed |\
awk '{print $4,$NF/500,$10}' > ${references}/${SPECIES}/Tcar_0.5kb_tiles.counts

### step 6-3: add counts from small RNA libraries
for lib_sRNA in Tetragonula_carbonaria_ovaries_sRNA_rep1 Tetragonula_carbonaria_ovaries_sRNA_rep2 Tetragonula_carbonaria_eggs_sRNA_rep1 Tetragonula_carbonaria_eggs_sRNA_rep2; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Tcarb.CSA7p2.gt22nt.0.5kbtiles.counts >> ${references}/${SPECIES}/Tcar_0.5kb_tiles.counts
done
### step 6-4: tile_analysis_plots.R was used to make scatter plots

### post-processing of the tile coverage data for small RNA libraries --- END ---


### measuring ping-pong linkage across 0.5kb tiles from the cluster R2_1 --- START ---

### step 7-1: extracting 0.5kb tiles that have greater than 1000 cpm from the plus- and 10 cpm from the minus-strands from the two ovary libraries
for lib_sRNA in Tetragonula_carbonaria_ovaries_sRNA_rep1 Tetragonula_carbonaria_ovaries_sRNA_rep2; do
awk -v LIB=${lib_sRNA} '{split($1,a,":"); split(a[2],b,"_"); split(b[1],c,"-"); if($3==LIB && a[1]=="R2_1" && c[1]>809500 && c[1]<989500) print}' ${references}/${SPECIES}/Tcar_0.5kb_tiles.counts |\
awk '{split($1,a,"_"); if(a[3]=="-" && $2>1000) print a[1]"_"a[2]}' >>  ${references}/${SPECIES}/Tcarb_0.5kb_abundant_tiles.R2_1.list

awk -v LIB=${lib_sRNA} '{split($1,a,":"); split(a[2],b,"_"); split(b[1],c,"-"); if($3==LIB && a[1]=="R2_1" && c[1]>809500 && c[1]<989500) print}' ${references}/${SPECIES}/Tcar_0.5kb_tiles.counts |\
awk '{split($1,a,"_"); if(a[3]=="+" && $2>10) print a[1]"_"a[2]}' >>  ${references}/${SPECIES}/Tcarb_0.5kb_abundant_tiles.R2_1.list
done
cat ${references}/${SPECIES}/Tcarb_0.5kb_abundant_tiles.R2_1.list | sort | uniq -c |\
awk '{if($1==4) print $2}' > ${references}/${SPECIES}/Tcarb_0.5kb_abundant_tiles.R2_1.A_AS.list

### step 7-2: measuring the ping-pong score per tile
mkdir -p ${references}/${SPECIES}/End.plots/
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-unique-mappers.Tcarb.CSA7p2.bed | awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count}' )`
mkdir -p ${analysis}/${lib_sRNA}/End.plots.cluster

### step 7-2-1: printing out piRNA read at every coordinate
awk '{split($1,a,":"); split(a[2],b,"-"); print a[1],b[1],b[2]
}' ${references}/${SPECIES}/Tcarb_0.5kb_abundant_tiles.R2_1.A_AS.list | while read CHR CHR_START CHR_END; do
printf "${CHR} ${CHR_START} ${CHR_END}""\n" |\
awk '{for(i=0; i<=$3-$2-1; i++) print $1,$2+i,$2+i+1,"AAAAAAAAAAAAAAAAAAAAAAA@1 . +""\n"$1,$2+i,$2+i+1,"AAAAAAAAAAAAAAAAAAAAAAA@1 . -"
}' | tr ' ' '\t' > ${references}/${SPECIES}/End.plots/${CHR}_${CHR_START}_${CHR_END}_pseudo.piRNAs.txt

### step 7-2-2: intersect reads that overlap more than half of their length.
printf "${CHR} ${CHR_START} ${CHR_END}""\n" | tr ' ' '\t' |\
bedtools intersect -f 0.5 -wa -a ${analysis}/${lib_sRNA}/${lib_sRNA}_genome-unique-mappers.Tcarb.CSA7p2.bed -b - > ${analysis}/${lib_sRNA}/End.plots.cluster/${lib_sRNA}_${CHR}_${CHR_START}_${CHR_END}.bed

### step 7-2-3: count 5end and 3end of plus and minus mapping reads (>22nt)
cat ${analysis}/${lib_sRNA}/End.plots.cluster/${lib_sRNA}_${CHR}_${CHR_START}_${CHR_END}.bed ${analysis}/Tcarb/End.plots/${CHR}_${CHR_START}_${CHR_END}_pseudo.piRNAs.txt |\
awk -v TOTAL=${TOTAL} -v CHR_START=${CHR_START} -v CHR_END=${CHR_END} '{split($4,a,"@"); if($6=="+" && length(a[1])>22 && $2>=CHR_START && $3<=CHR_END) {PLUS5[$2]+=a[2]; PLUS3[$3-1]+=a[2]
} else if($6=="-" && length(a[1])>22 && $2>=CHR_START && $3<=CHR_END) {MINUS5[$3-1]+=a[2]; MINUS3[$2]+=a[2]
}} END {for(var in PLUS5) print var,(PLUS5[var]-1)/TOTAL*1000000,"plus_5end""\n"var,(PLUS3[var]-1)/TOTAL*1000000,"plus_3end""\n"var,(MINUS5[var]-1)/TOTAL*1000000,"minus_5end""\n"var,(MINUS3[var]-1)/TOTAL*1000000,"minus_3end"
}' > ${analysis}/${lib_sRNA}/End.plots.cluster/${lib_sRNA}_${CHR}_${CHR_START}_${CHR_END}_genome-unique.piRNAs.counts

DIRECTORY="${analysis}/${lib_sRNA}/End.plots.cluster"
LENGTH=$(echo "${CHR_END} - ${CHR_START} + 1" | bc)
### step 7-2-4: plotting 5' and 3' ends of sense and antisense mapping piRNAs
Rscript ${scripts}/plotting_piRNAs_genome.unique.R DIRECTORY=${DIRECTORY} LENGTH=${LENGTH} LIB=${lib_sRNA} CHR=${CHR} CHR_START=${CHR_START} CHR_END=${CHR_END}
### step 7-2-5: measure pingpong signature
Rscript ${scripts}/pingpong_cluster.R DIRECTORY=${DIRECTORY} LIB=${lib_sRNA} CHR=${CHR} CHR_START=${CHR_START} CHR_END=${CHR_END}
done

### measuring ping-pong linkage across 0.5kb tiles from the cluster R2_1 --- END ---
