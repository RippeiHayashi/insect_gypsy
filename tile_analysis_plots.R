library(ggplot2)
library(reshape2)

pdf.options(useDingbats = FALSE)

setwd("${references}/${SPECIES}")

### Acheta domesticus kb tile counts
table=read.table("Adom_0.5kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0

### rep1
table_d2 <- subset(table_d, table_d$`Acheta_domesticus_ovaries_sRNA_rep1` > 1 & table_d$`Acheta_domesticus_eggs_sRNA_rep1` > 1)
p<-ggplot(table_d2)+
  geom_point(aes(x=`Acheta_domesticus_ovaries_sRNA_rep1`,
                 y=`Acheta_domesticus_eggs_sRNA_rep1`),size=2,shape=1,alpha=0.5,colour="gray")+
  geom_point(data=subset(table_d2, cluster_8446_72 > 0.5),
             aes(x=`Acheta_domesticus_ovaries_sRNA_rep1`,
                 y=`Acheta_domesticus_eggs_sRNA_rep1`),size=2,shape=1,alpha=0.6,colour="red")+
  geom_point(data=subset(table_d2, cluster_8446_26_left > 0.5 | cluster_8446_26_right > 0.5),
             aes(x=`Acheta_domesticus_ovaries_sRNA_rep1`,
                 y=`Acheta_domesticus_eggs_sRNA_rep1`),size=2,shape=1,alpha=0.6,colour="purple")+
  labs(title="Adom piRNAs, ovaries vs eggs per 0.5kb tiles, cluster_8446_72(red), cluster_8446_26(purple)",
       x="Acheta_domesticus_ovaries_sRNA_rep1 / RPKM", y="Acheta_domesticus_eggs_sRNA_rep1 / RPKM")+
  scale_x_log10(limits=c(1,10000))+
  scale_y_log10(limits=c(1,10000))+
  coord_fixed()+
  theme_bw()

### print out the scatter plot in pdf
pdf(file="Acheta_domesticus_ovary-vs-egg_kb-tile_rep1.pdf",width=8,height=8)
p
dev.off()

### writing the table in a comma separated file
write.csv(table_d, file="Acheta_domesticus_ovary-vs-egg_tile_coverage.csv", quote=F, row.names=F)

### sessionInfo
> sessionInfo()
R version 4.4.3 (2025-02-28)
Platform: x86_64-apple-darwin20
Running under: macOS Ventura 13.7.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

time zone: Australia/Sydney
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] reshape2_1.4.4 ggplot2_3.5.1 

loaded via a namespace (and not attached):
 [1] R6_2.6.1         magrittr_2.0.3   gtable_0.3.6     glue_1.8.0      
 [5] stringr_1.5.1    tibble_3.2.1     pkgconfig_2.0.3  lifecycle_1.0.4 
 [9] cli_3.6.4        scales_1.3.0     grid_4.4.3       vctrs_0.6.5     
[13] withr_3.0.2      compiler_4.4.3   plyr_1.8.9       tools_4.4.3     
[17] munsell_0.5.1    pillar_1.10.1    Rcpp_1.0.14      colorspace_2.1-1
[21] rlang_1.1.5      stringi_1.8.7   

