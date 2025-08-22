#get in variables from bash
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

### use together with ${scripts}/Tetragonula_carbonaria.sh

library(reshape2)
pdf.options(useDingbats = FALSE)

v=c(DIRECTORY,"/",LIB,"_",CHR,"_",CHR_START,"_",CHR_END,"_genome-unique.piRNAs.counts")
vname=paste(v,collapse="")
table=read.table(vname,header=F)
table_d <- dcast(table,V1~V3, value.var="V2")

### barplot
#scale=as.numeric(scale)
scale=max(table[,2])
widthw=as.numeric(LENGTH)/1000*40
#norm=5

p=c(DIRECTORY,"/",LIB,"_",CHR,"_",CHR_START,"_",CHR_END,"_genome-unique.piRNAs.pdf")
pname=paste(p,collapse="")
m=c(LIB,", ",CHR," ",CHR_START," ",CHR_END," genome-unique.piRNAs, black: plus_5end, gray: plus_3end, red: minus_5end, pink: minus_3end")
mains=paste(m,collapse="")

pdf(file=pname,width=widthw,height=7)
barplot(table_d$plus_5end,
        main=mains,
        ylim=c(-scale,scale),
        col="black",
        border=NA)
abline(h=0)
par(new=TRUE)
barplot(table_d$plus_3end,
        ylim=c(-scale,scale),
        col="gray",
        border=NA)
par(new=TRUE)
barplot(-table_d$minus_5end,
        ylim=c(-scale,scale),
        col="red",
        border=NA)
par(new=TRUE)
barplot(-table_d$minus_3end,
        ylim=c(-scale,scale),
        col="pink",
        border=NA)
dev.off()
