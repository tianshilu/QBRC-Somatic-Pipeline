########  read data  ##########

# cns_file="~/software/cancer/somatic/example/cnv/tumor.cns"
# refFlat="/home2/twang6/data/genomes/hg38/hs38d1.fa_cnvkit.txt"
# somatic="/project/bioinformatics/Xiao_lab/shared/neoantigen/data/somatic/exome_seq/005T/somatic_mutations_hg38.txt" 
# output="~/software/cancer/somatic/example/cnv/cnv_results.txt"

args = commandArgs(trailingOnly=TRUE)
cns_file=args[1]
refFlat=args[2]
somatic=args[3]
output=args[4]

cns=read.table(cns_file,stringsAsFactors = F,header=T)
refFlat=read.table(refFlat,stringsAsFactors = F)
all_genes=unique(refFlat[,1])

########  assemble CNV data  ##############

cnv_gene=data.frame(gene=all_genes,cnv=0,count=0)
for (i in 1:dim(cns)[1])
{
  keep=cnv_gene$gene %in% strsplit(cns$gene[i],",")[[1]]
  cnv_gene$cnv[keep]=cnv_gene$cnv[keep]+cns$log2[i]
  cnv_gene$count[keep]=cnv_gene$count[keep]+1
}

########  normalize  ###############

cnv_gene$cnv=cnv_gene$cnv/cnv_gene$count
cnv_gene$cnv=cnv_gene$cnv-median(cnv_gene$cnv,na.rm = T)
cnv_gene$cnv[is.na(cnv_gene$cnv)]=0
cnv_gene$cnv=2^(cnv_gene$cnv)*2
cnv_gene=cnv_gene[order(cnv_gene$gene),]

######  adjust by CNV  ##########

if (somatic!="1")
{
  somatic=read.table(somatic,stringsAsFactors = F,header=T,sep="\t",colClasses=c("Ref"="character","Alt"="character"))
  tumor_fraction=quantile(somatic$Tumor_alt/(somatic$Tumor_alt+somatic$Tumor_ref),0.8)
  cnv_gene$cnv=(cnv_gene$cnv-(1-tumor_fraction)*2)/tumor_fraction
}

#######  write data  #############

# final adjustment
cnv_gene$cnv=round(cnv_gene$cnv,d=3)
cnv_gene$cnv[cnv_gene$cnv>4]=4
cnv_gene$cnv[cnv_gene$cnv<0]=0

# write to file
write.table(cnv_gene[,c("gene","cnv")],file=output,quote=F,row.names = F,sep="\t")

######  genome coverage of exome-seq data  ########

cnn_normal=read.table(sub("tumor.cns","normal.targetcoverage.cnn",cns_file),stringsAsFactors = F,header=T)
cnn_tumor=read.table(sub("tumor.cns","tumor.targetcoverage.cnn",cns_file),stringsAsFactors = F,header=T)
cnn_normal=cnn_normal[cnn_normal$gene!="-",]
cnn_tumor=cnn_tumor[cnn_tumor$gene!="-",]

coverage=matrix(paste(round(c(mean(cnn_normal$depth),mean(cnn_tumor$depth)),d=1),"x",sep=""),nrow = 1)
colnames(coverage)=c("Normal","Tumor")
write.table(coverage,file=sub("CNV_gene_level.txt","coverage.txt",output),quote=F,row.names = F,
            sep="\t")
