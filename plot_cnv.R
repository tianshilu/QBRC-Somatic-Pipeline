# plotting results for human CNV calling
# you may modify the following codes to suit some special needs

# the first argument is the path to the CNV design file 
# it is a tab-delimited file, with two columns: sample_id, file
# file is the path to the CNV calling file (CNV_gene_level.txt)

# the second argument is the output pdf file

# the third argument is the folder to the reference genome files

# Rscript plot_cnv.R ./example/cnv_design.txt ./example/cnv.pdf /home2/twang6/data/genomes/hg38/

#########  prepare  ################
args = commandArgs(trailingOnly = TRUE)
design=read.table(args[1],stringsAsFactors = F,header=T)

chrom_size_file=list.files(args[3],pattern="chrom.sizes.txt",full.names = T)
ch_size=read.table(chrom_size_file,stringsAsFactors = F,row.names = 1)
ch_size=ch_size[paste("chr",1:22,sep=""),,drop=F]
colnames(ch_size)="len"
ch_size$len=ch_size$len/1.0
ch_size$cumulative_len=ch_size$len[1]
for (i in 2:dim(ch_size)[1]) {ch_size$cumulative_len[i]=ch_size$cumulative_len[i-1]+ch_size$len[i]}
max_len=ch_size$cumulative_len[dim(ch_size)[1]]

gene_loc_file=list.files(args[3],pattern="_cnvkit.txt",full.names = T)
gene_loc=read.table(gene_loc_file,stringsAsFactors = F)
gene_loc=gene_loc[,c(1,3,5)]
gene_loc=gene_loc[!duplicated(gene_loc$V1),]
rownames(gene_loc)=gene_loc$V1

#########  plot  ###################

pdf(args[2],height=10,width=6)
par(mfrow=c(5,1))

for (j in 1:dim(design)[1])
{
  # read data
  cnv_result=read.table(design$file[j],stringsAsFactors = F,header=T)
  cnv_result=cnv_result[cnv_result$gene %in% gene_loc$V1,]
  gene_loc1=gene_loc[cnv_result$gene,]
  gene_loc1$V5=(gene_loc1$V5+ch_size[gene_loc1$V3,"cumulative_len"]-ch_size[gene_loc1$V3,"len"])/max_len
  
  # set up plotting area
  plot(1:2,1:2,xlim=c(0,1),ylim=c(0,4),type="n",ylab="CNV",xlab="",xaxt="n",main=design$sample_id[j])
  abline(v=0)
  for (i in 1:dim(ch_size)[1]) 
  {
    abline(v=ch_size$cumulative_len[i]/max_len)
    text((ch_size$cumulative_len[i]-ch_size$len[i]/2)/max_len, cex=0.7,
         par("usr")[3] - 0.2, labels = paste("chr",i,sep=""), srt = 45, pos = 1, xpd = TRUE)
  }
  
  # plot CNVs
  cnv_result$cnv[cnv_result$cnv>4]=4
  points(x=gene_loc1$V5,y=cnv_result$cnv,pch=19,cex=0.1,col=c("red","blue")[1+(cnv_result$cnv<2)])
  abline(h=2)
}

dev.off()
