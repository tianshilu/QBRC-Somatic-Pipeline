# summarizing results for human CNV calling and quality check metrics
# you may modify the following codes to suit some special needs

# the first argument is the path to the CNV design file 
# it is a tab-delimited file, with two columns: sample_id, folder
# folder is the path to the CNV calling folder 

# the second argument is the output folder

# the third argument is the folder to the reference genome files

# Rscript plot_cnv.R ./example/cnv_design.txt ./example/cnv_summarization 
#   /home2/twang6/data/genomes/hg38/

#########  prepare  ################

args = commandArgs(trailingOnly=TRUE)

design=read.table(args[1],stringsAsFactors = F,header=T)
if (any(colnames(design)!=c("sample_id","folder"))) {stop("Did you forgot the header line?")}

if (!file.exists(args[2])) {dir.create(args[2])}

chrom_size_file=list.files(args[3],pattern="chrom.sizes.txt",full.names = T)
ch_size=read.table(chrom_size_file,stringsAsFactors = F,row.names = 1)
#ch_size=ch_size[paste("chr",c(1:22,"X","Y"),sep=""),,drop=F]
ch_size=ch_size[paste("chr",c(1:22),sep=""),,drop=F]
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

#########  CNV  ###################

cnv_result_all=c()

# plotting function
plot_cnv<-function(title,ch_size,gene_loc,cnv_result)
{
  cnv_result=cnv_result[cnv_result$gene %in% gene_loc$V1,]
  gene_loc1=gene_loc[cnv_result$gene,]
  gene_loc1$V5=(gene_loc1$V5+ch_size[gene_loc1$V3,"cumulative_len"]-ch_size[gene_loc1$V3,"len"])/max_len
  
  plot(1:2,1:2,xlim=c(0,1),ylim=c(0,4),type="n",ylab="CNV",xlab="",xaxt="n",main=title)
  abline(v=0)
  for (i in 1:dim(ch_size)[1]) 
  {
    abline(v=ch_size$cumulative_len[i]/max_len)
    x1=ch_size$cumulative_len[i]-ch_size$len[i]
    x2=ch_size$cumulative_len[i]
    col=c('lightgrey','gray95')[1+(round(i/2)==i/2)]
    rect(x1/max_len,-0.12,x2/max_len,4.12,col=col,border=col)
    text((x1+x2)/2/max_len,cex=0.7,par("usr")[3]-0.2,labels=rownames(ch_size)[i],srt=45,pos=1,xpd=TRUE)
  }
  
  # plot CNVs
  for (i in 1:dim(gene_loc1)[1])
  {
    segments(x0=gene_loc1$V5[i],x1=gene_loc1$V5[i],y0=2,y1=cnv_result$cnv[i],lwd=0.1,
             c("red","blue")[1+(cnv_result$cnv[i]<2)])
  }
  abline(h=2)
}

pdf(paste(args[2],"/cnv_plot.pdf",sep=""),height=10,width=6)
par(mfrow=c(5,1))

for (j in 1:dim(design)[1])
{
  # read data
  cnv_result=read.table(paste(design$folder[j],"/CNV_gene_level.txt",sep=""),
                        stringsAsFactors = F,header=T)
  cnv_result_all=cbind(cnv_result_all,cnv_result$cnv)
  rownames(cnv_result_all)=cnv_result$gene
  
  # set up plotting area
  plot_cnv(design$sample_id[j],ch_size,gene_loc,cnv_result)
}

# generate a summary plot
cnv_result[,2]=apply(cnv_result_all,1,mean)
plot_cnv("Summary",ch_size,gene_loc,cnv_result)
dev.off()

colnames(cnv_result_all)=design$sample_id
write.csv(cnv_result_all,file=paste(args[2],"/cnv_results.csv",sep=""))

#########  summarize quality metrics  #######################

quality=c()
fields=c("Basic Statistics","Per base sequence quality","Per tile sequence quality",
         "Per sequence quality scores","Per base sequence content",
         "Per sequence GC content","Per base N content",
         "Sequence Length Distribution","Sequence Duplication Levels","Overrepresented sequences",
         "Adapter Content","Kmer Content" )

for (j in 1:dim(design)[1])
{
  # coverage
  coverage=read.table(paste(design$folder[j],"/coverage.txt",sep=""),stringsAsFactors = F,header=T)

  # tumor quality
  R12=list.files(paste(design$folder[j],"/tumor/fastqc",sep=""),full=T)
  R12=R12[!grepl("html",R12)]
  R12=R12[!grepl("zip",R12)]
  R12_quality_tumor=c(read.table(paste(R12[1],"/summary.txt",sep=""),sep="\t",stringsAsFactors=F)[,1],
                      read.table(paste(R12[2],"/summary.txt",sep=""),sep="\t",stringsAsFactors=F)[,1])
  
  tmp=read.table(paste(R12[1],"/summary.txt",sep=""),sep="\t",stringsAsFactors=F)[,2]
  names(R12_quality_tumor)=c(paste("T1_",tmp,sep=""),paste("T2_",tmp,sep=""))
  R12_quality_tumor=R12_quality_tumor[c(paste("T1_",fields,sep=""),paste("T2_",fields,sep=""))]
  
  # normal quality
  R12=list.files(paste(design$folder[j],"/normal/fastqc",sep=""),full=T)
  R12=R12[!grepl("html",R12)]
  R12=R12[!grepl("zip",R12)]
  R12_quality_normal=c(read.table(paste(R12[1],"/summary.txt",sep=""),sep="\t",stringsAsFactors=F)[,1],
                       read.table(paste(R12[2],"/summary.txt",sep=""),sep="\t",stringsAsFactors=F)[,1])

  tmp=read.table(paste(R12[1],"/summary.txt",sep=""),sep="\t",stringsAsFactors=F)[,2]
  names(R12_quality_normal)=c(paste("N1_",tmp,sep=""),paste("N2_",tmp,sep=""))
  R12_quality_normal=R12_quality_normal[c(paste("N1_",fields,sep=""),paste("N2_",fields,sep=""))]
  
  # combine
  quality=rbind(quality,c(coverage,R12_quality_normal,R12_quality_tumor))
}

rownames(quality)=design$sample_id
write.csv(quality,file=paste(args[2],"/quality.csv",sep=""))
