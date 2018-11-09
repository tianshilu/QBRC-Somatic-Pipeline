# you may modify the following codes to suit some special needs
# need the maftools, gplots and corrplot R packages
# best use R>=3.4.0

# exclude TG/PDX samples from this analysis
# strongly recommend NOT to add un-matched mutation calling results to this analysis.

# the first argument is the path to the filter design file 
# it is a tab-delimited file, with three columns: sample_id, patient_id, folder
# folder is the path to the somatic mutation calling folder
# if patient id is not available, set to NA for all samples

# the second argument is the output folder to place all filtering results

# the third argument is the reference genome build, hg38, hg19, etc

# the fourth argument is the path to the reference genome file

# the fifth argument is the minimum VAF of the mutations in the tumor sample (recommended: 0.01-0.05)

# the sixth argument is whether to filter out extremely long genes. TRUE/FALSE. Default is FALSE
# see below for list of genes. These genes usually turn out to have somatic mutations in 
# any cohort of patients

# Rscript filter.R ./example/filter.txt ./example/filter hg38 \
#   /home2/twang6/data/genomes/hg38/hs38d1.fa 0.01 FALSE

######  setting up  ##########

name_pkg <- c("corrplot","NMF","maftools","gplots")
bool_nopkg <- !name_pkg %in% rownames(installed.packages())
if(sum(bool_nopkg) > 0){
  install.packages(name_pkg[bool_nopkg],repos="https://cloud.r-project.org")
}
invisible(lapply(name_pkg, library, character.only = T)) # load multiple packages

options(scipen=999)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {stop("Error: Not the correct number of arguments!")}
design=args[1]
path=args[2]
refBuild=args[3]
ref_genome=args[4]
min_tumor_vaf=as.numeric(args[5])
filter_long=args[6]=="TRUE"
long_genes=c("TTN","KCNQ1OT1","MUC16","ANKRD20A9P","TSIX","SYNE1","ZBTB20","OBSCN",
  "SH3TC2","NEB","MUC19","MUC4","NEAT1","SYNE2","CCDC168","AAK1","HYDIN","RNF213",      
  "LOC100131257","FSIP2")  

design=read.table(design,stringsAsFactors = F,header=T)
if (colnames(design)[1]!="sample_id") 
  {stop("Error: Did you forgot the column headers for the design file?")}
design=design[order(design$patient_id,design$sample_id),]
if (!file.exists(path)) {dir.create(path)}
if (!file.exists(paste(path,"/each",sep=""))) {dir.create(paste(path,"/each",sep=""))}

########  check parental origin (for sample mislabeling)  ##################

if (any(is.na(design$patient_id))) 
{
  design$patient_id=design$sample_id
}else
{
  # read germline mutations
  mutations=list()
  for (i in 1:dim(design)[1])
  {
    file=paste(design$folder[i],"/germline_mutations_",refBuild,".txt",sep="")
    tmp=read.table(file, stringsAsFactors = F,sep="\t",header = T)
    mutations[[design$sample_id[i]]]=paste(tmp$Chr,tmp$Start,tmp$Ref,tmp$Alt)
  }
  
  # find overlap between germline mutations
  germline_overlap=matrix(NA,nrow=dim(design)[1],ncol=dim(design)[1],
    dimnames=list(design$sample_id,design$sample_id))
  for (i in 1:dim(design)[1])
  {
    for (j in i:dim(design)[1])
    {
      a=mutations[[i]]
      b=mutations[[j]]
      germline_overlap[i,j]=sum(a %in% b)/(length(a)+length(b))*2
      germline_overlap[j,i]=germline_overlap[i,j]
    }
  }
  
  # plot overlap matrix
  pdf(file=paste(path,"/germline_overlap.pdf",sep=""),width=6,height=6)
  heatmap.2(germline_overlap,Rowv=F,Colv=F,dendrogram="none",symm=T,trace="none",srtCol=45,
            density.info="none",key.title=NA)
  dev.off()
}

#######  read and filter mutations  ##########

mutations=c()
for (i in 1:dim(design)[1])
{
  # read data
  file=paste(design$folder[i],"/somatic_mutations_",refBuild,".txt",sep="")
  tmp=read.table(file,stringsAsFactors = F,sep="\t",header = T)
  tmp$sample_id=design$sample_id[i]
  tmp$patient_id=design$patient_id[i]
  
  # filter
  tmp=tmp[tmp$Tumor_alt/(tmp$Tumor_alt+tmp$Tumor_ref)>=min_tumor_vaf,]
  tmp=tmp[tmp$Func.refGene %in% c("exonic","exonic;splicing","splicing;exonic","splicing"),] # UTR or coding regions
  tmp=tmp[tmp$ExonicFunc.refGene!="synonymous SNV",] # non-S mutations
  tmp=tmp[!(tmp$SIFT_pred=="T" & tmp$Polyphen2_HVAR_pred=="B"),] # damaging missense mutations
  genes=table(tmp$Gene.refGene) # too many mutations on the same gene
  tmp=tmp[tmp$Gene.refGene %in% names(genes)[genes<=4],]
  if (filter_long) {tmp=tmp[!tmp$Gene.refGene %in% long_genes,]}
  mutations=rbind(mutations,tmp)
}

# further filter possible artefacts (the exact mutation in too many samples)
mutations$mutation=paste(mutations$Chr,mutations$Start,mutations$Ref,mutations$Alt)
tmp=mutations[,c("patient_id","mutation")]
tmp=tmp[!duplicated(tmp),]
tmp=table(tmp$mutation)
artefact=tmp[tmp>max(length(unique(mutations$patient_id))*0.2,1)]
cat(paste("Filtering ",round(sum(mutations$mutation %in% names(artefact))/dim(mutations)[1]*100),
          "% of mutations due to being exactly the same\n"))
mutations=mutations[!mutations$mutation %in% names(artefact),]

###########  write results  ###############

# all mutations
for (i in 1:dim(mutations)[1]) # newer versions of annovar do not label fs mutations completely
{
  annotations=strsplit(mutations$AAChange.refGene[i],",")[[1]]
  mutations$AAChange.refGene[i]=paste(sapply(strsplit(annotations,":"),function(x) {
    if (substr(x[length(x)],1,2)!="p.") {x=c(x,"p.X0X")}
    paste(x,collapse=":")
  }),collapse=",")
}
write.csv(mutations,file=paste(path,"/all_mutations.csv",sep=""),row.names = F)

tmp=mutations[,c("Chr","Start","End","Ref","Alt","Gene.refGene","ExonicFunc.refGene",
                 "AAChange.refGene","sample_id","Func.refGene")]
colnames(tmp)[9]="Tumor_Sample_Barcode"
tmp$GeneDetail.refGene=NA

for (i in 1:dim(tmp)[1]) # newer versions of annovar do not label fs mutations properly
{
  if (!grepl("frameshift",tmp$ExonicFunc.refGene[i])) {next}
  if (nchar(tmp$Ref[i])>nchar(tmp$Alt[i]))
  {
    tmp$ExonicFunc.refGene[i]=gsub("frameshift substitution","frameshift deletion",tmp$ExonicFunc.refGene[i])
  }else
  {
    tmp$ExonicFunc.refGene[i]=gsub("frameshift substitution","frameshift insertion",tmp$ExonicFunc.refGene[i])
  }
}
write.table(tmp[,c("Chr","Start","End","Ref","Alt","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene",
  "AAChange.refGene","Tumor_Sample_Barcode","Func.refGene")],file=paste(path,"/all_mutations.txt",sep=""),row.names = F,sep="\t",quote=F)

# summary of mutations
genes=unique(strsplit(paste(mutations$Gene.refGene,collapse=";"),";")[[1]])
samples=unique(mutations$sample_id)
sum_mut=matrix("",nrow=length(genes),ncol=length(samples))
rownames(sum_mut)=genes
colnames(sum_mut)=samples
for (i in 1:dim(mutations)[1])
{
  x=strsplit(mutations$Gene.refGene[i],";")[[1]]
  y=mutations$sample_id[i]
  z=paste(mutations$Func.refGene[i],mutations$ExonicFunc.refGene[i],sep=" ")
  sum_mut[x,y]=sub("^;","",paste(sum_mut[x,y],z,sep=";"),perl=T)
}
write.csv(sum_mut,file=paste(path,"/summary_mutations_details.csv",sep=""))
write.csv(1*(sum_mut!=""),file=paste(path,"/summary_mutations.csv",sep=""))

pdf(file=paste(path,"/mutations_heatmap.pdf",sep=""),width=8,height=20)
tmp=1-1*(sum_mut!="")
tmp=tmp[apply(tmp,1,mean)<0.9,]
tmp=tmp[order(apply(tmp,1,mean)),]
heatmap.2(tmp,Rowv=F,Colv=F,dendrogram="none",trace="none",srtCol=45,density.info="none")
dev.off()

for (sample in unique(mutations$sample_id))
{
  mutations_sample=mutations[mutations$sample_id==sample,]
  write.csv(mutations_sample,file=paste(path,"/each/",sample,".csv",sep=""),row.names = F)
  write.table(mutations_sample[,!colnames(mutations_sample) %in% c("sample_id","patient_id","mutation")],
    file=paste(path,"/each/",sample,".txt",sep=""),row.names = F,sep="\t",quote=F)
  
  mutations_sample$"#CHROM"=mutations_sample$Chr
  mutations_sample$"POS"=mutations_sample$Start
  mutations_sample$"ID"="."
  mutations_sample$"REF"=mutations_sample$Ref
  mutations_sample$"ALT"=mutations_sample$Alt
  mutations_sample$"QUAL"="."
  mutations_sample$"FILTER"="PASS"
  mutations_sample$"INFO"="."
  mutations_sample$"FORMAT"="RD:AD"
  mutations_sample$"NORMAL"=paste(ceiling(mutations_sample$Normal_ref),
                                  ceiling(mutations_sample$Normal_alt),sep=":")
  mutations_sample$"TUMOR"=paste(ceiling(mutations_sample$Tumor_ref),
                                 ceiling(mutations_sample$Tumor_alt),sep=":")
  write.table(mutations_sample[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",
    "FORMAT","NORMAL","TUMOR")],file=paste(path,"/each/",sample,".vcf",sep=""),row.names = F,
    sep="\t",quote=F)
}

##########  plotting  ###################

# get data into maf format
laml=annovarToMaf(annovar=paste(path,"/all_mutations.txt",sep=""),refBuild)
write.table(laml,file=paste(path,"/all_mutations.maf",sep=""),quote=F,sep="\t",row.names = F)
laml = read.maf(maf = paste(path,"/all_mutations.maf",sep=""), useAll = TRUE)

# summary
pdf(file=paste(path,"/summary.pdf",sep=""),width=8,height=6)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,top=20)
dev.off()

# vaf
pdf(file=paste(path,"/vaf.pdf",sep=""),width=6,height=4)
mutations$vaf=mutations$Tumor_alt/(mutations$Tumor_alt+mutations$Tumor_ref)
plot(density(mutations$vaf),xlab="VAF",main="Variant allele frequencies",lwd=2,ylim=c(0,10))
for (sample in unique(mutations$sample_id))
{
  if (sum(mutations$sample_id==sample)<2) {next}
  lines(density(mutations$vaf[mutations$sample_id==sample]),lwd=1,col="orange")
}
dev.off()

# oncoplot
pdf(file=paste(path,"/oncoplot.pdf",sep=""),width=10,height=10)
oncoplot(maf = laml, top = 50, removeNonMutated = TRUE)
dev.off()

# lollipop plot
# disabled for now, AAchange annotations are not correct
#dir=paste(path,"/lollipop",sep="")
#if (!file.exists(dir)) {dir.create(file.path(dir))}
#for (gene in getGeneSummary(laml)$Hugo_Symbol[1:20])
#{
#  pdf(file=paste(path,"/lollipop/",gene,".pdf",sep=""), width=10, height=3)
#  tryCatch({lollipopPlot(maf = laml, gene = gene, AACol = 'AAChange',labelPos="all")},
#    error=function(e) 1)
#  dev.off()
#}

#plot for gene cloud
pdf(file=paste(path,"/genecloud.pdf",sep=""), width=10, height=10)
geneCloud(input = laml, minMut = 5)
dev.off()

# identify significantly driver genes by maftools
# results maybe problematic, as it used AAChange for now.
#laml.sig = oncodrive(maf = laml, AACol = 'AAChange', minMut = 5, pvalMethod = 'zscore')
#pdf(file=paste(path,"/driver.pdf",sep=""), width=4, height=4)
#plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
#dev.off()

#somatic signature analysis
laml.tnm = trinucleotideMatrix(maf = laml, ref_genome,  ignoreChr=NULL, useSyn = TRUE)

laml.sign = extractSignatures(mat = laml.tnm, n=4,nTry = 10, plotBestFitRes = FALSE)
write.csv(laml.sign$contributions,paste(path,"/mut_sig.csv",sep=""),row.names=TRUE)

pdf(file=paste(path,"/mut_sig.pdf",sep=""), width=5, height=5)
plotSignatures(laml.sign)
corrplot::corrplot(corr = laml.sign$coSineSimMat, col = RColorBrewer::brewer.pal(n = 9, name = 'Oranges'), is.corr = FALSE, tl.cex = 0.6, tl.col = 'black', cl.cex = 0.6)
dev.off()
