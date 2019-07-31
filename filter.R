# you may modify the following codes to suit some special needs
# need the corrplot, NMF, maftools, gplots, ape, and phangorn R packages
# best use R>=3.4.0. Can be executed from Rstudio or Rscript, but not in plain R console

# exclude TG/PDX samples from this analysis
# strongly recommend NOT to add un-matched mutation calling results to this analysis.

# the first argument is the path to the filter design file 
# it is a tab-delimited file, with three columns: sample_id, patient_id, folder
# folder is the path to the somatic mutation calling folder
# if patient id is not available, set to NA for all samples

# the second argument is the output folder to place all filtering results

# the third argument is the reference genome build, hg38, hg19, mm10

# the fourth argument is the path to the reference genome file

# the fifth argument is the minimum VAF of the mutations in the tumor sample (recommended: 0.01-0.05)

# the sixth argument is whether to filter out extremely long genes. TRUE/FALSE. Default is FALSE
# see below for list of genes. These genes usually turn out to have somatic mutations in 
# any cohort of patients

# Rscript filter.R ./example/filter.txt ./example/filter hg38 \
#   /home2/twang6/data/genomes/hg38/hs38d1.fa 0.01 FALSE

######  setting up  ##########

library(corrplot)
library(NMF)
library(maftools)
library(gplots)
library(phangorn)
library(ape) 
options(scipen=999)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {stop("Error: Not the correct number of arguments!")}
design=args[1]
path=args[2]
refBuild=args[3]
ref_genome=args[4]
min_tumor_vaf=as.numeric(args[5])
filter_long=args[6]=="TRUE"

if (Sys.getenv("RSTUDIO") == "1")
{
  library(rstudioapi)
  scriptPath=dirname(rstudioapi::getSourceEditorContext()$path)
}else
{
  args=commandArgs(trailingOnly = F)
  scriptPath=normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
}
source(paste(scriptPath,"/somatic_script/filter_functions.R",sep=""))

cosmic_genes=read.csv(paste(scriptPath,"/somatic_script/cancer_gene_census.csv",sep=""),row.names=1,
                      stringsAsFactors = F)

long_genes=c("TTN","KCNQ1OT1","MUC16","ANKRD20A9P","TSIX","SYNE1","ZBTB20","OBSCN",
  "SH3TC2","NEB","MUC19","MUC4","NEAT1","SYNE2","CCDC168","AAK1","HYDIN","RNF213",      
  "LOC100131257","FSIP2","MUC5B")  

design=read.table(design,stringsAsFactors = F,header=T)
design=design[!duplicated(design),]
if (colnames(design)[1]!="sample_id") 
  {stop("Error: Did you forgot the column headers for the design file?")}
design=design[order(design$patient_id,design$sample_id),]
if (!file.exists(path)) {dir.create(path)}

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
    tmp=read.table(file, stringsAsFactors = F,sep="\t",header = T,
                   colClasses=c("Ref"="character","Alt"="character"))
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
  tmp=read.table(file,stringsAsFactors = F,sep="\t",header = T,
                 colClasses=c("Ref"="character","Alt"="character"))
  tmp$sample_id=design$sample_id[i]
  tmp$patient_id=design$patient_id[i]
  
  # filter
  tmp=tmp[tmp$Tumor_alt/(tmp$Tumor_alt+tmp$Tumor_ref)>=min_tumor_vaf,]
  tmp=tmp[tmp$Func.refGene %in% c("exonic","exonic;splicing","splicing;exonic","splicing"),] # UTR or coding regions
  tmp=tmp[tmp$ExonicFunc.refGene!="synonymous SNV",] # non-S mutations
  if (dim(tmp)[1]==0) {next}
  if (!"SIFT_pred" %in% colnames(tmp)) {tmp$SIFT_pred=tmp$Polyphen2_HVAR_pred="."} # mouse
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
artefact=tmp[tmp>max(length(unique(mutations$patient_id))*0.2,2)]
cat(paste("Filtering ",round(sum(mutations$mutation %in% names(artefact))/dim(mutations)[1]*100),
          "% of mutations due to being exactly the same\n"))
mutations=mutations[!mutations$mutation %in% names(artefact),]

###########  write results  ###############

# all mutations
cosmic_role=cosmic_genes$Role.in.Cancer[match(tolower(mutations$Gene.refGene),tolower(rownames(cosmic_genes)))]
write.csv(cbind(mutations,cosmic_role),file=paste(path,"/all_mutations.csv",sep=""),row.names = F)

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

cosmic_role=cosmic_genes$Role.in.Cancer[match(tolower(rownames(sum_mut)),tolower(rownames(cosmic_genes)))]
write.csv(cbind(cosmic_role,sum_mut),file=paste(path,"/summary_mutations_details.csv",sep=""))
write.csv(cbind(cosmic_role,1*(sum_mut!="")),file=paste(path,"/summary_mutations.csv",sep=""))

# plot heatmap
tryCatch({
  tmp <-sum_mut
  tmp=tmp[apply(tmp!="",1,sum)/dim(tmp)[2]>0.05,]
  tmp=tmp[tolower(rownames(tmp)) %in% tolower(rownames(cosmic_genes)),]
  tmp=tmp[order(-apply(tmp!="",1,sum)),]
  percentages=paste(" (",round(apply(tmp!="",1,sum)/dim(tmp)[2]*100,d=1),"%)",sep="")
  rownames(tmp)=paste(rownames(tmp),percentages)
  
  pdf(file=paste(path,"/mutations_heatmap_cosmic.pdf",sep=""),
      width=max(dim(tmp)[2]/6+4,8),height=max(dim(tmp)[1]/10+5,8))
  pat_id=design$patient_id[match(colnames(sum_mut),design$sample_id)]
  
  tmp[tmp==""]<-0
  tmp[tmp=="exonic nonsynonymous SNV"]<-1
  tmp[tmp=="exonic frameshift substitution"]<-2
  tmp[tmp=="exonic stopgain"]<-3
  tmp[tmp=="exonic nonframeshift substitution"]<-4
  tmp[tmp=="splicing ."]<-5
  tmp[tmp=="exonic unknown"]<-6
  tmp[tmp=="exonic stoploss"]<-7
  tmp[!(tmp <= 8 )]<-8  #please change
  tmp=as.data.frame(tmp,stringsAsFactors=F)
  for (i in 1:dim(tmp)[2]) {tmp[,i]=as.numeric(tmp[,i])}
  tmp=as.matrix(tmp)
  
  breaks = c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5)
  col= c("gray", "red", "blue", "green", "cyan", "gold", "orchid", "plum","black")
  heatmap.2(tmp,breaks=breaks,col=col,lmat=rbind( c(4,5,2), c(6,1,3) ), key=F, extrafun=myplot, 
            Rowv=F,Colv=F,dendrogram="none",trace="none",srtCol=45,density.info="none",  
            colsep=which(pat_id[-1]!=pat_id[-length(pat_id)]),sepwidth=c(0.03,0.3),
            sepcolor="white", lhei=c(1,dim(tmp)[1]/10), lwid=c(1,6,1), keysize=1, 
            key.par = list(cex=0.5), cexRow= 0.8,cexCol = 0.8)
  dev.off()
},error=function(e) {print("Not plotting heatmap!")})

##########  plotting  ###################

# phylo tree
# pdf(file=paste(path,"/phylo_tree.pdf",sep=""),width=6,height=6)
# tmp=table(design$patient_id)
# phylo_tree_pats=names(tmp[tmp>=3])
# 
# for (phylo_tree_pat in phylo_tree_pats)
# {
#   # extract mutation data
#   mutations_pat=mutations[mutations$patient_id==phylo_tree_pat,]
#   mutations_pat$vaf=mutations_pat$Tumor_alt/(mutations_pat$Tumor_alt+mutations_pat$Tumor_ref)
#   mutations_pat$vaf=mutations_pat$vaf/quantile(mutations_pat$vaf,0.8)
#   mutations_pat$vaf[mutations_pat$vaf>1]=1
#   mutations_mat=matrix(0,ncol=length(unique(mutations_pat$mutation)),
#                        nrow=length(unique(mutations_pat$sample_id)))
#   colnames(mutations_mat)=unique(mutations_pat$mutation)
#   rownames(mutations_mat)=unique(mutations_pat$sample_id)
#   for (i in 1:dim(mutations_pat)[1])
#     {mutations_mat[mutations_pat$sample_id[i],mutations_pat$mutation[i]]=mutations_pat$vaf[i]}
# 
#   # transform
#   mutations_mat=rbind(mutations_mat[1,],mutations_mat)
#   rownames(mutations_mat)[1]="Normal"
#   mutations_mat[1,]=0
#   mutations_mat=t(mutations_mat)
# 
#   # plot
#   vaf=mutations_mat
#   thr = 0.05
#   vaf_bin <- vaf
#   vaf_bin[vaf>=thr] <- 1
#   vaf_bin[vaf<thr] <- 0
#   vaf_bin <- as.data.frame(vaf_bin)
#   phydat <- phyDat(vaf_bin,type="USER",levels=c(0,1))
#   partree <- pratchet(phydat,trace = F)
#   partree <- acctran(partree,phydat)
#   tree <- as.phylo(partree)
#   plot.phylo(tree,main=phylo_tree_pat,
#              type = "unrooted",direction = "rightwards",edge.width = 3,cex = 1.2)
# }
# dev.off()

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
dir=paste(path,"/lollipop",sep="")
if (!file.exists(dir)) {dir.create(file.path(dir))}
for (gene in getGeneSummary(laml)$Hugo_Symbol[1:30])
{
 pdf(file=paste(path,"/lollipop/",gene,".pdf",sep=""), width=10, height=3)
 tryCatch({lollipopPlot(maf = laml, gene = gene, AACol = 'AAChange',labelPos="all",repel=T,
  labPosAngle=45,domainLabelSize=1.5,printCount=T)},
   error=function(e) 1)
 dev.off()
}

#somatic signature analysis
laml.tnm = trinucleotideMatrix(maf = laml, ref_genome,  ignoreChr=NULL, useSyn = TRUE)

laml.sign = extractSignatures(mat = laml.tnm, n=4,nTry = 10, plotBestFitRes = FALSE)
write.csv(laml.sign$contributions,paste(path,"/mut_sig.csv",sep=""),row.names=TRUE)

pdf(file=paste(path,"/mut_sig.pdf",sep=""), width=5, height=5)
plotSignatures(laml.sign)
corrplot::corrplot(corr = laml.sign$coSineSimMat, col = RColorBrewer::brewer.pal(n = 9, name = 'Oranges'), is.corr = FALSE, tl.cex = 0.6, tl.col = 'black', cl.cex = 0.6)
dev.off()
