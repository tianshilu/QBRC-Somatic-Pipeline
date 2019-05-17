args = commandArgs(trailingOnly=TRUE)
# args=c("/home2/twang6/software/cancer/somatic/example/evolution/design.txt",
#   "/home2/twang6/software/cancer/somatic/example/evolution/results","hg38")
design=read.table(args[1],stringsAsFactors = F,sep="\t")
output=args[2]
build=args[3]

dir.create(paste(output,"/pyclone",sep=""))
in_files=""
tumour_contents=""

for(i in 1:dim(design)[1])
{
  # prepare mutation file 
  file=paste(design[i,3],"/somatic_mutations_",build,".txt",sep="")
  tmp=read.table(file,stringsAsFactors = F,sep="\t",header = T,colClasses=c("Ref"="character","Alt"="character"))
  tmp$vaf=tmp$Tumor_alt/(tmp$Tumor_alt+tmp$Tumor_ref)
  tmp=tmp[tmp$vaf>=0.01,]
  tmp$vaf[tmp$vaf>quantile(tmp$vaf,0.95)]=quantile(tmp$vaf,0.95)
  tmp=tmp[tmp$Func.refGene %in% c("exonic","exonic;splicing","splicing;exonic","splicing"),] # UTR or coding regions
  tmp$mutation_id=paste(tmp$Chr,tmp$Start,tmp$Ref,tmp$Alt)
  tmp=tmp[sub("chr","",tmp$Chr) %in% as.character(1:100),]
  tumour_contents=paste(tumour_contents,min(1,2*max(tmp$vaf)))
  pyclone=tmp
  
  # prepare CNV file
  file=paste(design[i,2],"/CNV_gene_level.txt",sep="")
  tmp=read.table(file,stringsAsFactors = F,sep="\t",header = T)
  rownames(tmp)=tmp$gene
  pyclone=pyclone[pyclone$Gene.refGene %in% rownames(tmp),]
  pyclone$cnv=tmp[pyclone$Gene.refGene,2]
    
  # prepare input for pyclone
  mutation=data.frame(mutation_id=pyclone$mutation_id,ref_counts=floor(pyclone$Tumor_ref),
    var_counts=floor(pyclone$Tumor_alt),normal_cn=2,minor_cn=0,major_cn=round(pyclone$cnv))
  mutation$major_cn[mutation$major_cn==0]=1
  write.table(mutation,file=paste(output,"/pyclone/",design[i,1],'.tsv',sep=""),
      quote=F,col.names = T,row.names = F,sep='\t')
  in_files=paste(in_files,paste(design[i,1],'.tsv',sep=""))
}

writeLines(c(in_files,tumour_contents),paste(output,"/pyclone/parameters.txt",sep=""))
