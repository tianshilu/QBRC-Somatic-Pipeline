options(scipen=999)
args=commandArgs(trailingOnly = TRUE)
path=args[1]
build=args[2]
type=args[3]

setwd(path)
base=paste(type,"_mutations_",build,".txt",sep="")

mutations_tmp=tryCatch({read.table(paste(base,"_tmp.txt",sep=""),stringsAsFactors = F,sep="\t",header=F)},
                       error=function(e) {q()}) 
mutations_tmp$V1=as.numeric(sub("line","",mutations_tmp$V1))
mutations_tmp$V3=sub(",$","",mutations_tmp$V3,perl=T)
  
mutations=read.table(base,stringsAsFactors = F,sep="\t",header=T,colClasses=c("Ref"="character","Alt"="character"))  
mutations[mutations_tmp$V1-1,"AAChange.refGene"]=mutations_tmp$V3
write.table(mutations,file=base,sep="\t",row.names=F,quote=F)

