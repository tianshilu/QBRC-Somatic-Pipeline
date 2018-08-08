options(scipen=999)
args=commandArgs(trailingOnly = TRUE)
mutations=tryCatch(read.table(args[1],stringsAsFactors = F,sep="\t"),
  error=function(e) {
    cat(paste("Warning: No lines in ",args[1],"!\n",sep=""))
    q()
  }
)

colnames(mutations)=c("Chr","Start","End","Ref","Alt","Caller","Normal_ref","Normal_alt",
                      "Tumor_ref","Tumor_alt")
annotations=read.table(args[2],stringsAsFactors = F,sep="\t",header=T)

as_numeric<-function(x)
{
  x=as.numeric(x)
  x[is.na(x)]=0
  x
}

annotations=cbind(mutations,annotations[,-c(1:5)])
suppressWarnings({annotations$cosmic70=as_numeric(annotations$cosmic70)})
suppressWarnings({annotations$esp6500siv2_all=as_numeric(annotations$esp6500siv2_all)})
suppressWarnings({annotations$ExAC_ALL=as_numeric(annotations$ExAC_ALL)})
suppressWarnings({annotations$X1000g2015aug_all=as_numeric(annotations$X1000g2015aug_all)})

#annotations=annotations[annotations$cosmic70<0.01,]
annotations=annotations[annotations$esp6500siv2_all<0.01,]
annotations=annotations[annotations$ExAC_ALL<0.01,]
annotations=annotations[annotations$X1000g2015aug_all<0.01,]

annotations=annotations[,c("Chr","Start","End","Ref","Alt","Caller","Normal_ref","Normal_alt",
 "Tumor_ref","Tumor_alt","Func.refGene","Gene.refGene","ExonicFunc.refGene",
 "AAChange.refGene","SIFT_pred","Polyphen2_HVAR_pred")]
write.table(annotations,file=args[3],sep="\t",row.names=F,quote=F)
