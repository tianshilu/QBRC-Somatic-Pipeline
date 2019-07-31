######  read vcfs  ##############

name_pkg <- c("statmod")
bool_nopkg <- !name_pkg %in% rownames(installed.packages())
if(sum(bool_nopkg) > 0){
  install.packages(name_pkg[bool_nopkg])
}
invisible(lapply(name_pkg, library, character.only = T)) # load multiple packages

options(scipen=999)
args=commandArgs(trailingOnly = TRUE)
wd=args[1]
normal=args[2]
setwd(wd)

######  define functions  ###########

extract_count<-function(fields,info,field)
{
  if (length(info)==0) {return(numeric(0))}
  split_info=strsplit(info,":")
  as.numeric(sapply(1:length(info),
    function(i) sub(",.*","",split_info[[i]][fields[[i]]==field])))
}

filter_vcf<-function(vcf,caller,type="somatic")
{
  cat(paste("Filtering",caller,type,"\n"))
  if (dim(vcf)[1]==0) {return(vcf)}
  vcf$V8=caller
  vcf=vcf[!grepl(",",vcf$V5),]
  vcf=vcf[vcf$V7=="PASS",]
  vcf$V3=vcf$V2+nchar(vcf$V4)-1
  vcf=vcf[vcf$V10!="." & vcf$V11!=".",]
  fields=strsplit(vcf$V9,":")
  
  # extract read count
  if (caller=="mutect")
  {
    if (dim(vcf)[1]==0)
    {
      vcf$normal_ref=vcf$normal_alt=vcf$tumor_ref=vcf$tumor_alt=numeric(0)
      return(vcf)
    }
    
    normal_ct=strsplit(sapply(strsplit(vcf$V10,":"),function(x) x[2]),",")
    vcf$normal_ref=as.numeric(sapply(normal_ct,function(x) x[1]))
    vcf$normal_alt=as.numeric(sapply(normal_ct,function(x) x[2]))
    
    tumor_ct=strsplit(sapply(strsplit(vcf$V11,":"),function(x) x[2]),",")
    vcf$tumor_ref=as.numeric(sapply(tumor_ct,function(x) x[1]))
    vcf$tumor_alt=as.numeric(sapply(tumor_ct,function(x) x[2]))
  }else if (caller %in% c("speedseq","shimmer","varscan","strelka"))
  {
    if (caller %in% c("speedseq","shimmer"))
    {
      RO="RO";AO="AO"
    }else if (caller %in% c("varscan"))
    {
      RO="RD";AO="AD"
    }else if (caller %in% c("strelka"))
    {
      RO="TAR";AO="TIR"
    }
    
    vcf$normal_ref=extract_count(fields,vcf$V10,RO)
    vcf$normal_alt=extract_count(fields,vcf$V10,AO)
    vcf$tumor_ref=extract_count(fields,vcf$V11,RO)
    vcf$tumor_alt=extract_count(fields,vcf$V11,AO)
  }else if (caller=="lofreq") # lofreq, only for unmatched samples
  {
    tmp=sapply(fields,function(x) strsplit(x,";")[[1]],simplify=F)
    total_count=sapply(tmp,function(x) as.numeric(sub("DP=","",grep("DP=",x,value=T))))
    var_count=round(total_count*
      sapply(tmp,function(x) as.numeric(sub("AF=","",grep("AF=",x,value=T)))))
    vcf$normal_alt=vcf$tumor_alt=var_count
    vcf$normal_ref=vcf$tumor_ref=total_count-var_count
  }else if (caller=="strelka_snp") 
  {
    if (dim(vcf)[1]>0)
    {
      vcf$normal_alt=vcf$tumor_alt=vcf$normal_ref=vcf$tumor_ref=NA
      for (i in 1:dim(vcf)[1])
      {
        tmp=strsplit(vcf$V10[i],":")[[1]]
        vcf$normal_ref[i]=as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V4[i],"U",sep="")]))
        vcf$normal_alt[i]=as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V5[i],"U",sep="")]))
        tmp=strsplit(vcf$V11[i],":")[[1]]
        vcf$tumor_ref[i]=as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V4[i],"U",sep="")]))
        vcf$tumor_alt[i]=as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V5[i],"U",sep="")]))
      }
    }else
    {
      vcf$normal_ref=vcf$normal_alt=vcf$tumor_ref=vcf$tumor_alt=numeric(0)
    }
  }else if (caller=="strelka_germline")
  {
    split_info=strsplit(vcf$V10,":")
    tmp=t(sapply(1:length(split_info),
      function(i) as.numeric(strsplit(split_info[[i]][fields[[i]]=="AD"],",")[[1]])))
    vcf$normal_alt=vcf$tumor_alt=tmp[,2]
    vcf$normal_ref=vcf$tumor_ref=tmp[,1]
  }
  
  # filter by read count and allele frequency
  if (caller!="strelka_germline") 
  {
    vcf=vcf[vcf$normal_ref+vcf$normal_alt>=7,]
    vcf=vcf[vcf$tumor_alt>=3,]
    if (type=="somatic")
    {
      vcf=vcf[vcf$normal_alt/(vcf$normal_ref+vcf$normal_alt)<
                vcf$tumor_alt/(vcf$tumor_ref+vcf$tumor_alt)/2,]
      vcf=vcf[vcf$normal_alt/(vcf$normal_ref+vcf$normal_alt)<0.05,]
    }else
    {
      vcf=vcf[vcf$normal_alt>=3,]
    } 
  }else # for tumor-only calling, make the calling super sensitive
  {
    vcf=vcf[vcf$normal_ref+vcf$normal_alt>=3,]
    vcf=vcf[vcf$normal_alt>=1,]
  }
  
  vcf
}

read_vcf<-function(file)
{
  cat(paste("Reading",file,"\n"))
  x=tryCatch({
    read.table(file,stringsAsFactors = F)
  }, error = function(e) { # give an empty VCF file
    cat(paste("Warning: failed to read",file,
              ". Maybe this caller didn't find any variants\n"))
    data.frame(V1=character(0),V2=numeric(0),V4=character(0),V5=character(0),V7=character(0),
      V8=character(0),normal_ref=numeric(0),normal_alt=numeric(0),
      tumor_ref=numeric(0),tumor_alt=numeric(0),stringsAsFactors = F)
  })
  unlink(file)
  x
}

#########  process tumor-only case first  ############

if (normal=="NA") 
{
  strelka=read_vcf("variants.vcf")
  strelka$V11=strelka$V10
  strelka_germline=filter_vcf(strelka,"strelka_germline","germline")

  write.table(strelka_germline[,c("V1","V2","V3","V4","V5","V8","normal_ref","normal_alt","tumor_ref",
    "tumor_alt")],file="germline_mutations.txt",col.names=F,row.names=F,sep="\t",quote=F)
  write.table(strelka_germline[,c("V1","V2","V3","V4","V5","V8","normal_ref","normal_alt","tumor_ref",
    "tumor_alt")],file="somatic_mutations.txt",col.names=F,row.names=F,sep="\t",quote=F)
  q()
}

#########  read somatic vcfs  ##################

speedseq=read_vcf("left_speedseq2.vcf")
speedseq=filter_vcf(speedseq,"speedseq")

mutect=read_vcf("left_mutect.vcf")
mutect=filter_vcf(mutect,"mutect")

shimmer=read_vcf("left_somatic_diffs.readct.vcf")
if (dim(shimmer)[1]>0) {shimmer$V7="PASS"}
shimmer=filter_vcf(shimmer,"shimmer")

varscan_indel=read_vcf("left_varscan.indel.Somatic.hc.vcf")
varscan_indel=filter_vcf(varscan_indel,"varscan")
varscan_snp=read_vcf("left_varscan.snp.Somatic.hc.vcf")
varscan_snp=filter_vcf(varscan_snp,"varscan")
varscan=rbind(varscan_indel,varscan_snp)

strelka_snp=read_vcf("left_passed.somatic.snvs.vcf")
strelka_snp=filter_vcf(strelka_snp,"strelka_snp")
if (dim(strelka_snp)[1]>0) {strelka_snp$V8="strelka"}
strelka_indel=read_vcf("left_passed.somatic.indels.vcf")
if (dim(strelka_indel)[1]>0) {strelka_indel$V7="PASS"}
strelka_indel=filter_vcf(strelka_indel,"strelka")
strelka=rbind(strelka_snp,strelka_indel)

lofreq_t=read_vcf("left_lofreq_t.vcf")
lofreq_t$V9=lofreq_t$V10=lofreq_t$V11=lofreq_t$V8
lofreq_t=filter_vcf(lofreq_t,"lofreq","germline")
lofreq_t$mutation=paste(lofreq_t$V1,lofreq_t$V2,lofreq_t$V4,lofreq_t$V5)

lofreq_n=read_vcf("left_lofreq_n.vcf")
lofreq_n$V9=lofreq_n$V10=lofreq_n$V11=lofreq_n$V8
lofreq_n=filter_vcf(lofreq_n,"lofreq","germline")
lofreq_n$mutation=paste(lofreq_n$V1,lofreq_n$V2,lofreq_n$V4,lofreq_n$V5)

lofreq_t$normal_ref=lofreq_n$normal_ref[match(lofreq_t$mutation,lofreq_n$mutation)]
lofreq_t$normal_alt=lofreq_n$normal_alt[match(lofreq_t$mutation,lofreq_n$mutation)]
lofreq=lofreq_t[,colnames(lofreq_t)!="mutation"]
lofreq_germline=lofreq # for germline filtering
lofreq=lofreq[is.na(lofreq$normal_alt) | (lofreq$normal_alt/(lofreq$normal_ref+lofreq$normal_alt)<
                lofreq$tumor_alt/(lofreq$tumor_ref+lofreq$tumor_alt)/2),]
lofreq=lofreq[is.na(lofreq$normal_alt) | 
                (lofreq$normal_alt/(lofreq$normal_ref+lofreq$normal_alt)<0.05),]
lofreq=lofreq[lofreq$V6>=10,]

#########  read and process germline vcfs  ###########

# varscan results
varscan_germline_indel=read_vcf("varscan.indel.Germline.vcf")
varscan_germline_snp=read_vcf("varscan.snp.Germline.vcf")
varscan_LOH_indel=read_vcf("varscan.indel.LOH.vcf")
varscan_LOH_snp=read_vcf("varscan.snp.LOH.vcf")
varscan_germline=rbind(varscan_germline_indel,varscan_germline_snp,varscan_LOH_indel,varscan_LOH_snp)
varscan_germline=filter_vcf(varscan_germline,"varscan","germline")

# lofreq results
keep=paste(varscan_germline$V2,varscan_germline$V3,varscan_germline$V4,varscan_germline$V5) %in%
  paste(lofreq_germline$V2,lofreq_germline$V3,lofreq_germline$V4,lofreq_germline$V5)
varscan_germline$V8[keep]="varscan,lofreq"

# output
write.table(varscan_germline[,c("V1","V2","V3","V4","V5","V8","normal_ref","normal_alt","tumor_ref",
  "tumor_alt")],file="germline_mutations.txt",col.names=F,row.names=F,sep="\t",quote=F)

#########  combine somatic vcfs  #############

vcf=rbind(lofreq,mutect,shimmer,speedseq,varscan,strelka)[,c("V1","V2","V4","V5","V8","normal_ref",
  "normal_alt","tumor_ref","tumor_alt")]
colnames(vcf)[1:5]=c("chr","pos","ref","alt","caller")
vcf$variant=paste(vcf$chr,vcf$pos,vcf$ref,vcf$alt)
if (dim(vcf)[1]==0) {stop("No valid mutations left in VCF file!\n")}

tmp1=aggregate(vcf$caller,by=list(vcf$variant),function(x) paste(x,collapse=","))
tmp2=aggregate(vcf[,c("normal_ref","normal_alt","tumor_ref","tumor_alt")],
               by=list(vcf$variant),function(x) mean(x,na.rm=T))
vcf=cbind(tmp1,tmp2[,-1])
vcf$chr=sapply(strsplit(vcf$Group.1," "),function(x) x[1])
vcf$pos=as.numeric(sapply(strsplit(vcf$Group.1," "),function(x) x[2]))
vcf$ref=sapply(strsplit(vcf$Group.1," "),function(x) x[3])
vcf$alt=sapply(strsplit(vcf$Group.1," "),function(x) x[4])
vcf$pos2=vcf$pos+nchar(vcf$ref)-1

# important! A variant must have been found by >=3 caller
vcf=vcf[grepl(",.*,",vcf$x),] 

vcf=vcf[,c("chr","pos","pos2","ref","alt","x","normal_ref","normal_alt","tumor_ref","tumor_alt")]
write.table(vcf,file="somatic_mutations.txt",col.names = F,row.names = F,sep="\t",quote=F)
