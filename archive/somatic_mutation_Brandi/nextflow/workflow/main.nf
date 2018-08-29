#!/usr/bin/env nextflow

// Default parameter values to run tests
params.bams="$baseDir/../test_data/*.bam"
params.design="$baseDir/../test_data/design.txt"
params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.targetpanel="$params.genome/UTSWV2.bed"

dbsnp="$params.genome/dbSnp.vcf.gz"
cosmic="$params.genome/cosmic.vcf.gz"

design_file = file(params.design)
bams=file(params.bams)

reffa=file("$params.genome/genome.fa")
index_path = file(params.genome)

target_panel = file(params.targetpanel)
dbsnp=file(dbsnp)

strelkaconfig="/cm/shared/apps/strelka/1.0.15/etc/strelka_config_bwa_default.ini"
confstrelka=file(strelkaconfig)

snpeff_vers = 'GRCh38.82';
if (params.genome == '/project/shared/bicf_workflow_ref/GRCm38') {
   snpeff_vers = 'GRCm38.82';
}
if (params.genome == '/project/shared/bicf_workflow_ref/GRCh37') {
   snpeff_vers = 'GRCh37.75';
}

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def prefix = []

new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    tidx = header.findIndexOf{it == 'TumorID'};
    nidx = header.findIndexOf{it == 'NormalID'};
    oneidx = header.findIndexOf{it == 'TumorBAM'};
    twoidx = header.findIndexOf{it == 'NormalBAM'};
    if (twoidx == -1) {
       twoidx = oneidx
       }      
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      prefix << tuple(row[tidx],row[nidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }
	  
} 
}

if( ! prefix) { error "Didn't match any input files with entries in the design file" }

process indexbams {
  input:
  set tid,nid,file(tumor),file(normal) from prefix
  output:
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into mutectbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into ssbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into shimmerbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into vscanbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into virmidbam
//  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into strelkabam
  script:
  """
  module load speedseq/20160506 samtools/intel/1.3
  sambamba index -t 30 ${tumor}
  sambamba index -t 30 ${normal}
  """
}

process sstumor {

  publishDir "$baseDir/output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from ssbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.sspanel.vcf.gz") into ssvcf
  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 speedseq/20160506 bcftools/intel/1.3 vcftools/0.1.11
  speedseq somatic -q 10 -w ${target_panel} -t 30 -o ${tid}.sssom ${reffa} ${normal} ${tumor}
  vcf-annotate -H -n --fill-type ${tid}.sssom.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar filter --pass '((QUAL >= 10) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bgzip > ${tid}_${nid}.sspanel.vcf.gz
  """
}
process mutect {

  publishDir "$baseDir/output", mode: 'copy'

  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from mutectbam

  output:
  set val("${tid}_${nid}"),file("${tid}_${nid}.pmutect.vcf.gz") into mutectvcf
  script:
  """
  module load python/2.7.x-anaconda gatk/3.5  bcftools/intel/1.3 bedtools/2.25.0 snpeff/4.2 vcftools/0.1.11
  cut -f 1 ${index_path}/genomefile.chr.txt | xargs -I {} -n 1 -P 10 sh -c "java -Xmx10g -jar \$GATK_JAR -R ${reffa} -D ${dbsnp} -T MuTect2 -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -I:tumor ${tumor} -I:normal ${normal} --cosmic ${cosmic} -o ${tid}.{}.mutect.vcf -L {}"
  vcf-concat ${tid}*.vcf | vcf-sort | vcf-annotate -n --fill-type | java -jar \$SNPEFF_HOME/SnpSift.jar filter -p '((FS <= 60) & GEN[*].DP >= 10)' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bgzip > ${tid}_${nid}.pmutect.vcf.gz
  """
}
process varscan {

  publishDir "$baseDir/output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from vscanbam
  output:
  set val("${tid}_${nid}"),file("${tid}_${nid}.varscan.vcf.gz") into varscanvcf
  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3 VarScan/2.4.2 speedseq/20160506 vcftools/0.1.11
  sambamba mpileup -L ${target_panel} -t 30 ${tumor} --samtools "-C 50 -f ${reffa}"  > t.mpileup
  sambamba mpileup -L ${target_panel} -t 30 ${normal} --samtools "-C 50 -f ${reffa}"  > n.mpileup
  VarScan somatic n.mpileup t.mpileup ${tid}.vscan --output-vcf 1
  VarScan copynumber n.mpileup t.mpileup ${tid}.vscancnv 
  vcf-concat ${tid}.vscan*.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar \$SNPEFF_HOME/SnpSift.jar filter '((exists SOMATIC) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bedtools intersect -header -a stdin -b ${target_panel} |bgzip > ${tid}_${nid}.varscan.vcf.gz
  """
}
process shimmer {

  publishDir "$baseDir/output", mode: 'copy'

  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from shimmerbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.shimmer.vcf.gz") into shimmervcf
  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3  shimmer/0.1.1 vcftools/0.1.11
  shimmer.pl --minqual 25 --ref ${reffa} ${normal} ${tumor} --outdir shimmer 2> shimmer.err
  perl $baseDir/scripts/add_readct_shimmer.pl
  vcf-annotate -n --fill-type shimmer/somatic_diffs.readct.vcf | java -jar \$SNPEFF_HOME/SnpSift.jar filter '(GEN[*].DP >= 10)' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' | bedtools intersect -header -a stdin -b ${target_panel} | bgzip > ${tid}_${nid}.shimmer.vcf.gz

  """
}

// process strelka {

//   publishDir "$baseDir/output", mode: 'copy'
//   cpus 32
//   input:
//   set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from strelkabam
//   output:
//   file("${tid}.strelka.vcf.gz") into strelkavcf
//   script:
//   """
//   module load bedtools/2.25.0 snpeff/4.2 speedseq/20160506 vcftools/0.1.11 strelka/1.0.15
//   cp ${confstrelka} config.ini
//   configureStrelkaWorkflow.pl --normal=${normal} --tumor=${tumor} --ref ${reffa} --config=./config.ini --output=strelka
//   make -j 32 -C strelka/
//   vcf-concat strelka/results/passed.somatic.*.vcf | vcf-sort | vcf-annotate -n --fill-type -n | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bedtools intersect -header -a stdin -b ${target_panel} |bgzip > ${tid}.strelka.vcf.gz
//   """
// }

process virmid {

  publishDir "$baseDir/output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from virmidbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.virmid.vcf.gz") into virmidvcf
  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 virmid/1.2 vcftools/0.1.11
  virmid -R ${reffa} -D ${tumor} -N ${normal} -s $cosmic t 30 -M 2000 -c1 10 -c2 10
  perl $baseDir/scripts/addgt_virmid.pl ${tumor}.virmid.som.passed.vcf
  perl $baseDir/scripts/addgt_virmid.pl ${tumor}.virmid.loh.passed.vcf
  vcf-concat *gt.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar \$SNPEFF_HOME/SnpSift.jar filter '((NDP >= 10) & (DDP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bedtools intersect -header -a stdin -b ${target_panel} |bgzip > ${tid}_${nid}.virmid.vcf.gz
  """
}

ssvcf .phase(mutectvcf)
      .map {p,q -> [p[0],p[1],q[1]]}
      .set { twovcf }
twovcf .phase(varscanvcf)
      .map {p,q -> [p[0],p[1],p[2],q[1]]}
      .set { threevcf }
threevcf .phase(virmidvcf)
      .map {p,q -> [p[0],p[1],p[2],p[3],q[1]]}
      .set { fourvcf }
fourvcf .phase(shimmervcf)
	.map {p,q -> [p[0],p[1],p[2],p[3],p[4],q[1]]}
      	.set { vcflist }

process integrate {
  publishDir "$baseDir/output", mode: 'copy'

  input:
  set fname,file(ss),file(mutect),file(vscan),file(virmid),file(shimmer) from vcflist
  
  output:
  set fname,file("${fname}.union.vcf.gz") into union
  script:
  """
  module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3
  module load vcftools/0.1.14
  perl $baseDir/scripts/unionize_vcf.pl -r ${index_path} ${ss} ${mutect} ${shimmer} ${vscan} ${virmid}
  sh integrate.sh
  perl $baseDir/scripts/uniform_integrated_vcf.pl ${fname}.temp.vcf
  bgzip ${fname}.union.vcf
  """
}

process annot {
  publishDir "$baseDir/output", mode: 'copy'

  input:
  set fname,unionvcf from union
  
  output:
  file("${fname}.annot.vcf.gz") into annot
  file("${fname}.stats.txt") into stats
  file("${fname}.statplot*") into plotstats

  script:
  if (params.genome == '/project/shared/bicf_workflow_ref/GRCh38')
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3
  tabix ${unionvcf}
  bcftools annotate -Oz -a ${index_path}/ExAC.vcf.gz -o ${fname}.exac.vcf.gz --columns CHROM,POS,AC_Het,AC_Hom,AC_Hemi,AC_Adj,AN_Adj,AC_POPMAX,AN_POPMAX,POPMAX ${unionvcf}
  tabix ${fname}.exac.vcf.gz 
  bcftools annotate -Oz -a ${index_path}/dbSnp.vcf.gz -o ${fname}.dbsnp.vcf.gz --columns CHROM,POS,ID,RS ${fname}.exac.vcf.gz
  tabix ${fname}.dbsnp.vcf.gz
  bcftools annotate -Oz -a ${index_path}/clinvar.vcf.gz -o ${fname}.clinvar.vcf.gz --columns CHROM,POS,CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNREVSTAT,CLNACC ${fname}.dbsnp.vcf.gz
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} ${fname}.clinvar.vcf.gz | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/cosmic.vcf.gz - |java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -v -db ${index_path}/dbNSFP.txt.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar gwasCat -db ${index_path}/gwas_catalog.tsv - |bgzip > ${fname}.annot.vcf.gz
  tabix ${fname}.annot.vcf.gz
  bcftools stats ${fname}.annot.vcf.gz > ${fname}.stats.txt
  plot-vcfstats -s -p ${fname}.statplot ${fname}.stats.txt
  """
  else
  """
  module load snpeff/4.2
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} ${unionvcf} |bgzip > ${fname}.annot.vcf.gz
  tabix ${fname}.annot.vcf.gz
  bcftools stats ${fname}.annot.vcf.gz > ${fname}.stats.txt
  plot-vcfstats -s -p ${fname}.statplot ${fname}.stats.txt
  """
}
