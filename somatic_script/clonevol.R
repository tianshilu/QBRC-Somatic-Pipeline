#######  set up  ############

# https://github.com/hdng/clonevol
# ClonEvol can plot both node-based tree (each clone is a node), or branch-based tree (each
# branch represents the evolution of a clone from its parental clone, and each node represents a
# point where the clone is established/founded)

library(clonevol)
args = commandArgs(trailingOnly=TRUE)
path=args[1] # /home2/twang6/software/cancer/somatic/example/evolution/results

#######3  read pyclone results  ############

cluster=read.table(paste(path,"/pyclone/tables/loci.tsv",sep=""),stringsAsFactors=F,sep="\t",header=T)
tmp=cluster[,c("mutation_id","cluster_id")]
tmp=tmp[!duplicated(tmp),]
rownames(tmp)=tmp$mutation_id

cluster_clonevol=matrix(NA,nrow=length(unique(cluster$mutation_id)),ncol=length(unique(cluster$sample_id)))
rownames(cluster_clonevol)=unique(cluster$mutation_id)
colnames(cluster_clonevol)=unique(cluster$sample_id)
for (i in 1:dim(cluster)[1])
  {cluster_clonevol[cluster$mutation_id[i],cluster$sample_id[i]]=cluster$cellular_prevalence[i]*100}
cluster_clonevol=data.frame(cluster=tmp[rownames(cluster_clonevol),"cluster_id"],
  cluster_clonevol,stringsAsFactors=F)
cluster_clonevol=cluster_clonevol[order(cluster_clonevol$cluster),]

##########  clonevol  ################

y = infer.clonal.models(variants = cluster_clonevol, cluster.col.name = 'cluster',
  ccf.col.names = colnames( cluster_clonevol)[-1],cancer.initiation.model='polyclonal',
  subclonal.test = 'bootstrap',subclonal.test.model = 'non-parametric',
  num.boots = 1000,cluster.center = 'mean',ignore.clusters = NULL,
  clone.colors = NULL,min.cluster.vaf = 0.01,sum.p = 0.05,alpha = 0.05,random.seed=1)

y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

setwd(path)
plot.clonal.models(y,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = 'output',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 1.5*(dim(cluster_clonevol)[2]-1),
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,2))
