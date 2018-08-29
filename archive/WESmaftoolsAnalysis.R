suppressPackageStartupMessages(require(maftools))
#read maf file for LAML
laml.maf="WESALL5callerQualityFilteredAF5.txt"
laml = read.maf(maf = laml.maf, useAll = TRUE)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)
write.mafSummary(maf = laml, basename = 'laml')
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,top=20, titvRaw = FALSE)
oncoplot(maf = laml, top = 50, removeNonMutated = TRUE)

#plot CNV data for one example AC2-200N
seg<-"scg2.txt"
pdf(file="WES59segPlotAC2New.pdf", width=8, height=1.2)
plotCBSsegments(cbsFile = seg, maf = laml, labelAll = TRUE)
dev.off()

#Changing colors for oncoplot (You can use any colors, here in this example we will use a color palette from RColorBrewer)
col = c("orange","red","purple","darkgreen","black","blue","yellow","cyan")
names(col) = c('In_Frame_Del','Nonsense_Mutation', 'In_Frame_Ins', 'Missense_Mutation', 'Multi_Hit', "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site")
oncoplot(maf = laml, top = 50, removeNonMutated = TRUE, colors=col,showTumorSampleBarcodes = TRUE)
oncoplot(maf = laml, top = 25, removeNonMutated = TRUE, colors=col)


#We will plot same top 30 mutated genes with fibrosis classification as annotation and using above defined colors.
laml.fab.anno = "Annotation.txt"
laml.fab.anno = read.delim(laml.fab.anno, sep = '\t')
head(laml.fab.anno)
oncoplot(maf = laml, top = 30, annotation = laml.fab.anno, removeNonMutated = TRUE, colors = col)

# add annotation
oncoplot(maf = laml, top = 10, annotation = laml.fab.anno, removeNonMutated = TRUE)
#plot for oncoscript for MLL pathway
pdf(file="oncostrip.pdf", width=5, height=3)
oncostrip(maf = laml, genes = c('KMT2B','KMT2C', 'KMT2D'), removeNonMutated = TRUE, showTumorSampleBarcodes = FALSE)
dev.off()

# rainfall plots
pdf(file="rainfall.pdf", width=5, height=3)
rf = rainfallPlot(maf = laml, detectChangePoints = TRUE, fontSize = 12, pointSize = 0.6)
dev.off()


#analyze transition transversion
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
pdf(file="TiTv.pdf", width=10, height=10)
plotTiTv(res = laml.titv)
dev.off()

#Lets plot lollipop plot for KMT2D
pdf(file="KMT2D.pdf", width=10, height=3)
KMT2D.lpop = lollipopPlot(maf = laml, gene = 'KMT2D', AACol = 'HGVSp')
dev.off()

#plot for mutation load comparison
pdf(file="ALLcompare.pdf", width=10, height=10)
laml.mutload = tcgaCompare(maf = laml, cohortName = 'Cirrhosis')
dev.off()

#plot for gene cloud
pdf(file="ALLGeneCloud.pdf", width=10, height=10)
geneCloud(input = laml, minMut = 5)
dev.off()

#plot for mutation allele frequency distribution
pdf(file="Tvaf.pdf", width=10, height=10)
vafPlot = plotVaf(maf = laml, vafCol = 't_AF', top=20,flip = TRUE)
dev.off()

#identify significantly driver genes by maftools
laml.sig = oncodrive(maf = laml, AACol = 'HGVSp', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)
write.csv(laml.sig,"DriverGenesMaftools5callerWED.csv",row.names=FALSE)
pdf(file="driver.pdf", width=5, height=5)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
dev.off()

#correct gene symbol for mutsigCV
laml.mutsig.corrected = prepareMutSig(maf = laml)
write.csv(laml.mutsig.corrected, "WESAllAF5forMutsig.corrected.csv",row.names=FALSE)

#plot for pfam protein domain
laml.pfam = pfamDomains(maf = laml, AACol = 'HGVSp', top = 10)
laml.pfam$proteinSummary[,1:7, with = FALSE]

##comparing two cohorts (MAFs)
cir.maf="WESALL5callerQualityFilteredAF5.txt"
cir= read.maf(maf = cir.maf,  useAll = TRUE)
LIHC.maf="data_mutations_extended.txt"
LIHC = read.maf(maf = LIHC.maf,  useAll = TRUE)
pt.vs.rt <- mafCompare(m1 = cir, m2 = LIHC, m1Name = 'Cirrhosis', m2Name = 'LIHC', minMut = 5)
pt.vs.rt1 <- mafCompare(m1 = LIHC, m2 = cir, m1Name = 'Cirrhosis', m2Name = 'LIHC', minMut = 5)

print(pt.vs.rt)
#write.csv(pt.vs.rt$results,"CirrhosisvsLIHC.csv",row.names=FALSE)
apl.pt.vs.rt.fp = forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, show = 'stat', color = c('royalblue', 'maroon'))

genes = c("TP53","CTNNB1","ALB","ARID1A","AXIN1","KMT2D","KMT2C","RB1","KEAP1")
coOncoplot(m1 = LIHC, m2 = cir, m1Name = 'HCC', m2Name = 'Cirrhosis', genes = genes, removeNonMutated = TRUE)

genes = c("ADPRHL1","MADCAM1","ZNF814","CGREF1","MAGEC1", "KMT2D" ,"AIM1L", "KRTAP5-10", "LOR","NEB","SHANK1","SLC9B1")
coOncoplot(m1 = cir, m2 = LIHC, m1Name = 'cirrhosis', m2Name = 'HCC', genes = genes, removeNonMutated = TRUE)


#check sample heterogeneity
HS37N.het = inferHeterogeneity(maf = laml, tsb = 'HS37N')
print(HS37N.het$clusterMeans)
plotClusters(clusters = HS37N.het)

#somatic signature analysis
laml.tnm = trinucleotideMatrix(maf = laml, ref_genome ='C:/Users/xluo4/Desktop/Projects/Zhu Hao Lab/genome.fa',ignoreChr=NULL, useSyn = TRUE, fn=NULL)
#check Apobec levels
plotApobecDiff(tnm = laml.tnm, maf = laml)
# Create data for the graph.
x <- c(5,81)
labels <- c("APOBEC enriched", "non-APOBEC enriched")
# Give the chart file a name.
png(file = "apobec.jpg")
# Plot the chart.
pie(x,labels,col=rainbow(2),radius=0.5)
library(plotrix)
pie3D(x,labels=labels,explode = 0.1)
# Save the file.
dev.off()

require('NMF')
laml.sign = extractSignatures(mat = laml.tnm, n=4,nTry = 10, plotBestFitRes = FALSE)
A<-laml.sign$contributions
write.csv(A,"signature4.csv",row.names=TRUE)

plotSignatures(laml.sign)
require('corrplot')
corrplot::corrplot(corr = laml.sign$coSineSimMat, col = RColorBrewer::brewer.pal(n = 9, name = 'Oranges'), is.corr = FALSE, tl.cex = 0.6, tl.col = 'black', cl.cex = 0.6)
