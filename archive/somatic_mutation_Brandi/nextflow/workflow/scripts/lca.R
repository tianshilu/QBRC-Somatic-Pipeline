#!/usr/bin/Rscript

library(R2jags)
library(getopt)
library(reshape2)

discard = runif(1) # fixes this error: 'Error in set.seed() : argument "seed" is missing, with no default'

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
'verbose', 'v', 2, "logical",
'countsFile', 'c', 1, "character",
'statsOutFile', 's', 1, "character",
'help' , 'h', 0, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

if ( is.null(opt$countsFile ) ) { 
  cat("please specify a --countsFile\n")
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

if ( is.null(opt$statsOutFile ) ) { 
  cat("please specify a --statsOutFile\n")
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

########################################
########################################

# This function generates a JAGS file, to be written out and used by R2Jags
# The JAGS file changes depending on the number of SNP calls
generateJagsFile <- function(numberOfSnpCallSets, outJagsFile) {

	 # numberOfSnpCallSets is the number of VCF files the user gives us

jagCode = 
"model {

  theta ~ dbeta(1,2);

  for (i in 1:N) {
    fp[i] ~ dbeta(1,2);
    fn[i] ~ dbeta(1,2);
    tp[i] <- 1-fn[i];
    tn[i] <- 1-fp[i];
  }

";

  jagCode = paste(jagCode, "for (j in 1:M) {\n\tmargprobs[j] <- theta * \n" )
  for (i in 1:numberOfSnpCallSets){ 
     if ( i == numberOfSnpCallSets ){
       jagCode = paste(jagCode, "\t(fn[", i, "]^(1-x[j,", i, "])) * (tp[", i, "]^x[j,", i, "]) \n", sep="") 
     }
     else {
       jagCode = paste(jagCode, "\t(fn[", i, "]^(1-x[j,", i, "])) * (tp[", i, "]^x[j,", i, "]) *\n", sep="") 
     }
  }
  jagCode = paste(jagCode, "\t+ (1-theta) * \n")
  for (i in 1:numberOfSnpCallSets){ 
     if ( i == numberOfSnpCallSets ){
          jagCode = paste(jagCode, "\t(tn[", i, "]^(1-x[j,", i, "])) * (fp[", i, "]^x[j,", i, "]);\n\n", sep="")
     }
     else {
          jagCode = paste(jagCode, "\t(tn[", i, "]^(1-x[j,", i, "])) * (fp[", i, "]^x[j,", i, "]) *\n", sep="")
     }
  }

  jagCode = paste(jagCode, "\tpostprobs[j] <- theta * \n", sep="")

  for (i in 1:numberOfSnpCallSets){ 
     if ( i == numberOfSnpCallSets ){
          jagCode = paste(jagCode, "\t(fn[", i, "]^(1-x[j,", i, "])) * (tp[", i, "]^x[j,", i, "])\n\t/ margprobs[j];\n", sep="")
     }
     else {
          jagCode = paste(jagCode, "\t(fn[", i, "]^(1-x[j,", i, "])) * (tp[", i, "]^x[j,", i, "]) *\n", sep="")
     }
  }
  jagCode = paste(jagCode, "}

counts ~ dmulti(margprobs, total);
}", sep="")

  fileConn<-file(outJagsFile)
  writeLines(jagCode, fileConn)
  close(fileConn)

}

########################################
########################################

# read in counts file
raw.data = read.table(opt$countsFile, header=F, skip=1 )

# counts <- c(2942808473, 17491655, 21576, 23189, 339805, 89159,
#             168214, 76044, 43138288, 530963, 22682, 22169,
#             462052, 129472, 2804257, 3454104);

counts = raw.data[,length(raw.data)]
M <- length(counts)
N <- log(M, base=2)
total <- sum(counts);

# make matrix of all possible 0/1 values for each SNP caller
cmdstring =  paste( "expand.grid(", paste( rep("c(0,1)", N), collapse = ", " ), ")", sep="")
x <- eval( parse(text = cmdstring ) )

jags.data <- c("N", "M", "x", "counts", "total");
jags.params <- c("theta", "fp", "fn", "margprobs", "postprobs");

jags.inits <- NULL;

# write out lca.jags file 
tempJagFile = tempfile();
generateJagsFile(N, tempJagFile);

jagsfit1 <- jags(data=jags.data, inits=jags.inits, param=jags.params, DIC=FALSE,
                                n.chains=1, n.iter=120000, n.thin=10, n.burnin=10000,
                                model.file=tempJagFile)

mcmc <- as.mcmc(jagsfit1);

# below is what Aaron was doing before, but I'm just going to write out the whole summary statistics file
# print(params <- summary(mcmc)$statistics[c(1:8,41),1, drop=F])
# print(pp <- summary(mcmc)$statistics[25:40,1, drop=F])
# print(cbind(x, pp, pp >= 0.5))

# write out summary statistics
summaryStats = summary(mcmc)$statistics

idx <- 2*N + 1 # skip the fp/fn parameters
data <- cbind(x, cts = counts, margprobs = summary(mcmc)$statistics[idx:(idx+M-1), 1])

bvrs <- matrix(NA, ncol=N, nrow=N);

for (i in 2:N) {
  for (j in 1:(i-1)) {
    obs <- dcast(data[,i]+data[,j] ~ ., sum, value.var="cts", data=data)[,3]
    exp <- dcast(data[,i]+data[,j] ~ ., sum, value.var="margprobs", data=data)[,3]
    bvrs[i,j] <- chisq.test(obs, p=exp, rescale.p=T)$statistic / 3 # BVR == Chi-square/d.f.

  }
}

sink(opt$statsOutFile)
cat("==== Summary Statistics ====\n")
summaryStats
cat("==== Bivariate residuals ====\n")
bvrs
sink()
