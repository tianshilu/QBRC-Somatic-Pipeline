myplot <- function() {
  colorLst <- c("red", "blue", "green", "cyan", "gold", "orchid", "plum","black")
  mut_type <- c("Exonic nonsynonymous SNV", "Exonic frameshift substitution", "Exonic stopgain", "Exonic nonframeshift substitution", "Splicing", "Exonic unknown", "Exonic stoploss","Multi-hit")
  
  plot(0, 9.5, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[1])
  text(0.1, 9.5, pos=4, mut_type[1]) 
  
  par(new=T)
  plot(0, 8, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[2])
  text(0.1, 8, pos=4, mut_type[2]) 
  
  par(new=T)
  plot(0, 6.5, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[3])
  text(0.1, 6.5, pos=4, mut_type[3]) 
  
  par(new=T)
  plot(0, 5, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[4])
  text(0.1, 5, pos=4, mut_type[4]) 
  
  par(new=T)
  plot(5, 9.5, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[5])
  text(5.3, 9.5, pos=4, mut_type[5]) 
  
  par(new=T)
  plot(5, 8, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[6])
  text(5.3, 8, pos=4, mut_type[6]) 
  
  par(new=T)
  plot(5, 6.5, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[7])
  text(5.3, 6.5, pos=4, mut_type[7]) 
  
  par(new=T)
  plot(5, 5, cex=3, xlim=c(0, 10), ylim=c(3, 10), axes=F,
       main="", xlab="", ylab="", pch=15, col=colorLst[8])
  text(5.3, 5, pos=4, mut_type[8]) 
}
