##' @export
plot.eqtl <- function(marker.info, probepos, phenonames,
                      markers, chr1, pos1,
                      main="",
                      plot.chr=c(1:19, "X"), add=FALSE, col="black"){

  chr.info <- get.chr.info(marker.info)
  minpos <- chr.info$start
  maxpos <- chr.info$end
  names(minpos) <- names(maxpos) <- chr.info$chr
  L <- maxpos - minpos
  csL <- cumsum(c(0, L))

  chrs <- c(1:19, "X")
  csL.move <- csL[-21]
  names(csL.move) <- chrs

  chr <- probepos[phenonames, "chr"]
  pos <- probepos[phenonames, "cM"]
  names(pos) <- phenonames
  cspos <- pos
  
  for(i.chr in c(1:19, "X")){
    ind <- which(chr==i.chr)
    cspos[ind] <- pos[ind] - minpos[i.chr] + csL.move[i.chr]
  }
  
  ## markers could be specified by mname or chr & pos.
  if(missing(chr1)){
    maxlod.cspos <- marker.info[markers, "pos"]
    maxlod.chr <- marker.info[markers, "chr"]
  }else if(missing(markers)){
    maxlod.cspos <- pos1
    maxlod.chr <- chr1
  }
  
  for(i.chr in chrs){
    ind <- which(maxlod.chr==i.chr)
    maxlod.cspos[ind] <- maxlod.cspos[ind] - minpos[i.chr] + csL.move[i.chr]
  }

  ## ##############################################################
  bgcolor <- "gray85"
  cex.lab <- 1.2
  col.lab <- "blue"

  if(!add){
    par(mar=c(2.6, 2.6, 1.6, 1.6), pty="s")
    par(col.lab=col.lab, cex.lab=cex.lab)

    plot(0, 0, type="n", xlim=c(-1, max(csL)+1), ylim=c(-1, max(csL)+1), 
         xaxt="n", yaxt="n", xlab="QTL position (cM)", 
         ylab="Probe position (cM)", main=main,
         xaxs="i", yaxs="i", mgp=c(1.1, 0, 0))
    
    for(i in 1:20)
        for(j in 1:20)
            rect(csL[i], csL[j], csL[i+1], csL[j+1], border=NA, 
                 col=c("gray80", "gray100")[((i + j) %% 2)+1])

    u <- par("usr")
    rect(u[1], u[3], u[2], u[4])

    midchr <- (csL[-1] + csL[-21])/2
    y <- u[3] - 0.3/par("mar")[3]*par("mai")[3]/par("fin")[2]*diff(u[3:4])
    x <- u[1] - 0.2/par("mar")[1]*par("mai")[1]/par("fin")[1]*diff(u[1:2])
    text(midchr, rep(y, 20), labels=c(1:19, "X"), xpd=TRUE, adj=c(0.5, 1), cex=0.9)
    text(rep(x, 20), midchr, labels=c(1:19, "X"), xpd=TRUE, adj=c(1, 0.5), cex=0.9)
  }
  
  for(i in plot.chr){
    ind <- which(maxlod.chr==i)
    points(maxlod.cspos[ind], cspos[phenonames[ind]], cex=0.5, pch=16, col=col)
  }

}
