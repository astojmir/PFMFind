#psfilename <- "test2.ps"
#file <- "tt2.txt"

readhist <- function(file, skip0)
  {
    title0 <- scan(file=file, what = list(x=""), nlines=1, skip=skip0, 
                  flush=TRUE, quiet = TRUE)
    points <- scan(file=file, what = list(x="", y=0, z=0), nlines=1, 
                   skip=skip0+1, flush=TRUE, quiet = TRUE)	
    values <- scan(file=file, what = list(x="", y=0), nlines=1, 
                   skip=skip0+2, flush=TRUE, quiet = TRUE)
    header <- scan(file=file, what = list(x="", y=""), nlines=1, 
                   skip=skip0+3, flush=TRUE, quiet = TRUE)
    hist0 <- scan(file=file, what = list(val=0, freq=0), nlines=values$y, 
                   skip=skip0+4, flush=TRUE, quiet = TRUE)
    counts0 <- hist0$freq
    breaks0 <- c(hist0$val[1] - 0.5, hist0$val+0.5)
    intensities0 <- counts0/points$y 	
    mids0 <- hist0$val

    hist1<- list(counts=counts0, breaks=breaks0, intensities=intensities0, 
	         mids=mids0, xname="Value", equidist=TRUE)
    title0 <- paste(title0, " \n\n Points:", points$y)
    plot.histogram(hist1, ylab="Frequency", main=title0)
    skip1 <- skip0 + 4 + values$y + 1	
    hlist <- list(mids=mids0, counts=counts0, skip=skip1)
  }
size <-  scan(file=file, what = list(x="", y="", z="", w=0), nlines=1, 
	skip=0, flush=TRUE, quiet = TRUE)
pages <- size$w

main0 <- paste("Distribution of distances from random fragments")

postscript(psfilename, horizontal=FALSE, onefile=TRUE, pointsize=10)
par(mfrow=c(4,1), oma=c(5,2,10,2))


skip0 <- 1

for (i in 1:pages)
  {
    rand.pt <- scan(file=file, what = list(x="", y="", z=""), 
	nlines=1, skip=skip0, flush=TRUE, quiet = TRUE)
    skip0 <- skip0+1
    main <- paste(main0, "\n\nFragment: ", rand.pt$z)

    hlist<- readhist(file, skip0)
    skip0 <- hlist$skip
    qdl.x <- hlist$mids
    qdl.y <- hlist$counts
    hlist<- readhist(file, skip0)
    skip0 <- hlist$skip
    qdr.x <- hlist$mids
    qdr.y <- hlist$counts
    hlist<- readhist(file, skip0)
    skip0 <- hlist$skip
    d.x <- hlist$mids
    d.y <- hlist$counts

plot(c(0,110), c(0,400), type="n", ylab="Frequency", xlab="Distance")
lines(d.x, d.y, type="o", col=1, pch=1)
lines(qdr.x, qdr.y, type="o", col=2, pch=5)
lines(qdl.x, qdl.y, type="o", col=4, pch=4)
legend(0, 300, c("Metric", "Left quasi-metric", "Right quasi-metric"), col = c(1,4,2), pch=c(1,4,5), lty=c(1,1,1))
    mtext(main, line=1, outer=TRUE)
  }

dev.off()
