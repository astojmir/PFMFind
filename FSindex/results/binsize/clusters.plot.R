psfilename <- "test2.ps"
file <- "sprot_8_6.txt"
Tfile <- "sprot_8_6_T.txt"

plottheoretical <- function(file)
  {
    A = read.table(file, header = TRUE, skip = 15, nrows=99)    
    plot(A[1:100,1], A[1:100,2], log="y", type="l", xlab="Bin Size",
         ylab="Probability")
    plot(A[1:100,1], A[1:100,4], add=TRUE, log="y", type="l", xlab="Bin Size",
         ylab="Probability")
    title("Theoretical Bin Size Distribution")
  }



readhist1 <- function(file, Tfile, skip0, msize, total)
  {
    par(mfrow=c(1,1), oma=c(5,2,10,2))
    size <- 0
    i <- 1
    breaks0 <- -0.5
    counts0 <- 1
    intensities0 <- 1
    mids0 <- 1
    totfreq <- 0
    
    while (size < msize)
      {
        hist0 <- scan(file=file, what = list(val=0, freq=0), nlines=1, 
                   skip=skip0+1+i, flush=TRUE, quiet = TRUE)        
        size <- hist0$val
        totfreq <- hist0$freq + totfreq
        
        if (totfreq/total < 0.999999)
          {          
            counts0 <- c(counts0, hist0$freq)
            breaks0 <- c(breaks0, size + 0.5)
            intensities0 <- c(intensities0,
                              ((hist0$freq/total)/(breaks0[i+1]-breaks0[i])))
            mids0 <- c(mids0, (breaks0[i+1]+breaks0[i])/2)
          }
        i <- i + 1
      }
    counts0 <- counts0[2:length(counts0)]
    intensities0 <- intensities0[2:length(intensities0)]
    mids0 <- mids0[2:length(mids0)]

    title0 <- "Bin Size Distribution"
    A = read.table(Tfile, header = TRUE, skip = 15, nrows=99)    
    plot(A[1:60,1], A[1:60,4], col="red", type="l", xlab="Bin Size",
         ylab="Relative Frequency")
    lines(mids0[1:60], intensities0[1:60])
    lines(A[1:60,1], A[1:60,2], col="blue")


    title(title0)
    skip1 <- skip0+3+i 
  }



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


    size <- 0
    breaks0 <- -0.5
    counts0 <- 1
    intensities0 <- 1
    mids0 <- 1
    totfreq <- 0
    
    for (i in 1:values$y)  
      {
        hist0 <- scan(file=file, what = list(val=0, freq=0), nlines=1, 
                   skip=skip0+3+i, flush=TRUE, quiet = TRUE)        
        size <- hist0$val
        totfreq <- hist0$freq + totfreq
        
        if (totfreq/points$y < 0.999999)
          {
            if (hist0$freq > 0)
              {
                counts0 <- c(counts0, hist0$freq)
                breaks0 <- c(breaks0, size + 0.5)
                intensities0 <- c(intensities0,
                                  (hist0$freq/(breaks0[i+1]-breaks0[i])))
                mids0 <- c(mids0, (breaks0[i+1]+breaks0[i])/2)
              }
          }
      }
    counts0 <- counts0[2:length(counts0)]
    intensities0 <- intensities0[2:length(intensities0)]
    mids0 <- mids0[2:length(mids0)]

#    hist1<- list(counts=counts0, breaks=breaks0, intensities=intensities0, 
#	         mids=mids0, xname="Bin Size", equidist=FALSE)
    title0 <- " "
    plot(mids0, intensities0, type="l", log="y") 
#    plot.histogram(hist1, ylab="Frequency", main=title0)
    skip1 <- skip0 + 4 + values$y + 1 


    
    
    
#    hist0 <- scan(file=file, what = list(val=0, freq=0), nlines=values$y, 
#                   skip=skip0+4, flush=TRUE, quiet = TRUE)
#    counts0 <- hist0$freq
#    breaks0 <- c(hist0$val[1] - 0.5, hist0$val+0.5)
#    intensities0 <- counts0/points$y 	
#    mids0 <- hist0$val
#
#    hist1<- list(counts=counts0, breaks=breaks0, intensities=intensities0, 
#	         mids=mids0, xname="Value", equidist=TRUE)
#    title0 <- paste(title0, " \n\n Points:", points$y,"/", points$z)
#    title0 <- " "
#    plot.histogram(hist1, ylab="Frequency", main=title0)
#    skip1 <- skip0 + 4 + values$y + 1	
  }

#db <- scan(file=file, what = list(x="", y=""), nlines=1, skip=2, 
#                  flush=TRUE, quiet = TRUE)
#matrix <- scan(file=file, what = list(x="", y=""), nlines=1, skip=2, 
#                  flush=TRUE, quiet = TRUE)
#alpha <- scan(file=file, what = list(x="", y=0), nlines=1, skip=3, 
#                  flush=TRUE, quiet = TRUE)
#beta <- scan(file=file, what = list(x="", y=0), nlines=1, skip=4, 
#                  flush=TRUE, quiet = TRUE)
#
#size <-  scan(file=file, what = list(x="", y=0, z=0), nlines=1, skip=5, 
#                  flush=TRUE, quiet = TRUE)
#pages <- size$y*size$z
#
#main0 <- paste("BLAST Seed Statistics - ", "Database: ", db$y,
#              "Scoring: ", matrix$y, "(", alpha$y, ",",
#              beta$y, ")\n\n")

total <-  scan(file=file, what = list(x1="", x2="", x3="", x4="", x5="", y=0),
               nlines=1, skip=15, flush=TRUE, quiet = TRUE)
msize <-  scan(file=file, what = list(x1="", x2="", x3="", x4="", x5="", y=0),
               nlines=1, skip=17, flush=TRUE, quiet = TRUE)

postscript(psfilename, horizontal=FALSE, onefile=TRUE, pointsize=10)
par(mfrow=c(3,2), oma=c(5,2,10,2))


skip0 <- 20
skip0 <- readhist1(file, Tfile, skip0, msize$y, total$y)
plottheoretical(Tfile)
skip0 <- skip0 + 4
#skip0 <- readhist(file, skip0)
#skip0 <- readhist(file, skip0)
#skip0 <- readhist(file, skip0)
#skip0 <- readhist(file, skip0)


#for (i in 1:pages)
#  {
#    expect <- scan(file=file, what = list(x="", y=0), nlines=1, skip=skip0, 
#                  flush=TRUE, quiet = TRUE)
#    skip0 <- skip0+1
#    length <- scan(file=file, what = list(x="", y=0), nlines=1, skip=skip0, 
#                  flush=TRUE, quiet = TRUE)
#    skip0 <- skip0+2
#    main <- paste(main0, "Expect: ", expect$y, " Length: ", length$y)
#
#    for (i in 1:8)
#      {
#    	skip0 <- readhist(file, skip0)
#      }
#    mtext(main, line=1, outer=TRUE)
#  }

dev.off()
