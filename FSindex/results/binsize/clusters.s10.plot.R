psfilename <- "test3.ps"
file <- "sprot_10_4_s10_F.txt"
Tfile <- "sprot_10_4_s10_F_T.txt"

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
    library(stepfun)
    title0 <- "Bin Size Distribution"
    A = read.table(Tfile, header = TRUE, skip = 15, nrows=99)    
#    plot(stepfun(A[0:30,1], c(0, A[0:30,4])), col.hor="red", col.vert="red", do.points=FALSE)
    plot(stepfun(A[1:30,1], c(0, A[1:30,2])), col.hor="blue", col.vert="blue", col.points="blue", pch="*")
    plot(stepfun(mids0[1:30], c(0, intensities0[1:30])), add=TRUE)



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

    title0 <- " "
    plot(mids0, intensities0, type="l", log="y") 
    skip1 <- skip0 + 4 + values$y + 1 
  }


total <-  scan(file=file, what = list(x1="", x2="", x3="", x4="", x5="", y=0),
               nlines=1, skip=16, flush=TRUE, quiet = TRUE)
msize <-  scan(file=file, what = list(x1="", x2="", x3="", x4="", x5="", y=0),
               nlines=1, skip=18, flush=TRUE, quiet = TRUE)

postscript(psfilename, horizontal=FALSE, onefile=TRUE, pointsize=10)
par(mfrow=c(3,2), oma=c(5,2,10,2))


skip0 <- 21
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
