

readhist1 <- function(file, Tfile, title0, range)
  {
    total1 <-  scan(file=file, what = list(x1="", x2="", x3="", x4="",
                   x5="", y=0), nlines=1, skip=16,
                   flush=TRUE, quiet = TRUE) 
    msize1 <-  scan(file=file,
                   what = list(x1="", x2="", x3="", x4="", x5="", y=0),
                   nlines=1, skip=18, flush=TRUE, quiet = TRUE)
    total <- total1$y
    msize <- msize1$y
    skip0 <- 21    
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

    A = read.table(Tfile, header = TRUE, skip = 15, nrows=range)    

    MA = max(A[1:range,2])
    Mm = max(intensities0[1:range])

    if (MA > Mm)
      {    
        plot(A[1:range,1], A[1:range,2], col="red", type="p", xlab="Bin Size",
             ylab="Relative Frequency")
        points(mids0[1:range], intensities0[1:range])
      }
    else
      {
        plot(mids0[1:range], intensities0[1:range], type="p", xlab="Bin Size",
             ylab="Relative Frequency")
        points(A[1:range,1], A[1:range,2], col="red")
      }
    
    title(title0)
    skip1 <- skip0+3+i 
  }





psfilename <- "test2.ps"
postscript(psfilename, horizontal=FALSE, onefile=TRUE, pointsize=10)
par(mfrow=c(3,2), oma=c(5,2,10,2))

file <- "sprot_8_6.txt"
Tfile <- "sprot_8_6_T.txt"

skip0 <- 20
skip0 <- readhist1("sprot_8_6.txt", "sprot_8_6_T.txt", "SwissProt, m=8, K=6", 100)
skip0 <- readhist1("sprot_10_4.txt", "sprot_10_4_T.txt", "SwissProt, m=10, K=4", 100)
skip0 <- readhist1("sprot_10_5.txt", "sprot_10_5_T.txt", "SwissProt, m=10, K=5", 99)
skip0 <- readhist1("sprot_10_4_s10.txt", "sprot_10_4_s10_T.txt", "SwissProt, m=10, K=4, skip=10", 40)
skip0 <- readhist1("sprot_10_4_s10_F.txt", "sprot_10_4_s10_F_T.txt", "SwissProt_F, m=10, K=4, skip=10", 40)
skip0 <- readhist1("nr_10_5.txt", "nr_10_5_T.txt", "nr, m=10, K=5", 199)



dev.off()

