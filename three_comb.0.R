require("plotrix")	# for multhist
require(ggplot2)	# for ggplot
require(reshape2)	# for melt
require(Rcmdr)		# for Hist

calc.intervals <- function(var, len) {
    # compute intervals of width len for an array of values
    m <- min(var)
    if (m > 0) m <- 0
    M <- max(var) + len - 1
    if (m == M) {
        return(c(m,M))
    }
    inter <- seq(from=0, to=M, by=len)
    if (m < 0) inter <- c(m, inter)
    return(inter)
}


for (pfu in c('0.1', '001', '010', '100')) {
    print(paste("PFU =", pfu))
    basename <- paste('three_matches/', pfu, 'pfu_exp1-3_R1-2_001', sep='')
    
    # create an empty data frame with the corresponding columns
    # to hold all the data for this growth level
    distances <- data.frame(NAM=character(),
                            SEG=character(),
                            SIZ1=integer(),
                            SIZ2=integer(),
                            LEFT1=integer(),
                            LEFT2=integer(),
                            RIGHT1=integer(),
                            RIGHT2=integer(),
                            REAR=character(),
                            stringsAsFactors=FALSE)
    
    for (exp in c('1', '2', '3')) {
        print(paste('    EXP =', exp))
        for (dir in c('1', '2')) {
	    print(paste('        R =', dir))
	    file <- paste('three_matches/', 
	                  pfu, 'pfu_exp', exp, '_R', dir, '_001.csv',
	                  sep='')
	    data <- read.table(file, header=T, sep='\t', stringsAsFactors=T)
	    # add this data to 'distances'
	    distances <- rbind(distances, data)
	}
    }
        	
    # save sorted distances in a dist file
    distfile <- paste(basename, '.csv', sep='')
    print(paste("    writing", distfile))
    sorted.distances <- distances[order(distances$SEG, distances$SIZ1),]
    #sorted.distances <- distances
    write.table(sorted.distances, file=distfile, 
                row.names=TRUE, col.names=TRUE,
                sep='\t', quote=FALSE)

    distances <- sorted.distances
    
    # remove anomalies
    distances <- subset(distances, REAR == 'N')

    # Now make all summary plots
    
    # interval length
    int.len = 25
    
    # split by segment
    dist.seg.A <- subset(distances, SEG == "A" )
    dist.seg.B <- subset(distances, SEG == "B" )

    # Plot the distribution of gap sizes for both segments in one figure
    # join all (SIZ1 and SIZ2) in one vector
    size <- rbind(distances$SIZ1, distances$SIZ2)
    size.m <- 0
    size.M <- max(size[ ! is.na(size) ])
    size.M <- max(size.M, size.M + int.len - 1)
    print(paste("    MAX SIZ A+B =", size.M))
    # separate by segment
    size.A <- rbind(dist.seg.A$SIZ1, dist.seg.A$SIZ2)
    size.B <- rbind(dist.seg.B$SIZ1, dist.seg.B$SIZ2)

    size.hist <- paste(basename, '.hist', int.len, 'AB.png', sep='')
    png(size.hist, width=1000, height=1000)
#	hist( size.A, 
#            breaks = calc.intervals(size, int.len),
#            border = "blue", 
#            main = basename,
#            sub = "Segment A")
#	hist( size.B, 
#            breaks = calc.intervals(size, int.len),
#            border = "red", 
#            add=T)
#	axis(1, at=int.len, labels=int.len)
        par(mfrow=c(2,1))
        hist(size.A, breaks = calc.intervals(size, int.len), border="blue", sub="Segment A")
	hist(size.B, breaks = calc.intervals(size, int.len), border="red", sub="Segment B")
        par(mfrow=c(1,1))
    dev.off()   

    # Plot the distribution of start positions
    # join all (SIZ1 and SIZ2) in one vector
    start <- rbind(distances$LEFT1, distances$LEFT2)
    start.m <- 0
    start.M <- max(start[ ! is.na(start) ])
    start.M <- max(start.M, start.M + int.len - 1)
    print(paste("    MAX LEFT A+B =", start.M))
    # separate by segment
    start.A <- rbind(dist.seg.A$LEFT1, dist.seg.A$LEFT2)
    start.B <- rbind(dist.seg.B$LEFT1, dist.seg.B$LEFT2)

    start.hist <- paste(basename, '.start', int.len, 'AB.png', sep='')
    png(start.hist, width=1000, height=1000)
#          multhist(list(start.A, start.B),
#              breaks = calc.intervals(start, int.len),
#              main = basename
#              )
        par(mfrow=c(2,1))
        hist(start.A, breaks = calc.intervals(start, int.len), border="blue", sub="Segment A")
	hist(start.B, breaks = calc.intervals(start, int.len), border="red", sub="Segment B")
        par(mfrow=c(1,1))

    dev.off()


    # Plot the distribution of end positions
    # join all (SIZ1 and SIZ2) in one vector
    end <- rbind(distances$RIGHT1, distances$RIGHT2)
    end.m <- 0
    end.M <- max(end[ ! is.na(end) ])
    end.M <- max(end.M, end.M + int.len - 1)
    print(paste("    MAX RIGHT A+B =", end.M))
    # separate by segment
    end.A <- rbind(dist.seg.A$RIGHT1, dist.seg.A$RIGHT2)
    end.B <- rbind(dist.seg.B$RIGHT1, dist.seg.B$RIGHT2)

    end.hist <- paste(basename, '.end', int.len, 'AB.png', sep='')
    png(end.hist, width=1000, height=1000)
#          multhist(list(end.A, end.B),
#              breaks = calc.intervals(end, int.len),
#              main = basename
#              )
        par(mfrow=c(2,1))
        hist(end.A, breaks = calc.intervals(end, int.len), border="blue", sub="Segment A")
	hist(end.B, breaks = calc.intervals(end, int.len), border="red", sub="Segment B")
        par(mfrow=c(1,1))

    dev.off()
}
