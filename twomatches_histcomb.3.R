require("plotrix")	# for multhist
require(ggplot2)	# for ggplot
require(reshape2)	# for melt

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


for (pfu in c("0.1", "001", "010", "100")) {
    print(paste("PFU =", pfu))

    for (dir in c(1, 2)) {
        print(paste("    DIR =", dir))
        basename <- paste("two_matches/", pfu, "pfu_exp1-3_R", dir, sep="")
        exp1 <- paste("two_matches/", pfu, "pfu_exp1_R", dir, "_001.csv", sep="")
        exp2 <- paste("two_matches/", pfu, "pfu_exp2_R", dir, "_001.csv", sep="")
        exp3 <- paste("two_matches/", pfu, "pfu_exp3_R", dir, "_001.csv", sep="")
        
        d1 <- read.table(exp1, header=T, sep='\t', stringsAsFactors=T)
        d2 <- read.table(exp2, header=T, sep='\t', stringsAsFactors=T)
        d3 <- read.table(exp3, header=T, sep='\t', stringsAsFactors=T)
        distances <- rbind(d1, d2, d3)
        print(paste("        DIST:", dim(distances)[1], dim(distances)[2]))

        # We now have all the sorted data for one direction of a given growth

        # save sorted distances in a dist file
        distfile <- paste(basename, '.csv', sep='')
        print(paste("        writing", distfile))
        sorted.distances <- distances[order(distances$SEG, distances$SIZ),]
        #sorted.distances <- distances
        write.table(sorted.distances, file=distfile, 
                    row.names=TRUE, col.names=TRUE,
                    sep='\t', quote=FALSE)

        distances <- sorted.distances
        # interval length
        int.len = 25

        # seggregate distances by segment
        dist.seg.A <- subset(distances, SEG == "A")
        dist.seg.B <- subset(distances, SEG == "B")

        # Make global plots for segment A
        print("        Plotting gap lengths for segment A")
        # plot gap length of each read
        plot.gaps <- paste(basename, '.gaps.A.png', sep='')
          png(plot.gaps, width=1000, height=1000)
        plot(dist.seg.A$SIZ, 
             type = "l", 
             main = basename, 
             sub= " Segment A",
             xlab = "Read (sorted by gap size)", 
             ylab = "Gap length")
        dev.off()
        
        # plot histogram of lengths at specified intervals
        print("        Frequecies of gap lengths for segment A")
        m <- min(dist.seg.A$SIZ[ ! is.na(dist.seg.A$SIZ)])
        if (m > 0) m <- 0
        M <- max(dist.seg.A$SIZ[ ! is.na(dist.seg.A$SIZ)])
        M <- max(M, M + int.len - 1)
        intervals <- seq(
                     from=0,
                     to=M,
                     by=int.len
                     )
        if (m < 0) intervals <- c(m, intervals)
        x.interval.labels <- seq(
                     from=0,
                     to=M,
                     by=100
                     )
        if (m < 0) x.interval.labels <- c(m, x.interval.labels)

        hist.gaps <- paste(basename, '.hist', int.len, '.A.png', sep='') 
        png(hist.gaps, width=1000, height=1000)
          hist(dist.seg.A$SIZ, 
               breaks = intervals, 
               col = "pink", 
               main = basename,
               sub = "Segment A")
        axis(1, at=x.interval.labels, labels=x.interval.labels)
        dev.off()

        # make global plots for segment B
        print("        Plotting gap lengths for segment B")
        plot.gaps <- paste(basename, '.gaps.B.png', sep='')
        png(plot.gaps, width=1000, height=1000)
          plot(dist.seg.B$SIZ, 
               type = "l", 
               main = basename, 
               sub = "Segment B",
               xlab="Read (sorted by gap size)", 
               ylab="Gap length")
        dev.off()

        print("        Frequecies of gap lengths for segment B")
        m <- min(dist.seg.B$SIZ[ ! is.na(dist.seg.B$SIZ)])
        if (m > 0) m <- 0
        M <- max(dist.seg.B$SIZ[ ! is.na(dist.seg.B$SIZ)])
        M <- max(M, M + int.len - 1)
        intervals <- seq(
                     from=0,
                     to=M + int.len - 1,
                     by=int.len
                     )
        if (m < 0) intervals <- c(m, intervals)
        x.interval.labels <- seq(
                     from=0,
                     to=M + int.len - 1,
                     by=100
                     )
        if (m < 0) x.interval.labels <- c(m, x.interval.labels)
        hist.gaps <- paste(basename, '.hist.', int.len, '.B.png', sep='') 
        png(hist.gaps, width=1000, height=1000)
          hist(dist.seg.B$SIZ, 
               breaks = intervals, 
               col = "pink", 
               main = basename,
               sub = "Segment B")
        axis(1, at=x.interval.labels, labels=x.interval.labels)
        dev.off()

        # Histograms showing location of gap start and end points
        #SEGMENT A 
        # Plot forward start and end positions for segment A
        print("        Start-end location (over)")
        rd.start <- dist.seg.A$LEFT
        rd.end <- dist.seg.A$RIGHT
        se.hist <- paste(basename, '.se.', int.len, '.A.png', sep='') 
        png(se.hist, width=1000, height=1000)
          hist( rd.start, 
              breaks = calc.intervals(rbind(rd.start, rd.end), int.len),
              border = "blue", 
              main = basename,
              sub = "Segment A")
          hist( rd.end, 
              breaks = calc.intervals(rbind(rd.start, rd.end), int.len),
              border = "red", 
              add=T)
          axis(1, at= int.len, labels=int.len)
        dev.off()   
        print("        Start-end location (side)")
        se.hist <- paste(basename, '.SE.', int.len, '.A.png', sep='') 
        png(se.hist, width=1000, height=1000)
          multhist(list(rd.start, rd.end),
              breaks = calc.intervals(c(rd.start, rd.end), int.len),
              main = basename,
              sub = "Segment A"
              )
        dev.off()
        
        gg.hist <- paste(basename, ".ES.", int.len, '.A.png', sep='')
        print(paste("        ggplotting to", gg.hist))
        #png(gg.hist, width=1000, height=1000)
          es <- data.frame(Start = rd.start, End = rd.end)
          ggplot(melt(es), 
                 aes(value, fill = variable)) + 
          geom_histogram(position = "dodge", bins=30) +
          #geom_histogram(position = "dodge", bins=(M-m)/int.len)+
          theme(text = element_text(size=6)) +
          labs(title=basename, sub="Segment A")
          ggsave(gg.hist, 
                 width=10, 
                 height=10, 
                 units="cm", 
                 device="png", 
                 dpi=300)
        #dev.off()


        # Plot forward start and end positions for segment B
        rd.start <- dist.seg.B$LEFT
        rd.end <- dist.seg.B$RIGHT
        se.hist <- paste(basename, '.se.', int.len, '.B.png', sep='') 
        png(se.hist, width=1000, height=1000)
          hist( rd.start, 
              breaks = calc.intervals(rbind(rd.start, rd.end), int.len),
              border = "blue", 
              main = basename,
              sub = "Segment B")
          hist( rd.end, 
              breaks = calc.intervals(rbind(rd.start, rd.end), int.len),
              border = "red", 
              add=T)
          axis(1, at= int.len, labels=int.len)
        dev.off()   
        se.hist <- paste(basename, '.SE.', int.len, '.B.png', sep='') 
        png(se.hist, width=1000, height=1000)
          multhist(list(rd.start, rd.end),
              breaks = calc.intervals(c(rd.start, rd.end), int.len),
              main = basename,
              sub = "Segment B"
              )
        dev.off()
        gg.hist <- paste(basename, ".ES.", int.len, '.B.png', sep='')
        print(paste("        ggplotting to", gg.hist))
        #png(gg.hist, width=1000, height=1000)
          es <- data.frame(Start = rd.start, End = rd.end)
          ggplot(melt(es), 
                 aes(value, fill = variable)) + 
          geom_histogram(position = "dodge", bins=30) +
          #geom_histogram(position = "dodge", bins=(M-m)/int.len)+
          theme(text = element_text(size=6)) +
          labs(title=basename, sub="Segment B") 
          ggsave(gg.hist, 
                 width=10, 
                 height=10, 
                 units="cm", 
                 device="png", 
                 dpi=300)
        #dev.off()
        
        # save read data in separate datasets
	if (dir == 1) {
	    r1.dist <- data.frame(distances, READ="D")
	} else {
	    r2.dist <- data.frame(distances, READ="R")
	}
    }

    # Now plot both reads together so we can compare their behaviour
    basename <- paste("two_matches/", pfu, "pfu_exp1-3_R1-2", sep="")
    dist <- rbind(r1.dist, r2.dist)

    # save sorted distances in a dist file
    distfile <- paste(basename, '.csv', sep='')
    print(paste("        writing", distfile))
    sorted.distances <- dist[order(dist$SEG, dist$SIZ),]
    #sorted.distances <- distances
    write.table(sorted.distances, file=distfile, 
                row.names=TRUE, col.names=TRUE,
                sep='\t', quote=FALSE)

    dist <- sorted.distances
    # interval length
    int.len = 25


    dist.seg.A <- subset(dist, SEG == "A")
    dist.seg.B <- subset(dist, SEG == "B")
    rd.start <- dist$LEFT
    rd.end <- dist$RIGHT
    
    # Plots for segment A
    png(paste(basename, ".sizes.A.png", sep=''), width=1000, height=1000)
    with(dist.seg.A, 
         Hist(SIZ, groups=READ, scale="percent", breaks=75, 
	      col="darkgray", main=paste(pfu, "pfu: size A"))
	 )
    dev.off()
    png(paste(basename, ".left.A.png", sep=''), width=1000, height=1000)
    with(dist.seg.A, 
         Hist(LEFT, groups=READ, scale="percent", breaks=75, 
	      col="darkgray", main=paste(pfu, "pfu: left A"))
	 )
    dev.off()
    png(paste(basename, ".right.A.png", sep=''), width=1000, height=1000)
    with(dist.seg.A, 
         Hist(RIGHT, groups=READ, scale="percent", breaks=75, 
	      col="darkgray", main=paste(pfu, "pfu: right A"))
	 )
    dev.off()

    png(paste(basename, ".sizes.B.png", sep=''), width=1000, height=1000)
    with(dist.seg.B, 
         Hist(SIZ, groups=READ, scale="percent", breaks=75, 
	      col="darkgray", main=paste(pfu, "pfu: size B"))
	 )
    dev.off()
    png(paste(basename, ".left.B.png", sep=''), width=1000, height=1000)
    with(dist.seg.B, 
         Hist(LEFT, groups=READ, scale="percent", breaks=75, 
	      col="darkgray", main=paste(pfu, "pfu: left B"))
	 )
    dev.off()
    png(paste(basename, ".right.B.png", sep=''), width=1000, height=1000)
    with(dist.seg.B, 
         Hist(RIGHT, groups=READ, scale="percent", breaks=75, 
	      col="darkgray", main=paste(pfu, "pfu: right B"))
	 )
    dev.off()
}
