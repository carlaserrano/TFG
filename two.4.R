require("plotrix")

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

b7files <- list.files("two_matches", pattern=".*b7", full.names=TRUE)

for (i in 1:length(b7files)) {
    # get name of B7 file
    b7file <- b7files[i]
    print(paste("Processing file", i, b7file))
    # get base name
    basename <- gsub('.{5}$', '', b7file)
    # create namefile name
    namefile <- paste(basename, '.names.txt', sep='')
    print(paste("    will use", namefile, "for names"))    
    
    # open the two files (blast and read-names)
    if (file.exists(b7file)) {
        if (file.info(b7file)$size == 0) {
            print(paste("   ERROR:", b7file, "is empty"));
            next
        }
        b7 <- read.table(b7file, stringsAsFactors=FALSE)
    } else {
        print(paste("    ERROR:", b7file, "doesn't exist!"))
        next
    }
    if (file.exists(namefile)) {
        if (file.info(namefile)$size== 0) {
            print(paste("   ERROR:", namefile, "is empty"));
            next
        }
        names <- read.table(namefile)
    } else {
        print(paste("    ERROR:", b7file, "doesn't exist!"))
        next
    }
  
    #distances <- matrix(c(rep(" ", dim(names)[1]), 
    #                      rep(0,   dim(names)[1])),
    #             ncol=2, byrow=FALSE)
    # Let us use a data frame instead of a matrix
    distances <- data.frame(NAM=rep(" ", dim(names)[1]),
                            SEG=rep(" ", dim(names)[1]),
                            SIZ=rep(0,   dim(names)[1]),
                            LEFT=rep(0,   dim(names)[1]),
                            RIGHT=rep(0,   dim(names)[1]),
                            FW=rep(NA,   dim(names)[1]),
                            REAR=rep(FALSE, dim(names)[1]),
                            stringsAsFactors=FALSE)
    
    # for each read name, locate its blast matches in the b7 file
    # and do the calculations
    i = 0
    for ( name in names[ ,1] ) {
       i = i + 1
       distances$NAM[i] = name
       print(distances$NAM[i])
       #if (i == 5) { stop() }
       #print(i)
       #buscar en b7 las filas cuya primera columna sea igual a "name"
       matches <- b7[b7[ ,1] == name, ]
       print(matches) 

       #matches[7, ] is q.start
       #matches[8, ] is q.end
       #matches[9, ] is s.start
       #matches[10, ] is s.end

       #SI estamos con two-matches
       s.start.1 <- matches[1, 9] #es s.start del primer match
       s.start.2 <- matches[2, 9] #es s.start del segundo match
       s.end.1 <- matches[1, 10]
       s.end.2 <- matches[2, 10]
       seg1 <- matches[1, 2]
       seg2 <- matches[2, 2]

                             
#            s1    s2    s3    s4
#                                      S4     S3  S2    S1     
#       ----------------------------------------------------------------
#            +----->     +----->
#            R1    R2    R3    R4
#                                      <-----+    <-----+
#                                      R4    R3   R2    R1
#            .
#       Inverted order:
#
#            +----->     +----->
#            R3    R4    R1    R2
#                                      <-----+    <-----+
#                                      R2    R1   R4    R3
#            
       
        # convenience variables to access matches elements by name
        segment <- 2
        read.1 <- 1; read.2 <- 2
	r.start <- 7; r.end <- 8
        g.start <- 9; g.end <- 10

	# detect location of each read match in chromosomal segments
        
	if (matches[read.1, segment] != matches[read.2, segment]) {
            # this is an inter-chromosomal (inter-segment) rearrangement
            # report it
            print(paste(b7file, "KO", i, "Inter-segment", 
                        matches[read.1, segment],
                        "-",
                        matches[read.2, segment],
                        "rearrangement detected"))
            print(matches)
            distances$SEG[i]  <- paste(matches[read.1, segment], matches[read.2, segment])
            distances$SIZ[i]  <- Inf
            distances$REAR[i] <- TRUE
            distances$FW[i]   <- NA
            next
        }
        distances$SEG[i] <- matches[read.1, segment]


	# detect read orientation
        if (matches[read.1, r.start] > matches[read.1, r.end]) {
            # match #1 is being reported in reverse _read order_
            # This shouldn't happen, so we'll notify it for inspection
            print(paste("WARNING:", i, "inverted first match, please check this"))
            print(matches)
            distances$SIZ[i] <- NA
            distances$FW[i]  <- NA
            next
        }
        if (matches[read.2, r.start] > matches[read.2, r.end]) {
            # match #2 is being reported in reverse _read_ order
            # This shouldn't happen, so we'll notify it for inspection
            print(paste("WARNING:", i, " inverted second match, please check this"))
            print(matches)
            distances$SIZ[i] <- NA
            distances$FW[i]  <- NA
            next
        }

	# Both matches are correctly oriented and match the same
        # genome segment:
        # check their relative positions and place them in the
        # correct genome order

        if (matches[read.1, r.start] < matches[read.2, r.start]) {
            # the first match is before the second _in the read_
            read.1.start <- matches[read.1, r.start]
            read.1.end <- matches[read.1, r.end]
            read.2.start <- matches[read.2, r.start]
            read.2.end <- matches[read.2, r.end]
            genome.1.start <- matches[read.1, g.start]
            genome.1.end <- matches[read.1, g.end]
            genome.2.start <- matches[read.2, g.start]
            genome.2.end <- matches[read.2, g.end]
        } else {
            # the second match is before the first _in the read_
            read.1.start <- matches[read.2, r.start]
            read.1.end <- matches[read.2, r.end]
            read.2.start <- matches[read.1, r.start]
            read.2.end <- matches[read.1, r.end]
            genome.1.start <- matches[read.2, g.start]
            genome.1.end <- matches[read.2, g.end]
            genome.2.start <- matches[read.1, g.start]
            genome.2.end <- matches[read.1, g.end]
        }        

	# At this point, read.1 is the first fragment of the read
        # and read.2 is the second fragment; genome.1 is where the
        # first read fragment matches and genome.2 where the second
        # read fragment matches in the genome
        
        # since the read fragments are in read order, we can now
        # interpret the genome fragments they represent and
        # calculate distance between matched regions in the genome
        
        
        # we now detect major rearrangements (each fragment goes in a
        # different direction
        if ((genome.1.end < genome.1.start) && (genome.2.end < genome.2.start)) {
            # we are matching the reverse strand in both fragments
            print(paste("OK", i, "matching genome in reverse strand for 1 and 2"))
            
        } else if (genome.1.end < genome.1.start) {
            # we are matching the reverse strand in fragment 1 but not in
            # fragment 2 (because otherwise the first check would have 
            # succeeded)
            #   <====|
            #          |---->
            
            # this implies a genomic rearrangement
	    print(paste(b7file, "KO", i, "matching genome in reverse for 1 and direct for 2"))
            distances$SIZ[i] <- NA
            distances$REAR[i] <- TRUE
            distances$FW[i] <- NA
            next
            
        } else if (genome.2.end < genome.2.start) {
            # we are matching the reverse strand in fragment 2 but not in
            # fragment 1 (because otherwise the first check would have 
            # succeeded)
            #    |====>
            #            <----|
            
            # this implies a genomic rearrangement
	    print(paste(b7file, "KO", i, "matching genome in direct for 1 and reverse for 2"))
            distances$SIZ[i] <- NA
            distances$REAR[i] <- TRUE
            distances$FW[i] <- NA
            next
            
	} else {
            # we are matching the direct strand in both fragments (because
            # all other checks have failed)
            
            # this implies a genomic rearrangement
            print(paste("OK", i, "matching genome in direct strand for 1 and 2"))

        }
        
        # There are no major rearrangements, but there may still be overlaps
        if (genome.1.start < genome.2.start) {
            distances$SIZ[i] <- genome.2.start - genome.1.end
            distances$LEFT[i] = genome.1.end
            distances$RIGHT[i] = genome.2.start
            distances$FW[i] = TRUE
            print(paste("SIZE =", distances$SIZ[i]))
            if (genome.1.end < genome.2.start) {
                # they are separated by a gap
                # |====>   |---->
                print("    |====>")
                print("               |---->")
                # size will be positive
            } else {
                # they overlap
                # |====>.......>
                #     |---->
                if (genome.1.end >= genome.2.end) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    |=======>")
                    print("      |--->")
                    # this could happen if tthere are tandemly repeated regions
                    # size will be negative to indicate overlap
                    distances$REAR[i] <- TRUE
                    # we'll consider overlaps to be the result of
                    # rearrangements
                } else {
                    # it is a simple overlap
                    # |====>
                    #     |--->
                    print("    |====>")
                    print("       |---->")
                    # size will be negative to indicate overlap
                    distances$REAR[i] <- TRUE
                }
            }
        } else {
            distances$SIZ[i]   <- genome.1.end - genome.2.start
            distances$LEFT[i]  <-  genome.2.start
            distances$RIGHT[i] <- genome.1.end
            distances$FW[i]    <- FALSE
            # genome.2.start goes before genome.1.start
            if (genome.1.end > genome.2.start) {
                # we have a gap
                #   <----|   <====|
                print("              <====|")
                print("    <----|")
            } else {
                # they overlap
                #  <......<====|
                #     <-----|
                if (genome.2.end < genome.1.end) {
                    # genome.2 is contained in genome.1
                    #  <==========|
                    #     <-----|
                    print("    <========|")
                    print("      <----|")
		    # this might be a blast glitch
                    # size will be negative to indicate overlap
                    distances$REAR[i] <- TRUE
                } else {
                    # it is a simple overlap
                    #     <=====|
                    # <-----|
                    print("       <====|")
                    print("    <----|")
                    distances$REAR[i] <- TRUE
                }
            }
        }
    }
       
    # save sorted distances in a dist file
    distfile <- paste(basename, '.csv', sep='')
    print(paste("    writing", distfile))
    sorted.distances <- distances[order(distances$SEG, distances$SIZ),]
    #sorted.distances <- distances
    write.table(sorted.distances, file=distfile, 
                row.names=TRUE, col.names=TRUE,
                sep='\t', quote=FALSE)
                
    # seggregate subsets
    distances.A <- subset(distances, SEG == "A")
    distances.B <- subset(distances, SEG == "B")

    # define interval length for histograms
    int.len = 25

    # Plot gap lengths for segment A
    print("Plotting gaps for segment A")
    plotgaps <- paste(basename, '.gaps.A.png', sep='')
    png(plotgaps, width=1000, height=1000)
    plot(distances.A$SIZ, type = "l")
    dev.off()

    print("Frequecies of gap lengths for segment A")
    m <- min(distances.A$SIZ[ ! is.na(distances.A$SIZ)])
    if (m > 0) m <- 0
    M <- max(distances.A$SIZ[ ! is.na(distances.A$SIZ)])
    M <- max(M, M + int.len - 1)
    intervals <- seq(
                 from=0,
                 to=M,
                 by=int.len
                 )
    if (m < 0) intervals <- c(m, intervals)
    x.intervals <- seq(
                 from=0,
                 to=M,
                 by=100
                 )
    if (m < 0) x.intervals <- c(m, x.intervals)
    
    histgaps <- paste(basename, '.hist', int.len, '.A.png', sep='') 
    png(histgaps, width=1000, height=1000)
    hist(distances.A$SIZ, breaks = intervals, col = "pink", main = basename)
    axis(1, at=x.intervals, labels=x.intervals)
    dev.off()

    # Plot gap lengths for segment B
    print("Plotting gaps for segment B")
    plotgaps <- paste(basename, '.gaps.B.png', sep='')
    png(plotgaps, width=1000, height=1000)
    plot(distances.B$SIZ, type = "l")
    dev.off()

    print("Frequecies of gap lengths for segment B")
    m <- min(distances.B$SIZ[ ! is.na(distances.B$SIZ)])
    if (m > 0) m <- 0
    M <- max(distances.B$SIZ[ ! is.na(distances.B$SIZ)])
    M <- max(M, M + int.len - 1)
    intervals <- seq(
                 from=0,
                 to=M + int.len - 1,
                 by=int.len
                 )
    if (m < 0) intervals <- c(m, intervals)
    x.intervals <- seq(
                 from=0,
                 to=M + int.len - 1,
                 by=100
                 )
    if (m < 0) x.intervals <- c(m, x.intervals)
    histgaps <- paste(basename, '.hist.', int.len, '.B.png', sep='') 
    png(histgaps, width=1000, height=1000)
    hist(distances.B$SIZ, breaks = intervals, col = "pink", main = basename)
    axis(1, at=x.intervals, labels=x.intervals)
    dev.off()
    
    
    # Histograms of gap start and end point locations
        
    int.len=50

    # Plot forward start and end positions for segment A
    print("Start-end location (over)")
    rd.start <- distances.A$LEFT
    rd.end <- distances.A$RIGHT
    se.hist <- paste(basename, '.se.', int.len, '.A.png', sep='') 
    png(se.hist, width=1000, height=1000)
      hist( rd.start, 
          breaks = calc.intervals(rbind(rd.start, rd.end), int.len),
          border = "blue", 
          main = basename)
      hist( rd.end, 
          breaks = calc.intervals(rbind(rd.start, rd.end), int.len),
          border = "red", 
          add=T)
      axis(1, at= int.len, labels=int.len)
    dev.off()   
    print("Start-end location (side)")
    se.hist <- paste(basename, '.SE.', int.len, '.A.png', sep='') 
    png(se.hist, width=1000, height=1000)
      multhist(list(rd.start, rd.end),
          breaks = calc.intervals(c(rd.start, rd.end), int.len),
          title=basename,
          subtitle="Segment A"
          )
    dev.off()
  
    # Plot forward start and end positions for segment B
    rd.start <- distances.B$LEFT
    rd.end <- distances.B$RIGHT
    se.hist <- paste(basename, '.se.', int.len, '.B.png', sep='') 
    png(se.hist, width=1000, height=1000)
      hist( rd.start, 
          breaks = calc.intervals(rbind(rd.start, rd.end), int.len),
          border = "blue", 
          main = basename)
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
          title=basename,
          subtitle="Segment B"
          )
    dev.off()
}
