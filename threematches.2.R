b7files <- list.files("three_matches", pattern=".*b7", full.names=TRUE)

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
    }
    
    #distances <- matrix(c(rep(" ", dim(names)[1]), 
    #                      rep(0,   dim(names)[1])),
    #             ncol=2, byrow=FALSE)
    # Let us use a data frame instead of a matrix
    distances <- data.frame(NAM=rep(" ", dim(names)[1]),
                            SEG=rep(" ", dim(names)[1]),
                            SIZ1=rep(0,   dim(names)[1]),
                            SIZ2=rep(0,   dim(names)[1]),
                            LEFT1=rep(0,   dim(names)[1]),
                            LEFT2=rep(0,   dim(names)[1]),
                            RIGHT1=rep(0,   dim(names)[1]),
                            RIGHT2=rep(0,   dim(names)[1]),
                            REAR=rep("N", dim(names)[1]),
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
       print(matches) }

       #matches[7, ] is q.start
       #matches[8, ] is q.end
       #matches[9, ] is s.start
       #matches[10, ] is s.end

       #SI estamos con two-matches
       s.start.1 <- matches[1, 9] #es s.start del primer match
       s.start.2 <- matches[2, 9] #es s.start del segundo match
       s.start.3 <- matches[3, 9]
       s.end.1 <- matches[1, 10]
       s.end.2 <- matches[2, 10]
       s.end.3 <- matches[3, 10]
       seg1 <- matches[1, 2]
       seg2 <- matches[2, 2]
       seg3 <- matches[3, 2]


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
        read.1 <- 1; read.2 <- 2; read.3 <- 3
	r.start <- 7; r.end <- 8
        g.start <- 9; g.end <- 10

	# detect location of each read match in chromosomal segments

	if (matches[read.1, segment] != matches[read.2, segment] | matches[read.1, segment] != matches[read.3, segment] | matches[read.2, segment] != matches[read.3, segment]) {
            # this is an inter-chromosomal (inter-segment) rearrangement
            # report it
            print(paste("KO", i, "Inter-segment", 
                        matches[read.1, segment],
                        "-",
                        matches[read.2, segment],
                        "-",
                        matches[read.3, segment],
                        "rearrangement detected"))
            print(matches)
            distances$SEG[i] = paste(matches[read.1, segment], matches[read.2, segment], matches[read.3, segment])
            distances$SIZ1[i] = Inf
            distances$REAR[i] = 'Y'
            next
        }
        distances$SEG[i] = matches[read.1, segment]


	# detect read orientation
        if (matches[read.1, r.start] > matches[read.1, r.end]) {
            # match #1 is being reported in anti-sense _for the read_
            # This shouldn't happen, so we'll notify it for inspection
            print(paste(i, ": inverted first match, please check this"))
            print(matches)
            distances$SIZ1[i] = NA
            next
        }
        if (matches[read.2, r.start] > matches[read.2, r.end]) {
            # match #2 is being reported in anti-sense for _the read_
            # This shouldn't happen, so we'll notify it for inspection
            print(paste(i, ": inverted second match, please check this"))
            print(matches)
            distances$SIZ1[i] = NA
            next
        }
        if (matches[read.3, r.start] > matches[read.3, r.end]) {
            # match #2 is being reported in anti-sense for _the read_
            # This shouldn't happen, so we'll notify it for inspection
            print(paste(i, ": inverted second match, please check this"))
            print(matches)
            distances$SIZ1[i] = NA
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
            read.3.start <- matches[read.3, r.start]
            read.3.end <- matches[read.3, r.end]
            genome.1.start <- matches[read.1, g.start]
            genome.1.end <- matches[read.1, g.end]
            genome.2.start <- matches[read.2, g.start]
            genome.2.end <- matches[read.2, g.end]
            genome.3.start <- matches[read.3, g.start]
            genome.3.end <- matches[read.3, g.end]
        } else {
            # the second match is before the first _in the read_
            read.1.start <- matches[read.2, r.start]
            read.1.end <- matches[read.2, r.end]
            read.2.start <- matches[read.1, r.start]
            read.2.end <- matches[read.1, r.end]
            read.3.start <- matches[read.3, r.start]
            read.3.end <- matches[read.3, r.end]
            genome.1.start <- matches[read.2, g.start]
            genome.1.end <- matches[read.2, g.end]
            genome.2.start <- matches[read.1, g.start]
            genome.2.end <- matches[read.1, g.end]
            genome.3.start <- matches[read.3, g.start]
            genome.3.end <- matches[read.3, g.end]
        }

	# At this point, read.1 is the first fragment of the read
        # and read.2 is the second fragment; genome.1 is where the
        # first read fragment matches and genome.2 where the second
        # read fragment matches in the genome

        # since the read fragments are in read order, we can now
        # interpret the genome fragments they represent and
        # calculate distance between matched regions in the genome

        if ((genome.1.end < genome.1.start) && (genome.2.end < genome.2.start) && (genome.3.end < genome.3.start)) {
            # we are matching the reverse strand in all fragments
            print(paste("OK", i, "matching genome in reverse strand for 1 and 2"))

        } else if (genome.1.end < genome.1.start && genome.2.end < genome.3.start) {
            # we are matching the reverse strand in fragment 1 but not in
            # fragment 2 (because otherwise the first check would have 
            # succeeded)
            #   <====|
            #          |---->
            #                 |=-=->
            # this implies a genomic rearrangement
	    print(paste("KO", i, "matching genome in reverse for 1 and direct for 2"))
            distances$SIZ1[i] = NA
            distances$REAR[i] = 'Y'
            next
        } else if (genome.2.end < genome.2.start && genome.1.end < genome.1.start) {
            # we are matching the reverse strand in fragment 1 but not in
            # fragment 2 (because otherwise the first check would have 
            # succeeded)
            #   <====|
            #          <----|
            #                 |=-=->
            # this implies a genomic rearrangement
	    print(paste("KO", i, "matching genome in reverse for 1 and direct for 2"))
            distances$SIZ1[i] = NA
            distances$REAR[i] = 'Y'
            next

        } else if (genome.2.end < genome.2.start && genome.3.end < genome.3.start) {
            # we are matching the reverse strand in fragment 2 but not in
            # fragment 1 (because otherwise the first check would have 
            # succeeded)
            #    |====>
            #            <----|
            #                    <=-=-|
            # this implies a genomic rearrangement
	    print(paste("KO", i, "matching genome in direct for 1 and reverse for 2"))
            distances$SIZ1[i] = NA
            distances$REAR[i] = 'Y'
            next
        } else if (genome.3.end < genome.3.start && genome.1.end < genome.2.start) {
            # we are matching the reverse strand in fragment 2 but not in
            # fragment 1 (because otherwise the first check would have 
            # succeeded)
            #     |====>
            #            |---->
            #                    <=-=-|
            # this implies a genomic rearrangement
	    print(paste("KO", i, "matching genome in direct for 1 and reverse for 2"))
            distances$SIZ1[i] = NA
            distances$REAR[i] = 'Y'
            next
         } else if (genome.2.end < genome.2.start && genome.2.end < genome.3.start) {
            # we are matching the reverse strand in fragment 2 but not in
            # fragment 1 (because otherwise the first check would have 
            # succeeded)
            #    |====>
            #            <----|
            #                   |=-=->
            # this implies a genomic rearrangement
	    print(paste("KO", i, "matching genome in direct for 1 and reverse for 2"))
            distances$SIZ1[i] = NA
            distances$REAR[i] = 'Y'
            next
         } else if (genome.1.end < genome.1.start && genome.3.end < genome.3.start) {
            # we are matching the reverse strand in fragment 2 but not in
            # fragment 1 (because otherwise the first check would have 
            # succeeded)
            #    <====|
            #           |--->
            #                   <=-=-|
            # this implies a genomic rearrangement
	    print(paste("KO", i, "matching genome in direct for 1 and reverse for 2"))
            distances$SIZ1[i] = NA
            distances$REAR[i] = 'Y'
            next

	} else {
            # we are matching the direct strand in both fragments (because
            # all other checks have failed)

            # this implies a genomic rearrangement
            print(paste("OK", i, "matching genome in direct strand for 1 and 2"))

          }
        if (genome.1.start < genome.2.start) {
            if (genome.1.end < genome.2.start && genome.2.end < genome.3.start) {
                # they are separated by a gap
                # |====>   |---->
                print("    |====>")
                print("               |---->")
                print("                         |=-=->")
		distances$SIZ1[i] = genome.2.start - genome.1.end
                distances$LEFT1[i] = genome.1.end
                distances$RIGHT1[i] = genome.2.start
                distances$SIZ2[i] = genome.3.start - genome.2.end
                distances$LEFT2[i] = genome.2.end
                distances$RIGHT2[i] = genome.3.start

            } else {
                # they overlap
                # |====>.......>
                #     |---->

                if (genome.1.end > genome.2.end && genome.1.end < genome.3.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    |=======>")
                    print("      |--->")
                    print("               |=-=->")
                    distances$SIZ1[i] = 0

                } else if (genome.2.end > genome.3.end && genome.1.end < genome.2.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    |=======>")
                    print("              |-------->")
                    print("               |=-=->")
                    distances$SIZ1[i] = 0

               } else if (genome.1.end > genome.2.end && genome.2.end > genome.3.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    |=======>")
                    print("      |--->")
                    print("         |=-=->")
                    distances$SIZ1[i] = 0

               } else if (genome.1.end > genome.2.end && genome.1.end > genome.3.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    |=======>")
                    print("     |-->")
                    print("         |=-=-=->")
                    distances$SIZ1[i] = 0

               } else if (genome.1.end > genome.2.start && genome.2.end > genome.3.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    |=======>")
                    print("           |-->")
                    print("             |=-=-=->")
                    distances$SIZ1[i] = 0
                }
            }
        } else {
            # genome.2.start goes before genome.1.start
            if (genome.1.end > genome.2.start && genome.2.end > genome.3.start) {
                # we have a gap
                #   <----|   <====|
                print("                  <====|")
                print("         <----|")
                print(" <-=-=|")
                distances$SIZ1[i] = genome.1.end - genome.2.start
                distances$LEFT1[i] =  genome.2.start
                distances$RIGHT1[i] = genome.1.end
                distances$SIZ2[i] = genome.2.end - genome.3.start
                distances$LEFT2[i] =  genome.3.start
                distances$RIGHT2[i] = genome.2.end

            } else {
                # they overlap
                #  <......<====|
                #     <-----|
                if (genome.2.end < genome.1.end && genome.1.end < genome.3.start) {
                    # genome.2 is contained in genome.1
                    #  <==========|
                    #     <-----|
                    print("          <========|")
                    print("            <----|")
                    print(" <-=-=|")
                    distances$SIZ1[i] = 0

               } else if (genome.2.end < genome.3.end && genome.1.end < genome.2.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    <======|")
                    print("              <--------|")
                    print("                <-=-=-|")
                    distances$SIZ1[i] = 0

               } else if (genome.1.end < genome.2.end && genome.2.end < genome.3.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    <=======|")
                    print("      <---|")
                    print("         <=-=-|")
                    distances$SIZ1[i] = 0

               } else if (genome.1.end < genome.2.end && genome.1.end < genome.3.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    <=======|")
                    print("     <--|")
                    print("         <=-=-=-|")
                    distances$SIZ1[i] = 0

               } else if (genome.1.end < genome.2.start && genome.2.end < genome.3.start) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    <=======|")
                    print("           <--|")
                    print("             <=-=-=-|")
                    distances$SIZ1[i] = 0
                }
           }


    # save sorted distances in a dist file
    distfile <- paste(basename, '.csv', sep='')
    print(paste("    writing", distfile))
    sorted.distances <- distances[order(distances$SEG, distances$SIZ1, distances$SIZ2),]
    #sorted.distances <- distances
    write.table(sorted.distances, file=distfile, 
                row.names=FALSE, col.names=TRUE,
                sep='\t', quote=FALSE)

    # seleccionar todas las columnas de distances, pero solo las filas
    # en las que el segmento sea 'A'
    l = 100


    distances.A <- subset(distances, SEG == "A")
    #m <- min(distances.A$SIZ1[ ! is.na(distances.A$SIZ1)])
    m <- 0
    #M <- max(distances.A$SIZ2[ ! is.na(distances.A$SIZ2)])
    noNA1 <- c(distances.A$SIZ1)
    noNA2 <- c(distances.A$SIZ2)
    noNA <- c(noNA1, noNA2)
    print(noNA)
    try(M <- max(distances.A$SIZ1[ ! is.na(distances.A$SIZ1)], distances.A$SIZ2))
    l = 100
    M <- max(M, M + l - 1)
    print(paste('MAX A =', M))
    intervals <- seq(
                 from=m,
                 to=M,
                 by=l
                 )
    x.intervals <- seq(
                 from=m,
                 to=M,
                 by=100
                 )

    plotgaps <- paste(basename, '.gaps.A.png', sep='')
    png(plotgaps, width=1000, height=1000)
    plot(noNA)
    dev.off()

    histgaps <- paste(basename, '.hist', l, '.A.png', sep='') 
    png(histgaps, width=1000, height=1000)
    try(hist(noNA, breaks = intervals, col = "pink"))
    try(axis(1, at=x.intervals, labels=x.intervals))
    dev.off()

    distances.B <- subset(distances, SEG == "B")
    #m <- min(distances.B$SIZ1[ ! is.na(distances.B$SIZ1)])
    m <- 0
    #M <- max(distances.B$SIZ2[ ! is.na(distances.B$SIZ2)])
    noNA1 <- c(distances.B$SIZ1)
    noNA2 <- c(distances.B$SIZ2)
    noNA <- c(noNA1, noNA2)
    print(noNA)
    try(M <- max(distances.B$SIZ1[ ! is.na(distances.B$SIZ1)], distances.B$SIZ2))
    l = 100
    M <- max(M, M + l - 1)
    print(paste('MAX B =', M))
    intervals <- seq(
                 from=m,
                 to=M,
                 by=l
                 )
    x.intervals <- seq(
                 from=m,
                 to=M,
                 by=100
                 )
    plotgaps <- paste(basename, '.gaps.B.png', sep='')
    png(plotgaps, width=1000, height=1000)
    plot(noNA)
    dev.off()

    histgaps <- paste(basename, '.hist.', l, '.B.png', sep='') 
    png(histgaps, width=1000, height=1000)
    try(hist(noNA, breaks = intervals, col = "pink"))
    try(axis(1, at=x.intervals, labels=x.intervals))
    dev.off()

}


