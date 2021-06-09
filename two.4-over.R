
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
            print(paste("KO", i, "Inter-segment", 
                        matches[read.1, segment],
                        "-",
                        matches[read.2, segment],
                        "rearrangement detected"))
            print(matches)
            distances$SEG[i] = paste(matches[read.1, segment], matches[read.2, segment])
            distances$SIZ[i] = Inf
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
            distances$SIZ[i] = NA
            next
        }
        if (matches[read.2, r.start] > matches[read.2, r.end]) {
            # match #2 is being reported in anti-sense for _the read_
            # This shouldn't happen, so we'll notify it for inspection
            print(paste(i, ": inverted second match, please check this"))
            print(matches)
            distances$SIZ[i] = NA
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
	    print(paste("KO", i, "matching genome in reverse for 1 and direct for 2"))
            distances$SIZ[i] = NA
            distances$REAR[i] = 'Y'
            next
            
        } else if (genome.2.end < genome.2.start) {
            # we are matching the reverse strand in fragment 2 but not in
            # fragment 1 (because otherwise the first check would have 
            # succeeded)
            #    |====>
            #            <----|
            
            # this implies a genomic rearrangement
	    print(paste("KO", i, "matching genome in direct for 1 and reverse for 2"))
            distances$SIZ[i] = NA
            distances$REAR[i] = 'Y'
            next
            
	} else {
            # we are matching the direct strand in both fragments (because
            # all other checks have failed)
            
            # this implies a genomic rearrangement
            print(paste("OK", i, "matching genome in direct strand for 1 and 2"))

        }
        
        
        if (genome.1.start < genome.2.start) {
            if (genome.1.end < genome.2.start) {
                # they are separated by a gap
                # |====>   |---->
                print("    |====>")
                print("               |---->")
		distances$SIZ[i] = genome.2.start - genome.1.end
            
            } else {
                # they overlap
                # |====>.......>
                #     |---->
                if (genome.1.end > genome.2.end) {
                    # genome.2 is contained within genome.1
                    # |=========>
                    #    |--->
		    print("    |=======>")
                    print("      |--->")
                    distances$SIZ[i] = 0
                    distances$LEFT = genome.1.start
                    distances$RIGHT = genome.1.end
                } else {
                    # it is a simple overlap
                    # |====>
                    #     |--->
                    print("    |====>")
                    print("       |---->")
                    distances$SIZ[i] = 0
                    distances$LEFT = genome.1.start
                    distances$RIGHT = genome.2.end
                  }
            }
        } else {
            # genome.2.start goes before genome.1.start
            if (genome.1.end > genome.2.start) {
                # we have a gap
                #   <----|   <====|
                print("              <====|")
                print("    <----|")
                distances$SIZ[i] = genome.1.end - genome.2.start
       
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
                    distances$SIZ[i] = 0
                    distances$LEFT = genome.1.end
                    distances$RIGHT = genome.2.start
                } else {
                    # it is a simple overlap
                    #     <=====|
                    # <-----|
                    print("       <====|")
                    print("    <----|")
                    distances$SIZ[i] = 0
                    distances$LEFT = genome.2.end
                    distances$RIGHT = genome.1.start
                }
            }
        }
    }
       
    
    #histogramas para los puntos de inicio y fin de los gaps en toda la longitud del genoma
require("plotrix")	# for multhist
require(ggplot2)	# for ggplot
require(reshape2)	# for melt


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

 
        # interval length
        int.len = 25

        # seggregate distances by segment
        dist.seg.A <- subset(distances, SEG == "A")
        dist.seg.B <- subset(distances, SEG == "B")


        # Histograms showing location of gap start and end points
        #SEGMENT A 
        # Plot forward start and end positions for segment A
        print("        Start-end location (over)")
        rd.start <- dist.seg.A$LEFT
        rd.end <- dist.seg.A$RIGHT
        se.hist <- paste(basename, '.se-overlap.', int.len, '.A.png', sep='') 
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
        se.hist <- paste(basename, '.SE-overlap.', int.len, '.A.png', sep='') 
        png(se.hist, width=1000, height=1000)
          multhist(list(rd.start, rd.end),
              breaks = calc.intervals(c(rd.start, rd.end), int.len),
              main = basename,
              sub = "Segment A"
              )
        dev.off()
        
        gg.hist <- paste(basename, ".ES-overlap.", int.len, '.A.png', sep='')
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
        se.hist <- paste(basename, '.se-overlap.', int.len, '.B.png', sep='') 
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
        se.hist <- paste(basename, '.SE-overlap.', int.len, '.B.png', sep='') 
        png(se.hist, width=1000, height=1000)
          multhist(list(rd.start, rd.end),
              breaks = calc.intervals(c(rd.start, rd.end), int.len),
              main = basename,
              sub = "Segment B"
              )
        dev.off()
        gg.hist <- paste(basename, ".ES-overlap.", int.len, '.B.png', sep='')
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
    }
   }
