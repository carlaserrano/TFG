# we will store allexp the depth data in a data.frame callexped "allexp"

rm(allexp)
# the first time 'allexp' will not exist yet
for (f in infiles) { 
    print(f) ; 
    # get the name of the experiment
    expnam <- gsub("_coverage.depth","", f); 
    print(expnam) ; 
    depth <- read.table(f) ; 
    colnames(depth) <- c('Chr', 'nt', expnam) ; 
    print(head(depth)) ; 
    
    # normalize the data from 0 to 1
    M <- max(depth[expnam])
    depth[expnam] = depth[expnam]/M

    if (!exists('allexp')) 
    	# the first time 'allexp' does not exist, we assign everything
        allexp <- depth 
    else
        # the rest of times, we only add the experiment depth column
	# NOTE: this will fail if the number of depth entries is not
        # the same in allexp files!!!
        #allexp <- data.frame(allexp, depth[expnam]) }
        # this should work instead
        allexp <- merge(allexp, depth, by=c(1,2), all=TRUE)

}
