for i in *IBDV_CNB.bam ;
do
  echo $i ;
    samtools depth $i \
    > `basename $i IBDV_CNB.bam`coverage.depth
done

for i in *coverage.depth ;
do
  echo $i ;
    R --vanilla <<END
    t <- read.table("$i")
    png("`basename $i coverage.depth`coverage_A.png", width=1000, height=1000)
    d <- t[t$V1=='A',3]
    plot(d)
    dev.off()
    png("`basename $i coverage.depth`coverage_B.png", width=1000, height=1000)
    d <- t[t$V1=='B',3]
    plot(d)
    dev.off()
END
done
  
for i in *IBDV_CNB.bam ;
do
  echo $i ;
    echo -n "Generating table of unmapped reads"
    # -f 4 means select the read if it is not mapped, even if its mate is
    # -f 12 means select if both the read and its mate are not mapped
    # -b means output as BAM file
    echo -n "."
    samtools view -b -f 4 $i\
        > `basename $i IBDV_CNB.bam`unmapped.bam
    echo -n "."
    samtools sort -n $i \                          
        -o `basename $i unmapped.bam`unmapped_sorted.bam
    echo -n "."
    bamToFastq -i $i\
        -fq `basename $i unmapped_sorted.bam`unmaped_R1.fq \
        -fq2 `basename $i unmapped_sorted.bam`unmaped_R2.fq
    echo ""
done             
