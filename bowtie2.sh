
reads_prefix=$1
echo "Assembling reads from project $reads_prefix"

R1=${reads_prefix}_R1_001.fastq
R2=${reads_prefix}_R2_001.fastq


#mkdir mapping
#cd mapping
 
#copiar los datos de secuencia de referencia y las lecturas en esta carpeta , bien directamente o utilizando el comando "cp"

if [ ! -e ibdv.btindex.1.bt2 ] ; then
    echo "building index for idbv"
    bowtie2-build ibdv.fasta ibdv.btindex
fi

#ibdv, es el nombre de nuestra secuencia de renferencia en este caso.

if [ ! -e ${reads_prefix}_mapped.sam ] ; then
    echo "Running Bowtie2 on $R1 $R2"
    echo "bowtie2 -x ibdv.btindex -1 $R1 -2 $R2 -S ${reads_prefix}_mapped.sam"
    bowtie2 -x ibdv.btindex -1 $R1 -2 $R2 -S ${reads_prefix}_mapped.sam
fi

if [ ! -e ${reads_prefix}_mapped_sorted.bam ] ; then
    echo "Sorting ${reads_prefix}_mapped.sam"
    echo "samtools sort ${reads_prefix}_mapped.sam -o ${reads_prefix}_mapped_sorted.bam"
    samtools sort ${reads_prefix}_mapped.sam -o ${reads_prefix}_mapped_sorted.bam
fi

if [ ! -e ${reads_prefix}_mapped_sorted.bam.bai ] ; then
    echo "indexing ${reads_prefix}_mapped_sorted.bam"
    echo "samtools index ${reads_prefix}_mapped_sorted.bam"
    samtools index ${reads_prefix}_mapped_sorted.bam
fi


if [ ! -e ${reads_prefix}_coverage.depth ] ; then
    samtools depth ${reads_prefix}_mapped_sorted.bam \
    > ${reads_prefix}_coverage.depth
fi

if [ ! -e ${reads_prefix}_coverage_A.png ] ; then
    R --vanilla <<END
    t <- read.table("$OUT_PFX.sorted.depth")
    # t$V1 contains the chromosome ID
    # t$V2 contains the position
    # t$V3 contains the depth at that position
    # get the list of chromosomes
    chromosomes <- levels(as.factor(t\$V1))
    # for each chromosome, plot its sequencing depth
    for (c in chromosomes) {
        png(paste("$OUT_PFX.sorted.depth", c, "png", sep='.'), 
            width=1000, height=1000)
        d <- t[t\$V1==c,3]
        plot(d)
        dev.off()
    }
END
fi


if [ ! -e ${reads_prefix}_unmapped_R1.fq ] ; then
    echo -n "Generating table of unmapped reads"
    # -f 4 means select the read if it is not mapped, even if its mate is
    # -f 12 means select if both the read and its mate are not mapped
    # -b means output as BAM file
    if [ ! -s $OUT_PFX.unmapped.bam ] ; then
    echo -n ":"
    samtools view -b -f 4 ${reads_prefix}_mapped_sorted.bam \
        > ${reads_prefix}_unmapped.bam
    fi
    if [ ! -s $$OUT_PFX.unmapped.sorted.bam ] ; then
        echo -n "."
        samtools sort -n ${reads_prefix}_unmapped.bam \
            -o ${reads_prefix}_unmapped_sorted.bam
    fi
    echo -n "_"
    bamToFastq -i ${reads_prefix}_unmapped_sorted.bam \
        -fq ${reads_prefix}_unmapped_R1.fq \
        -fq2 ${reads_prefix}_unmapped_R2.fq
    echo ""
fi
