source /home/ruofany/my_python/3.7.12/bin/activate
    for f in *LMNB*bam
    do
        BEDPREFIX=$(echo $f | sed  's/_.*//') #change according to naming pattern
        input=$(find . -name "$BEDPREFIX*IGG*.bam")
        echo $BEDPREFIX
        echo $input
        bsub -M 50000 -n 16 -o bamC.out bamCompare -b1 $f -b2 $input -of bedgraph  -o  "$f.bedgraph" -p 16 -bs 5000 --smooth 15000 --effectiveGenomeSize 2913022398  --exactScaling --normalizeUsing RPKM  --scaleFactorsMethod None -bl /home/ruofany/bowtie2Index38/hg38-blacklist.v2.bed
    done
