step=2000000

metadata=$1
outpath=$2


for chr in {2L,2R,3L,3R,X}; do
    # get length of chromosome
    header=`tabix -H data/vcf/ag1000g.agam.n1470.filtered_variants.$chr\.ann.vcf.gz | grep "##contig=<ID=$chr,length="`
    length=`echo $header | awk '{sub(/.*=/,"");sub(/>/,"");print}'`

    # subset vcf and run genome scan
    endwindow=$step
    for startwindow in `seq 1 $step $length`; do
        tabix -h data/vcf/ag1000g.agam.n1470.filtered_variants.$chr\.ann.vcf.gz $chr\:$startwindow\-$endwindow > data/$outpath\-tmp.vcf
        python scripts/spatial_genome_scan.py \
        --vcf data/$outpath\-tmp.vcf \
        --sample_data $metadata \
        --out out/$outpath\_$chr\_$startwindow\-$endwindow

        endwindow=$((endwindow+step))
        rm data/$outpath\-tmp.vcf
   done
done 
