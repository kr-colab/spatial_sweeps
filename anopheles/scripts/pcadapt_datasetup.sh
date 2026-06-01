#!/bin/bash
# set up Anopheles data for pcadapt
vcfpath=/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/phased_vcf/gamb/
metadata=data/admixture_k1_metadata.txt
bcfmeta=data/admixture_k1_metadata.bcftools.txt

# get only individuals from analysis
awk '{ print $3 }' $metadata > $bcfmeta; sed -i '1d' $bcfmeta
for chr in 2L 2R 3L 3R; do
    echo "filtering samples from chromosome $chr..."
    bcftools view -S $bcfmeta $vcfpath\gamb.$chr\.phased.n1470.derived.vcf.gz > /sietch_colab/crehmann/filtered_ag3/$chr.vcf
done

vcfpath="/sietch_colab/crehmann/filtered_ag3/"
# filter out inversions from 2L and 2R
echo "Filtering out 2La..."
bgzip $vcfpath\2L.vcf
tabix -p vcf $vcfpath\2L.vcf.gz
echo "Filtering out 2La..."
header=$(tabix -H $vcfpath\2L.vcf.gz | grep "##contig=<ID=2L,length=")
length=$(echo $header | awk '{sub(/.*=/,"");sub(/>/,"");print}')
echo -e "2L\t1\t20524058\n2L\t421655321\t$length" > data/2La_mask.bed
bcftools view -R data/2La_mask.bed $vcfpath\2L.vcf.gz > $vcfpath\2L.masked.vcf

echo "Filtering out 2Rb..."
bgzip $vcfpath\2R.vcf
tabix -p vcf $vcfpath\2R.vcf.gz
header=$(tabix -H $vcfpath\2R.vcf | grep "##contig=<ID=2R,length=")
length=$(echo $header | awk '{sub(/.*=/,"");sub(/>/,"");print}') 
echo -e "2R\t1\t18575300\n2R\t26767588\t$length" > data/2Rb_mask.bed
bcftools view -R data/2Rb_mask.bed $vcfpath\2R.vcf.gz > $vcfpath\2R.2Rbmasked.vcf
bgzip $vcfpath\2R.2Rbmasked.vcf
tabix -p vcf $vcfpath\2R.2Rbmasked.vcf.gz
echo "Filtering out 2Rc..."
echo -e "2R\t1\t26750000\n2R\t31473000\t$length" > data/2Rc_mask.bed
bcftools view -R data/2Rc_mask.bed $vcfpath\2R.2Rbmasked.vcf.gz > $vcfpath\2R.2Rbmasked.2Rcmasked.vcf

# convert to plink files
echo "converting to PLINK..."
plink --vcf $vcfpath\2L.masked.vcf --make-bed --out $vcfpath\2L --allow-extra-chr
plink --vcf $vcfpath\2R.2Rbmasked.2Rcmasked.vcf --make-bed --out $vcfpath\2R --allow-extra-chr
for chr in 3L 3R; do
    plink --vcf $vcfpath$chr.vcf --make-bed --out $vcfpath$chr --allow-extra-chr
done

cat /sietch_colab/crehmann/filtered_ag3/2L.masked.vcf | sed '/^#/d' | awk '{print $2}' > /sietch_colab/crehmann/filtered_ag3/2L.positions
cat /sietch_colab/crehmann/filtered_ag3/2R.2Rbmasked.2Rcmasked.vcf | sed '/^#/d' | awk '{print $2}' > /sietch_colab/crehmann/filtered_ag3/2R.positions
cat /sietch_colab/crehmann/filtered_ag3/3L.vcf | sed '/^#/d' | awk '{print $2}' > /sietch_colab/crehmann/filtered_ag3/3L.positions
cat /sietch_colab/crehmann/filtered_ag3/3R.vcf | sed '/^#/d' | awk '{print $2}' > /sietch_colab/crehmann/filtered_ag3/3R.positions