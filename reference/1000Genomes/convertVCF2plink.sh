## set parameters
# parent directory with 1000Genomes data
dir=/path/to/parent/1000Genomes/directory

## convert vcf to plink format
for chr in `seq 1 22`; do
    plink2 --vcf $dir/vcf/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.vcf.gz \
        --make-bed \
        --out $dir/plink/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes
done

## Filter European samples
#CEU, Utah residents with Northern and Western European ancestry
#FIN, Finnish in Finland
#GBR, British in England and Scotland
#IBS, Iberian populations in Spain
#TSI, Toscani in Italy
awk -F, '$3 == "CEU" || $3 == "FIN" || $3 == "GBR" || $3 == "IBS" || $3 =="TSI" {print $1," ",$2}' $dir/20130606_sample_info_worksheetSampleInfo.csv > $dir/1000Genomes_EUR_samples.txt

## extract European samples
for chr in `seq 1 22`; do
    plink2 --bfile $dir/plink/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes \
        --keep $dir/1000Genomes_EUR_samples.txt \
        --make-bed \
        --out $dir/plink/EUR.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes
done

