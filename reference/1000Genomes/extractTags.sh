## set parameters
# Data directory
dir=/path/to/parent/1000Genomes/directory
# r2 threshold for tag selection
r2=0.1
# kb window for tag selection
kb=500

for chr in `seq 1 22`; do
    # compute all tags
    plink --bfile $dir/plink/EUR.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes \
          --show-tags all \
          --tag-r2 $r2 --tag-kb $kb \
          --out $dir/plink/EUR.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.${kb}kb.r${r2}
    # filter tags without identifier in first column
    awk '$1 != "."' $dir/plink/EUR.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.${kb}kb.r${r2}.tags.list \
    > $dir/plink/EUR.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.${kb}kb.r${r2}.tags.list.rsIDs
done
