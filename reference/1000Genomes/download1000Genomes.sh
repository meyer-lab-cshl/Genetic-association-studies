## set parameters
# ftp server for 1000Genomes data
ftp=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp
# parent directory to download 1000Genomes data into
dir=/path/to/parent/1000Genomes/directory

## move into directory
cd $dir

## download sample information
wget $ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx
## manually open xlsx file and save first worksheet with sample IDs and
## ethnicities to $dir/20130606_sample_info_worksheetSampleInfo.csv

## download vcf files
for chr in `seq 1 22`; do
    wget $ftp/release/20130502/supporting/vcf_with_sample_level_annotation/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.vcf.gz
done
