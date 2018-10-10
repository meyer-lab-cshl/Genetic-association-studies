sig=$1
tags=$2
out=$3

head -n 1 $tags > $out
awk 'NR==FNR {snp[$1]=$0; next} $1 in snp {print $0}' $sig $tags >> $out
