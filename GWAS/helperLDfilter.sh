sig=$1
tags=$2
out=$3

awk 'NR==FNR {snp[$1]=$0; next} $1 in snp {print $0}' $sig $tags > $out
