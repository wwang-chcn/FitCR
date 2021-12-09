#! /bin/bash

# Nov-1-2018

# USAGE: $0 <fragments.bed> <fragments_length>

function print_help {
        echo "USAGE: $0 <fragments.bed>"
}


awk -v extsize=${2} '
{
        mid=($2+$3)/2
	if(mid-extsize/4<0) printf $1"\t0\t%d\t"$4"\t"$5"\t"$6"\n", mid+extsize/4; else printf $1"\t%d\t%d\t"$4"\t"$5"\t"$6"\n", mid-extsize/4, mid+extsize/4
}' $1 | sort -S 1% -k1,1 -k2,2g > ${1::(${#1}-4)}_shift.bed
