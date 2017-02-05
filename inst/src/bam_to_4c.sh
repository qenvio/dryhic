#!/usr/bin/env bash

# help message

show_help() {

	cat << EOF

	Usage: ${0##*/} [options] <in.bam> <chr:start-end>

	Get a virsual 4C contact matrix from a hic-bam file for the bait region <chr:start-end>. Only cis contacts are reported.

            -h          display this help and exit
            -f INT      only include reads with all FLAG bits set in INT [0]
            -F INT      only include reads with none of the FLAG bits set in INT [783]
            -w INT      size of the window in bp (resolution) [10000]
            -s INT      size of flanking region interrogated [2000000]

EOF

}


# get arguments

filterin=0
filterex=783
resolution=10000
window=2000000

while getopts "hf:F:w:s:" opt
do

    case "$opt" in
        h)
            show_help
            exit 0
            ;;
        f)  filterin=$OPTARG
            ;;
        F)  filterex=$OPTARG
            ;;
        w)  resolution=$OPTARG
            ;;
        s)  window=$OPTARG
            ;;
        '?')
            show_help >&2
            exit 1
            ;;

    esac

done

shift "$((OPTIND-1))"

if [[ $# -eq 0 ]]
then
	
    show_help >&2
    exit 1
	
fi	

inbam=$1
region=$2

# get coordinates

info=($(echo $region | tr ":-" " "))

chromosome=${info[0]}

start=$(expr ${info[1]} - $window)
end=$(expr ${info[2]} + $window)


# ensure cis-contacts

if [[ $filterex -lt 1024 ]]
then

	filterex=$((filterex+1024))

fi

# make contact matrix

samtools view -f $filterin -F $filterex $inbam $region | \
    awk -v w=$resolution -v OFS="\t" -v s=$start -v e=$end \
        '$8 >= s && $8 <= e{j = int($8 / w) * w; a[j] += 1}END{for(k in a) print k, a[k]}'
