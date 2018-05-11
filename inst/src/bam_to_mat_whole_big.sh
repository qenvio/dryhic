#!/usr/bin/env bash

# help message

show_help() {

	cat << EOF

	Usage: ${0##*/} [options] <in.bam> [<chr:start-end>|<chr>]

	Get the whole genome contact matrix from a hic-bam file.

            -h          display this help and exit
            -f INT      only include reads with all FLAG bits set in INT [0]
            -F INT      only include reads with none of the FLAG bits set in INT [783]
            -w INT      size of the window in bp (resolution) [100000]

EOF

}


# get arguments

filterin=0
filterex=783
resolution=100000

while getopts "hf:F:w:r:" opt
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

# make contact matrix

samtools view -f $filterin -F $filterex $inbam $region | \
	awk -v w=$resolution -v OFS="\t" -v W=1000000 \
    '{if($7 == "=") $7 = $3; i = int($4 / w) * w; j = int($8 / w) * w; a[$3 OFS i OFS $7 OFS j] += 1; if((NR % W) == 0){for(k in a) print k, a[k]; split("", a)}}END{for(k in a) print k, a[k]}'
