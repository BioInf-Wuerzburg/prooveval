#!/bin/bash

# Author: Thomas Hackl, thomas.hackl@uni-wuerzburg.de

## TODO

shopt -s extglob

# Variable defaults
DEBUG=1;

# reliable detect binary folder
pushd `dirname $0` > /dev/null
BD=`pwd`
popd > /dev/null

SEQCHUNKER=SeqChunker;
THREAD_NUM=1;
OUT="prooveval";

show_usage()
{
	echo "Usage: prooveval.sh [OPTIONS] -g <GMAPSRY> -- <prooveval passthrough parameter> > STATS"
	echo " e.g.: prooveval.sh -t 20 -g ecoli.fa.sry -- -r EcGen.fa -q ecoli.fa > ecoli.fa.stats"
}

show_help()
{
  show_usage
  cat <<EOF

Parallelization wrapper for prooveval using xargs. Input is splitted on the fly
using SeqChunker, output is merged using prooveval-merge.pl.

  -g|--gmap-sry               Gmap summary file, required.
  -t|--threads                Number of parallel processes [$THREAD_NUM]
  -o|--out                    Out file prefix [$OUT]
  --seqchunker-bin            Path to SeqChunker script [$SEQCHUNKER]
  --debug
  -h|--help

EOF
}

# Execute getopt
ARGS=`getopt --name "prooveval.sh" --alternative\
    --options "t:g:h" \
    --longoptions "gmap-sry:,threads:,seqchunker-bin:help,debug" \
    -- "$@"`

# no bad args - are passed trough to blat
#Bad arguments
#[ $? -ne 0 ] && exit 1;

# A little magic
eval set -- "$ARGS"


# Now go through all the options
while true; do
    case "$1" in
        -t|--threads)
            [ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
            if [[ "$2" =~ [^0-9] ]]; then
                echo "$1: '$2' not an INT" 1>&2;
                exit 1;
            fi;
            THREAD_NUM=$2;
            shift 2;;

		-g|--gmap-sry)
			[ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
			GMAP_SRY=$2
			shift 2;;
        	      	
      	--seqchunker-bin)
			[ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
			SEQCHUNKER=$2
			shift 2;;
    	
    	--debug)
        	REFINE_ALL=1
        	shift;;        	           	    	        	        
        	           	    	        	                   	    	        	                   	    	        	        
	    -h|--help)
	       show_help && exit 0;;
	       
        --)
            shift
            break;;
            
        *)
	    
  esac
done

if [ $# -lt 3 ]; then
	echo "REF QRY and GMAP_SRY required" 1>&2;
	show_usage;
	exit 1;
fi;


CMDS=$(for CC in $(seq 1 $THREAD_NUM); do
	# prepend stream with additional '>', it is used up by blat for format checking
	CMD="$SEQCHUNKER"' -n '"$THREAD_NUM"' -f '"$CC"' -l '"$CC"' '"$GMAP_SRY"' | prooveval '"$@"' --out '"$OUT.$CC"'.stats --gmap-sry -'
	echo "'$CMD'";
done;)
[ $DEBUG -gt 1 ] && echo "$CMDS";
echo "$CMDS" | xargs -n1 -P $THREAD_NUM bash -c 

# merge results
perl "$BD"/prooveval-merge.pl "$OUT".*.stats;

rm "$OUT".*.stats;
