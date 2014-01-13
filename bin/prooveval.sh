#!/bin/bash

# Author: Thomas Hackl, thomas.hackl@uni-wuerzburg.de

## TODO

shopt -s extglob

# Variable defaults
DEBUG=1;

SEQCHUNKER=SeqChunker;
THREAD_NUM=1;

show_usage()
{
  echo "Usage: prooveval.sh [--threads <INT>] [--seqchunker-bin /path/to/SeqChunker] <prooveval parameter>"
}

show_help()
{
  show_usage
  cat <<EOF

Run multiple instances of prooveval using xargs. Input is splitted on the fly
using SeqChunker, output is merged here.

  -g|--gmap-sry               Gmap summary file, required.
  -t/--threads                Number of parallel processes [$THREAD_NUM]
  --seqchunker-bin            Path to SeqChunker script [$SEQCHUNKER]
  --debug
  -h|--help

EOF
}

# Execute getopt
ARGS=`getopt --name "prooveval.sh" --alternative\
    --options "t:g:" \
    --longoptions "gmap-sry:,threads:,seqchunker-bin:" \
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

		--gmap-sry)
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
	CMD="$SEQCHUNKER"' -n '"$THREAD_NUM"' -f '"$CC"' -l '"$CC"' '"$GMAP_SRY"' | prooveval '"$@"' --out '"$CC"' --gmap-sry -'
	echo "'$CMD'";
done;)
[ $DEBUG -gt 1 ] && echo "$CMDS";
echo "$CMDS" | xargs -n1 -P $THREAD_NUM bash -c 




