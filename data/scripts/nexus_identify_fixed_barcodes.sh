#!bin/bash

# Melanie Weilert
# January 2020
# Purpose: In a sampled subset of a .fastq.gz file, determine the chipnexus barcode frequencies.
# This is useful when you are trying to determine which set of ChIP-nexus fixed barcodes were used in a previous sample.

#Set defaults
nsamples=250000
fixstart=6
fixend=9
helpmessage="Usage:\n\
determine_nexus_barcodes.sh \n\
-h  Display this help message.\n\
-i  Input .fastq.gz filepath.\n\
-o  Output .txt file with ordered barcode frequencies.\n\
-n  Integer listing number of reads to sample. [default = 250000]\n\
-s  Sequence position that fixed barcode starts. [default = 6]\n\
-e  Sequence position that fixed barcode ends. [default = 9]"

#Option parsing
while getopts "i:o:n:s:e:h" opt; do
  case ${opt} in
    i )
      input=$OPTARG
      ;;
    o )
      output=$OPTARG
      ;;
    n )
      nsamples=$OPTARG
      ;;
    s )
      fixstart=$OPTARG
      ;;
    e )
      fixend=$OPTARG
      ;;
    h )
      echo -e $helpmessage
      exit 1
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument. Try -h" 1>&2
      ;;
  esac
done
#shift $((OPTIND -1)) #parse through options
if [ -z "$*" ]; then echo "No args, please type: bash determine_nexus_barcodes.sh -h"; exit 0; fi #create error message.


#Define variables
sampledfastq=$(basename ${input})\.sampled.tmp.fastq
trimmedfastq=$(basename ${input})\.temp.fastq

#Subset values
nlines=$(($nsamples * 4))
gzip -cd $input | head -n $nlines | fastx_trimmer -f $fixstart -l $fixend -o $trimmedfastq -i -


#Parse, sort by fragment order, find unique values, order by most common
paste - - - - <$trimmedfastq | cut -f2 | sort | uniq -c | sort -k 1nr > $output

paste - - - - <$trimmedfastq | cut -f2 | sort | uniq -c | sort -k 1nr | head -n 4| awk '{print $2}' | sort | paste -s -d, > $output.str.txt

rm $trimmedfastq

echo "Finished!"
