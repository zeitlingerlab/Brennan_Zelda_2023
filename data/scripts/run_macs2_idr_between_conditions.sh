#!bin/bash

# Melanie Weilert
# October 2019
# Purpose: Given a string representing a factor, run between-condition pairwise (non-repeating) combinations of IDR.

# Notes: This assumes that there IS an "oracle" peak set to compare to.
# Notes: Make sure your within-condition samples are compariable.
# Notes: Running this script repeatedly can corrupt the .log into making weird statements. Clear idr file environment before rerunning.

# Based on Snakemake file setup.

MACS2_FILE_A=$1 # input sample name
IDR_OUTPUT_DIR=$2 # output idr directory
CONDITION=$3 #prefix string indicating desired condition to grep for
MACS2_COMBINED_FILE=$4 #MACS2 combined between-replicates file

# Define variables (subset MACS2 more strictly than "across conditions")
MACS2_FILEPATH=$(echo ${MACS2_FILE_A} | sed 's!\(.*\)/.*!\1!')
MACS2_FILES_B=$(ls ${MACS2_FILEPATH}/${CONDITION}*_peaks.narrowPeak) #get filepath, then select all that contain condition prefix
PREFIX_A=$(echo ${MACS2_FILE_A} | sed 's!.*/!!' | sed 's/_peaks.*//') #echo, basename, get everything before _peaks
echo "Searching for other files at ${MACS2_FILEPATH}..."
echo "Files found ${MACS2_FILES_B}..."
echo "Prefix assigned ${PREFIX_A}..."


#Make the IDR if it is not existing
mkdir ${IDR_OUTPUT_DIR}

echo "Running IDR for ${PREFIX_A}..."
#Run pairwise replicates of idr
for MACS2_FILE_B in $MACS2_FILES_B; do
  echo ${MACS2_FILE_A}
  echo ${MACS2_FILE_B}
  PREFIX_B=$(echo ${MACS2_FILE_B} | sed 's!.*/!!' | sed 's/_peaks.*//') #echo, basename, get everything before _peaks

  #Check for combination
  IDR_OUTPUT_A_B=${PREFIX_A}_vs_${PREFIX_B}_idr.txt
  IDR_OUTPUT_B_A=${PREFIX_B}_vs_${PREFIX_A}_idr.txt

  #Write to log

  echo -e "\n{${IDR_OUTPUT_A_B}}" >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log

  if [ ${MACS2_FILE_A} != ${MACS2_FILE_B} ]; then #don't run file with itself
    if [[ -f ${IDR_OUTPUT_DIR}/${IDR_OUTPUT_A_B} ]] || [[ -f ${IDR_OUTPUT_DIR}/${IDR_OUTPUT_B_A} ]]; then #if the file exists, do nothing
      echo -e "Skipped because file exists." >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log
      echo "File exists, skipping IDR."
    else
      echo -e "IDR ran on this combo." >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log
      idr --samples ${MACS2_FILE_A} ${MACS2_FILE_B} --idr-threshold 0.05 \
      --input-file-type narrowPeak --output-file ${IDR_OUTPUT_DIR}/${IDR_OUTPUT_A_B} --peak-list ${MACS2_COMBINED_FILE}
    fi
  else
    echo -e "Skipped because file is same." >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log
  fi

done


echo "Goodbye!"
