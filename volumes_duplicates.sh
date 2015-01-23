#!/bin/bash 

export FSLDIR="/usr/local/fsl5"
. /usr/local/fsl5/etc/fslconf/fsl.sh

usage()
{
cat << EOF
usage: $0 options

This script checks for imaging sequences that were run multiple times.

OPTIONS:
  -h  Show this message
  -N  Full path to your raw NIfTI directory (will be created by this script if it doesn't already exist)
  -D  Duplicates file 
  -O  Output file
   
EOF
}

#Read arguments entered by user
while getopts "hN:D:O:" OPTION; do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    N)
      NIFTIDIR=$OPTARG
      ;;
    D)
      DUPS=$OPTARG
      ;;
    O)
      OUTPUT=$OPTARG
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

printf "ID\tParadigm\tSeries\tXdim\tYdim\tZdim\tTdim\n" >> $OUTPUT

tail -n +2 "$DUPS" | while read ID RUNS SEQ; do 
  NIFTI=$NIFTIDIR/$ID/ACE_${ID}_${SEQ}_srs*.nii.gz
  X=$(/usr/local/fsl5/bin/fslinfo $NIFTI | grep -P "^dim1" | awk '{print $2}') 
  Y=$(/usr/local/fsl5/bin/fslinfo $NIFTI | grep -P "^dim2" | awk '{print $2}') 
  Z=$(/usr/local/fsl5/bin/fslinfo $NIFTI | grep -P "^dim3" | awk '{print $2}') 
  T=$(/usr/local/fsl5/bin/fslinfo $NIFTI | grep -P "^dim4" | awk '{print $2}') 
  SRS=$(echo $NIFTI | sed 's/.*_srs\([0-9]*\).nii.gz/\1/g')
  printf "${ID}\t${SEQ}\t${SRS}\t${X}\t${Y}\t${Z}\t${T}\n" >> $OUTPUT 
done 
