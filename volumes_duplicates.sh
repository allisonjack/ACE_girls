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
  echo $NIFTIDIR/$ID/ACE_${ID}_${SEQ}_srs*.nii.gz | sed 's/\s/\n/g' | sed 's/\//\t/g' | awk 'BEGIN {OFS = "\t"} {print $(NF-1), $NF}' >> $NIFTIDIR/nifti_tmp.txt 
done

cat $NIFTIDIR/nifti_tmp.txt | while read ID NIFTI; do 
  X=$(fslinfo $NIFTIDIR/$ID/$NIFTI | grep -P "^dim1" | awk '{print $2}') 
  Y=$(fslinfo $NIFTIDIR/$ID/$NIFTI | grep -P "^dim2" | awk '{print $2}') 
  Z=$(fslinfo $NIFTIDIR/$ID/$NIFTI | grep -P "^dim3" | awk '{print $2}') 
  T=$(fslinfo $NIFTIDIR/$ID/$NIFTI | grep -P "^dim4" | awk '{print $2}') 
  SRS=$(echo $NIFTI | sed 's/_/\t/g' | awk '{print $4}' | grep -o -P "\d*")
  SEQ=$(echo $NIFTI | sed 's/_/\t/g' | awk '{print $3}')
  printf "${ID}\t${SEQ}\t${SRS}\t${X}\t${Y}\t${Z}\t${T}\n" >> $OUTPUT 
done 

rm $NIFTIDIR/nifti_tmp.txt
