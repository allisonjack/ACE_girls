#!/bin/bash

export FSLDIR="/usr/local/fsl5"
. /usr/local/fsl5/etc/fslconf/fsl.sh

usage()
{
cat << EOF
usage: $0 options

This script uses fsl_motion_outliers to create a matrix of corrupt volumes

OPTIONS:
  -h  Show this message
  -i  Non-motion-corrected functional image
  -S  List of ID numbers, used to label output
  -o  Output directory

EOF
}

INPUT=
SUBJ=
OUTDIR=

while getopts "hi:S:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    i)
      INPUT=$OPTARG
      ;;
    S)
      SUBJ=$OPTARG
      ;;
    o)
      OUTDIR=$OPTARG
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

SUBJDIR=$OUTDIR/$SUBJ
mkdir -p $SUBJDIR
TMPDIR=$OUTDIR/countfiles_tmp
mkdir -p $TMPDIR

SUBJCT=$TMPDIR/motion_outliers_count_${SUBJ}.txt

#---Run outlier detection---#
OUTPUT=$SUBJDIR/motion_outliers_${SUBJ}.txt
DVARS=$SUBJDIR/dvars_values_${SUBJ}.txt

/usr/local/fsl5/bin/fsl_motion_outliers -i $INPUT -o $OUTPUT --dvars -s $DVARS

#---Create an empty outliers file if no outliers were detected---#
touch $OUTPUT

#---Count # of columns (thus, # outliers) in matrix fsl produced---#
MOCOUNT=$(awk '{print NF}' $OUTPUT | sort -nu | tail -n 1)
