#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

Prepare confound matrix by pasting standard motion parameters and motion outliers together

OPTIONS:
   -h   Show this message
   -S   List of ID numbers, used to label output
   -F   Prestatistics feat directory
   -o   Output directory
   
EOF
}

SUBJ=
FEATDIR=
OUTDIR=

while getopts "hi:S:F:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    S)
      SUBJ=$OPTARG
      ;;
    F)
      FEATDIR=$OPTARG
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
OUTLIERS=$SUBJDIR/motion_outliers_${SUBJ}.txt
PARFILE=$FEATDIR/mc/prefiltered_func_data_mcf.par

#--Paste together motion outliers matrix and motion parameters into one file---#
paste $OUTLIERS $PARFILE >> $SUBJDIR/confound_evs_${SUBJ}.txt
