#!/bin/bash

export FSLDIR="/usr/local/fsl5"
. /usr/local/fsl5/etc/fslconf/fsl.sh

usage()
{
cat << EOF
usage: $0 options

This script uses the FSL rmsdiff tool to calculate:
    1) absolute motion relative to refvol (in mm)
    2) motion relative to previous volume (in mm)

OPTIONS:
   -h  Show this message
   -F  The FEAT directory from which to calculate rms movement 
   -S  List of ID numbers, used to label output
   -L  One column text file listing timepoints 1 to N
   -M  One column text file listing timepoints 0 to N-1
   -v  Integer value: total number of volumes in the paradigm
   -o  Output directory
   
EOF
}

FEATDIR=
SUBJ=
NLIST=
NMINUSLIST=
VOLUMES=
OUTPUT=

while getopts "hF:S:L:M:v:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    F)
      FEATDIR=$OPTARG
      ;;
    S)
      SUBJ=$OPTARG
      ;;
    L)
      NLIST=$OPTARG
      ;;
    M)
      NMINUSLIST=$OPTARG
      ;;
    v)
      VOLUMES=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

#---Input----#
MATDIR=$FEATDIR/mc/prefiltered_func_data_mcf.mat
REGDIR=$FEATDIR/reg
REFVOL=$REGDIR/example_func.nii.gz

IDENT=/usr/local/fsl5/etc/flirtsch/ident.mat

#---Output---#

mkdir -p $OUTPUT

ABS=$OUTPUT/${SUBJ}_rmsdiff_abs.txt
REL=$OUTPUT/${SUBJ}_rmsdiff_rel.txt
RMSTMP=$OUTPUT/${SUBJ}_rms_tmp.txt
RMSTOT=$OUTPUT/${SUBJ}_rms_tot.txt

#---Calculate RMS movement---#

for ((i = 0; i < $VOLUMES; ++i)); do
        printf -v num '%04d' $i
        /usr/local/fsl5/bin/rmsdiff $MATDIR/MAT_${num} $IDENT $REFVOL >> $ABS
done

while IFS= read -r N && IFS= read -r NMINUS <&3 ; do
        /usr/local/fsl5/bin/rmsdiff $MATDIR/MAT_${N} $MATDIR/MAT_${NMINUS} $REFVOL >> $REL
done <$NLIST 3<$NMINUSLIST

#-----Paste together files and label for R----#

paste $ABS $REL >> $RMSTMP
awk '{print "'$SUBJ'", $0}' $RMSTMP >> $RMSTOT
sed -i '1iID\tABS\tREL' $RMSTOT

#-----Clean up-----#

rm $ABS
rm $REL
rm $RMSTMP

#----Creates text files that list N and N-1----#
#----(Only needs to be run once)---#

#       for ((i = 1; i < 154; ++i)); do
#               printf -v num '%04d' $i
#               echo $num >> $NLIST
#       done
# 
#       for ((i = 0; i < 153; ++i)); do
#               printf -v num '%04d' $i
#               echo $num >> $NMINUSLIST
#       done
