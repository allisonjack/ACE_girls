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
   -h   Show this message
   -F   Non motion corrected input file
   -I   List of ID numbers
   -V   Number of volumes 
   -O   Output directory
   -N   List of volume numbers starting at 001
   -M   List of volume numbers starting at 000
   
EOF
}

while getopts "hF:I:V:O:N:M:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    F)
      FX=$OPTARG
      ;;
    I)
      ID=$OPTARG
      ;;
    V)
      VOLS=$OPTARG
      ;;
    O)
      OUTPUT=$OPTARG
      ;;
    N)
      NLIST=$OPTARG
      ;;
    M)
      NMINUSLIST=$OPTARG
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

mkdir -p $DIFFDIR
mkdir -p $OUTPUT

#---Input----#

#Get name of MCFLIRT matrix directory
MATDIR=$(echo $FX | sed 's/.nii.gz/_mcf.mat/')

#Extract the middle volume from the raw functional image
HALFVOLS=$(($VOLS/2))
/usr/local/fsl5/bin/fslroi $FX $OUTPUT/example_func $HALFVOLS 1

REFVOL=$OUTPUT/example_func
IDENT=/usr/local/fsl5/etc/flirtsch/ident.mat

#---Output---#

ABS=$OUTPUT/${ID}_rmsdiff_abs.txt
REL=$OUTPUT/${ID}_rmsdiff_rel.txt
RMSTMP=$OUTPUT/${ID}_rms_tmp.txt
RMSTOT=$OUTPUT/${ID}_rms_tot.txt

#---Calculate RMS movement---#

#Absolute RMS movement
for ((i = 0; i < $VOLS; ++i)); do
        printf -v num '%04d' $i
        /usr/local/fsl5/bin/rmsdiff $MATDIR/MAT_${num} $IDENT $REFVOL >> $ABS
done

#Relative RMS movement
while IFS= read -r N && IFS= read -r NMINUS <&3 ; do
        /usr/local/fsl5/bin/rmsdiff $MATDIR/MAT_${N} $MATDIR/MAT_${NMINUS} $REFVOL >> $REL
done <$NLIST 3<$NMINUSLIST

#-----Paste together files and label for R----#

paste $ABS $REL >> $RMSTMP
awk '{print "'$ID'", $0}' $RMSTMP >> $RMSTOT
sed -i '1iID\tABS\tREL' $RMSTOT

#-----Clean up-----#

rm $ABS
rm $REL
rm $RMSTMP
