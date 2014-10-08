#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

OPTIONS:
   -h           Show this message
   -o           Output directory
   
EOF
}

OUTDIR=

while getopts "ho:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
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

TMPDIR=$OUTDIR/countfiles_tmp
NOW=$(date "+%Y-%m-%d_%H%M%S")
ALLCT=$OUTDIR/motion_outliers_count_${NOW}.txt

#Print header
printf "ID\tOutlierCount\n" >> $ALLCT

#Concatenate all subj-specific outlier count files
#Replace blank entries with zeroes
cd $TMPDIR
cat *.txt | awk '
        BEGIN { FS = OFS = "\t" }
        { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1
        ' >> $ALLCT
#Clean up
cd $OUTDIR
rm -rf countfiles_tmp
