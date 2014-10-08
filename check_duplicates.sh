#!/bin/bash 
usage()
{
cat << EOF
usage: $0 options

This script checks for imaging sequences that were run multiple times.

OPTIONS:
  -h  Show this message
  -N  Full path to your raw NIfTI directory (will be created by this script if it doesn't already exist)
  -D  Full path to DICOM directory (This is where you had IDAGet send your files)
   
EOF
}

NIFTIDIR=
DCMDIR=

#Read arguments entered by user
while getopts "hN:D:" OPTION; do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    N)
      NIFTIDIR=$OPTARG
      ;;
    D)
      DCMDIR=$OPTARG
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

#Create a timestamp
NOW=$(date +"%m-%d-%Y_%H%M")

#Extract list of subject IDs
cd $DCMDIR
ls | grep "IDAGet_[0-9]*.output" | while read -r line; do
  DCMINPUT=$(ls $line | head -1)
  echo $DCMINPUT | grep -o -P "(?<=ACE_).*(?=_MR)" | while read SUBJ; do
    cd $NIFTIDIR

    #Set up to create a file with ID #s
    printf "$SUBJ\n" >> IDs_main_tmp.txt

    #List duplicates and write to temporary file
    find . -name "*.nii.gz" | sed 's/\(.*\)_srs[0-9]*.nii.gz/\1/g' | sort | uniq -c -d >> "duplicates_${NOW}_tmp.txt"
    done
  done

cd $NIFTIDIR

#Extract number of duplications for each item
sed 's/^ *//g' "duplicates_${NOW}_tmp.txt" | awk 'BEGIN {print "No.Runs";} {print $1}' >> "dups_tmp.txt"

#Extract name of sequence duplicated
grep -o "[a-zA-Z]*_*[a-zA-Z0-9]*$" "duplicates_${NOW}_tmp.txt" | sed 's/^_//g' | sed '1i Sequence' >> "seq_tmp.txt"

#Extract ID number for that row
awk 'BEGIN {FS = "/"; print "ID"}; {print $2}' "duplicates_${NOW}_tmp.txt" >> "ID_tmp.txt"

#Print a file with info on duplicates
paste ID_tmp.txt dups_tmp.txt seq_tmp.txt >> "duplicates_${NOW}_tmp2.txt"
grep -v "localizer" "duplicates_${NOW}_tmp2.txt" >> "duplicates_${NOW}_tmp3.txt"

awk 'NR == 1; NR > 1 {print $0 | "sort -nr"}' "duplicates_${NOW}_tmp3.txt" | uniq >> "duplicates_${NOW}.txt"
#sort duplicates_${NOW}_tmp3.txt | uniq >> "duplicates_${NOW}.txt"

#Print a file listing the IDs that were processed
sort IDs_main_tmp.txt | uniq >> "IDs_${NOW}.txt"

#Clean up
rm *_tmp*.txt
