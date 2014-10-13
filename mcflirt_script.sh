#!/bin/bash

export FSLDIR="/usr/local/fsl5"
. /usr/local/fsl5/etc/fslconf/fsl.sh

usage()
{
cat << EOF
usage: $0 options

OPTIONS:
   -h   Show this message
   -i   Non-motion-corrected input file 


EOF
}

while getopts "hi:" OPTION
do
  case $OPTION in
   h)
    usage
    exit 1
    ;;
  i)
    INPUT=$OPTARG
    ;;
  ?)
    usage
    exit
    ;;
  esac
done

/usr/local/fsl5/bin/mcflirt -in $INPUT -mats -plots
