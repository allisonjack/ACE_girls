#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script updates your prestats fsf file.

OPTIONS:
   -h   Show this message
   -f   Input design.fsf file
   -u   Output design.fsf file
   -d   4D data
   -o   Output directory name
   -v   Total volumes
   -i   Initial structural image
   -m   Main structural image
   -s   Standard space
   
EOF
}

INPUTFSF=
UPDATEFSF=
DATA=
OUTPUT=
VOLUMES=
INITIALSX=
MAINSX=
STANDARD=

while getopts "hf:u:d:o:v:i:m:s:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    f)
      INPUTFSF=$OPTARG
      ;;
    u)
      UPDATEFSF=$OPTARG
      ;;
    d)
      DATA=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    v)
      VOLUMES=$OPTARG
      ;;
    i)
      INITIALSX=$OPTARG
      ;;
    m)
      MAINSX=$OPTARG
      ;;
    s)
      STANDARD=$OPTARG
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

#--Replaces keywords in template design file with information stored in variables--#

for i in $INPUTFSF; do
  sed -e 's@DATA@'$DATA'@g' \
  -e 's@OUTPUT@'$OUTPUT'@g' \
  -e 's@999@'$VOLUMES'@g' \
  -e 's@INITIALSX@'$INITIALSX'@g'\
  -e 's@MAINSX@'$MAINSX'@g' \
  -e 's@STANDARD@'$STANDARD'@g' <$i> $UPDATEFSF
done
