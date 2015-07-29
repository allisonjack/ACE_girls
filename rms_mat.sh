#!/bin/bash

#----Creates text files that list N and N-1----#
#----(Only needs to be run once)---#

NLIST=/ifshome/USER/ACE/TEMPLATES/mat_N.txt
NMINUSLIST=/ifshome/USER/ACE/TEMPLATES/mat_N-1.txt

#The number of volumes in your paradigm
VOLS=172
#That minus one
VOLSminus1=171

#Say you have 172 total volumes. The mat_N file will run from 1 to 171.
#The mat_N-1 file will run from 0 to 170. This is because FSL codes your
#first volume as 0 (NOT 1).

for ((i = 1; i < $VOLS; ++i)); do
  printf -v num '%04d' $i
  echo $num >> $NLIST
done

for ((i = 0; i < $VOLSminus1; ++i)); do
  printf -v num '%04d' $i
  echo $num >> $NMINUSLIST
done
