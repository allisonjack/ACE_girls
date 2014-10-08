#!/bin/bash

export FSLDIR="/usr/local/fsl5"
. /usr/local/fsl5/etc/fslconf/fsl.sh

export FREESURFER_DIR="/usr/local/freesurfer-5.0.0_64bit"
export FREESURFER_HOME="/usr/local/freesurfer-5.0.0_64bit"
. /usr/local/freesurfer-5.0.0_64bit/FreeSurferEnv.sh

usage()
{
cat << EOF
usage: $0 options

This script organizes and NIfTI-converts DICOM files grabbed by IDAGet.

OPTIONS:
  -h  Show this message
  -N  Full path to your raw NIfTI directory (will be created by this script if it doesn't already exist)
  -D  Full path to DICOM directory (This is where you had IDAGet send your files)
EOF
}

NIFTIDIR=
DCMDIR=

#Read arguments entered by user
while getopts "hN:D:" OPTION
do
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

#--Create a directory in which to put the nifti files we generate--#
mkdir -p $NIFTIDIR

#--Move into the directory with all images grabbed by IDAGet--#
cd $DCMDIR

#--Figure out DICOM directory name--#
IDADIR=$(echo $DCMDIR | grep -o "output-[0-9]*$")

#--Get name of the first DICOM image in the folder and store it in a variable--#
DCMINPUT=$(ls | head -1)

#--Get subject ID--#
SUBJ=$(echo $DCMINPUT | grep -o -P "(?<=ACE_).*(?=_MR)")

#--Specify/abbreviate names for newly-generated NIfTI files based on DICOM file name--#
NIFTI=$(ls | head -1 | sed 's/_*br_raw_[0-9]*_[0-9]*.dcm//g' \
| sed 's/_MR_BOLD_-//g' \
| sed 's/_*[0-9]*$//g' \
| sed 's/_MR_T2-AX__In-plane_Structural/_T2Ax/g' \
| sed 's/_MR_localizer/_localizer/g' \
| sed 's/_MR_trufi_localizer_3-plane/_localizer/g' \
| sed 's/_MR_MPRAGE_-_BWM/_MPRAGE/g')

#--Create subject subfolders w/in main NIfTI directory--#
SUBJDIR=$NIFTIDIR/$SUBJ
mkdir -p $SUBJDIR

IDAOUT=$NIFTIDIR/$SUBJ/$IDADIR
mkdir -p $IDAOUT

#--Convert to NIfTI files and save information about the conversion in a text file--#
/usr/local/freesurfer-5.0.0_64bit/bin/mri_convert -it siemens_dicom -ot nii $DCMINPUT $IDAOUT/${NIFTI}.nii.gz >> $IDAOUT/info_${NIFTI}.txt

#--Grab series number from info file--#
SRS=$(grep "SeriesNo" $IDAOUT/info_${NIFTI}.txt | awk '{print $2}')

#--Add series number to nifti filename in case there are two runs of the same sequence--#
NEWNIFTI=${NIFTI}_srs${SRS}.nii.gz
mv $IDAOUT/${NIFTI}.nii.gz $IDAOUT/$NEWNIFTI

NEWINFO=info_${NIFTI}_srs${SRS}.txt
mv $IDAOUT/info_${NIFTI}.txt $IDAOUT/$NEWINFO

#--Pull files up into main subject directory--#
mv $IDAOUT/$NEWNIFTI $SUBJDIR/$NEWNIFTI
mv $IDAOUT/$NEWINFO $SUBJDIR/$NEWINFO

#--Reorient images to match orientation of standard space template--#
#(This is not a registration step, only applies 90, 180 or 270 degree rotations about
#the different axes as necessary to get labels in same position as standard template)
/usr/local/fsl5/bin/fslreorient2std $SUBJDIR/$NEWNIFTI $SUBJDIR/$NEWNIFTI

echo $DCMDIR

#Clean up
rm -rf $IDAOUT
