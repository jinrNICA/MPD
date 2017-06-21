#!/bin/sh

FILE_LIST=$1
OUTDIR=$2
DCAFILE=$3

cd $OUTDIR

while read FILENAME; do
	AFTERSLASHNAME=${FILENAME##*/}
	BASENAME=${AFTERSLASHNAME%.*}
	sbatch /lustre/nyx/hades/user/parfenov/real-flow/example-resolutions.sh ${FILENAME} ${OUTDIR}/${BASENAME}_res_out.root ${DCAFILE}
done <$FILE_LIST
