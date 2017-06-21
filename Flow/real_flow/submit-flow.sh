#!/bin/sh

FILE_LIST=$1
OUTDIR=$2
RES_FIT_FILE=$3
DCAFILE=$4

cd $OUTDIR

while read FILENAME; do
	AFTERSLASHNAME=${FILENAME##*/}
	BASENAME=${AFTERSLASHNAME%.*}
	sbatch /lustre/nyx/hades/user/parfenov/real-flow/example-flow.sh ${FILENAME} ${OUTDIR}/${BASENAME}_flow_out.root ${RES_FIT_FILE} ${DCAFILE}
done <$FILE_LIST
