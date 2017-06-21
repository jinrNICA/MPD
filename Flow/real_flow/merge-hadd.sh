#!/bin/sh

INDIR=$1
NFILES=$2

cd ${INDIR}/MERGE/
FILE_LIST=${INDIR}/MERGE/file_list
find ${INDIR}/ -iname "*.root" > ${FILE_LIST}

split -l ${NFILES} -d ${FILE_LIST} hadd_files_
find $PWD -name "hadd_files_*" > ALLHADD
 
while read FILENAME; do
	TASKNAME=${FILENAME}.sh
	cp /lustre/nyx/hades/user/svintsov/mpdroot_batch/real-flow/hadd-example.sh ${TASKNAME}
	sbatch ${TASKNAME} ${FILENAME}
done <ALLHADD
