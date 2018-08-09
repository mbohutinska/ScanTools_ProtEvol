#!/bin/bash -e
#PBS -N snpEff
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=2:mem=24gb:scratch_local=1gb
#PBS -m abe
#PBS -j oe

module add jdk-8
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
DATADIR="/storage/plzen1/home/holcovam/ScanTools/all300"
cp -r /storage/plzen1/home/holcovam/programs/snpEff $SCRATCHDIR || exit 1
cp $DATADIR/All.NOfiltr.scaff* $SCRATCHDIR || exit 1
cd $SCRATCHDIR || exit 2
echo data loaded at `date`

#TODO
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff1.vcf.gz > All.NOfiltr.ann.scaff1.vcf.gz
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff2.vcf.gz > All.NOfiltr.ann.scaff2.vcf.gz
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff3.vcf.gz > All.NOfiltr.ann.scaff3.vcf.gz
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff4.vcf.gz > All.NOfiltr.ann.scaff4.vcf.gz
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff5.vcf.gz > All.NOfiltr.ann.scaff5.vcf.gz
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff6.vcf.gz > All.NOfiltr.ann.scaff6.vcf.gz
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff7.vcf.gz > All.NOfiltr.ann.scaff7.vcf.gz
java -Xmx24g -jar snpEff/snpEff.jar LyV2 All.NOfiltr.scaff8.vcf.gz > All.NOfiltr.ann.scaff8.vcf.gz

rm All.NOfiltr.scaff*
rm -r snpEff
cp $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
