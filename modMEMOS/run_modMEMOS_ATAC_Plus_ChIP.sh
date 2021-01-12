#!/bin/bash
#$ -cwd
#$ -S /bin/bash


while read -r TFBED TFNAME; do

echo $TFBED  $TFNAME
./modMEMOS_ATAC_Plus_ChIP.sh $TFBED $TFNAME

done < TFnames.txt
