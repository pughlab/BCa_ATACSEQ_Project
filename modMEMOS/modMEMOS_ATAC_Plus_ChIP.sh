#!/bin/bash

#SBATCH -t 72:00:00
#SBATCH --mem=20G
#SBATCH -J MEMOS
#SBATCH -p all
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out

OUTDIR=""
WDIR="" #working directory


BCADIR="" #full path to the ATAC Catalogue file
ATAC=$BCADIR/Catalogue_ATAC_BCa.bed #ATAC Catalogue


MUTVCF="" # full path to mutation file


#full path to fasta reference file
FASTA=/mnt/work1/data/genomes/human/hg19/iGenomes/Sequence/WholeGenomeFasta/genome.fa

# name of the TF bed ChIP-seq file
TFBED=$1
TFNAME=$2

###### Overlap ChIP-seq from TF of interest with ATAC-seq catalogue 

## Load module
module load bedtools/2.23.0

BASEATAC=$(basename "$ATAC")
BASETF=$(basename "$TFBED")
bedtools intersect \
-a $TFBED \
-b $ATAC \
-wa \
> "$OUTDIR/${BASETF%.*}-overlapped-${BASEATAC}"

###### convert bed files to fasta 
bedtools getfasta \
-fi $FASTA \
-bed "$OUTDIR/${BASETF%.*}-overlapped-${BASEATAC}" \
-fo "$OUTDIR/${BASETF%.*}-overlapped-${BASEATAC%.*}.fa"


###### Run modMEMOS-wrapper script
for j in {0,10,20,30,40,50,100,200,300,400,500,1000}; do

sbatch  $WDIR/modMEMOS-wrapper.sh $TFNAME "$OUTDIR/${BASETF%.*}-overlapped-${BASEATAC%.*}.fa" $MUTVCF 100 "$OUTDIR/${BASETF%.*}-bg-ATAC-${j}flank" \
0.0001 ${j} $ATAC

done
