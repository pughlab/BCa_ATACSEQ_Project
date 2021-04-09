#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --mem=20G
#SBATCH -J MEMOSv03
#SBATCH -p all
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out

# 1=TF
# 2=FASTA for regions to be scanned
# 3=Mutation file (vcf)
WDIR="/mnt/work1/users/pughlab/projects/KOMBAT/atacseq_analysis/Sorted_ATACSeq_hg19/MEMOS/MEMOS-wrapper"
TFDIR=/mnt/work1/users/pughlab/projects/KOMBAT/atacseq_analysis/Sorted_ATACSeq_hg19/MEMOS/TF-pfms
HG19=$WDIR/hg19.genome
###
module load python/2.7
module load R/3.2.2
module load bedtools/2.23.0

###
TF=$1
MUTATION=$3
NUM=$4
OUTDIR=$5
PVAL=$6
FLANK=$7
FASTA=$OUTDIR/tmp.fa
BLACKLIST='/mnt/work1/users/pughlab/projects/KOMBAT/atacseq_analysis/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed'
INCLUDE='/mnt/work1/users/pughlab/projects/KOMBAT/atacseq_analysis/Sorted_ATACSeq_hg19/Chip_ATAC/ATAC-only-ICGCEU-v04/bedsymlinks/ALLMotif_ATAC.bed'
#INCLUDE='/mnt/work1/users/pughlab/projects/KOMBAT/atacseq_analysis/Sorted_ATACSeq_hg19/Chip_ATAC/ENCODE-Remap_ChIP_TF/hg19/MCF7/Merged_Peaks/ALLChIP_Catalogue.bed'
###
mkdir -p $OUTDIR
###
echo "creating a bed file with required columns from the mutation vcf file ..."
grep -F '#' $MUTATION > $OUTDIR/vcfheader.tmp
grep -Fv '#' $MUTATION > $OUTDIR/vcfbody.tmp
awk -F "\t" '{print $1 "\t" $2 "\t" $2 "\t" $4 ">" $5}' $OUTDIR/vcfbody.tmp > $OUTDIR/mutations.bed
MUT=$OUTDIR/mutations.bed

echo "reformat the fasta file; creating tmp.fa ..."
cp $2 "${OUTDIR}/tmp.fa"
$WDIR/formatFA.sh $FASTA
#FASTA=/mnt/work1/users/lupiengroup/People/parisa/MOODSanalysis/20PCa_H3K27acpeaks_extended500bp_merged.fa
############################# MOODS #############################
echo "started running moods for $TF ..."
for pfm in $TFDIR/$TF/*.pfm
do
	pfmName="${pfm##*/}"
	echo "motif name = $pfmName ..."
	echo "running MOODS with pvalue = $PVAL ..."
	echo "parameters:  -m $pfm -s $FASTA -p $PVAL > $OUTDIR/${pfmName%.*}-${PVAL}.out "
	python /mnt/work1/software/MOODS/1.9.2/scripts/moods_dna.py -m $pfm -s $FASTA -p $PVAL > $OUTDIR/${pfmName%.*}-${PVAL}.out
	n="$(awk -F' ' '{print NF; exit}' $pfm)"
	n=$(($n-1))
	awk -v x=$n -F '[, ]' '{print $1"\t"$2+$5"\t"$2+$5+x"\t"$6 "\t" $7 "\t" $8}' $OUTDIR/${pfmName%.*}-${PVAL}.out > $OUTDIR/${pfmName%.*}-${PVAL}.bed

####################### Permutation Test #########################
	echo "started running permutation test for $TF ..."
	mkdir $OUTDIR/${TF}_shuffletest_${PVAL}
	cd $OUTDIR/${TF}_shuffletest_${PVAL}
	$WDIR/shuffleTest.sh  $OUTDIR/${pfmName%.*}-${PVAL}.bed  $FLANK  $MUT  $HG19  $NUM 1000  ${BLACKLIST}  ${INCLUDE} ${TF}
	echo "done with test!"
done
