###################
#### FUNCTIONS ####
###################
## Creates the exclusion bed
function createExclusionBed(){
  INBED=$1
  MAXFLANK=$2
  GENOME=$3
  BLACKLIST=$4
  OUTBED=$5
  
  echo "Creating an 'excl.bed' file using a ${MAXFLANK}bp flank"
  cut -f1-3 ${INBED} | bedtools slop -i - -g ${GENOME} -b ${MAXFLANK} > total_flank.bed
  echo "Adding Blacklisted regions to excluded bed..."
  cat total_flank.bed ${BLACKLIST} > total_flank.tmp
  mv total_flank.tmp total_flank.bed
  cut -f1-3 total_flank.bed | bedtools sort -i - | bedtools merge -i - > ${OUTBED}
  rm total_flank.bed
  
  echo "Exclusion file created to ${OUTBED}..."
}

## Gets number of rows equivalent to 80% of input bed file rows
function getSubsetN(){
  inBed=$1  
  frac=0.80
  trows=$(< "${inBed}" wc -l)
  nrows=$(awk -vp=${trows} -vq=${frac} 'BEGIN{printf "%.0f" ,p * q}')
  echo ${nrows}
}

## Shuffles a file, selects the number of rows inputted, and outputs the result
function shufAndSubset(){
  inBed=$1
  outBed=$2
  n=$3
  
  cut -f1-3 ${inBed} | shuf -  > shuf.tmp
  head -n ${n} shuf.tmp | bedtools sort -i - > ${outBed}
  rm shuf.tmp
}

###################
#### VARIABLES ####
###################



MOTIFBED=$1
FLANK=$2
MUTBED=$3
GENOME=$4
NUM=$5
MAXFLANK=$6
BLACKLIST=$7
DNase=$8
TF=$9


#### MAIN ####
##############
## Create a 100 line bed file with non-overlapping segments
rm ${TF}-${FLANK}_ResampCount.txt  ${TF}-${FLANK}_ShuffleCount.txt  ${TF}-${FLANK}_placeholder.bed ${TF}-${FLANK}_resampled.bed placeholder.bed

## Create the exclusion bed file
createExclusionBed ${MOTIFBED} ${MAXFLANK} ${GENOME} ${BLACKLIST} exclude.bed
createExclusionBed ${DNase} ${MAXFLANK} ${GENOME} ${DNase} include.bed

## Calculate subset size
nrows=$(getSubsetN ${MOTIFBED})
echo "Calculated the number of rows to subset to: ${nrows}"

for i in `seq 1 ${NUM}`; do
  printf "Input ${i}: shuffle, resample"
  ## Shuffle the bed file
  echo "creating ${TF}-${FLANK}_resampled.bed"
  shufAndSubset ${MOTIFBED} ${TF}-${FLANK}_resampled.bed ${nrows}
  printf ", slopping, sorting"
  bedtools slop -i ${TF}-${FLANK}_resampled.bed -g ${GENOME} -b ${FLANK} > resample_flank.tmp
  bedtools sort -i resample_flank.tmp | bedtools merge -i - > ${TF}-${FLANK}_resampled.bed
  printf ", counting"
  bedtools intersect -wa -a ${MUTBED} -b ${TF}-${FLANK}_resampled.bed | sort -u | wc -l >> ${TF}-${FLANK}_ResampCount.txt
  
  ## Randomly shuffle a butt load of intervals
  if [ ! -f "placeholder.bed" ]; then
    echo "Generating mega placeholder.bed"
    for i in $(seq 1 10); do
      cat ${TF}-${FLANK}_resampled.bed >> placeholder.bed
    done
    bedtools shuffle -i placeholder.bed -g ${GENOME} -excl exclude.bed -incl include.bed -noOverlapping > ${TF}-${FLANK}_placeholder.bed
  fi
  
  printf " || Randomly shuffling"
  shufrows=$(< "${TF}-${FLANK}_resampled.bed" wc -l)
  shufAndSubset ${TF}-${FLANK}_placeholder.bed ${TF}-${FLANK}_shuffled.bed ${shufrows}
  printf ", ${shufrows} rows, counting...\n"
  bedtools intersect -wa -a ${MUTBED} -b ${TF}-${FLANK}_shuffled.bed | sort -u | wc -l >> ${TF}-${FLANK}_ShuffleCount.txt
done

### Create csv file [Resample,Shuffle]
paste -d "," ${TF}-${FLANK}_ResampCount.txt  ${TF}-${FLANK}_ShuffleCount.txt > ${TF}-${FLANK}.csv
