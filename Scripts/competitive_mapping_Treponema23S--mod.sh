#!/bin/bash

# Usage : competitive_mapping_Treponema23S.sh <sampleid> <resultsdir> <reference.fasta> <read1> <read2> <sitelist.tsv>


if [ $# -ne 6 ]; then
    echo "Usage : $0 <sampleid> <resultsdir> <reference.fasta> <read1> <read2> <sitelist.tsv>"
    exit 1
fi


sampleid=$1
resultsdir=$2
reference=$3
read1=$4
read2=$5
sitelist=$6

threads=4
mem=8

module load bwa/0.7.17=pl5.22.0_2
module load samtools/1.6--h244ad75_4
module load picard/2.21.4
module load gatk/3.7.0
module load bcftools/1.5--3


samtools=samtools
bwa=bwa
PICARDPATH=picard
GATK=gatk
BCFTOOLS=bcftools

mkdir -pv $resultsdir/

# Map to reference and filter out unmapped and improperly paired reads
$bwa index $reference

# remove insistence on pairing
$bwa mem -t $threads -M $reference $read1 $read2 | $samtools view -@ $threads -hb -F4 -q 25 - | $samtools sort -@ $threads - > $resultsdir/$sampleid\.map.bam


# Mark duplicates (and remove) - note for amplicon sequencing this is inappropriate, so has been commented out below
$PICARDPATH AddOrReplaceReadGroups TMP_DIR=$resultsdir/tmpdir/ INPUT=$resultsdir/$sampleid\.map.bam OUTPUT=$resultsdir/$sampleid\.fixed.bam RGLB=SureSelect RGPL=Illumina RGPU=$sampleid RGSM=$sampleid RGCN=WTSI CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT
#$PICARDPATH MarkDuplicates TMP_DIR=$resultsdir/tmpdir/ INPUT=$resultsdir/$sampleid\.fixed.bam OUTPUT=$resultsdir/$sampleid\.marked.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=$resultsdir/$sampleid\.fixed.picard-metrics


# for deep amplicon seq, don't want to mark dups
cp $resultsdir/$sampleid\.fixed.bam $resultsdir/$sampleid\.marked.bam


$samtools index $resultsdir/$sampleid\.marked.bam

cp $reference $resultsdir/tmpref.fa
$samtools faidx $resultsdir/tmpref.fa

# Perform indel realignment
$PICARDPATH CreateSequenceDictionary R= $resultsdir/tmpref.fa O= $resultsdir/tmpref.dict

$GATK -I $resultsdir/$sampleid\.marked.bam -R $resultsdir/tmpref.fa -T RealignerTargetCreator -o $resultsdir/tmp.intervals
$GATK -I $resultsdir/$sampleid\.marked.bam  -R $resultsdir/tmpref.fa -T IndelRealigner --filter_bases_not_stored -targetIntervals $resultsdir/tmp.intervals -o $resultsdir/$sampleid\.indelr.bam

$samtools sort -@ $threads $resultsdir/$sampleid\.indelr.bam > $resultsdir/$sampleid\.indelr.sort.bam
$samtools index $resultsdir/$sampleid\.indelr.sort.bam

# cleanup
rm $resultsdir/$sampleid\.map.bam*
rm $resultsdir/$sampleid\.fixed.ba*
rm $resultsdir/$sampleid\.marked.ba*
rm $resultsdir/$sampleid\.indelr.ba*


# use mpileup
# -A, --count-orphans     do not discard anomalous read pairs (left this off, as don't want to count them here)
$samtools mpileup -t DP,DP4 -C 50 -L 1000 -d 1000 -m 5 -x -ugf $reference $resultsdir/$sampleid\.indelr.sort.bam > $resultsdir/$sampleid\.mpileup

# Call variants and all sites
echo "$sampleid 1" > $resultsdir/$sampleid\.ploidy

$BCFTOOLS call -P 0.001 -O b -A -M -S $resultsdir/$sampleid\.ploidy -c $resultsdir/$sampleid\.mpileup > $resultsdir/$sampleid\.bcf
$BCFTOOLS index $resultsdir/$sampleid\.bcf

$BCFTOOLS call -P 0.001 -O b -A -M -v -S $resultsdir/$sampleid\.ploidy -c $resultsdir/$sampleid\.mpileup > $resultsdir/$sampleid\.variant.bcf
$BCFTOOLS index $resultsdir/$sampleid\.variant.bcf 

# run analysis scripts
source /lustre/scratch118/infgen/team216/mb29/tools/miniconda3/etc/profile.d/conda.sh
conda activate python2.7
#~/scripts/bcf-to-minorvars_v1.2.py -r Nichols_23s_228816-237492 -i $resultsdir/$sampleid\.bcf -t 10
#~/scripts/Bam2_Read_Phased_Haplotypes_v0.4.py -p $resultsdir/ -b $sampleid\.indelr.sort.bam -w 6370 -o $resultsdir/$sampleid\.


# Extract key sites for resistance testing
myrefheader=$(grep '>' $reference | head -1 | perl -pe 's/\>//g')
~mb29/scripts/bcf-summarise-specific-sites.py -i $resultsdir/$sampleid\.bcf -v $sitelist -r $myrefheader -m 5 -p 5 -t 10

conda deactivate


# Cleanup some more
rm $resultsdir/tmpref*
rm -r $resultsdir/tmpdir/


