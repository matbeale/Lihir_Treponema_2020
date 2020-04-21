#!/bin/bash


# Set default options for flags (will be overwritten if specified)

suffix=fastq.gz
threads=4
queue=normal
memory=8
maxbatchsize=500
resultsbase=$PWD
mydate=$(date +%Y-%m-%d)
filetype=fastq

echo "    `basename $0` [ options ] <full-path-to-folder-with-fastqs>" 


input_path=$1

# Capture command line options and positional arguments

# Introduce date into results path
resultsdir=$resultsbase/23s-comp-mapping_$mydate
MEMORYfull=$(($memory * 1000))



mkdir $resultsdir/
mkdir $resultsdir/logs/

filelist=$(ls -1 $input_path/*_1.$suffix )
filecount=$(ls -1 $input_path/*_1.$suffix | wc -l)
ls -1 $input_path/*_1.$suffix | perl -pe 's/^.+\///g' | perl -pe "s/\_1\.$suffix//g" > $resultsdir/full.seqlist

# echo "reading in files : $filelist"
echo "There are $filecount input files"


# split file list into groups of 500 with prefix 'batch'
split -a 3 -d -l $maxbatchsize $resultsdir/full.seqlist $resultsdir/batch
# Rename splits on a 1-1000 scheme (not 0-999)
for f in $resultsdir/batch* ; do mv ${f} $(echo ${f} | perl -ne 'm/batch(0*)(\d+)/g; print "compmap.Batch.",$2+1,"\n";'); done
mv compmap.Batch.* $resultsdir/

# Make list of splits and get counts into a variable
grouplist=$(ls -1 $resultsdir/compmap.Batch*)
ls -1 $resultsdir/compmap.Batch* | perl -pe 's/^.+\///g' > $resultsdir/full.grouplist

groupcount=$(cat $resultsdir/full.grouplist | wc -l)
echo "Splitting analysis into $groupcount groups of $maxbatchsize genomes"


i=1 ; while [ $i -le $groupcount ]; do cat $resultsdir/compmap.Batch.$i | wc -l > $resultsdir/compmap.Batch.$i\.batchsize ; i=$(($i+1)) ; done


# Make script files (in batches)
# map-to-contigs.sh <sampleid> <resultsdir> <reference.fasta> <read1> <read2> <sitelist.tsv>

i=1; while [ $i -le $groupcount ]; do
	while read seq ; do echo "~mb29/bsub_scripts/competitive_mapping_Treponema23S.sh $seq $resultsdir/$seq/ ~mb29/references/Competitive-mapping/Treponema_23s_copy1_+_Strepdysgalactiae.16s+23s+5s.fa $input_path/$seq\_1.$suffix $input_path/$seq\_2.$suffix ~mb29/references/Competitive-mapping/23s_sites.Treponema-comp-mapping.tsv" >> $resultsdir/compmap.Batch.$i\.scriptlist.sh ; done < $resultsdir/compmap.Batch.$i ; i=$(($i+1)) ; done 



# submit jobs using James Hadfield's run_array script (makes each line of a script or list of scripts into a separate entry in the array)
i=1 ; while [ $i -le $groupcount ]; do
current_batchsize=$(cat $resultsdir/compmap.Batch.$i\.batchsize)
#echo "testing $current_batchsize"

compmap_job=$(bsub -J compmap.batch$i\[1-$current_batchsize]%50 -o $resultsdir/logs/compmap$i\.%J.%I.o -e $resultsdir/logs/compmap$i\.%J.%I.e -q $queue -n $threads -M$MEMORYfull -R "span [hosts=1] select[mem>$MEMORYfull] rusage[mem=$MEMORYfull]" ~mb29/scripts/working/run_array $resultsdir/compmap.Batch.$i\.scriptlist.sh )

compmap_job_jobID=$(echo $compmap_job | sed 's/[^0-9]*//g')
echo $compmap_job | sed 's/[^0-9]*//g' >> $resultsdir/submission.jobids.txt
echo "competitive mapping batch job submitted with ID $compmap_job_jobID"

i=$(($i+1))
done



# Merge data from loop (not quite running yet)
#jobcompletion=$(cat 23s-comp-mapping_2020-03-02/submission.jobids.txt | perl -pe 's/\n/\&/g' | perl -pe 's/\&$//g')
#echo "this is the jobcompletion output $jobcompletion "

#compmap_merge=$(bsub -w "ended($jobcompletion)" -J compmap_merge -o $resultsdir/logs/compmap_merge.%J.o -e $resultsdir/logs/compmap_merge.%J.e -q normal -n 1 -M$MEMORYfull -R "span [hosts=1] select[mem>$MEMORYfull] rusage[mem=$MEMORYfull]" "while read line ; do grep -v $resultsdir/$line/$line\.bcf.called-sites.tsv >> $resultsdir/23S.collated.called-sites.$mydate\.tsv ; done < $resultsdir/full.seqlist")

# capture and echo job id to console
#compmap_merge_jobID=$(echo $compmap_merge | sed 's/[^0-9]*//g')
# echo "comparative mapping merge job submitted with ID $compmap_merge_jobID (will wait for completion of batch)"



