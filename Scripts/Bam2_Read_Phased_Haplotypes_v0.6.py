#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Loop to do multiple samples AND multiple windows

import sys
import subprocess
import os
import pysam
import numpy as np
import scipy
import pandas as pd
from Bio import SeqIO
import re
from optparse import OptionParser
from optparse import OptionGroup



print "\nWritten by Mat Beale (mathew.beale@sanger.ac.uk), Wellcome Trust Sanger Institute, October 2017\n"
print "Takes a single or list of reference-mapped bam files and positional windows and extracts read-phased haplotypes, outputing haplotype counts and sequence in .fasta, and a list of reads supporting each haplotype. Haplotypes are assigned and counted by sample (SHap) and over all bam files as a group (GHap)\n"

#######
# Define options

usage = "Usage: %prog [options]"
parser = OptionParser(usage=usage)

group = OptionGroup(parser, "General Options")
group.add_option("-p", "--path",action="store", type="string", dest="inpath",help="Specify path of .bam files [%default]", default='./')
group.add_option("-o", "--outprefix",action="store", type="string", dest="outprefix",help="Specify prefix for naming output files [%default]", default='')
parser.add_option_group(group)

group = OptionGroup(parser, "Specify Bam files to Analyse","Must specify either a bam file or a list")
group.add_option("-b", "--bam",action="store", type="string", dest="singlebam",help="Specify filename for a single bam file")
group.add_option("-B", "--bamlist",action="store", type="string", dest="bamlist",help="Specify filename of line-separated list of bam files")
parser.add_option_group(group)

group = OptionGroup(parser, "Specify Windows to Analyse","Must specify either a single number or a list")
group.add_option("-w", "--window",action="store", type="string", dest="singlewindow",help="Specify start position of window (relative to reference) to analyse")
group.add_option("-W", "--windowlist",action="store", type="string", dest="windowlist",help="Specify filename of line separated list of windows to test (start position)")
group.add_option("-A","--windowauto",action="store", type="string", dest="windowauto",help="Specify automated sliding window generation along genome, incrementing start of each window by specified value (WARNING: may take a LONG time and use a lot of memory)")
parser.add_option_group(group)

group = OptionGroup(parser, "Haplotype Calling Parameters","Allows you to adjust the parameters used for assigning haplotypes")
group.add_option("-s","--windowsize",action="store", type="string",dest="windowsize",help="Specify size of haplotype windows to extract [%default]", default=75)
group.add_option("-c","--mincov",action="store", type="string",dest="mincov",help="Specify minimum count to include haplotype [%default]", default=2)
group.add_option("--baseq",action="store", type="string",dest="baseq",help="Specify minimum base quality to include read [%default]", default=25)
group.add_option("--mapq",action="store", type="string",dest="mapq",help="Specify minimum map quality to include read [%default]", default=25)
group.add_option("--readq",action="store", type="string",dest="readq",help="Specify minimum read quality to include read [%default]", default=25)
parser.add_option_group(group)

(options, args) = parser.parse_args()

########
# Check the mandatory options have been specified

if not (options.bamlist or options.singlebam):   # if bam files are not given
    parser.error('No bam files were specified - use  (-b or -B)')
if not (options.singlewindow or options.windowlist or options.windowauto):   # if windows are not given
    parser.error('No windows were specified - use  -w, -W or -WA')


########
# Capture Options and make sense of them

bampath = options.inpath

myseqlist = []
if options.bamlist:
    with open(options.bamlist) as f:
        for line in f:
            line = line.strip('\n')
            myseqlist.append(line)
else:
    myseqlist.append(options.singlebam)
samplelist = myseqlist

windowsearchlist = []
if options.windowlist:
    with open(options.windowlist) as f:
        for line in f:
            line = line.strip('\n')
            line = int(line)
            windowsearchlist.append(line)
elif options.windowauto:
    refs = []
    bam_file = bampath + myseqlist[0]
    pybam = pysam.Samfile(bam_file, "rb")
    for reference, length in zip(pybam.references, pybam.lengths):
        refs.append((length))
        pybam.close()
    windowsearchlist.extend(range(1,refs[0],int(options.windowauto)))
else:
    windowsearchlist.append(int(options.singlewindow))

windowsize = int(options.windowsize)
hap_read_threshold = int(options.mincov)
minimum_base_quality = int(options.baseq)
minimum_mapq = int(options.mapq)
minimum_read_quality = int(options.readq)

###################
print("windowsize is %i" % windowsize)


##################
# Set up loop parameters
window_unique_haplotypes = pd.DataFrame()
my_read_ids = pd.DataFrame()

my_read_ids = pd.DataFrame()


#################
# Primary Loop - iterate through all reads that span the window, filtering for poor quality, and extract the unique sequences (and read names)

for currentwindow in windowsearchlist:     
    print("Performing analysis for window starting at " + str(currentwindow) + " bp")
    #for currentwindow in bamsearchlist:
    for currentseq in samplelist:
    #for currentseq in currentseq:
        print("Analysing sample " + currentseq)    
        group_haplotypes = pd.DataFrame()
        bam_file = bampath + currentseq
        samfile = pysam.AlignmentFile(bam_file, "rb" )
        mainreference = samfile.references[0]       
        trimmed_reads = []
        read_names = []
#        trimmed_reads = set()
#        read_names = set()        
        endrem_collect = []      
        seq_my_read_ids = pd.DataFrame()
        
        for pileupcolumn in samfile.pileup(mainreference, currentwindow-1, currentwindow+windowsize):
            position = pileupcolumn.reference_pos
            
            for pileupread in pileupcolumn.pileups:
                read = pileupread.alignment
                querypos = pileupread.query_position
                mapstart = pileupread.alignment.query_alignment_start
                mapend = pileupread.alignment.query_alignment_end
                refstart = pileupread.alignment.reference_start
                refend = pileupread.alignment.reference_end
                startremove = currentwindow - refstart - 1
                endremove = refend - (currentwindow + windowsize)
                readlength = pileupread.alignment.query_length
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                    # skip optical duplicates (e.g. PCR duplicates reads)
                elif pileupread.alignment.is_duplicate or pileupread.alignment.is_qcfail:
                    continue                
                elif pileupread.alignment.is_qcfail:
                    continue
                    # skip reads with mapq below threshold
                elif pileupread.alignment.mapping_quality < minimum_mapq:
                    continue
                    # skip reads with no mapq specified
                elif read.mapping_quality == 255:
                    continue
                    # skip mean qscore of the read below threshold
                elif np.mean(read.query_qualities) < minimum_read_quality:
                    continue
                    # skip reads with a base quality below threshold
                elif read.query_qualities[pileupread.query_position] < minimum_base_quality:
                    continue
                    # skip reads with less than the minimum window length after soft clipping
                elif len(pileupread.alignment.query_alignment_sequence) < windowsize:
                    continue
                    # skip read if maps to reference after beginnning of window
                elif pileupread.alignment.reference_start > (currentwindow-1):
                    continue
                    # skip read if maps to reference after beginnning of window
                elif pileupread.alignment.reference_end < currentwindow + (windowsize):
                    continue
                    # skip if aligned sequence contains any 'N's (i.e. ambiguous base calls)
                elif bool(re.search('N', pileupread.alignment.query_alignment_sequence)):
                    continue
                trimmed_reads.append(pileupread.alignment.query_alignment_sequence[startremove:(startremove+windowsize)])                           
                read_names.append(pileupread.alignment.query_name)
        samfile.close()  

        # Loop can pull out each read multiple times - need to deduplicate
        all_reads = pd.DataFrame({'read_id':read_names,'seq':trimmed_reads},columns=['read_id','seq'])
        all_reads.drop_duplicates(subset='read_id',keep='first',inplace=True)
        trimmed_reads = list(all_reads['seq'])
        read_names = list(all_reads['read_id'])
        
        # Find Unique haplotypes in window
        unique_haplotypes = [[x,trimmed_reads.count(x)] for x in set(trimmed_reads)]
        # sort haplotypes
        unique_haplotypes.sort(key=lambda haplotype: haplotype[1],reverse=True)

        # split into different lists for downstream
        unique_haplotypes_df = pd.DataFrame(unique_haplotypes,columns=['read', 'count'])
        unique_haplotypes_df['Sample'] = currentseq
        unique_haplotypes_df['WindowStart'] = currentwindow
        unique_haplotypes_df['WindowEnd'] = currentwindow + windowsize -1
        
        # Filter out haplotypes below minimum threshold
        unique_haplotypes_df= unique_haplotypes_df[unique_haplotypes_df['count'] >= int(hap_read_threshold)]       
        
        # uniquely label each haplotype    
        Seq_haplotype = []
        for currentrow in range(0,len(unique_haplotypes_df['read'])): 
            Seq_haplotype.append('SHap-' + str(currentrow))
        unique_haplotypes_df['Seq_haplotype'] = Seq_haplotype              
        myhapnames = []

        window_unique_haplotypes = pd.concat([window_unique_haplotypes,unique_haplotypes_df])
        seq_my_read_ids['id'] = read_names
        seq_my_read_ids['reads'] = trimmed_reads
        seq_my_read_ids['window'] = currentwindow

        my_read_ids = my_read_ids.append(pd.DataFrame(data = seq_my_read_ids), ignore_index=True)

    # generates a table of all potential haplotypes for assignment
    group_haplotypes = pd.DataFrame(window_unique_haplotypes.groupby(['read','count']).sum()).reset_index()
    group_haplotypes = group_haplotypes.sort_values(by=['count'],ascending=False)


####################
# Find Group haplotypes by window (accross samples)

raw_haplotypes2 = pd.DataFrame()
for currentwindow in list(set(my_read_ids['window'])):
    rawgroup_haplotypes = my_read_ids[my_read_ids['window']==currentwindow]
    rawgroup_haplotypes = pd.DataFrame([[x,list(rawgroup_haplotypes['reads']).count(x)] for x in list(set(rawgroup_haplotypes['reads']))])
    rawgroup_haplotypes.columns = ['haplotype', 'Gcount']
    rawgroup_haplotypes['window'] = currentwindow 
    rawgroup_haplotypes = rawgroup_haplotypes.sort_values(by=['Gcount'],ascending=False).reset_index(drop=True)
    group_hap = []
    for currentrow in range(0,len(rawgroup_haplotypes)):       
           group_hap.append('GHap-' + str(currentrow))
    rawgroup_haplotypes['GHap'] = group_hap
    raw_haplotypes2 = pd.concat([raw_haplotypes2, rawgroup_haplotypes], axis=0)


group_unique_haplotypes = window_unique_haplotypes.merge(raw_haplotypes2,left_on='read', right_on='haplotype',how='inner', copy=False, sort=False)
myhapnames = []
for currentrow in range(0,len(group_unique_haplotypes)): # uniquely label each haplotype           
    currentseq = group_unique_haplotypes.loc[currentrow,'Sample']
    windowstart = group_unique_haplotypes.loc[currentrow,'WindowStart']
    windowend = group_unique_haplotypes.loc[currentrow,'WindowEnd']
    seqhap = group_unique_haplotypes.loc[currentrow,'Seq_haplotype']    
    grouphap = group_unique_haplotypes.loc[currentrow,'GHap']    
    mycount = group_unique_haplotypes.loc[currentrow,'count']    
    myhapnames.append(currentseq + '_Window-' + str(windowstart) + '-' + str(windowend) + '_' + seqhap + '_' + grouphap + '_Count-' + str(mycount))

group_unique_haplotypes['hapnames'] = myhapnames


group_label_reads = group_unique_haplotypes.merge(my_read_ids,left_on='read', right_on='reads',how='right', copy=False, sort=False)
group_label_reads = group_label_reads.loc[:,['id','hapnames']]


# Generage output files
group_unique_haplotypes.to_csv(options.outprefix+"bam2haplotype.tsv", sep='\t',index=False,header=True)
group_unique_haplotypes_df2 = group_unique_haplotypes.loc[:,['hapnames','read']]
group_unique_haplotypes_df2.to_csv(options.outprefix+"bam2haplotype.tab", sep='\t',index=False,header=False)
group_label_reads.to_csv(options.outprefix+"bam2haplotype.readlist.tsv", sep='\t',index=False,header=True)


count = SeqIO.convert(options.outprefix+"bam2haplotype.tab", "tab", options.outprefix+"bam2haplotype.fasta", "fasta")
print("Converted %i haplotypes to fasta" % count)



