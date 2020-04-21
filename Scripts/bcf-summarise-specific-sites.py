#!/usr/bin/env python


import sys
import subprocess
import os
import math
import pandas as pd
import numpy as np
#import ggplot as gg

# To generate new annotations for snpEff: -
#	Create folder in ~/programs/snpEff/data/ with name of genome
#	Create entry in ~/programs/snpEff/snpEff.config with details of genome (search for Treponema for examples)
# 	Deposit WGS fasta (labelled sequences.fa) and gff3 (labelled genes.gff) in folder.
# 	Ensure sequence headers inside file correspond to the one used by the target bcf/vcf files (can parse using one liner)
# 	Build database in snpEff, e.g. : java -jar ~/programs/snpEff/snpEff.jar build -gff3 -v Neisseria_gonorrhoeae_FA_1090


print "\nWritten by Mat Beale (mathew.beale@sanger.ac.uk), Wellcome Trust Sanger Institute, March 2020"
print "Takes a input bcf file, tsv of variant positions, and the name of a mapped reference, and checks the specified sites for consensus and minority variants\n"

#######

from optparse import OptionParser
#def get_options():
usage = "Usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", action="store", type="string", dest="vcfinput",help="Specify input FILE")
parser.add_option("-r", action="store", type="string", dest="reference",help="Specify name of primary reference in bcf file")
parser.add_option("-v",action="store", type="string", dest="variantlist",help="Specify tsv in format '<site label> <variant position>' variant positions to check")

parser.add_option("-m", action="store",default=5, type="int",dest="minreads", help="Minimum number of reads required to call minor variant [%default]")
parser.add_option("-p", action="store", type="int",dest="minperc", default=5, help="Minimum percentage minor variant to accept [%default]")
parser.add_option("-t", action="store", type="int",dest="covthreshold", default=20, help="Minimum number of reads at any site to be included in analysis (variant and non variant) [%default]")

(options, args) = parser.parse_args()


#	return parser.parse_args()
#	if len(args) < 1:
#        	parser.error("incorrect number of arguments - must specify input using '-i'\n")
#get_options()

####

vcfinput = options.vcfinput

minreads = options.minreads
minfreq = options.minperc
covthreshold = options.covthreshold
reference = options.reference
variantlist=options.variantlist
####


# input name of bcf file - will want to tweak this to take command line arguments and maybe work from a sam/bam aswell
vcfoutput = vcfinput + ".temp"
#print vcfoutput


# Use bcftools to extract the required columns into a tsv file
bcftools_path = 'bcftools'
snpeff_path = '/nfs/users/nfs_m/mb29/programs/snpEff/'

print '\nExtracting relevant columns from '+ vcfinput +' into '+ vcfoutput


def bcfcommand():
	bcfrun = [bcftools_path,"query","-f",'%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t%INFO/PV4\n',"-o",vcfoutput,vcfinput]
	subprocess.call(bcfrun)

# need to restore this to make the script work again
bcfcommand()


### Now pull the tsv into the python script
vcfdf = pd.read_table(vcfoutput,sep='\t', header=None,names=['CHROM','POS','REF','ALT','DP','DP4','PV4','ANN'])
### Get rid of temporary output
subprocess.call(["rm",vcfoutput])


print 'Analysing '+ vcfoutput +' using reference ' +reference
print 'Calling minor variants >=' + str(minfreq) + '% with >=' + str(minreads) + ' supporting reads at sites with >=' + str(covthreshold) + ' reads'


# Label positions where there is a potential variant
vcfdf['PossVariant'] = np.where(vcfdf['ALT']=='.', 'None', 'Variant')


# Identify INDELS
def isINDEL(myrows):
    INDEL = []
    for currentrow in myrows:
        reflength = len(vcfdf.loc[currentrow]['REF'])
        altlength = len(vcfdf.loc[currentrow]['ALT'])
        if reflength >1 or altlength > 1:
            INDEL.append('TRUE')
        else:
            INDEL.append('FALSE')
    return INDEL
vcfdf['INDEL'] = isINDEL(range(len(vcfdf)))


# Expand 'DP4' columnn into 4 separate columns, then add back in
DP4 = pd.DataFrame(vcfdf['DP4'].map(eval).tolist())
DP4.columns = ['RF','RR','AF','AR']
# Remerge with vcfdf
vcfdf = pd.concat([vcfdf[:], DP4[:]], axis=1)

# Split up the PV4 (tests of bias)
# Deal with columns without any values by generating '1' four times, else use PV4 values
vcfdf['PV4'] = np.where(vcfdf['PV4']=='.', '1,1,1,1', vcfdf['PV4'])

# Now split into columns and remerge (as previously)
PV4 = pd.DataFrame(vcfdf['PV4'].map(eval).tolist())
PV4.columns = ["StrandBias","BaseQBias","MapQBias","TailDistBias"]

vcfdf = pd.concat([vcfdf[:], PV4[:]], axis=1)


# Calculate supporting reads from each side
vcfdf['RefReads'] = vcfdf[['RF','RR']].sum(axis=1)
vcfdf['AltReads'] = vcfdf[['AF','AR']].sum(axis=1)

# Also mark variants where count is below read threshold and coverage threshold
vcfdf['MinReadSup'] = np.where(vcfdf[['RefReads','AltReads']].sum(axis=1) >= minreads, 'ok','low')
vcfdf['CoverageThreshold'] = np.where(vcfdf[['RefReads','AltReads']].sum(axis=1) >= covthreshold, 'ok','low')

# If count is below read threshold, blank down to 0
vcfdf['RefReads'] = np.where(vcfdf['RefReads']<= minreads, 0, vcfdf['RefReads'])
vcfdf['AltReads'] = np.where(vcfdf['AltReads']<= minreads, 0, vcfdf['AltReads'])

# get variant percentages
vcfdf['RefPerc'] = (vcfdf['RefReads'] / vcfdf[['RefReads','AltReads']].sum(axis=1))*100
vcfdf['AltPerc'] = (vcfdf['AltReads'] / vcfdf[['RefReads','AltReads']].sum(axis=1))*100

# Work out the minor variant allele percentage
vcfdf['MinorPerc'] = vcfdf[['RefPerc','AltPerc']].min(axis=1)

# Determine variants that fail the bias tests
vcfdf['TestBias'] = 'Fail'
vcfdf.loc[(vcfdf['StrandBias'] >= 0.05) & 
          (vcfdf['BaseQBias'] >=0.05) &
          (vcfdf['MapQBias'] >=0.05) &
          (vcfdf['TailDistBias'] >=0.05), 'TestBias'] = 'Pass'

# Determine Minor SNPs
vcfdf['MinorSNP'] = 'None'
vcfdf.loc[(vcfdf['MinorPerc'] >= minfreq) & 
          (vcfdf['MinReadSup'] =="ok") &
          (vcfdf['INDEL'] == 'FALSE') &
          (vcfdf['PossVariant'] == 'Variant') &
          (vcfdf['StrandBias'] >= 0.05) & 
          (vcfdf['BaseQBias'] >=0.05) &
          (vcfdf['MapQBias'] >=0.05) &
          (vcfdf['TailDistBias'] >=0.05), 'MinorSNP'] = 'Minor'


### Import list of sites to check
variantdf = pd.read_table(variantlist,sep='\t', header=None,names=['label','position'])

### Subset dataframe to only keep sites in variant list
vardatadf = vcfdf[vcfdf.POS.isin(variantdf['position'])]

vardatadf['VariantPresent'] = 'NA'
vardatadf.loc[(vcfdf['CoverageThreshold']=='ok'),'VariantPresent'] = 'No'

vardatadf.loc[(vardatadf['MinorSNP'] == 'Minor') &
              (vardatadf['CoverageThreshold'] == 'ok'),'VariantPresent'] = 'Hetero'

vardatadf.loc[(vardatadf['TestBias'] == 'Pass') &
              (vardatadf['INDEL'] == 'FALSE') &
              (vardatadf['PossVariant'] == 'Variant') &
              (vardatadf['AltPerc'] >= 100 - minfreq) &
              (vardatadf['CoverageThreshold'] == 'ok'), 'VariantPresent'] = 'Yes'

vardatadf['VariantLabel'] = variantdf['label'].tolist()

vardatadf2 = vardatadf[['VariantLabel','POS','VariantPresent','REF','ALT','DP','AltPerc']]
vardatadf2 = vardatadf2.set_index('VariantLabel')

######
# make df single line
vardataSingle = vardatadf2.unstack().to_frame().T

vardataSingle.columns = vardataSingle.columns.map('{0[0]}_{0[1]}'.format)


vardataSingle.insert(0, 'Sample', vcfinput)


#print 'test single line'
#print vardataSingle


################################
vcfdf.drop(['PV4','DP4'],axis=1,inplace=True)



# Extract minor variant positions
minorvar = vcfdf.loc[(vcfdf['MinorSNP']=='Minor') & (vcfdf['CoverageThreshold']=='ok') & (vcfdf['CHROM']==reference)]


# Print full dataframe and minor vars only to separate tab delimited files
vcfdf.to_csv(vcfinput + ".analysed.tsv",sep='\t', index=False)
#minorvar.to_csv(vcfinput + ".minorvars.tsv",sep='\t', index=False)

# Print single line output
vardataSingle.to_csv(vcfinput + ".called-sites.tsv",sep='\t', index=False)




