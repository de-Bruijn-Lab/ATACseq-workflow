#!/usr/bin/env python

### 15/01/18
### In house workflow for ATACseq data processing
# This workflow is for after peak-calling
# Of note - different sample sheets for intersections and DiffBind. This allows us to specify biological replicates for DiffBind, and technical replicates for intersecting!

#####################

import argparse, os, csv, glob, sys, re
import pandas as pd

"""
Argument Parsing
"""

parser = argparse.ArgumentParser(description='Activate this in Working Directory.')

parser.add_argument(
	'-genome',
	help = 'Genome used for mapping during previous pipeline processing',
	required = True
)
parser.add_argument(
	'-s',
	help = 'Required. Specify the Samples.txt file to be used. This is a tab deliminated file that specifies, in order, Sample group (Must be prefix of file name) peak file location, bam file location, and BAM READ COUNT.',
	required = True
)
parser.add_argument(
	'-sDB',
	help = 'Required. Specify the Samples.csv to be used. must conform to DiffBind format. Note - Factor column is used for batches by my script!.',
	required = True
)


args = parser.parse_args()

"""
PreProcessing
"""

### Should be more verbose to disclose all parameters used

os.system("mkdir -p ./reads_normalised")
os.system("mkdir -p ./peaks_coverage")
os.system("mkdir -p ./tmp")
os.system("mkdir -p ./GoldStandardAnnotation")
print('#######################################\n')
print('de Bruijn Lab ATACseq analysis pipeline\n Joe Harman - 16/11/2017\n')
print('#######################################\n\n')

bedList = []
peakList = []
bamList = []
bamNormList = []
groupList = []
goldStandardList = []
readCountList = []
cwd = os.getcwd()


"""
Convert to bed files and merge peaks
"""

with open(args.s,'r') as f:
	next(f) # skip headings
	reader=csv.reader(f,delimiter='\t')
	for group,fPeak,fBam, readCount in reader:
		groupList.append(group)
		peakList.append(fPeak)
		bamList.append(fBam)
		readCountList.append(readCount)	
        
		# Convert to .bed format
		print('Converting .narrowPeak to .bed format, and merging overlapping peaks\n\n')
		fbed = os.path.dirname(fPeak) + '/' + os.path.splitext(os.path.basename(fPeak))[0] + '.bed'
		#os.system('cat ' + fPeak + ' | cut -f 1-3 > ' + fbed)  # restore this if blacklisted regions should be kept!
		os.system('cat ' + fPeak + ' | cut -f 1-3 | intersectBed -v -a stdin -b /t1-data/user/jharman/Scripts/Full_blacklist_clip_FINAL.bed > ' + fbed) # Removal of blacklisted genomic regions!
        
		# Merge peaks
		fbedMerge = os.path.dirname(fbed) + '/' + os.path.splitext(os.path.basename(fbed))[0] + '_mergePeaks.bed'
		os.system('bedtools merge -d 150 -i ' + fbed + ' > ' + fbedMerge)
		bedList.append(fbedMerge)    
		
		
# Print samples

groupListSet = set(groupList)
minCount = min(readCountList)
x = [bedList, peakList, bamList, groupList]

print('\n\nSample Data:\n\n')
for v in zip(*x):
        print(v)
        

"""
Gold Standard Peaks
"""

print('Intersecting sample replicates! Peaks retained if found in 2 or more replicates - (gold standard peak sets)\n')

for group in groupListSet:
	# Extract files matching Group
	r = re.compile('.*/' + group + '.*')
	groupSub = filter(r.match, bedList)
	
	# setup variables
	groupSubTmp = []
	groupLen = len(groupSub)

	# setup files & directories
	newf = './peaks_intersect/' + group + '_intersect_all.bed'
	os.system("mkdir -p ./peaks_intersect")
	os.system('> ' + newf)
	
	# Perform intersect with every relevant combination (1 vs all others, to create all pairwise intersections)
	if groupLen > 1:
		for i in range(groupLen):
		
			groupSubTmp = list(groupSub)
			del groupSubTmp[i]
		
			###fbedIntersect = os.path.dirname(groupSub[i]) + '/' + os.path.splitext(os.path.basename(groupSub[i]))[0] + 'intersect.bed'
		
			# Piping used here to directly append to file
			os.system('bedtools intersect -wa -a ' + groupSub[i] + ' -b ' + (' '.join(groupSubTmp)) + ' >> ' + newf)
	else: 
		print('\tWarning, only one biological replicate for ' + group + '... No intersection performed!\n')
		os.system(groupSub + ' >> ' + newf)
		
	# Pipe intersected data through sorting, then merging
	os.system('sort -k1,1 -k2,2n ' + newf + ' | bedtools merge -i stdin > ' + './peaks_intersect/' + group + '_goldStandardPeaks.bed')
	os.system('''awk \'{printf \"'''+ group +'''_peak_%d\\t%s\\n\", NR, $0}\' ''' + './peaks_intersect/' + group + '_goldStandardPeaks.bed > ./peaks_intersect/' + group + '_goldStandardPeaks_tmp.bed')
	os.system('''awk -v OFS='\t' '{print $2,$3,$4,$1,"","."}' ''' + './peaks_intersect/' + group + '_goldStandardPeaks_tmp.bed' + ' > ./peaks_intersect/' + group + '_goldStandardPeaks_ID.bed')
	goldStandardList.append('./peaks_intersect/' + group + '_goldStandardPeaks_ID.bed')
	
print('\tgold standard peak sets from all sample groups combined, and overlapping peaks merged, to create a total peak set\n')
os.system('cat ./peaks_intersect/*_goldStandardPeaks_ID.bed | sort  -k1,1 -k2,2n | bedtools merge -c 4,5,6 -o distinct -i stdin > ./peaks_intersect/allSamples_goldStandardPeaks_ID.bed')
goldStandardList.append('./peaks_intersect/allSamples_goldStandardPeaks_ID.bed')

# Peak annotation
print('Annotating gold standard peak sets\n\n')
for f in glob.iglob("./peaks_intersect/*_goldStandardPeaks_ID.bed"):
	os.system("annotatePeaks.pl " + f + " " + args.genome + " -annStats ./GoldStandardAnnotation/Ann_Stats_" + os.path.splitext(os.path.basename(f))[0] + ".txt > ./GoldStandardAnnotation/" + os.path.splitext(os.path.basename(f))[0] + "_annotated.txt")


"""
Coverage
"""
print('Calculating coverage (reads per genomic region)\n\n')

os.system('echo "chr	start	end	PeakID" > ./peaks_coverage/allPeaksCoverage.txt')
os.system('cat ./peaks_intersect/allSamples_goldStandardPeaks_ID.bed | cut -f 1-4 >> ./peaks_coverage/allPeaksCoverage.txt')

os.system("cat ./peaks_intersect/allSamples_goldStandardPeaks_ID.bed | cut -f 1-4 | bedtools nuc -fi /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin | awk ' {print $6}' > ./tmp/tmp_GC.txt")
os.system('paste ./peaks_coverage/allPeaksCoverage.txt ./tmp/tmp_GC.txt > ./tmp/tmp2.txt' )
os.system('cat ./tmp/tmp2.txt > ./peaks_coverage/allPeaksCoverage.txt')

for bam in bamList:
	os.system('echo \"' + bam + '\" > ./tmp/tmp.txt')
	os.system('bedtools coverage -a ./peaks_intersect/allSamples_goldStandardPeaks_ID.bed -b ' + bam + "| awk ' {print $6}' >> ./tmp/tmp.txt")
	os.system('paste ./peaks_coverage/allPeaksCoverage.txt ./tmp/tmp.txt > ./tmp/tmp2.txt' )
	os.system('cat ./tmp/tmp2.txt > ./peaks_coverage/allPeaksCoverage.txt')
	
os.system('Rscript /t1-data/user/jharman/Scripts/ATACseq/CQN_norm.R')

"""
DiffBind
"""

print("\n\nProcessing data with DiffBind (R)")
os.system("Rscript /t1-data/user/jharman/Scripts/ATACseq/DiffBind.R " + args.sDB)
print("DiffBind processing complete!")

flags = []
flagStatus = []

with open("CompFlags.txt",'r') as f:
	next(f) # skip headings
	reader = csv.reader(f, delimiter = "\t")
	for A,B in reader:
		flags.append(A)
		flagStatus.append(B)
tmp = [x for (x, y) in zip(flags, flagStatus) if y == 'FALSE']
print("\t DiffBind Comparisons on: " + ", ".join(tmp))


for f in os.listdir("./DiffBind_Results"):
	os.system("annotatePeaks.pl ./DiffBind_Results/" + f + " " + args.genome + " -annStats ./DiffBind_Results/Ann_Stats_" + os.path.splitext(f)[0] + ".txt > ./DiffBind_Results/" + os.path.splitext(f)[0] + "_annotated.txt")
	
	a = pd.read_csv("./DiffBind_Results/" + os.path.splitext(f)[0] + "_annotated.txt", sep="\t")
	a = a.rename(columns = {a.columns.values[0]: "PeakID", "Peak Score":"width"})
	a["Gene_Peak"] = a["Gene Name"].map(str) + "_" + a["PeakID"].map(str)
	b = pd.read_csv("./DiffBind_Results/" + f, sep="\t")
	b.drop(b.columns[[0,1,2,4,5]], axis=1, inplace=True)
	
	merged = a.merge(b, on='PeakID')
	merged.to_csv("./DiffBind_Results/" + os.path.splitext(f)[0] + "_annotated_merged.csv", index=False)


"""
Bam normalisation
"""
print('Normalising bam files\n\n')
os.system("mkdir -p ./reads_normalised")

for group in groupListSet:
	bamFiles = [bam for g, bam in zip(groupList, bamList) if group in g]
	newBam = './reads_normalised/' + group + '_combined.bam'
	newBamSort = './reads_normalised/' + group + '_combined.sorted.bam'
	newBamBW = './reads_normalised/' + group + '_combined.sorted.bw'
	
	os.system('samtools merge ' + newBam + ' ' + ' '.join(bamFiles) )
	os.system('samtools sort ' + newBam + ' > ' + newBamSort)
	os.system('samtools index ' + newBamSort)
	
	os.system('bamCompare -b1 ' + newBamSort + ' -b2 /t1-data/user/jharman/Scripts/Background_ATAC_mouse_Merge_Sorted.bam -o ' + newBamBW + ' --scaleFactorsMethod readCount --ratio ratio --binSize 1 --ignoreForNormalization chrX chrM --bl /t1-data/user/jharman/Scripts/Full_blacklist_clip_FINAL.bed --skipNAs')
	
	

print('COMPLETE!\n\n')

