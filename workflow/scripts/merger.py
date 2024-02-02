'''
contact:    wlee9829@gmail.com
date:       2020_11_02
python:     python3.8
script:     merger_v2.py

This script is modified from the original merger.py to process the gap1filter.sam 
containing the "SA" tag, which triggers an error during CRSSANT assembly.

1. gaptypes.py in CRSSANT will generate 5 different types of sam files:
cont.sam;  gap1.sam;  gapm.sam;  trans/rri.sam;  homo.sam;  bad.sam
2. gap1filter.sam and trans/rri.sam are typically combined to produce "alignfile.sam" 
which is then used to assemble alignments to DGs and NGs using the crssant.py script.
3. gap1filter.sam contains some chimeric reads with "SA" tag, generated by two distant 
chimeric alignments. The "SA" tag triggers an error when performing CRSSANT analysis (line 163).
4. This script can be used to remove the "SA" tag from reads in gap1filter.sam and then 
merge gap1filter.sam and trans/rri.sam to one sam file.
5. Optionally, since the gap1filter.sam reads alone are enough for analyzing self structures and 
long range interactions of some genes (e.g., virus genomic RNA), merger_v2.py is modified to 
allow for removing the "SA" tag from only gap1filter.sam.
'''

# Import packages
import sys, argparse, os, re, time
from datetime import datetime
def timenow(): return str(datetime.now())[:-7]


gap1sam = open(str(snakemake.input.gap1filter),"r")
outsam = open(snakemake.output[0],"w")

####################################################################################################

print(timenow()+" Starting merger_v2.py...")
print(timenow()+" Processing gap1filter.sam for SA tag...")
for line in gap1sam:
	if line[0]=="@": outsam.write(line)
	elif "SA:Z:" not in line.split()[-1]: outsam.write(line)
	elif "SA:Z:" in line.split()[-1]: outsam.write('\t'.join(line.split()[:-2])+'\n')
gap1sam.close()

try: 
	trans = open(snakemake.input.trans,"r")
	print(timenow()+" Processing trans/rri.sam ...")
	for line in trans:
		if line[0]!="@": outsam.write(line)
	print(timenow()+" SA tag removal and merging completed successfully.")
	trans.close()

except:
	print(timenow()+" SA tag removal completed successfully.")
outsam.close()