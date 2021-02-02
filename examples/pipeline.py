
# PIPELINE For the construction for ancesteral sequences based on the tutorial provided by Rubin Studer:
# https://evosite3d.blogspot.com/2014/09/tutorial-on-ancestral-sequence.html

# Dependence all on path and aliased ()
# mafft - multiple sequence alignment
# phyml - large tree building using maximum likelihood
# codeml - (i.e. PAML) ancesteral state reconstruction


__author__ = "Joshua E Goldford"
__date__ = 'Feb. 4th 2017'


import subprocess
import os
from helpers import fasta2phylip,rst2fasta,PAMLparams

# the file ID is the prefix used before the .fasta identifier in the file name
FileID = 'K00384_renamed_trimmed';

if not os.path.exists(os.getcwd() + '/' + FileID +'.fasta'):
	print('error:  fasta file not found')
	exit();


# 1. Multiple sequence alignment using MAFFT
cmd = 'mafft-linsi %(fid)s.fasta > %(fid)s_MSA.fasta' % {'fid': FileID};
subprocess.call(cmd,shell=True);

# 2. Convert Fasta to Phylip
inFile = '%(fid)s_MSA.fasta' % {'fid': FileID};
outFile = '%(fid)s_MSA.phylip' % {'fid': FileID};

# only work woth phylip format for sequences
seqsFile = outFile;
# call helper function
fasta2phylip(inFile,outFile);

# 3. Call phyml (Alias is not working for some reason)
cmd = '/Library/PhyML-3.1/PhyML-3.1_macOS-MountainLion -i %(fid)s -d aa -m JTT -c 4 -a e -b 0' % {'fid': seqsFile};
subprocess.call(cmd,shell=True);
# move file to TreeName
treeFile = FileID+'.tree';
cmd = 'mv %(seqs)s_phyml_tree.txt %(tree)s' %  {'seqs': seqsFile,'tree': treeFile};
subprocess.call(cmd,shell=True);

PAMLoutputFile = FileID+'.mlc';

# 4 construct PAML parameters
params = PAMLparams(seqsFile,treeFile,PAMLoutputFile);
# note that you can simply change values by chainging values of attributes (e.g. params.verbose = 0)
text_file = open("control_file.ctl", "w")
text_file.write("%s" % params.toText())
text_file.close()

# 5 call PAML
cmd = 'codeml control_file.ctl';
subprocess.call(cmd,shell=True);

# 6 convert PAML outut to Fasta File
rst2fasta('rst','ancesteral_sequences.fasta');
