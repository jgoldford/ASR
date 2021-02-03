import sys
import os
from Bio import SeqIO
import numpy as np
import csv


class PAMLparams:
	def __init__(self, seqfile,treefile,outfile):
		self.seqfile = seqfile;
		self.treefile = treefile;
		self.outfile = outfile;
		self.noisy = 9;
		self.verbose = 2;
		self.runmode = 0;
		self.seqtype = 2;
		self.clock = 0;
		self.aaDist = 0;
		self.aaRatefile = '/Users/joshuagoldford/Documents/paml4.8/dat/wag.dat';
		self.model = 2;
		self.icode= 0;
		self.Mgene= 0;
		self.fix_alpha = 0;
		self.alpha = 0.5;
		self.Malpha = 1;
		self.ncatG = 4;
		self.getSE = 0;
		self.RateAncestor = 1;
		self.Small_Diff = 0.5e-6;
		self.cleandata = 0;
		self.method = 1;

	def toText(self):
		return '''
		  seqfile = %(seqfile)s    * sequence data filename
		 treefile = %(treefile)s   * tree structure file name
		  outfile = %(outfile)s    * main result file name

			noisy = %(noisy)s  * 0,1,2,3,9: how much rubbish on the screen
		  verbose = %(verbose)s  * 0: concise; 1: detailed, 2: too much
		  runmode = %(runmode)s   * 0: user tree;  1: semi-automatic;  2: automatic
					   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

		  seqtype = %(seqtype)s   * 1:codons; 2:AAs; 3:codons-->AAs

			clock = %(clock)s   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
		   aaDist = %(aaDist)s  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
	   aaRatefile = %(aaRatefile)s   * only used for aa seqs with model=empirical(_F)

					   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

			model = %(model)s
					   * models for codons:
						   * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
					   * models for AAs or codon-translated AAs:
						   * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
						   * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

			icode = %(icode)s  * 0:universal code; 1:mammalian mt; 2-10:see below
			Mgene = %(Mgene)s     * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff             
					   * AA: 0:rates, 1:separate

		fix_alpha = %(fix_alpha)s   * 0: estimate gamma shape parameter; 1: fix it at alpha
			alpha = %(alpha)s * initial or fixed alpha, 0:infinity (constant rate)
		   Malpha = %(Malpha)s   * different alphas for genes
			ncatG = %(ncatG)s   * # of categories in dG of NSsites models


			getSE = %(getSE)s  * 0: don't want them, 1: want S.E.s of estimates

	 RateAncestor = %(RateAncestor)s  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

	   Small_Diff = %(Small_Diff)s
		cleandata = %(cleandata)s  * remove sites with ambiguity data (1:yes, 0:no)?
		   method = %(method)s  * Optimization method 0: simultaneous; 1: one branch a time''' % self.__dict__
		

def fasta2phylip(infile,outfile):
	print("CONVERT FASTA TO PHYLIP")

	sequence_dict = {}
	print(infile)
	for record in SeqIO.parse(open(infile, "r"), "fasta") :
		tab = record.id.split(" ")
		#print record.title
		sequence = str(record.seq).replace(" ","")
		#print sequence, len(sequence)
		sequence_dict[tab[0]]= sequence
		if "U" in sequence:
			print(tab[0])
			sys.exit()
		
	print(len(sequence_dict))
	#sys.exit()

	# Test length of the alignment:
	alignment_length = 0
	for gene in sequence_dict:
		if (alignment_length != 0) and (len(sequence_dict[gene]) != alignment_length):
			print("Error in alignment length, exit on error !!!")
			sys.exit()
		else:
			alignment_length = len(sequence_dict[gene])

	number_of_seq = len(sequence_dict)
	print("Number of sequences:\t"+str(number_of_seq))
	print("Alignment length:\t"+str(alignment_length))
	print("Ratio =\t"+str(alignment_length/3))

	if alignment_length%3 != 0:
		print("Warning: Hum, your alignment didn't code for nucleotides")

	### Write PHYLIP file

	phyfile = open(outfile,"w")

	name_length = 50

	if len(sys.argv) > 3:
		name_length = int(sys.argv[3])

	phyfile.write(str(number_of_seq)+"\t"+str(alignment_length)+"\n")

	for gene in sequence_dict:
		if len(gene) > name_length:
			gene_name = gene[0:name_length].replace(" ","")
			if gene_name[-1] == "_":
				gene_name = gene_name[0:-1]
			##elif gene_name[-2] == "_":
	##            gene_name = gene_name[0:-2]    
		else:
			gene_name = gene
		phyfile.write(gene_name+"  "+sequence_dict[gene]+"\n")
	phyfile.close()


def rst2fasta(inFile,outFile):
	file_in = open(inFile,"r")
	file_out = open(outFile,'w');

	while 1:
		line = file_in.readline()
		if line == "":
			break
		line= line.rstrip()    # Remove end-of-line
		if line[0:4] == "node":
			tab = line.split("         ");
			line = ">"+tab[0].replace(" ","").replace("#","");
			file_out.write(line)
			file_out.write('\n');

			if len(tab)>1:
				line = tab[1].replace(" ","")
				file_out.write(line)
				file_out.write('\n');
			else:
				line = tab[1].replace(" ","")
				file_out.write(line)
				file_out.write('\n');

	file_in.close()
	file_out.close();


def getAncesteralTree(fileName,outfile):
	with open(fileName, "r") as fid:
		tree = fid.readlines()[14];
		tree = tree.rstrip()
	with open(outfile,"w") as fout:
		fout.write(tree);
	return tree;


def appendAttributesToNodes(fileIn,treeFile):
	with open(fileIn, 'rU') as csvfile:
		file = csv.reader(csvfile, delimiter=',');
		headers = file.next();
		lines = list(file)
		node_dict = {}
		for idx,line in enumerate(lines):
			if line[3] is not '-':
				attr = np.log(float(line[3])/float(line[4]));
			else:
				attr = 0; 
			node_dict[line[0]] = attr;
	tree_tab = []
	tree_file = open(treeFile,"r")
	#tree_file = open("FinalTree.tree","r")
	while 1:
		line = tree_file.readline()
		if line == "":
			break
		tab = line.split()
		tree_tab = tree_tab+tab
	tree_file.close()
	#print tree_tab
	new_tree_line = ""
	#print(node_dict.keys())
	for item in tree_tab:
		if 'node'+item in node_dict:
			print("node"+item, node_dict["node"+item])
			item = node_dict['node'+item]
			item = round(item,3)
			item = str(item)
		new_tree_line = new_tree_line+item
	print(new_tree_line)



