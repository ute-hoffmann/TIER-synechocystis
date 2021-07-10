"""
written for Python 3.5
author: Ute Hoffmann

input: fasta with RNA sequences
script uses RNAfold to calculate minimal free energy (mfe, delta G) in a sliding window of certain width (default: 50) and returns table with mfe values

usage:
RNAfold_inputFasta_outputMfeTable.py   [fasta_file] [name_output_file] [width of sliding window]
              
OR: use as function read_fasta_calc_mfe_write_table(fasta, mfe_file, width)

! does not check input parameters ! be careful when setting them !
"""

import sys
sys.path.append("/usr/local/lib/python3.7/site-packages/RNA")
import _RNA as RNA
from pysam import FastxFile
from random import sample

def read_fasta_calc_mfe_write_table(fasta, mfe_file, width=50):
	print("now opening " + fasta)
	# open fasta file
	sequences = FastxFile(fasta)

	# prepare dicts
	dict_for_mfe = {}
	dict_for_control = {}
	dict_for_mfe["names"] = []
	dict_for_control["names"] = []

	# initialize dict according to length of sequences
	j = 1
	for entry in sequences:
		# only consider 1st sequence
		if j==2:
			break
		j += 1
		
		# extract sequence
		seq = entry.sequence
		start=0
		end=width
		
		# iterate through length of seq and create index_centres
		for i in range(len(seq)-width+1):
			index_centre = (start+1+end)/2
			dict_for_mfe[index_centre] = []	
			dict_for_control[index_centre] = []	
			start += 1
			end += 1
	sequences.close()
	
	## Actual RNAfold execution
	# re-open sequences file to start from 1st sequence again
	sequences = FastxFile(fasta)
	for entry in sequences:
		# extract sequence
		seq = entry.sequence
		
		# extract name and add to "names"
		dict_for_mfe["names"].append(entry.name)
		dict_for_control["names"].append(entry.name)
		
		# create shuffled sequence as control
		seq_shuffle = "".join(sample(seq,len(seq)))

		# define width of sliding window
		start=0
		end=width

		# sliding window, return mfe
		for i in range(len(seq)-width+1):
			index_centre = (start+1+end)/2
			
			# extract sequence in sliding window from long sequences
			sub_seq = seq[start:end]
			sub_seq_shuffle=seq_shuffle[start:end]
			
			# calculate mfe and store numbers in dicts
			(ss,mfe) = RNA.fold(sub_seq)
			dict_for_mfe[index_centre].append(mfe)
			(ss,mfe) = RNA.fold(sub_seq_shuffle)
			dict_for_control[index_centre].append(mfe)
			
			# go 1nt further to the 'right'
			start += 1
			end += 1
	sequences.close()

	## Print results
	# first row: name of peak
	# first column: index of middle nt / window was centred around this position
	# each subsequent column: mfe at this position in one of the sequences - each column represents a sequence
	print("now writing " + mfe_file)
	with open(mfe_file, "w") as out:
		for i in dict_for_mfe.keys():
			out.write(str(i) + "\t")
			counter = 1
			for mfe in dict_for_mfe[i]:
				counter += 1
				if counter <= len(dict_for_mfe[i]):
					out.write(str(mfe) + "\t")
				else:
					out.write(str(mfe))
			out.write("\n")
			
	print("now writing " + mfe_file[:-4] + "_control.tsv")
	with open(mfe_file[:-4] + "_control.tsv", "w") as out:
		for i in dict_for_control.keys():
			out.write(str(i) + "\t")
			counter = 1
			for mfe in dict_for_control[i]:
				counter += 1
				if counter <= len(dict_for_control[i]):
					out.write(str(mfe) + "\t")
				else:
					out.write(str(mfe))
			out.write("\n")

	
if __name__ == "__main__":
	RNA.cvar.temperature=39.0 # set correct temperature so that folding energies are in correct range, scales with temperature, default 37.0C
	fasta_file=sys.argv[1]
	output_tsv=sys.argv[2]
	fold_width=int(sys.argv[3])
	read_fasta_calc_mfe_write_table(fasta_file, output_tsv, fold_width)

