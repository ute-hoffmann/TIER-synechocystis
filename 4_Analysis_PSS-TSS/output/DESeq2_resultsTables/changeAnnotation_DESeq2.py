"""
written for Python 3.7
author: Ute Hoffmann
add annotation to DESeq2 table (generally usable for everything, but dedicated for DESeq2 tables)
assumes that file with annotation information is in same directory, called "annotation_locusTags_stand13012021.csv"
assumes that column with locus tags is the first in table
usage:
changeAnnotation_DESeq2.py   [old_file] [new_annotated_file]
"""
import sys

def readInOverview(file_name):
	dict_overview = {}
	header=True
	with open(file_name) as overview:
		for line in overview:
			if header:
				header=False
				continue
			locus_tag = line.split("\t")[0].strip('"')
			gene_name = line.split("\t")[1].strip('"')
			gene_product = line.split("\t")[2].strip('"')
			dict_overview[locus_tag] = [gene_name, gene_product]			
	return dict_overview

def addAnnotation(old_file, new_file, annotation_dict, header=False, DESeq2_noHeader=False, sep="\t"):
	with open(old_file) as old:
		with open(new_file, "w") as new:
			if DESeq2_noHeader: # for DESeq2 files from usegalaxy.eu: they don't have headers
				new.write("locus_tag"+ sep + "gene_name" + sep + "gene_product" + sep + "mean_counts" + sep + "log2FC" + sep + "stderr-log2FC" + sep + "Wald" + sep + "p" + sep + "p.adj-BH\n")
				DESeq2_noHeader=False
			for line in old:
				if header:
					new_line = "locus_tag" + sep + "gene_name" + sep + "gene_product" + sep +  line.split(sep, 1)[1]
					new.write(new_line)
					header=False
					continue
				split_line = line.split(sep,1)
				locus_tag = split_line[0].strip('"')		
				try:
					new.write(locus_tag + sep + annotation_dict[locus_tag][0] + sep + annotation_dict[locus_tag][1] + sep + split_line[1])
				except KeyError:
					if "6803t" in locus_tag:
						new.write(locus_tag + sep + locus_tag + sep + "tRNA" + sep + split_line[1])
						continue
					try: 
						locus_noAppend = locus_tag.split("_")[0]
						append = locus_tag.split("_")[1]
						new.write(locus_tag + sep + annotation_dict[locus_noAppend][0] + sep + append + " of " + annotation_dict[locus_noAppend][1] + sep + split_line[1])
					except IndexError:
						new.write(locus_tag + sep + "NA" + sep + "NA" + sep + split_line[1])
						with open("logfile.txt", "a") as log:
							log.write(locus_tag + "\twas not identified. Printed NA\n")
					except KeyError:
						new.write(locus_tag + sep + "NA" + sep + "NA" + sep + split_line[1])
						with open("logfile.txt", "a") as log:
							log.write(locus_tag + "\twas not identified. Printed NA\n")
	return

if __name__ == "__main__":
	# reset log file
	with open("logfile.txt", "w") as log:
		log.write("")
	# read in annotation
	annotation_file = "annotation_locusTags_stand13012021.csv"
	locus_annotation = readInOverview(annotation_file)
	# do the actual work
	addAnnotation(sys.argv[1], sys.argv[2], locus_annotation, header=True)
	
