"""
written for Python 3.5
author: Ute Hoffmann
add annotation to DESeq2 table (generally usable for everything, but dedicated for DESeq2 tables), assumes \t is separator
usage:
TUs_add_info.py   [old_file] [new_annotated_file]
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
			TU = line.split("\t")[1].strip('"')
			location = line.split("\t")[0].strip('"')
			sense_tags = line.split("\t")[2].strip('"')
			annotation = line.split("\t")[3].strip('"')
			antisense_tags = line.split("\t")[4].strip('"').strip("\n")
			dict_overview[TU] = [location, sense_tags, annotation, antisense_tags]			
	return dict_overview

def addAnnotation(old_file, new_file, annotation_dict, column_locus_tags=0, header=False, DESeq2_noHeader=False):
	with open(old_file) as old:
		with open(new_file, "w") as new:
			if DESeq2_noHeader: # for DESeq2 files from usegalaxy.eu: they don't have headers
				new.write("locus_tag\tlocalization\tsense_tags\tannotation\tantisense_tags\tmean_counts\tlog2FC\tstderr-log2FC\tWald\tp\tp.adj-BH\n")
				DESeq2=False
			for line in old:
				if header:
					new_line = "locus_tag\tlocalization\tsense_tags\tannotation\tantisense_tags\t" + line.split("\t", 1)[1]
					new.write(new_line)
					header=False
					continue
				split_line = line.split("\t",1)
				locus_tag = split_line[0].strip('"')
				try:
					new.write(locus_tag + "\t" + annotation_dict[locus_tag][0] + "\t" + annotation_dict[locus_tag][1] + "\t" + annotation_dict[locus_tag][2] + "\t" + annotation_dict[locus_tag][3] + "\t" + split_line[1])
				except KeyError:
					new.write(locus_tag + "\tNA\tNA\t" + split_line[1])
					with open("logfile.txt", "a") as log:
						log.write(locus_tag + "\twas not identified. Printed NA\n")
	return

if __name__ == "__main__":
	# reset log file
	with open("logfile.txt", "w") as log:
		log.write("")
	# read in annotation
	annotation_file = "TUs_genes.csv"
	locus_annotation = readInOverview(annotation_file)
	# do the actual work
	addAnnotation(sys.argv[1], sys.argv[2], locus_annotation, header=True)
	
