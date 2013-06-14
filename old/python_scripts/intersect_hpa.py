#!/usr/bin/python

import argparse
import re



def hpa_line_get_gene(line):
	# compile the regular expression to match tissue lines
	p = re.compile('^(")(ENSG\d+)(".*)')
	m = p.match(line)
	if (m):
		return m.group(2)
	else:
		print "ERROR: ENSG pattern didn't match"
		
def hpa_line_get_APE_score(line):
	# compile the regular expression to match tissue lines
	# TODO this does not work yet:
	p = re.compile('.*"APE","([a-zA-Z ]*)"')
	m = p.match(line)
	if (m):
		return m.group(1)
	else:
		p = re.compile('.*"Staining","([a-zA-Z ]*)"')
		m = p.match(line)
		if (not m):
			print "ERROR: neither Staining nor APE pattern matched"
			print "for line:"
			print line

def hpa_gene_set(infile, min_score):
	# stores the result set
	genes = set()
	score2num = {"Very low":0,"Low":1,"Medium":2,"High":3}
	# discard first line
	infile.readline()
	
	# loop through all lines and add gene to set
	for line in iter(infile.readline, b''):
		gene = hpa_line_get_gene(line)
		score = hpa_line_get_APE_score(line)
		if (gene and score):
			if (score2num[score] >= score2num[min_score]):
				genes.add(gene)
	
	# reset file
	infile.seek(0)
	
	# return result
	return genes
	

# reads line for line from the input file
# and if the read gene is present in the set `genes`
# it writes the line to the output file
def hpa_filter_genes(infile, outfile, genes):
	# write first line directly to output
	outfile.write(infile.readline())
	
	# loop through all lines and add gene to set
	for line in iter(infile.readline, b''):
		gene = hpa_line_get_gene(line)
		if gene in genes:
			outfile.write(line)
	

def intersect_hpa(infile1,infile2,outfile1,outfile2,min_score,outproteins=None):
	# get the genes from both files
	print "Processing file 1..."
	genes1 = hpa_gene_set(infile1, min_score)
	print "Processing file 2..."
	genes2 = hpa_gene_set(infile2, min_score)
	
	# intersect the genes
	print "Intersecting genes"
	both = genes1.intersection(genes2)
	
	# output the common genes
	print "Writing output 1:"
	hpa_filter_genes(infile1, outfile1, both)
	print "Writing output 2:"
	hpa_filter_genes(infile2, outfile2, both)
	
	if (outproteins):
		for gene in both:
			outproteins.write(gene+"\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Intersects two HPA tissue data sets,\nand outputs only those genes that are scored with `APE`\nin both datasets.")
	parser.add_argument('inputfile1',type=argparse.FileType("r"))
	parser.add_argument('inputfile2',type=argparse.FileType("r"))
	parser.add_argument('outputfile1',type=argparse.FileType("w"))
	parser.add_argument('outputfile2',type=argparse.FileType("w"))
	parser.add_argument('-s','--min-score', dest='min_score',help='the minimum APE reliability score. genes that score less than this in either one of the sets, are not included in the intersection. (default=Medium)', default='Medium')
	parser.add_argument('-o','--output-genes', dest='out_genes',default=None,type=argparse.FileType('w'),help='outputs the genes that are in the intersection, separated by newline')
	args = parser.parse_args()
	
	intersect_hpa(args.inputfile1,args.inputfile2,args.outputfile1,args.outputfile2,args.min_score,args.out_genes)
	