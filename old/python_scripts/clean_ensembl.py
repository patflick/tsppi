#!/usr/bin/python

# get the ensembl ID matching from
# http://www.ensembl.org/biomart/martview/
# using homo sapiens version 70 database
# and output the attributes
#   Ensembl Gene ID
#   Ensembl Protein ID
# with NO filters

INPUT_FILE = './data/ensembl_ID_matching.csv'
OUTPUT_FILE = './data/ensemblID_matching_cleared.csv'

# import regular expressions
import re

def removeEmpty(infilename, outfilename):
	with open(infilename, 'r') as infile:
		with open(outfilename, 'w') as outfile:
			# read header line of input file
			outfile.write(infile.readline())
			for line in iter(infile.readline, b''):
				if (re.match("ENSG\d+,ENSP",line)):
					outfile.write(line)

def main():
	removeEmpty(INPUT_FILE, OUTPUT_FILE)
	#print re.match("ENSG\d+,ENSP", "ENSG00000262457,ENSP00000466098")
	
main()
	