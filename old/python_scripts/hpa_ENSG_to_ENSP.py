#!/usr/bin/python

INPUT_TISSUE = "./data/tissues/cerebral.cortex_neuronal.cells.csv"
OUTPUT_TISSUE = "./data/tissues/mapping/cerebral.cortex_neuronal.cells.csv"

ENSEMBL_MATCHING = "./data/ensemblID_matching_cleared.csv"

import re

# source of mdict: http://code.activestate.com/recipes/440502
class mdict(dict):
    def __setitem__(self, key, value):
        """add the given value to the list of values for this key"""
        self.setdefault(key, []).append(value)

def readGeneProteinMapping(filename):
	# create new mdict
	geneProteinMapping = mdict()
	# open file to read
	with open(filename, 'r') as infile:
		# read header, and forget about it
		infile.readline()
		# read all lines, and add genes with protein lists to dictionary
		for line in iter(infile.readline, b''):
			# strip whitespaces
			[gene, protein] = line.split(",")
			gene = gene.strip(' \t\n\r')
			protein = protein.strip(' \t\n\r')
			geneProteinMapping[gene] = protein
	return geneProteinMapping
	

def geneToProteinsHPA(tissuefilename, outfilename, proteinmapping):
	# open tissue file to read
	with open(tissuefilename, "r") as infile:
		with open(outfilename, "w") as outfile:
			# transfer header:
			outfile.write(infile.readline())
			# compile the regular expression to match tissue lines
			p = re.compile('^(")(ENSG\d+)(".*)')
			# counters
			matched = 0
			genes_matched = 0
			skipped = 0
			# interate through all lines
			for line in iter(infile.readline, b''):
				m = re.match('^"(ENSG\d+)".*', line)
				gene = m.group(1)
				if (proteinmapping.has_key(gene)):
					genes_matched += 1
					for i in proteinmapping[gene]:
						outfile.write(p.sub(r'\1' + i + r'\3',line))
						matched += 1
				else:
					skipped += 1
	print "Genes matched:", genes_matched
	print "Proteins inserted:", matched
	print "Skipped:", skipped
	
	
if __name__ == '__main__':
	mapping = readGeneProteinMapping(ENSEMBL_MATCHING)
	geneToProteinsHPA(INPUT_TISSUE,OUTPUT_TISSUE,mapping)
	

