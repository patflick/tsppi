#!/usr/bin/python

# filters the human ppi to only include proteins that are in the tissue file

import sys
import os
import re

def filterPPI(ppifilename, tissuefilename):
	# keep proteins in an empty set
	proteins = set()
	with open(tissuefilename, "r") as tissuefile:
		# ignore header:
		tissuefile.readline()
		# compile RE
		p = re.compile('^(")(ENSP\d+)(".*)')
		
		# read all lines
		for line in iter(tissuefile.readline, b''):
			m = p.match(line)
			protein = m.group(2)
			proteins.add(protein)
			
	sys.stderr.write("finished filling set\n")
	counter = 0
	
	# now filter all proteins in PPI network
	with open(ppifilename, "r") as ppifile:
		p = re.compile(r'^9606\.(ENSP\d+) 9606\.(ENSP\d+) \d+$')
		for line in iter(ppifile.readline, b''):
			m = p.match(line)
			p1 = m.group(1)
			p2 = m.group(2)
			counter += 1
			if (counter % 100000 == 0):
				sys.stderr.write("processed " + str(counter) + " lines...\n")
			if (p1 in proteins and p2 in proteins):
				sys.stdout.write(line)

	



if __name__ == "__main__":
	if len(sys.argv) < 3:
		sys.exit('Usage: %s <tissuefile> <ppi_file>', sys.argv[0])
		
	for i in range(2):
		if not os.path.exists(sys.argv[i]):
			sys.exit('ERROR: File %s was not found!' % sys.argv[i])
			
	tissuefilename = sys.argv[1]
	ppifilename = sys.argv[2]
	
	# call routines to filter and output to stdout
	filterPPI(ppifilename, tissuefilename)