'''
Created on 15 feb 2013

@author: flick
'''
import argparse
import re
import sys

def filterPPI(ppifile, genesfile, ppioutfile):
    # keep genes in an empty set
    genes = set()

    # compile RE
    p = re.compile('^(ENSG\d+)$')
    
    # read all lines
    for line in iter(genesfile.readline, b''):
        m = p.match(line)
        gene = m.group(1)
        genes.add(gene)
    
    genesfile.close()
    sys.stderr.write("finished filling set\n")
    counter = 0
    
    # now filter all proteins in PPI network
    ppioutfile.write("Gene1 Gene2 Score\n")
    
    p = re.compile(r'^9606\.(ENSG\d+) 9606\.(ENSG\d+) (\d+)$')
    for line in iter(ppifile.readline, b''):
        m = p.match(line)
        if (not m):
            print "wtf: " + line
        g1 = m.group(1)
        g2 = m.group(2)
        counter += 1
        if (counter % 100000 == 0):
            sys.stderr.write("processed " + str(counter) + " lines...\n")
        if (g1 in genes and g2 in genes):
            ppioutfile.write(g1 + " " + g2 + " " + m.group(3) + '\n')
            
    ppifile.close()

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intersects a PPI with a set of genes and outputs the resulting (smaller) PPI network that contains only the edges between the genes in the set.")
    parser.add_argument('ppi_file',type=argparse.FileType("r"), help='the PPI network to be filtered, this is assumed to contain ENSG <-> ENSG edges.')
    parser.add_argument('genes_file',type=argparse.FileType("r"), help='a file that contains the genes that are to be intersected with the ppi network. This is assumed to be a file in which each line contains only one ENSG gene ID.')
    parser.add_argument('ppi_output_file',type=argparse.FileType("w"), help='the output file for the filtered PPI network')
    args = parser.parse_args()
    
    filterPPI(args.ppi_file, args.genes_file, args.ppi_output_file)