'''
Created on 15 feb 2013

@author: flick
'''

import argparse
import re
import sys

# source of mdict: http://code.activestate.com/recipes/440502
class mdict(dict):
    def __setitem__(self, key, value):
        """add the given value to the list of values for this key"""
        self.setdefault(key, []).append(value)

def readGeneProteinMapping(infile):
    # create new mdict
    geneProteinMapping = mdict()
    # read header, and forget about it
    infile.readline()
    # read all lines, and add genes with protein lists to dictionary
    for line in iter(infile.readline, b''):
        # strip whitespaces
        [gene, protein] = line.split(",")
        gene = gene.strip(' \t\n\r')
        protein = protein.strip(' \t\n\r')
        geneProteinMapping[gene] = protein
    infile.close()
    return geneProteinMapping

def readProteinGeneMapping(infile):
    # create new mdict
    proteinGeneMapping = dict()
    # read header, and forget about it
    infile.readline()
    # read all lines, and add genes with protein lists to dictionary
    for line in iter(infile.readline, b''):
        # strip whitespaces
        both = line.split(",")
        if (len(both) < 2):
            print "wtf: " + line
        [gene, protein] = both
        gene = gene.strip(' \t\n\r')
        protein = protein.strip(' \t\n\r')
        if (proteinGeneMapping.has_key(protein)):
            print "wtf a protein has more than one gene?"
        proteinGeneMapping[protein] = gene
    infile.close()
    return proteinGeneMapping

def replace_P2G(ppifile,ppioutfile,P2G_mapping):
    p = re.compile(r'^9606\.(ENSP\d+) 9606\.(ENSP\d+) \d+$')
    gene_edges = set()
    counter = 0
    dupl = 0
    unmatched_protein_ctr = 0
    umatched_protein = set()
    

    line = ppifile.readline()
    while (line):

    # loop through all lines
    #for line in iter(ppifile.readline, b''):
        counter += 1
        m = p.match(line)
        p1 = m.group(1)
        p2 = m.group(2)
        if (not P2G_mapping.has_key(p1)):
            unmatched_protein_ctr += 1
            if (not p1 in umatched_protein):
                umatched_protein.add(p1)
            line = ppifile.readline()
            continue
        if (not P2G_mapping.has_key(p2)):
            unmatched_protein_ctr += 1
            if (not p2 in umatched_protein):
                umatched_protein.add(p2)
            line = ppifile.readline()
            continue
        g1 = P2G_mapping[p1]
        g2 = P2G_mapping[p2]
        # TODO figure out whether the PPI is directed or not
        gene_edge1 = g1 + g2
        gene_edge2 = g2 + g1
        if (not gene_edge1 in gene_edges and not gene_edge2 in gene_edges):
            gene_edges.add(gene_edge1)
            gene_edges.add(gene_edge2)
            line = re.sub(p1, g1, line)
            line = re.sub(p2, g2, line)
            ppioutfile.write(line)
        else:
            dupl += 1
            
        if (counter % 100000 == 0):
            sys.stderr.write("processed " + str(counter) + " lines...\n")
        line = ppifile.readline()
    sys.stderr.write("DONE!\n")
    sys.stderr.write(" -> processed a total of " + str(counter) + " lines\n")
    sys.stderr.write(" -> removed a total of " + str(dupl) + " duplicates \n")
    sys.stderr.write(" -> unmatched protein interactions: " + str(unmatched_protein_ctr) + " \n")
    sys.stderr.write(" -> unmatched unique proteins : " + str(len(umatched_protein)) + " \n")
    for p in umatched_protein:
        print p

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replaces proteins (ENSP) with genes (ENSG) in a PPI network, and removes all duplicate edges that result from the replacement.")
    parser.add_argument('ppi_file',type=argparse.FileType("r"), help='the string-db PPI network file to be filtered.')
    parser.add_argument('mapping_file',type=argparse.FileType("r"), help='the file that consists of the mapping between ESNG and ENSP identifiers (one gene maps to multiple proteins)')
    parser.add_argument('ppi_output_file',type=argparse.FileType("w"), help='the output file for the PPI network of genes (ENSG)')
    args = parser.parse_args()
    
    mapping = readProteinGeneMapping(args.mapping_file)
    replace_P2G(args.ppi_file, args.ppi_output_file, mapping)
    
