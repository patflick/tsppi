

# for combinations()
import itertools
# for matrix and mean
import numpy
# for (de)serialization of the mapping
import json
# for checking if file exists
import os.path

from pappi.go.fastdag import GODag
from pappi.go.prebuf_similarity import GoPreBufSimilarity
from pappi.go.fast_similarity import GoFastSimilarity


class GoGenePreBufSimilarity(GoPreBufSimilarity):
    def __init__(self, obo_file, sim_file, term_mapping_file,
                 gene_bpscore_file, gene_mapping_file, sql_conn, verbose=True):
        # initialize the super class
        GoPreBufSimilarity.__init__(self, obo_file, sim_file,
                                    term_mapping_file, sql_conn, verbose)

        # load the data file or create and fill it
        self._load_or_create_bpscore(gene_bpscore_file, gene_mapping_file)

    def _load_or_create_bpscore(self, bpscore_file, mapping_file):
        # both files (score matrix and term mapping) need to exist,
        # otherwise both have to be created from scratch
        if os.path.isfile(bpscore_file) and os.path.isfile(mapping_file):
            # load matrix and mapping file
            bp_score_arrays = numpy.load(bpscore_file)
            score_matrx = bp_score_arrays['bp_score']
            row_sums = bp_score_arrays['row_sums']
            with open(mapping_file, 'r') as f:
                mapping = json.loads(f.read())
        else:
            score_matrix, row_sums, mapping = self._fill_bpscore_matrix()
            # saving for future reuse:
            #  save matrix
            numpy.savez(bpscore_file, bp_score=score_matrix, row_sums=row_sums)
            with open(mapping_file, 'w') as f:
                f.write(json.dumps(mapping))

        # set as members
        self.bpscore_matrix = score_matrix
        self.gene_mapping = mapping
        self.mapped_genes = set(mapping.keys())
        self.row_sums = row_sums

    def _fill_bpscore_matrix(self, verbose=True):
        # get unique genes:
        all_genes = set(self.assoc.keys())
        annotated_genes = set()

        for gene, terms in self.assoc.items():
            inter = terms.intersection(self.go_dag.terms)
            if len(inter) > 0:
                annotated_genes.add(gene)

        if verbose:
            print("Number of annotated genes: " + str(len(annotated_genes))
                  + "/" + str(len(all_genes)))

        # precompute only for annotated genes, others get a score of 0.0 anyway
        genes = annotated_genes
        # number of genes
        nGenes = len(genes)

        # map terms to numberical range [0, nTerms-1]
        i = 0
        gene_mapping = dict()
        idx_2_gene = list()
        for t in genes:
            gene_mapping[t] = i
            idx_2_gene.append(t)
            i = i + 1

        if verbose:
            print("Pre-calculating BPScore between all " + str(nGenes)
                  + " genes...")
        bpscore_matrix = numpy.zeros((nGenes, nGenes))

        # fill matrix in row major
        for i in range(0, nGenes):
            g1 = idx_2_gene[i]
            if verbose:
                if i % 10 == 0:
                    print("filling row " + str(i) + "/" + str(nGenes))
            for j in range(i, nGenes):
                g2 = idx_2_gene[j]
                bpscore = GoFastSimilarity.gene_pairwise_score(self, g1, g2)
                bpscore_matrix[i][j] = bpscore
                bpscore_matrix[j][i] = bpscore

        # fill rows sums
        row_sums = numpy.zeros(nGenes)
        for i in range(0, nGenes):
            # set the row sum but ignore self
            row_sums[i] = numpy.sum(bpscore_matrix[i]) - bpscore_matrix[i][i]

        # return the score matrix but also the gene mapping
        return (bpscore_matrix, row_sums, gene_mapping)

    def gene_pairwise_score(self, gene1, gene2):
        # if either gene does not exist: return 0
        if not gene1 in self.gene_mapping or not gene2 in self.gene_mapping:
            return 0.0
        # get indeces into the score matrix
        i1 = self.gene_mapping[gene1]
        i2 = self.gene_mapping[gene2]
        score = self.bpscore_matrix[i1][i2]
        return score

    def gene_cluster_score(self, genes):
        """
        Score the cluster by calculating the avg of all internal combinations
        and additionally calculating the avg of all combinations where only
        one gene is in this set. The set of all possible genes is take to be
        the same set as the genes covered by the bpscore matrix.
        """
        # only score those genes that are mapped
        genes = genes.intersection(self.mapped_genes)

        # size of the cluster
        n = len(self.mapped_genes)
        k = len(genes)

        # get all the ids -> only map each gene to an integer once
        # use in sorted order for more cache efficient access in the
        # score matrix
        gene_ids = list(sorted(self.gene_mapping[g] for g in genes))

        # sum of all bpscores comparing this cluster with all other genes
        ext_cluster_sum = 0.0
        in_cluster_sum = 0.0
        # get sum of pairwise bp_scores
        for i1 in gene_ids:
            # get sum of scores in row that belong to same cluster
            row_cluster_sum = 0.0
            for i2 in gene_ids:
                # ignore self
                if i2 == i1:
                    continue
                row_cluster_sum = row_cluster_sum + self.bpscore_matrix[i1][i2]
            ext_cluster_sum = ext_cluster_sum + (self.row_sums[i1]
                                                 - row_cluster_sum)
            in_cluster_sum = in_cluster_sum + row_cluster_sum

        # get the average scores
        avg_in_cluster = in_cluster_sum / (k*(k-1))
        avg_ext_cluster = ext_cluster_sum / (k*(n-k))
        return (avg_in_cluster, avg_ext_cluster)




                

        pass
