
# TODO remove this abomination :)
import sys
sys.path.append("/home/patrick/dev/bio/fastSemSim-0.8.1")

# for combinations()
import itertools
import numpy

import fastSemSim.Ontology.ontologies
from fastSemSim.Ontology.AnnotationCorpus import AnnotationCorpus
from fastSemSim.SemSim.SetSemSim import SetSemSim


class GoSimilarity:
    # the biological process root node
    #BP_root = "GO:0008150"
    BP_root = 8150

    def __init__(self, obo_file, assoc_file, sql_conn, verbose=True):
        # load gene ontology with fastSemSim
        if verbose:
            print("loading gene ontology")
        self.go = fastSemSim.Ontology.ontologies.load(obo_file)

        # load annotation class with fastSemSim
        if verbose:
            print("creating annotation corpus")
        self.ac = AnnotationCorpus(self.go)

        # set parameters for accessions
        ac_params = {}
        # only keep go term and accession name in memory
        ac_params['simplify'] = True
        if verbose:
            print("parsing association file")
        self.ac.parse(assoc_file, "GOA")

        # check for consistency
        if not self.ac.isConsistent():
            raise Exception

        if verbose:
            print("loading gene associations")
        # TODO: might do this jointly with fastSemSim's associations
        # this load all HGNC gene -> GO term associations into memory
        self.gene_assoc = self.load_go_associations(sql_conn)

        if verbose:
            print("creating BPscore object")
        # load the semantic similarity class for SimRel (Schlicker et al. 2006)
        # using the `max` function as mixer, resulting in an implementation
        # for the BPScore metric from Schlicker et al. 2006
        self.semsim_class = SetSemSim(self.go, self.ac, TSS="SimRel",
                                      MSS="max")

        #print(self.semsim_class.util.lineage)

    #################################################
    #  load GO associations for genes in given set  #
    #################################################
    def load_go_associations(self, sql_conn, only_genes=None):
        cur = sql_conn.cursor()
        cur.execute("SELECT * FROM go_gene_assoc")
        assoc = {}
        for row in cur.fetchall():
            gene = row[0]
            go_term = row[1]
            if only_genes and gene not in only_genes:
                continue
            if type(go_term) is str:
                go_term = self.name2id(go_term)
            if gene in assoc:
                assoc[gene].add(go_term)
            else:
                assoc[gene] = set([go_term])
        return assoc

    def name2id(self, name):
        return int(name.split(":")[1])

    def pairwise_bp_score(self, gene1, gene2):
        # get similarity score for BP subtree
        score = self.semsim_class.SemSim(self.gene_assoc[gene1],
                                         self.gene_assoc[gene2])
                                        # root=self.BP_root)
        # in case one gene is not annotated with a BP term, None is returned
        # in that case the genes have 0 similarity
        if score is None:
            score = 0.0
        return score

    def avg_bp_score(self, gene_set):
        """
        Calculates the average BPscore for a set of genes by first determining
        the BPscore between all pairs of genes in the set and then averaging
        over all BPscores.
        """
        genes = set()
        for g in gene_set:
            if g in self.gene_assoc:
                genes.add(g)
        scores = []
        for g1, g2 in itertools.combinations(genes, 2):
            scores.append(self.pairwise_bp_score(g1, g2))
        # get mean
        avg = numpy.mean(scores)
        return avg
