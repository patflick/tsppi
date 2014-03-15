

# for combinations()
import itertools
import numpy

# TODO remove this abomination :)
import sys
sys.path.append("/home/patrick/dev/bio/fastSemSim-0.8.1")
# import needed fastSemSim modules and classes
import fastSemSim.Ontology.ontologies
from fastSemSim.Ontology.AnnotationCorpus import AnnotationCorpus
from fastSemSim.SemSim.SetSemSim import SetSemSim

# import base class
from pappi.go_similarity import GoSimilarity


class GoFastSemSimSimilarity(GoSimilarity):
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
        self.assoc = self.load_go_associations(sql_conn)

        if verbose:
            print("creating BPscore object")
        # load the semantic similarity class for SimRel (Schlicker et al. 2006)
        # using the `max` function as mixer, resulting in an implementation
        # for the BPScore metric from Schlicker et al. 2006
        self.semsim_class = SetSemSim(self.go, self.ac, TSS="SimRel",
                                      MSS="max")

        #print(self.semsim_class.util.lineage)


    def gene_pairwise_score(self, gene1, gene2):
        # get similarity score for BP subtree
        score = self.semsim_class.SemSim(self.assoc[gene1],
                                         self.assoc[gene2],
                                         root=self.BP_root)
        # in case one gene is not annotated with a BP term, None is returned
        # in that case the genes have 0 similarity
        if score is None:
            score = 0.0
        return score
