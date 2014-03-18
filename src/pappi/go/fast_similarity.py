

# for combinations()
import itertools
import numpy

from pappi.go.fastdag import GODag
from pappi.go.similarity import GoSimilarity


# TODO abstract interface into shared base class
class GoFastSimilarity(GoSimilarity):
    # the biological process root node
    BP_root = 8150

    def __init__(self, obo_file, sql_conn, verbose=True):
        # load gene ontology with fastSemSim
        if verbose:
            print("loading gene ontology")
        # import only the BP namespace:
        self.go_dag = GODag(obo_file, only_namespace="biological_process")

        # load annotation class with fastSemSim
        if verbose:
            print("loading associations")
        self.assoc = self.load_go_associations(sql_conn)

        if verbose:
            print("initializing GO Dag probabilities and IC")
        self.go_dag.term_probability(self.assoc)
        self.go_dag.term_IC()

    def _simRel_score(self, term1, term2):
        # get the MICA (common ancestor with maximal IC)
        max_ic_anc = self.go_dag.max_ic_anc(term1, term2)
        # get IC and p() of that term
        lca_IC = self.go_dag.IC[max_ic_anc]
        lca_p = self.go_dag.p[max_ic_anc]
        # calc the denominator for the score
        denom = self.go_dag.IC[term1] + self.go_dag.IC[term2]
        if denom != 0:
            score = 2*lca_IC / denom * (1 - lca_p)
        else:
            score = 0.0
        if (score > 1.0 or score < 0.0):
            raise Exception("SimRel score is invalid")
        return score

    def term_pairwise_score(self, term1, term2):
        return self._simRel_score(term1, term2)

    def gene_pairwise_score(self, gene1, gene2):
        # get associated GO-Terms:
        terms1 = self.assoc[gene1]
        terms2 = self.assoc[gene2]

        # get scores:
        scores = self.terms_sets_scores(terms1, terms2)

        # get maximum score
        if len(scores) == 0:
            # set to zero in case there are no scores
            score = 0.0
        else:
            score = max(scores)
        if (score > 1.0 or score < 0.0):
            raise Exception
        return score
