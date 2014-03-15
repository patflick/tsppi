

# for combinations()
import itertools
import numpy

from pappi.go_fastdag import GODag
from pappi.go_similarity import GoSimilarity

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
        #print("terms: " + term1 + ", " + term2)
        #term1 = name2id(term1)
        #term2 = name2id(term2)
        LCAs = self.go_dag.get_lca(term1, term2)
        # FIXME: take p from actual max lca
        lca_IC = max(self.go_dag.IC[t] for t in LCAs)
        lca_p = min(self.go_dag.p[t] for t in LCAs)
        denom = self.go_dag.IC[term1] + self.go_dag.IC[term2]
        if denom != 0:
            score = 2*lca_IC / denom * (1 - lca_p)
        else:
            score = 0.0
        if (score > 1.0 or score < 0.0):
            raise Exception("SimRel score is invalid")
        return score


    def gene_pairwise_score(self, gene1, gene2):
        terms1 = self.assoc[gene1]
        terms2 = self.assoc[gene2]
        scores = []

        # fixme: filter by set operations
        for t1 in terms1:
            if not self.go_dag.has_term(t1):
                continue
            for t2 in terms2:
                if not self.go_dag.has_term(t2):
                    continue
                score = self._simRel_score(t1, t2)
                scores.append(score)
        if len(scores) == 0:
            score = 0.0
        else:
            score = max(scores)
        if (score > 1.0 or score < 0.0):
            raise Exception
        return score
