

# for combinations()
import itertools
# for mean
import numpy
# the shared association loading function
from pappi.go.utils import load_go_associations_sql


class GoSimilarity:

    def __init__(self):
        pass

    def load_go_associations(self, sql_conn, only_genes=None):
        return load_go_associations_sql(sql_conn, only_genes)

    def name2id(self, name):
        return int(name.split(":")[1])

    # defines the similarity measure for two GO terms
    def term_pairwise_score(self, term1, term2):
        raise NotImplementedError

    def terms_sets_scores(self, terms1, terms2):
        scores = []

        # filter to only have terms from the Go DAG
        terms1 = terms1.intersection(self.go_dag.terms)
        terms2 = terms2.intersection(self.go_dag.terms)
        for t1 in terms1:
            for t2 in terms2:
                score = self.term_pairwise_score(t1, t2)
                scores.append(score)
        return scores

    # defines how scores are combined (i.e. max for BPscore)
    def gene_pairwise_score(self, gene1, gene2):
        raise NotImplementedError

    def gene_set_score(self, gene_set):
        """
        Calculates the average BPscore for a set of genes by first determining
        the BPscore between all pairs of genes in the set and then averaging
        over all BPscores.
        """
        genes = set()
        for g in gene_set:
            if g in self.assoc:
                genes.add(g)
        scores = []
        for g1, g2 in itertools.combinations(genes, 2):
            scores.append(self.gene_pairwise_score(g1, g2))
        # get mean
        avg = numpy.mean(scores)
        return avg
