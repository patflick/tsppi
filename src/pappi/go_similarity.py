

# for combinations()
import itertools
# for mean
import numpy

class GoSimilarity:

    def __init__(self):
        pass

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
            if g in self.gene_assoc:
                genes.add(g)
        scores = []
        for g1, g2 in itertools.combinations(genes, 2):
            scores.append(self.gene_pairwise_score(g1, g2))
        # get mean
        avg = numpy.mean(scores)
        return avg
