

# for combinations()
import itertools
# for matrix and mean
import numpy
# for (de)serialization of the mapping
import json
# for checking if file exists
import os.path

from pappi.go_fastdag import GODag
from pappi.go_fast_similarity import GoFastSimilarity


class GoPreBufSimilarity(GoFastSimilarity):
    def __init__(self, obo_file, sim_file, mapping_file, sql_conn, verbose=True):
        # initialize the super class
        GoFastSimilarity.__init__(self, obo_file, sql_conn, verbose)

        # load the data file or create and fill it
        self._load_or_create(sim_file, mapping_file)

    def _load_or_create(self, sim_file, mapping_file):
        # both files (score matrix and term mapping) need to exist,
        # otherwise both have to be created from scratch
        if os.path.isfile(sim_file) and os.path.isfile(mapping_file):
            # load matrix and mapping file
            score_matrix = numpy.load(sim_file)
            with open(mapping_file, 'r') as f:
                mapping_str = json.loads(f.read())
            # JSON dict is {string -> int}, but we need {int -> int}
            mapping = dict()
            for k, v in mapping_str.items():
                mapping[int(k)] = v
        else:
            score_matrix, mapping = self._fill_sim_matrix()
            # saving for future reuse:
            #  save matrix
            numpy.save(sim_file, score_matrix)
            #  save mapping, note this saves {string -> int} rather than
            #  {int -> int}, we'll need to take care of that when loading
            with open(mapping_file, 'w') as f:
                f.write(json.dumps(mapping))

        # set as members
        self.score_matrix = score_matrix
        self.term_mapping = mapping


    def _fill_sim_matrix(self, verbose=True):
        # get the unique go_terms
        terms = set()
        for term_set in self.assoc.values():
            terms = terms.union(term_set)
        # intersect with terms in GO Dag
        terms = terms.intersection(self.go_dag.terms)

        # number of terms
        nTerms = len(terms)

        # map terms to numberical range [0, nTerms-1]
        i = 0
        term_mapping = dict()
        idx_2_term = list()
        for t in terms:
            term_mapping[t] = i
            idx_2_term.append(t)
            i = i + 1

        if verbose:
            print("Pre-calculating SemSim between all " + str(len(terms)) + " terms...")
        sim_vals = numpy.zeros((nTerms, nTerms))

        # fill matrix in row major
        for i in range(0, nTerms):
            t1 = idx_2_term[i]
            if verbose:
                if i % 10 == 0:
                    print("filling row " + str(i) + "/" + str(nTerms))
            for j in range(i+1, nTerms):
                t2 = idx_2_term[j]
                sim = self._simRel_score(t1, t2)
                sim_vals[i][j] = sim
                sim_vals[j][i] = sim

        # fill diagonal with (1-p)
        if verbose:
            print("filling diagonal...")
        for i in range(0, nTerms):
            t = idx_2_term[i]
            # SimRel score for identical terms is (1-p)
            sim_vals[i][i] = 1.0 - self.go_dag.p[t]

        # return the score matrix but also the term mapping
        return (sim_vals, term_mapping)


    def term_pairwise_score(self, term1, term2):
        # if either term does not exist: return 0
        if not term1 in self.term_mapping or not term2 in self.term_mapping:
            return 0.0
        # get indeces into the score matrix
        i1 = self.term_mapping[term1]
        i2 = self.term_mapping[term2]
        score = self.score_matrix[i1][i2]
        return score

    # TODO: maybe implement prebuffering for pairwise gene rather than pairwise
    #       term
