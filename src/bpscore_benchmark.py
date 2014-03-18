#!/usr/bin/python3

# for timing
import time

# for the data connection
import pappi.sql
from pappi.data_config import *

# import the GO association loading function
from pappi.go.utils import load_go_associations_sql

# import similarity scorer to be benchmarked
from pappi.go.fast_similarity import GoFastSimilarity
from pappi.go.fastSemSim_similarity import GoFastSemSimSimilarity
from pappi.go.prebuf_similarity import GoPreBufSimilarity
from pappi.go.gene_prebuf_similarity import GoGenePreBufSimilarity


class BPScore_Benchmarker:
    def __init__(self):
        # get database connection
        self.con = pappi.sql.get_conn(DATABASE)
        self.genes = self.get_benchmark_genes(self.con)
        self.scorers = []
        self.init_time = []
        self.run_times = dict()  # dict {number of genes -> list of run times}

    def init_scorers(self):
        # initialize all the scorers (and save the initalization time)
        start = time.time()
        self.scorers.append(GoFastSimilarity(GO_OBO_FILE, self.con, True))
        self.init_time.append(time.time() - start)

        start = time.time()
        self.scorers.append(GoFastSemSimSimilarity(GO_OBO_FILE, GO_ASSOC_FILE,
                                                   self.con))
        self.init_time.append(time.time() - start)

        start = time.time()
        self.scorers.append(GoPreBufSimilarity(GO_OBO_FILE, GO_SCORE_FILE,
                                               GO_SCORE_MAP_FILE, self.con, True))
        self.init_time.append(time.time() - start)

        start = time.time()
        self.scorers.append(GoGenePreBufSimilarity(GO_OBO_FILE, GO_SCORE_FILE,
                                                   GO_SCORE_MAP_FILE,
                                                   GO_BPSCORE_FILE,
                                                   GO_BPSCORE_MAP_FILE, self.con,
                                                   True))
        self.init_time.append(time.time() - start)

    def benchmark_scorers(self, nGenes):
        # get a set of genes with the given size
        benchmark_genes = set(self.genes[0:nGenes])
        # score the gene set with all scorers
        score_time = []
        for scorer in self.scorers:
            start = time.time()
            score = scorer.gene_set_score(benchmark_genes)
            score_time.append(time.time() - start)
        # save run time to class table
        self.run_times[nGenes] = score_time

    def get_benchmark_genes(self, sql_conn):
        # load Gene->GO-Term associations to get a set of genes to be used in
        # the benchmark
        assoc = load_go_associations_sql(sql_conn)
        # use a list for fast/efficient range access
        genes = list(assoc.keys())
        return genes

    def run_benchmark(self):
        for n in range(10, 1001, 10):
            print("benchmarking for n = " + str(n) + " genes...")
            self.benchmark_scorers(n)

    def print_timings(self):
        print("scored by " + str(len(self.scorers)) + " scorers")
        print()
        print("n\t" + "\t".join(self.scorers[i].__class__.__name__
              for i in range(0, len(self.scorers))))
        print("init\t" + "\t".join(str(self.init_time[i])
              for i in range(0, len(self.scorers))))
        for n in sorted(self.run_times.keys()):
            score_time = self.run_times[n]
            print(str(n) + "\t" + "\t".join(str(s) for s in score_time))

# the main benchmark:
if __name__ == '__main__':
    print("loading benchmarking class...")
    benchmarker = BPScore_Benchmarker()
    print("benchmark init times...")
    benchmarker.init_scorers()
    print("benchmark scoring...")
    benchmarker.run_benchmark()
    # print the actual timing results
    benchmarker.print_timings()
