
#include "ModularityPerCluster.h"

std::vector<double> ModularityPerCluster::get_modularities(const NetworKit::Partition& zeta, const NetworKit::Graph G)
{
    typedef NetworKit::node node;
    typedef NetworKit::index index;
    typedef NetworKit::edgeweight edgeweight;

    //double expCov; // term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
    //double modularity;  // mod = coverage - expected coverage
    double totalEdgeWeight = 0.0;

    // get total edge weight
    totalEdgeWeight = G.totalEdgeWeight(); // compute total edge weight in G
    if (totalEdgeWeight == 0.0) {
        throw std::invalid_argument("Modularity is undefined for graphs without edges (including self-loops).");
    }

    // compute volume of each cluster
    std::vector<double> incidentWeightSum(zeta.upperBound(), 0.0);  //!< cluster -> sum of the weights of incident edges for all nodes
    G.parallelForNodes([&](node v) {
        // add to cluster weight
        index c = zeta[v];
        assert (zeta.lowerBound() <= c);
        assert (c < zeta.upperBound());
#pragma omp atomic
        incidentWeightSum[c] += G.weightedDegree(v) + G.weight(v,v); // account for self-loops a second time
    });

    // compute intra-cluster edge weights per cluster
    std::vector<double> intraEdgeWeight(zeta.upperBound(), 0.0); // cluster -> weight of its internal edges
    // TODO: Make parallel, protect intraEdgeWeight[c]
    G.forEdges([&](node u, node v, edgeweight ew)
    {
        assert (u < zeta.numberOfElements());
        assert (v < zeta.numberOfElements());
        index c = zeta[u];
        index d = zeta[v];
        if (c == d) {
            assert ((zeta.lowerBound()) <= c && (c < zeta.upperBound()));
            intraEdgeWeight[c] += ew;
        } // else ignore edge
    });

    // calculate the modularity per cluster
    std::vector<double> modularity(zeta.upperBound(), 0.0);
    #pragma omp parallel for
    for (index c = zeta.lowerBound(); c < zeta.upperBound(); ++c)
    {
        double cov = intraEdgeWeight[c] / totalEdgeWeight;
        double incidentFrac = incidentWeightSum[c] / totalEdgeWeight / 2.0;
        double expCov = incidentFrac * incidentFrac;
        modularity[c] = cov - expCov;
    }

    return modularity;
}

