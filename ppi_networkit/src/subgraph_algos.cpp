#include "subgraph_algos.h"

// for clustering coefficients
#include <properties/ClusteringCoefficient.h>
//#include <properties/GraphProperties.h>
// betweenness
#include <centrality/Betweenness.h>
// for PLP clustering
#include <community/PLP.h>

#include "graph_algos.h"


typedef NetworKit::count count_t;

namespace tsppi {
namespace algo {

/*********************************************************************
 *                              Degrees                              *
 *********************************************************************/

std::vector<std::vector<NetworKit::count> > subgraph_degrees(const Subgraphs& subgraphs)
{
    std::vector< std::vector< count_t > > result(subgraphs.numberOfSubgraphs(), std::vector<count_t>(subgraphs.graph.numberOfNodes()));

    // use OpenMP parallelism on the most outer loop
    subgraphs.graph.parallelForNodes([&](int u){
        subgraphs.graph.forNeighborsOf(u, [&](int v)
        {
            // create the edge label by boolean AND of the two node labels
            boost::dynamic_bitset<> edge_label = subgraphs.node_exists[u] & subgraphs.node_exists[v];

            // add bitset to vector
            // for each tissue
            for (unsigned int t = 0; t < subgraphs.numberOfSubgraphs(); ++t)
            {
                result[t][u] += edge_label[t];
            }
        });
    });

    return result;
}

std::vector<NetworKit::count> subgraph_max_degrees(const Subgraphs& subgraphs)
{
    // init the output:
    std::vector<count_t> result(subgraphs.graph.numberOfNodes());

    // use OpenMP parallelism on the most outer loop
    subgraphs.graph.parallelForNodes([&](int u){
        std::vector<count_t> degree_per_subgraph(subgraphs.numberOfSubgraphs());

        subgraphs.graph.forNeighborsOf(u, [&](int v)
        {
            // create the edge label by boolean AND of the two node labels
            boost::dynamic_bitset<> edge_label = subgraphs.node_exists[u] & subgraphs.node_exists[v];

            // add bitset to vector
            // for each tissue
            for (unsigned int t = 0; t < subgraphs.numberOfSubgraphs(); ++t)
            {
               degree_per_subgraph[t] += edge_label[t];
            }
        });

        result[u] = *std::max_element(degree_per_subgraph.begin(), degree_per_subgraph.end());
    });

    return result;
}

// ~ Co-Expression Degree (used in Bossi&Lehner)
std::vector<NetworKit::count> subgraph_edge_exist_degrees(const Subgraphs& subgraphs)
{
    // init the output:
    std::vector<count_t> result(subgraphs.graph.numberOfNodes());

    // use OpenMP parallelism on the most outer loop
    subgraphs.graph.parallelForNodes([&](int u){
        subgraphs.graph.forNeighborsOf(u, [&](int v)
        {
            // create the edge label by boolean AND of the two node labels
            boost::dynamic_bitset<> edge_label = subgraphs.node_exists[u] & subgraphs.node_exists[v];

            // count the coexpressed neighbors, by checking if the edge is
            // expressed in any tissue
            if (edge_label.any())
            {
                result[u]++;
            }
        });
    });

    return result;
}



/*********************************************************************
 *     maximum and minumum exists overlap with neighboring nodes     *
 *********************************************************************/

std::vector<NetworKit::count> subgraph_neighbor_min_exists_count(const Subgraphs& subgraphs)
{
    // init the output:
    std::vector<count_t> result(subgraphs.graph.numberOfNodes(), subgraphs.numberOfSubgraphs());

    // use OpenMP parallelism for the outer most loop
    subgraphs.graph.parallelForNodes([&](int u){
        subgraphs.graph.forNeighborsOf(u, [&](int v)
        {
            // create the edge label by boolean AND of the two node labels
            boost::dynamic_bitset<> edge_label = subgraphs.node_exists[u] & subgraphs.node_exists[v];

            // only check neighbors which are coexpressed in at least one tissue
            if (edge_label.any())
            {
                result[u] = std::min<count_t>(result[u], subgraphs.node_exists[v].count());
            }
        });
    });

    return result;
}

std::vector<NetworKit::count> subgraph_neighbor_max_exists_count(const Subgraphs& subgraphs)
{
    // init the output:
    std::vector<count_t> result(subgraphs.graph.numberOfNodes(), 0);

    // use OpenMP parallelism for the outer most loop
    subgraphs.graph.parallelForNodes([&](int u){
        subgraphs.graph.forNeighborsOf(u, [&](int v)
        {
            // create the edge label by boolean AND of the two node labels
            boost::dynamic_bitset<> edge_label = subgraphs.node_exists[u] & subgraphs.node_exists[v];

            // only check neighbors which are coexpressed in at least one tissue
            if (edge_label.any())
            {
                result[u] = std::max<count_t>(result[u], subgraphs.node_exists[v].count());
            }
        });
    });

    return result;
}


/*********************************************************************
 *                         Aggregate graphs                          *
 *********************************************************************/

NetworKit::Graph subgraph_edge_coexist_count_graph(const Subgraphs& subgraphs)
{
  // copy the graph and then set all the edge weights

  // create new graph with same amount of nodes, but set it as weighted
  count_t n = subgraphs.graph.numberOfNodes();
  // beware, this creates a n*n*sizeof(double) matrix
  // which for STRING (n ~= 13000) is approx: 1.3 GiB
  NetworKit::Graph result(n, true);

  subgraphs.graph.forEdges([&](NetworKit::node u, NetworKit::node v){
    boost::dynamic_bitset<> expr_u = subgraphs.node_exists[u];
    boost::dynamic_bitset<> expr_v = subgraphs.node_exists[v];
    // use bitwise AND to get the edge expression vector
    boost::dynamic_bitset<> edge_expr = expr_u & expr_v;

    // get the edge expression score as fraction of tissues in which the
    // edge is expressed
    //double edge_score = edge_expr.count() * 1.0 / edge_expr.size();

    // get the number of expressed tissues per protein of the edge
    int c_x = expr_u.count();
    int c_y = expr_v.count();
    int c_xy = edge_expr.count();

    // calculate the coexpression score
    double nom = c_xy * 1.0;
    int max_c = std::max(c_x, c_y);
    double coexpr_score;
    if (max_c == 0)
    {
      coexpr_score = 0.0;
    }
    else
    {
      coexpr_score = nom / (double)max_c;
    }

    // add the edge with the calculated edge score
    result.addEdge(u, v, coexpr_score);
  });

  return result;
}

NetworKit::Graph subgraph_edge_correlation_graph(const Subgraphs& subgraphs)
{
  // copy the graph and then set all the edge weights

  // create new graph with same amount of nodes, but set it as weighted
  count_t n = subgraphs.graph.numberOfNodes();
  // beware, this creates a n*n*sizeof(double) matrix
  // which for STRING (n ~= 13000) is approx: 1.3 GiB
  NetworKit::Graph result(n, true);

  subgraphs.graph.forEdges([&](NetworKit::node u, NetworKit::node v){
    boost::dynamic_bitset<> expr_u = subgraphs.node_exists[u];
    boost::dynamic_bitset<> expr_v = subgraphs.node_exists[v];
    // use bitwise AND to get the edge expression vector
    boost::dynamic_bitset<> edge_expr = expr_u & expr_v;

    // get the edge expression score as fraction of tissues in which the
    // edge is expressed
    //double edge_score = edge_expr.count() * 1.0 / edge_expr.size();

    // calculate the shifted correlation
    int c_x = expr_u.count();
    int c_y = expr_v.count();
    int c_xy = edge_expr.count();
    int n = edge_expr.size();


    double corr_xy; // the correlation between x and y
    if (c_x == 0 || c_x == n || c_y == 0 || c_y == n)
    {
      corr_xy = 0.0;
    }
    else
    {
      int nom = n*c_xy - c_x*c_y;
      int denom = (n*c_x - c_x*c_x)*(n*c_y - c_y*c_y);
      corr_xy = 1.0*nom / sqrt(1.0*denom);
    }

    // corr(x,y) is in [-1, 1], scale to range [0,1]
    double scaled_corr = 0.5 + 0.5*corr_xy;

    // add the edge with the calculated edge score
    result.addEdge(u, v, scaled_corr);
  });

  return result;
}


/*********************************************************************
 *               Clustering Coefficients for Subgraphs               *
 *********************************************************************/

/*
 * Simplest implementation: create subgraphs and run NetworKit's CC algo
 */
std::vector<std::vector< double > > subgraph_cc(const Subgraphs& subgraphs)
{
  std::vector<std::vector< double > > result(subgraphs.numberOfSubgraphs());

  // call the NetworKit clustering coefficient code for each subgraph
  subgraphs.foreach_subgraph([&](/* const */ NetworKit::Graph& g, count_t id){
    // call NetworKit's implementation for CC
    NetworKit::ClusteringCoefficient cc;
    // TODO const qualified graph into NetworKit implementation
    result[id] = cc.exactLocal(g);
  });

  return result;
}


/*
 * Create subgraphs and run our CC implementation (neighbor combinations)
 * for each one.
 */
std::vector<std::vector< double > > subgraph_cc_neighbor_comb(const Subgraphs& subgraphs)
{
  std::vector<std::vector< double > > result(subgraphs.numberOfSubgraphs());

  // call our implementation (neighbor combinations) for the clustering
  // coefficients for each subgraph
  subgraphs.foreach_subgraph([&](const NetworKit::Graph& g, count_t id){
    result[id] = cc_neighbor_comb(g);
  });

  return result;
}


/*
 * `Vectorized` version of the neighbor combination CC algorithm, that
 * calculates the clustering coefficients for all nodes in all subgraphs
 * at the same time, without having to create a graph instance for each
 * subgraph.
 */
std::vector<std::vector< double > > subgraph_cc_neighbor_comb_vec(const Subgraphs& subgraphs)
{
  // initialize results with zero
  std::vector<std::vector< double > > result(subgraphs.numberOfSubgraphs(),
                    std::vector<double>(subgraphs.graph.numberOfNodes(),0.0));

  // use NetworKit types
  typedef NetworKit::node node;
  typedef NetworKit::count count;

    subgraphs.graph.balancedParallelForNodes([&](node u) {
        // get neighbors
        std::vector<node> neighbors;
        subgraphs.graph.forNeighborsOf(u, [&](node v){
          if (v != u) // ignore self loops
          {
            neighbors.push_back(v);
          }
        });
        count d = neighbors.size();

        if (d >= 2) // otherwise just leave at zero!
        {
          std::vector<unsigned int> triangles(subgraphs.numberOfSubgraphs(), 0);
          std::vector<unsigned int> degrees(subgraphs.numberOfSubgraphs(), 0);

          // get degrees
          for (std::size_t i = 0; i < neighbors.size(); ++i)
          {
            boost::dynamic_bitset<> edge_label = subgraphs.node_exists[u] & subgraphs.node_exists[neighbors[i]];
            for (unsigned int t = 0; t < subgraphs.numberOfSubgraphs(); ++t)
            {
              degrees[t] += edge_label[t];
            }

          }

          // for all pairs of neighbors:
          for (std::size_t i = 0; i < neighbors.size()-1; ++i)
          {
            for (std::size_t j = i+1; j < neighbors.size(); ++j)
            {
              // does the global graph have this edge?
              if (subgraphs.graph.hasEdge(neighbors[i], neighbors[j]))
              {
                // the triangle only exists if all three nodes are expressed
                boost::dynamic_bitset<> triangle_label = subgraphs.node_exists[u] & subgraphs.node_exists[neighbors[i]] & subgraphs.node_exists[neighbors[j]];
                for (unsigned int t = 0; t < subgraphs.numberOfSubgraphs(); ++t)
                {
                  triangles[t] += triangle_label[t];
                }
              }
            }
          }

          for (unsigned int t = 0; t < subgraphs.numberOfSubgraphs(); ++t)
          {
            unsigned int deg = degrees[t];
            if (deg >= 2)
              result[t][u] = (double)triangles[t] * 2.0 / (double)(deg * (deg - 1));
          }
        }
    });

  return result;
}



/*********************************************************************
 *                     Betweenness for Subgraphs                     *
 *********************************************************************/

/*
 * Simplest implementation: create subgraphs and run NetworKit's betweenness
 * algorithm.
 */
std::vector<std::vector< double > > subgraph_betweenness(const Subgraphs& subgraphs)
{
  std::vector<std::vector< double > > result(subgraphs.numberOfSubgraphs());

  // "inefficient" method: create a subgraph for each tissue and then get the
  // betweenness centrality within that subgraph
  subgraphs.parallel_foreach_subgraph([&](const NetworKit::Graph& g, count_t id) {
    NetworKit::Betweenness betweeness(g);
    betweeness.run();
    result[id] = betweeness.scores();
  });

  return result;
}

/*
 * This implements the betweenness calculation for the custom
 * TSPPI graph datastructure, so that it is not necessary to copy
 * and create subgraphs for every tissue.
 */
std::vector<std::vector<double> > subgraph_betweenness_fast(const Subgraphs& subgraphs)
{
  typedef NetworKit::node node;
  typedef NetworKit::count count;
  static const count inf = std::numeric_limits<count>::max();
  count z = subgraphs.graph.numberOfNodes();

  // allocate result
  std::vector<std::vector< double > > result(subgraphs.numberOfSubgraphs());


#ifdef BENCHMARK_USE_OMP
#pragma omp parallel for schedule(runtime)
#endif
  for (unsigned int t = 0; t < subgraphs.numberOfSubgraphs(); ++t)
  {
    result[t] = std::vector<double>(z, 0.0);

    subgraphs.graph.forNodes([&] (node s) {
        /* first check if this node is expressed in the current tissue
         * In case it is not expressed, then the betweenness is equal to zero
         * ( due to degree zero)
         */
        if(!subgraphs.node_exists[s][t])
          return;

        /* nodes in order of increasing distance from s. */
        //std::stack<node> increasing;
        std::vector<node> increasing(z);
        count incr_count = 0;


        /* determine the shortest path tree from s via bfs. */
        std::vector<std::vector<node>> parents(z);          /* parents in dag. */
        /* number and length of shortest paths (= inf means unvisited). */
        std::vector<count> nshort(z), lshort(z, inf);
        std::queue<node> bfs_queue;                    /* working queue. */

        /* bfs that also computes the shortest path dag. */
        bfs_queue.push(s);
        lshort[s] = 0;
        nshort[s] = 1;
        while (!bfs_queue.empty()) {
          node v = bfs_queue.front();
          bfs_queue.pop();
          increasing[incr_count++] = v;
          //increasing.push(v);

          subgraphs.graph.forNeighborsOf(v, [&] (node w) {
            // check if the neighbor is expressed <=> the edge is expressed
            if (subgraphs.node_exists[w][t])
            {
              /* w found for first time -> enqueue */
              if (lshort[w] == inf) {
                bfs_queue.push(w);
                lshort[w] = lshort[v] + 1;
              }

              /* another shortest path to w via v. */
              if (lshort[w] == lshort[v] + 1) {
                nshort[w] = nshort[w] + nshort[v];
                parents[w].push_back(v);
              }
            }
          });
        }

        /* now compute the dependencies in order of decreasing distance. */
        std::vector<double> dependency(z, 0);
        //while (!increasing.empty()) {
        for (count i = incr_count; i > 0; --i)
        {
          //node w = increasing.top();
          node w = increasing[i-1];
          //increasing.pop();

          for (node v: parents[w]) {
            /* recursive formula: see lecture. */
            dependency[v] += double(nshort[v])/nshort[w] * (1 + dependency[w]);
          }
          if (w != s) {
            result[t][w] += dependency[w];
          }
        }
    });

  }

  return result;
}


/*********************************************************************
 *                   PLP Clustering for Subgraphs                    *
 *********************************************************************/

/*
 * Simplest implementation: create subgraphs and run NetworKit's betweenness
 * algorithm.
 */
std::vector<NetworKit::Partition> subgraph_PLP(const Subgraphs& subgraphs)
{
  std::vector<NetworKit::Partition> result(subgraphs.numberOfSubgraphs());

  // create a graph instance for each graph and runt he NetworKit algo
  subgraphs.parallel_foreach_subgraph([&](const NetworKit::Graph& g, count_t id) {
    NetworKit::PLP clusterer;
    result[id] = clusterer.run(g);
  });

  return result;
}

/*
 * Modified version of the PLP algorithm that runs directly on the
 * Subgraphs datastructure by checking if an edge exists for each traversal.
 */
std::vector<NetworKit::Partition> subgraph_PLP_vec(const Subgraphs& subgraphs, NetworKit::count updateThreshold)
{
  // TODO: unify all these in our namespace (not every single function)
  typedef NetworKit::index index;
  typedef NetworKit::node node;
  typedef NetworKit::count count;
  using Partition = NetworKit::Partition;
  typedef index label; // a label is the same as a cluster id
  constexpr index none = NetworKit::none;



  count n = subgraphs.graph.numberOfNodes();
  // assume the nodes ids are contiguous
  index z = n;
  count nTissues = subgraphs.numberOfSubgraphs();


  // update threshold heuristic
  if (updateThreshold == none) {
    updateThreshold = (count) (n / 1e5);
  }

  // set unique label for each node
  std::vector<Partition> labels(nTissues, Partition(z));
  subgraphs.graph.parallelForNodes([&](node v) {
    for (count t = 0; t < nTissues; ++t)
      labels[t][v] = v;
  });

  // number of nodes which have been updated in last iteration
  std::vector<count> nUpdated(nTissues, n);

  count nIterations = 0; // number of iterations

  /**
   * == Dealing with isolated nodes ==
   *
   * The pseudocode published does not deal with isolated nodes (and therefore does not terminate if they are present).
   * Isolated nodes stay singletons. They can be ignored in the while loop, but the loop condition must
   * compare to the number of non-isolated nodes instead of n.
   *
   * == Termination criterion ==
   *
   * The published termination criterion is: All nodes have got the label of the majority of their neighbors.
   * In general this does not work. It was changed to: No label was changed in last iteration.
   */

  std::vector<boost::dynamic_bitset<> > activeNodes(z, boost::dynamic_bitset<>(nTissues)); // record if node must be processed
  for (boost::dynamic_bitset<>& bs : activeNodes)
  {
    // set all to 1
    bs.set();
  }
  std::vector<boost::dynamic_bitset<> > hasUpdated(z, boost::dynamic_bitset<>(nTissues)); // record if node must be processed

  //Aux::Timer runtime;

  // propagate labels
  while (*std::max_element(nUpdated.begin(), nUpdated.end()) > updateThreshold) { // as long as a label has changed...
    //runtime.start();
    nIterations += 1;
    //INFO("[BEGIN] LabelPropagation: iteration #" , nIterations);

    // reset updated
    nUpdated.assign(nTissues, 0);

    subgraphs.graph.balancedParallelForNodes([&](node v){
      // check if this node is still active
      if (! activeNodes[v].any())
        return;
      if (! subgraphs.node_exists[v].any())
      {
        activeNodes[v].reset();
        return;
      }

      // get count of co-expressed neighboirs
      count tissue_deg = 0;
      subgraphs.graph.forNeighborsOf(v, [&](node w){
          boost::dynamic_bitset<> edge_expr = subgraphs.node_exists[v] & subgraphs.node_exists[w];
          tissue_deg += (count) edge_expr.any();
      });
      if (tissue_deg == 0)
      {
        activeNodes[v].reset();
        return;
      }

      std::vector<std::map<label, count> > labelWeights(nTissues); // neighborLabelCounts maps label -> frequency in the neighbors
      //std::vector<std::vector<count> > labelWeights(nTissues, std::vector<count>(z,0));

      // weigh the labels in the neighborhood of v
      subgraphs.graph.forNeighborsOf(v, [&](node w) {
        boost::dynamic_bitset<> edge_expr = subgraphs.node_exists[v] & subgraphs.node_exists[w];
        for (count t = 0; t < nTissues; ++t)
        {
          labelWeights[t][labels[t][w]] += (count)edge_expr[t]; // add weight of edge {v, w}
        }
      });

      // get heaviest label
      for (count t = 0; t < nTissues; ++t)
      {
        if (subgraphs.node_exists[v][t])
        {
          label heaviest = std::max_element(labelWeights[t].begin(),
                  labelWeights[t].end(),
                  [](const std::pair<label, count>& p1, const std::pair<label, count>& p2) {
                    return p1.second < p2.second;})->first;
          //label heaviest = std::max_element(labelWeights[t].begin(), labelWeights[t].end()) - labelWeights[t].begin();

          if (labels[t][v] != heaviest) { // UPDATE
            labels[t][v] = heaviest; //labels[v] = heaviest;
            nUpdated[t] += 1; // TODO: atomic update?
            hasUpdated[v][t] = true;
          }
        }
      }

      activeNodes[v] = subgraphs.node_exists[v] & hasUpdated[v];

      if (hasUpdated[v].any())
      {
        subgraphs.graph.forNeighborsOf(v, [&](node u) {
          boost::dynamic_bitset<> edge_expr = subgraphs.node_exists[v] & subgraphs.node_exists[u];
          activeNodes[u] |= edge_expr & hasUpdated[v];
        });
      }

      hasUpdated[v].reset();
    });

    /*
    std::cout << " nUpdated: " << std::endl;
    for (count t = 0; t < nTissues; ++t)
    {
      std::cout << nUpdated[t] << ", ";
    }
    std::cout << std::endl;
    */

    // for each while loop iteration...

    //runtime.stop();
    //DEBUG("[DONE] LabelPropagation: iteration #" , nIterations , " - updated " , nUpdated , " labels, time spent: " , runtime.elapsedTag());


  } // end while

  //std::cout << "Fast - number of iteratons: " << nIterations << std::endl;

  return labels;

}

} // namespace algo
} // namespace tsppi


