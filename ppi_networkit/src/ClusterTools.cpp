#include "ClusterTools.h"


void clusterSizeHist(const NetworKit::Partition& clustering)
{
  std::vector<count> cluster_sizes = clustering.subsetSizes();
  std::map<count, unsigned int> size_hist;
  for (auto size : cluster_sizes)
  {
    // if the element is available
    if (size_hist.count(size) > 0)
    {
      size_hist[size] += 1;
    }
    else
    {
      size_hist[size] = 1;
    }
  }

  // print out hist of cluster sizes
  for (auto size_count : size_hist)
  {
    std::cout << size_count.first << ": " << size_count.second << std::endl;
  }
}

void clusterSizeHist(const std::vector<NetworKit::Partition>& clusterings)
{
  std::map<count, unsigned int> size_hist;
  for (const NetworKit::Partition& clustering : clusterings)
  {
      std::vector<count> cluster_sizes = clustering.subsetSizes();
      for (auto size : cluster_sizes)
      {
          size_hist[size] += 1;
      }
  }

  // print out hist of cluster sizes
  for (auto size_count : size_hist)
  {
    std::cout << size_count.first << ": " << size_count.second << std::endl;
  }
}
