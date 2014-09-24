
#ifndef CLUSTER_TOOLS_H
#define CLUSTER_TOOLS_H

#include <structures/Partition.h>

typedef NetworKit::count count;

void clusterSizeHist(const NetworKit::Partition& clustering);
void clusterSizeHist(const std::vector<NetworKit::Partition>& clusterings);

#endif // CLUSTER_TOOLS_H
