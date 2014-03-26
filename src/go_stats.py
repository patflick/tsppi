

import pappi.go.fastdag as fastdag
from pappi.data_config import *


def dag_graph_size(dag):
    # get number of nodes
    n = len(dag.terms)

    # get number of edges
    m = sum(len(p) for p in dag.parents.values())

    return (n, m)

if __name__ == "__main__":
    namespaces = ["biological_process", "molecular_function",
                  "cellular_component"]
    # first import all:
    dag = fastdag.GODag(GO_OBO_FILE)
    print("Whole graph size: " + str(dag_graph_size(dag)))

    for namespace in namespaces:
        new_dag = fastdag.GODag(GO_OBO_FILE, only_namespace=namespace)
        print("namespace '" + namespace + "' graph size: "
              + str(dag_graph_size(new_dag)))
