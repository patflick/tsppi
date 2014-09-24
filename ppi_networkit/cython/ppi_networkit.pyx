# This cython module extends the NetworKit python module
# The NetworKit cython/ directory has to be in the setup include path for
# this to compile correctly.
# NOTE that this creates it's own python module including the NetworKit
#	  wrappers and thus is not directly compatible with the original
#	  NetworKit python module.

# using include enables us to use all NetworKit wrappers,
# which is especially helpful to be able to use the NetworKit::Graph wrapper
include "_NetworKit.pyx"

# import C types not yet imported in NetworKit's code
from libc.stdint cimport uint32_t
from cython.operator cimport dereference as deref

#######################
#  Subgraphs wrapper  #
#######################

cdef extern from "Subgraphs.h":
	cdef cppclass _Subgraphs "tsppi::Subgraphs":
		_Graph graph
		_Graph getSubgraph(uint32_t)
		uint64_t numberOfNodes()
		uint64_t numberOfEdges()
		uint64_t numberOfSubgraphs()


#######################
#  PPI Graph wrapper  #
#######################

cdef extern from "PpiGraph.h":
	cdef cppclass _PpiGraph "tsppi::PpiGraph":
		_Graph graph
		_PpiGraph()
		_PpiGraph(const _PpiGraph&)
		string getPpiName() except +
		string getGeneName(uint32_t nodeid) except +
		vector[string] getAllGenes() except +

cdef extern from "PpiGraph.h":
	cdef cppclass _TsPpiGraph "tsppi::TsPpiGraph":
		_Subgraphs subgraphs
		_TsPpiGraph()
		_TsPpiGraph(const _TsPpiGraph&)
		string ppi_name
		uint64_t numberOfTissues() except +
		string getGeneName(uint32_t nodeid) except +
		string getTissueName(uint32_t tissue_id) except +

cdef extern from "graph_algos.h" namespace "tsppi::algo":
	vector[uint64_t] graph_degrees(const _Graph&)

cdef extern from "subgraph_algos.h" namespace "tsppi::algo":
	# aggregate graphs
	_Graph subgraph_edge_coexist_count_graph(const _Subgraphs&)
	_Graph subgraph_edge_correlation_graph(const _Subgraphs&)
	# degrees
	vector[vector[uint64_t]] subgraph_degrees(const _Subgraphs& subgraphs)
	vector[uint64_t] subgraph_max_degrees(const _Subgraphs& subgraphs)
	vector[uint64_t] subgraph_edge_exist_degrees(const _Subgraphs& subgraphs)
	# neighbor exist counts
	vector[uint64_t] subgraph_neighbor_min_exists_count(const _Subgraphs& subgraphs)
	vector[uint64_t] subgraph_neighbor_max_exists_count(const _Subgraphs& subgraphs)
	# CC (use the fastest variant)
	vector[vector[double]] subgraph_cc_neighbor_comb_vec(const _Subgraphs& subgraphs)
	# BW (use faster variant)
	vector[vector[double]] subgraph_betweenness_fast(const _Subgraphs& subgraphs)
	# PLP
	vector[_Partition] subgraph_PLP(const _Subgraphs& subgraphs)



###################################
#  Subgraph python wrapper class  #
###################################


cdef class Subgraphs:
	cdef _Subgraphs* _this

	cdef _set_this(self, _Subgraphs* instance):
		self._this = instance

	def numberOfNodes(self):
		return self._this.numberOfNodes()

	def numberOfEdges(self):
		return self._this.numberOfEdges()

	def numberOfSubgraphs(self):
		return self._this.numberOfSubgraphs()

	def getSubgraph(self, index):
		cdef _Graph _g = self._this.getSubgraph(index)
		g = Graph()
		g.setThis(new _Graph(_g))
		return g

	def getFullGraph(self):
		g = Graph()
		# Needs copy, since Graph() implements __dealloc__
		g.setThis(new _Graph(self._this.graph))
		return g

	# algorithms for subgraphs:
	# degrees
	def calcDegrees(self):
		return subgraph_degrees(deref(self._this))

	def calcGlobalDegrees(self):
		return graph_degrees(self._this.graph)

	def calcMaxDegrees(self):
		return subgraph_max_degrees(deref(self._this))

	def calcEdgeExistDegrees(self):
		return subgraph_edge_exist_degrees(deref(self._this))

	def calcNeighborMinExistsCount(self):
		return subgraph_neighbor_min_exists_count(deref(self._this))

	def calcNeighborMaxExistsCount(self):
		return subgraph_neighbor_max_exists_count(deref(self._this))

	def calcCC(self):
		# use the fastest CC implementation
		return subgraph_cc_neighbor_comb_vec(deref(self._this))

	def calcBW(self):
		# use the faster bw implementation
		return subgraph_betweenness_fast(deref(self._this))


#######################################
#  Tissue specific PPI graph wrapper  #
#######################################

cdef class PpiGraph:
	cdef _PpiGraph* _this

#	def __cinit__(self):
#		# otherwise create the graph
#		self._this = new _PpiGraph()

	cdef _set_this(self, _PpiGraph* instance):
		self._this = instance

	def getPpiName(self):
		return pystring(self._this.getPpiName())

	def getGraph(self):
		g = Graph()
		# Needs copy, since Graph() implements __dealloc__
		g.setThis(new _Graph(self._this.graph))
		return g

	def getGeneName(self, gene_id):
		return pystring(self._this.getGeneName(gene_id))

	def getAllGenes(self):
		vec = self._this.getAllGenes()
		result = []
		for i in vec:
			result.append(pystring(i))
		return result

	def getDegrees(self):
		return graph_degrees(self._this.graph)

cdef class TsPpiGraph:
	cdef _TsPpiGraph* _this

	cdef _set_this(self, _TsPpiGraph* instance):
		self._this = instance

	def getPpiName(self):
		return pystring(self._this.ppi_name)

	def getSubgraphs(self):
		sg = Subgraphs()
		sg._set_this(&self._this.subgraphs)
		return sg

	def getNumberOfTissues(self):
		return self._this.subgraphs.numberOfSubgraphs()

	def getGraph(self):
		g = Graph()
		# Needs copy, since Graph() implements __dealloc__
		g.setThis(new _Graph(self._this.subgraphs.graph))
		return g

	def getTsGraph(self, tissue):
		cdef _Graph _g = self._this.subgraphs.getSubgraph(tissue)
		g = Graph()
		# Needs copy, since Graph() implements __dealloc__
		g.setThis(new _Graph(_g))
		return g

	def getEdgeCorrelationGraph(self):
		cdef _Graph _g = subgraph_edge_correlation_graph(self._this.subgraphs)
		g = Graph()
		g.setThis(new _Graph(_g))
		return g

	def getEdgeCoexprCountGraph(self):
		cdef _Graph _g = subgraph_edge_coexist_count_graph(self._this.subgraphs)
		g = Graph()
		g.setThis(new _Graph(_g))
		return g

	def getGeneName(self, gene_id):
		return pystring(self._this.getGeneName(gene_id))

	def getTissueName(self, tissue_id):
		return pystring(self._this.getTissueName(tissue_id))

	def getAllTissues(self):
		result = []
		for i in range(self.getNumberOfTissues()):
			result.append(self.getTissueName(i))
		return result


	def getAllGenes(self):
		result = []
		for i in range(self._this.subgraphs.numberOfNodes()):
			result.append(self.getGeneName(i))
		return result

	def getDegrees(self):
		return graph_degrees(self._this.subgraphs.graph)

##########################################################
#  SQLite IO wrapper (loads tissue specific networks)	#
##########################################################

cdef extern from "SQLiteIO.h":
	cdef cppclass _SQLiteIO "tsppi::SQLiteIO":
		_SQLiteIO(string db) except +
		_PpiGraph load_ppi_graph(string ppi_name) except +
		_TsPpiGraph load_tsppi_graph(string ppi_name, string expr_name) except +

cdef class SQLiteIO:
	"""An undirected, optionally weighted graph"""
	cdef _SQLiteIO* _this

	def __cinit__(self, database):
		self._this = new _SQLiteIO(stdstring(database))

	def load_tsppi_graph(self, ppi_name, expr_name):
		p = stdstring(ppi_name)
		e = stdstring(expr_name)
		tsppi_graph = TsPpiGraph()
		cdef _TsPpiGraph tg = self._this.load_tsppi_graph(p, e)
		# FIXME: remove unnecessary copy
		tsppi_graph._set_this(new _TsPpiGraph(tg))
		return tsppi_graph

	def load_ppi_graph(self, ppi_name):
		p = stdstring(ppi_name)
		ppi_graph = PpiGraph()
		cdef _PpiGraph pg = self._this.load_ppi_graph(p)
		# FIXME: remove unnecessary copy
		ppi_graph._set_this(new _PpiGraph(pg))
		return ppi_graph



#######################################################################
#                        Modularity per cluster                       #
#######################################################################

cdef extern from "ModularityPerCluster.h":
	cdef cppclass _ModularityPerCluster "ModularityPerCluster":
		_ModularityPerCluster() except +
		vector[double] get_modularities(_Partition _zeta, _Graph _G) except +

cdef class ModularityPerCluster:
	cdef _ModularityPerCluster _this

	def getModularities(self, Partition zeta, Graph G):
		return self._this.get_modularities(zeta._this, deref(G._this))

