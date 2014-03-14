#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import math
#from exceptions import EOFError

typedef_tag, term_tag = "[Typedef]", "[Term]"

# TODO: this is duplicate from go_similarity
def name2id(name):
    return int(name.split(":")[1])

def after_colon(line):
    # macro for getting anything after the :
    return line.split(":", 1)[1].strip()

def read_until(handle, start):
    # read each line until it has a certain start, and then puts
    # the start tag back
    while 1:
        pos = handle.tell()
        line = handle.readline()
        if not line:
            break
        if line.startswith(start):
            handle.seek(pos)
            return
    raise EOFError("%s tag cannot be found" % start)


class OBOReader:
    """
    parse obo file, usually the most updated can be downloaded from
    http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

    >>> reader = OBOReader()
    >>> for rec in reader:
            print rec

    """

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        #try:
        self._handle = open(obo_file)
        #except:
        #    print("download obo file first\n "
        #                         "[http://geneontology.org/ontology/"
        #                         "obo_format_1_2/gene_ontology.1_2.obo]",
        #                         file=sys.stderr)
        #    sys.exit(1)

    def __iter__(self):

        line = self._handle.readline()
        if not line.startswith(term_tag):
            read_until(self._handle, term_tag)
        while 1:
            yield self.next()

    def next(self):

        lines = []
        line = self._handle.readline()
        if not line or line.startswith(typedef_tag):
            raise StopIteration

        # read until the next tag and save everything in between
        while 1:
            pos = self._handle.tell()   # save current postion for roll-back
            line = self._handle.readline()
            if not line or (line.startswith(typedef_tag)
                            or line.startswith(term_tag)):
                self._handle.seek(pos)  # roll-back
                break
            lines.append(line)

        rec = GOTerm()
        for line in lines:
            if line.startswith("id:"):
                rec.id = after_colon(line)
            if line.startswith("alt_id:"):
                rec.alt_ids.append(after_colon(line))
            elif line.startswith("name:"):
                rec.name = after_colon(line)
            elif line.startswith("namespace:"):
                rec.namespace = after_colon(line)
            elif line.startswith("is_a:"):
                rec._parents.append(after_colon(line).split()[0])
            elif (line.startswith("is_obsolete:") and
                  after_colon(line) == "true"):
                rec.is_obsolete = True

        return rec


class GOTerm:
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(self):
        self.id = ""                # GO:xxxxxx
        self.int_id = 0
        self.name = ""              # description
        self.namespace = ""         # BP, CC, MF
        self._parents = []          # is_a basestring of parents
        self.parents = []           # parent records
        self.children = []          # children records
        self.level = -1             # distance from root node
        self.is_obsolete = False    # is_obsolete
        self.alt_ids = []           # alternative identifiers

    def __str__(self):
        obsolete = "obsolete" if self.is_obsolete else ""
        return "%s\tlevel-%02d\t%s [%s] %s" % (self.id, self.level, self.name,
                                               self.namespace, obsolete)

    def __repr__(self):
        return "GOTerm('%s')" % (self.id)

    def has_parent(self, term):
        for p in self.parents:
            if p.id == term or p.has_parent(term):
                return True
        return False

    def has_child(self, term):
        for p in self.children:
            if p.id == term or p.has_child(term):
                return True
        return False

    def get_all_parents(self):
        all_parents = set()
        for p in self.parents:
            all_parents.add(p.id)
            all_parents |= p.get_all_parents()
        return all_parents

    def get_all_children(self):
        all_children = set()
        for p in self.children:
            all_children.add(p.id)
            all_children |= p.get_all_children()
        return all_children


class GODag:

    # a set of integer IDs of all GO terms present in the DAG
    terms = set()
    # a dict of ID->GoTerm node (full struct)
    nodes = dict()
    # a dict of ID->level
    level = dict()
    # a dict of ID->[set of ID]
    parents = dict()
    # a dict of ID->[set of ID]
    children = dict()
    # all ancestors
    ancestors = dict()
    # all decendents
    decendents = dict()
    # a dict of ID-> Information content
    IC = dict()

    def __init__(self, obo_file="gene_ontology.1_2.obo", only_namespace=None, load_obsolete=False):
        self.load_obo_file(obo_file, only_namespace)

    def load_obo_file(self, obo_file, only_namespace=None, load_obsolete=False):

        print("load obo file %s" % obo_file, file=sys.stderr)
        obo_reader = OBOReader(obo_file)
        # parse all file entries:
        for rec in obo_reader:
            # filter by namespace (in case that option is set)
            if not only_namespace is None and rec.namespace != only_namespace:
                continue
            # filter out obsolete (done by default)
            if not load_obsolete and rec.is_obsolete:
                continue
            # get and initialize integer id
            rec.int_id = name2id(rec.id)
            # set data records
            self.nodes[rec.int_id] = rec
            # and add the integer id to the set of all IDs
            self.terms.add(rec.int_id)
            # let the alternative IDs refer to the real (current) ones
            # TODO: do I need to add them to the set of all terms?
            for alt in rec.alt_ids:
                self.nodes[name2id(alt)] = rec

        # initialize the internal data structure for the DAG
        self.populate_terms()
        print(str(len(self.terms)) + " nodes imported", file=sys.stderr)

    def populate_terms(self):

        # recusively find the level of a term
        def depth(term):
            if not term in self.level:
                level = 0
                # if there are no parents -> root at level 0
                if len(self.parents[term]) != 0:
                    level = min(depth(t) for t in self.parents[term]) + 1
                # save the level
                self.level[term] = level
            return self.level[term]

        # recursively find all decendents for a term
        def rec_dec(term):
            if not term in self.decendents:
                dec = set()
                for c in self.children[term]:
                    dec = dec.union(rec_dec(c))
                self.decendents[term] = dec
            else:
                dec = self.decendents[term]
            # add self before returning
            # create a copy first
            dec = set(dec)
            dec.add(term)
            return dec

        # recursively find all ancesotrs for a term
        def rec_anc(term):
            if not term in self.ancestors:
                anc = set()
                for p in self.parents[term]:
                    anc = anc.union(rec_anc(p))
                self.ancestors[term] = anc
            else:
                anc = self.ancestors[term]
            # add self before returning
            anc = set(anc)  # create a copy first
            anc.add(term)
            return anc

        # initialize the parents data structure (a dict from id->set of ids)
        self.roots = set()
        for rec in self.nodes.values():
            parents = set([name2id(x) for x in rec._parents])
            # check if this is a root node, if yes then add to roots
            if len(parents) == 0 and not rec.is_obsolete:
                self.roots.add(rec.int_id)
            # save the parents
            self.parents[rec.int_id] = parents

        # initialize children
        for term in self.terms:
            self.children[term] = set()

        # populate the `children` data
        for term in self.terms:
            # populate `children`
            for p in self.parents[term]:
                self.children[p].add(term)

        # initialize the ancestors and decendents
        for term in self.terms:
            rec_anc(term)
            rec_dec(term)

        # initialize the `1evel` data
        for term in self.terms:
            # populate `level`
            if not term in self.level:
                depth(term)


    def _to_id(self, term):
        if type(term) is str:
            return name2id(term)
        else:
            return term

    def has_term(self, term):
        term = self._to_id(term)
        return term in self.terms

    def paths_to_top(self, term, verbose=False):
        """ Returns all possible paths to the root node

            Each path includes the term given. The order of the path is
            top -> bottom, i.e. it starts with the root and ends with the
            given term (inclusively).

            Parameters:
            -----------
            - term:
                the id of the GO term, where the paths begin (i.e. the
                accession 'GO:0003682')

            Returns:
            --------
            - a list of lists of GO Terms
        """
        # error handling consistent with original authors
        term = self._to_id(term)
        if not term in self:
            print("Term %s not found!" % term, file=sys.stderr)
            return

        def _paths_to_top_recursive(term):
            if self.level[term] == 0:
                return [[term]]
            paths = []
            for parent in self.parents[term]:
                top_paths = _paths_to_top_recursive(parent)
                for top_path in top_paths:
                    top_path.append(term)
                    paths.append(top_path)
            return paths

        return _paths_to_top_recursive(term)


    def term_probability(self, associations):
        # initialize annotation counting dict
        anno = dict()
        for t in self.terms:
            anno[t] = 0

        # count number of annotations per term
        for key, terms in associations.items():
            for t in terms:
                t = self._to_id(t)
                if not t in self.terms:
                    print("term not found: " + str(t))
                    continue
                anno[t] = anno[t] + 1

        # determine cummulative frequency
        self.freq = dict()
        for t in self.terms:
            child_freq = sum(anno[s] for s in self.decendents[t])
            self.freq[t] = anno[t] + child_freq

        # from frequency determine probability
        self.p = dict()
        total_freq = sum(self.freq[r] for r in self.roots)
        for term, f in self.freq.items():
            self.p[term] = f / total_freq


    def term_IC(self, terms=None):
        # check that the probabilities (from the frequencies)
        # have already been determined
        if len(self.p) == 0:
            raise Exception

        # get the set of terms to use
        if terms is None:
            terms = self.terms
        else:
            # convert to integer ids
            if type(terms[0]) is str():
                terms = set([name2id(t) for t in terms])
            terms = self.terms.intersect(terms)

        # initialize and fill the information content
        self.IC = dict()
        for i in terms:
            if not self.p[i] > 0.0:
                self.IC[i] = 0.0
            else:
                self.IC[i] = - math.log10(self.p[i])


    def get_lca_option1(self, term1, term2):
        """
        Finds the lowest common ancestor for the two given terms.
        """
        # search upwards level by level, till an overlap is found
        term1 = self._to_id(term1)
        term2 = self._to_id(term2)

        # get the difference of levels of the two terms
        # and define `left` as the term with biggest level (furthest from root)
        left = set()
        right = set()
        l1 = self.level[term1]
        l2 = self.level[term2]
        if l1 < l2:
            diff_levels = l2 - l1
            left.add(term2)
            right.add(term1)
        else:
            diff_levels = l1 - l2
            left.add(term1)
            right.add(term2)

        # extend smaller set
        for i in range(diff_levels):
            new_level = set()
            # union all parents from all current nodes
            for t in left:
                new_level = new_level.union(self.parents[t])
            left = new_level


        # while there is no overlap
        LCAs = left.intersection(right)
        while len(LCAs) == 0:
            new_level = set()
            # union all parents from all current nodes
            for t in left:
                new_level = new_level.union(self.parents[t])
            left = new_level
            new_level = set()
            # union all parents from all current nodes
            for t in right:
                new_level = new_level.union(self.parents[t])
            right = new_level
            # get new intersection
            LCAs = left.intersection(right)

        return LCAs

    def get_lca(self, term1, term2):
        # TODO: try different versions of finding the LCA
        return self.get_lca_option1(term1, term2)
