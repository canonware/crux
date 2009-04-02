from Crux.Taxa cimport Taxon
from Crux.Tree cimport Tree, Node, Edge, Ring
from Crux.Tree.Bipart cimport Vec, Bipart

cdef class Trprob:
    def __init__(self, unsigned nruns):
        self._tree = None

        self._nobs = [0] * nruns
        self._trees = []

    cdef void _observe(self, unsigned run, Tree tree) except *:
        self._nobs[run] += 1
        self._trees.append(tree)

    cdef uint64_t getNobs(self):
        return len(self._trees)
    property nobs:
        """
            Total number of observations across all runs.
        """
        def __get__(self):
            return self.getNobs()

    cdef double getMse(self):
        cdef double mse, mean, diff
        cdef unsigned i, nruns
        cdef uint64_t nobs

        mse = 0.0
        nobs = len(self._trees)
        nruns = len(self._nobs)
        mean = <double>nobs / <double>nruns
        for 0 <= i < nruns:
            diff = <double>self._nobs[i] - mean
            mse += (diff * diff) / <double>nobs
        return mse
    property mse:
        """
            Mean squared error (MSE) of observation frequencies among runs.
        """
        def __get__(self):
            return self.getMse()

    cdef Tree getTree(self):
        cdef Edge edge, edgeI
        cdef Tree tree, treeI
        cdef Bipart bipart, bipartI
        cdef dict t2e
        cdef Vec vec, vecI
        cdef uint64_t i, ntrees
        cdef unsigned j
        cdef list nodes
        cdef Node node

        if self._tree is None:
            # Create tree with 0-length branches.
            self._tree = self._trees[0].dup()
            for edge in self._tree.getEdges():
                edge.length = 0.0
            bipart = self._tree.getBipart()

            # Create a lookup table of taxon-->edge for leaf edges.
            t2e = {}
            nodes = self._tree.getNodes()
            for 0 <= j < len(nodes):
                node = <Node>nodes[j]
                if node.getDegree() == 1:
                    t2e[node.taxon] = node.ring.edge

            # Iterate over sampled trees.
            ntrees = len(self._trees)
            for 0 <= i < ntrees:
                treeI = <Tree>self._trees[i]
                bipartI = treeI.getBipart()

                # Add internal branch lengths.
                for 0 <= j < len(bipart.edgeVecs):
                    vec = <Vec>bipart.edgeVecs[j]
                    vecI = <Vec>bipartI.edgeVecs[j]
                    vec.edge.length += vecI.edge.length / <double>ntrees

                # Add leaf branch lengths.
                nodes = treeI.getNodes()
                for 0 <= j < len(nodes):
                    node = <Node>nodes[j]
                    if node.getDegree() == 1:
                        edgeI = node.ring.edge
                        edge = <Edge>t2e[node.taxon]
                        edge.length += edgeI.length / <double>ntrees

        return self._tree
    property tree:
        """
            Tree with mean branch lengths.
        """
        def __get__(self):
            return self.getTree()

    cdef unsigned getNruns(self):
        cdef unsigned nruns

        nruns = 0
        for 0 <= i < len(self._nobs):
            if self._nobs[i] > 0:
                nruns += 1
        return nruns
    property nruns:
        """
            Number of runs in which tree was observed at least once.
        """
        def __get__(self):
            return self.getNruns()

cdef class Part:
    def __init__(self, Vec vec, unsigned nruns):
        self.vec = vec
        self._nobs = [0] * nruns
        self._edges = []

    cdef void _observe(self, unsigned run, Edge edge) except *:
        self._nobs[run] += 1
        self._edges.append(edge)

    cdef uint64_t getNobs(self):
        return len(self._edges)
    property nobs:
        """
            Total number of observations across all runs.
        """
        def __get__(self):
            return self.getNobs()

    cdef double getMse(self):
        cdef double mse, mean, diff
        cdef unsigned i, nruns
        cdef uint64_t nobs

        mse = 0.0
        nobs = len(self._edges)
        nruns = len(self._nobs)
        mean = <double>nobs / <double>nruns
        for 0 <= i < nruns:
            diff = <double>self._nobs[i] - mean
            mse += (diff * diff) / <double>nobs
        return mse
    property mse:
        """
            Mean squared error (MSE) of observation frequencies among runs.
        """
        def __get__(self):
            return self.getMse()

    cdef double getMean(self):
        cdef double mean
        cdef uint64_t i, nobs
        cdef Edge edge

        mean = 0.0
        nobs = len(self._edges)
        for 0 <= i < nobs:
            edge = <Edge>self._edges[i]
            mean += edge.length / <double>nobs
        return mean
    property mean:
        """
            Branch length mean.
        """
        def __get__(self):
            return self.getMean()

    cdef double getVar(self):
        cdef double var, mean, diff
        cdef uint64_t i, nobs
        cdef Edge edge

        var = 0.0
        mean = self.getMean()
        nobs = len(self._edges)
        for 0 <= i < nobs:
            edge = <Edge>self._edges[i]
            diff = (edge.length - mean)
            var += (diff * diff) / <double>nobs
        return var
    property var:
        """
            Branch length variance.
        """
        def __get__(self):
            return self.getVar()

    cdef unsigned getNruns(self):
        cdef unsigned nruns

        nruns = 0
        for 0 <= i < len(self._nobs):
            if self._nobs[i] > 0:
                nruns += 1
        return nruns
    property nruns:
        """
            Number of runs in which partition was observed at least once.
        """
        def __get__(self):
            return self.getNruns()

cdef class Sumt:
    def __init__(self, list treeLists):
        self._treeLists = treeLists
        self._trprobs = None
        self._parts = None

    cdef void _summarizeTrprobs(self) except *:
        cdef dict d
        cdef unsigned nruns, i
        cdef uint64_t j
        cdef list trees, decos
        cdef Tree tree
        cdef Bipart bipart
        cdef Trprob trprob
        cdef object deco

        d = {}
        nruns = len(self._treeLists)
        for 0 <= i < nruns:
            trees = self._treeLists[i]
            for 0 <= j < len(trees):
                tree = <Tree>trees[j]
                bipart = tree.getBipart()
                if bipart in d:
                    trprob = <Trprob>d[bipart]
                else:
                    trprob = Trprob(nruns)
                    d[bipart] = trprob
                trprob._observe(i, tree)

        decos = [(len((<Trprob>d[bipart])._trees), bipart, <Trprob>d[bipart]) \
          for bipart in d]
        decos.sort(reverse=True)
        self._trprobs = [deco[2] for deco in decos]
    cdef list getTrprobs(self):
        if self._trprobs is None:
            self._summarizeTrprobs()
            return self._trprobs
    property trprobs:
        """
            List of all observed trees, ordered from most- to least-represented.
            Each list element is a Trprob instance.
        """
        def __get__(self):
            return self.getTrprobs()

    cdef void _summarizeParts(self) except *:
        cdef dict d
        cdef unsigned nruns, i, k
        cdef uint64_t j
        cdef list trees, vecs, decos
        cdef Tree tree
        cdef Bipart bipart
        cdef Vec vec
        cdef Part part
        cdef object deco

        d = {}
        nruns = len(self._treeLists)
        for 0 <= i < nruns:
            trees = self._treeLists[i]
            for 0 <= j < len(trees):
                tree = <Tree>trees[j]
                bipart = tree.getBipart()
                vecs = bipart.edgeVecs
                for 0 <= k < len(vecs):
                    vec = <Vec>vecs[k]
                    if vec in d:
                        part = <Part>d[vec]
                        assert vec.cmp(part.vec) == 0
                    else:
                        part = Part(vec, nruns)
                        d[vec] = part
                    part._observe(i, vec.edge)

        decos = [(<Part>d[vec].nobs, vec, <Part>d[vec]) for vec in d]
        decos.sort(reverse=True)
        self._parts = [deco[2] for deco in decos]
    cdef list getParts(self):
        if self._parts is None:
            self._summarizeParts()
        return self._parts
    property parts:
        """
            List of all observed partitions, ordered from most- to
            least-represented.  Each list element is a Part instance.
        """
        def __get__(self):
            return self.getParts()

    cdef bint _incorpPart(self, Tree tree, unsigned ntaxa, Node node, \
      dict e2v, Vec incorp, double support) except *:
        cdef Ring ring
        cdef Edge edge
        cdef Vec vec, v
        cdef list compat, incompat
        cdef unsigned i
        cdef Node nodeCompat, sib

        # Create a temporary vector to be used for set operations.
        v = Vec(None, ntaxa)

        # Partition edges into two sets, according to whether their vecs are
        # subsets of vec.
        compat = []
        incompat = []
        for ring in node.ring:
            edge = ring.edge
            vec = <Vec>e2v[edge]
            v.reset()
            v.merge(incorp)
            v.merge(vec)
            if v.cmp(incorp) == 0:
                compat.append(edge)
            else:
                incompat.append(edge)

        # Merge bipartitions in the compat set to check whether they completely
        # comprise the incorp bipartition.
        v.reset()
        for 0 <= i < len(compat):
            edge = <Edge>compat[i]
            vec = <Vec>e2v[edge]
            v.merge(vec)
        if v.cmp(incorp) != 0:
            return False
        assert len(compat) > 0
        assert len(incompat) > 0

        # Incorporate.
        nodeCompat = Node(tree)
        for 0 <= i < len(compat):
            edge = <Edge>compat[i]
            sib = edge.ring.node
            if sib is node:
                sib = edge.ring.other.node
            edge.detach()
            edge.attach(nodeCompat, sib)
        edge = Edge(tree)
        edge.aux = support
        edge.attach(node, nodeCompat)
        return True

    cpdef Tree getConTree(self, minSupport=0.0):
        """
            Generate a consensus tree.  Iteratively process bipartitions with
            at least minSupport posterior support from most- to
            least-supported and greedily incorporate them if they are
            compatible with previously incorporated bipartitions.

            To get a strict consensus tree, use minSupport=1.0.  To get a
            majority-rule consensus tree, use minSupport=0.5.

            Compute each branch length by averaging the lengths across all
            trees that include the branch.

            Store bipartition support values (commonly referred to as "node
            support values") in the edges' aux attributes.
        """
        cdef Tree tree
        cdef Node base, node
        cdef Edge edge
        cdef uint64_t N
        cdef list taxa, nodes, trees
        cdef unsigned ntaxa, nruns, i, j
        cdef uint64_t k
        cdef Taxon taxon
        cdef Part part
        cdef double support
        cdef dict e2v, v2e, v2nobs
        cdef Bipart bipart
        cdef Vec vec

        # Create a star tree with 0-length branches.
        tree = Tree(None, None, False)
        base = Node(tree)
        tree.setBase(base)
        taxa = (<Tree>(<list>self._treeLists[0])[0]).getTaxa()
        ntaxa = len(taxa)
        for 0 <= i < ntaxa:
            taxon = <Taxon>taxa[i]
            edge = Edge(tree)
            edge.aux = 1.0
            node = Node(tree)
            node.setTaxon(taxon)
            edge.attach(base, node)

        N = 0
        for 0 <= i < len(self._treeLists):
            N += len(<list>self._treeLists[i])

        if self._parts is None:
            self._summarizeParts()
        # Iteratively incorporate bipartitions.
        for 0 <= i < len(self._parts):
            part = <Part>self._parts[i]
            support = <double>part.getNobs() / <double>N
            if support < minSupport:
                # This and all subsequent bipartitions lack adequate support
                # for incorporation.
                break

            # Create a lookup table of edge-->vec.
            bipart = Bipart(tree, True)
            e2v = {}
            for 0 <= j < len(bipart.edgeVecs):
                vec = <Vec>bipart.edgeVecs[j]
                e2v[vec.edge] = vec

            # Iterate over incompletely resolved nodes and try to incorporate
            # bipartition.
            nodes = tree.getNodes()
            for 0 <= j < len(nodes):
                node = <Node>nodes[j]
                if node.getDegree() > 3:
                    if self._incorpPart(tree, ntaxa, node, e2v, part.vec, \
                      support):
                        # Success.
                        break

            if len(tree.getEdges()) == (2*ntaxa)-3:
                # Tree is fully resolved.
                break


        # Create a lookup tabel of vec-->edge for the final tree topology.
        bipart = Bipart(tree, True)
        v2e = {}
        for 0 <= i < len(bipart.edgeVecs):
            vec = <Vec>bipart.edgeVecs[i]
            v2e[vec] = vec.edge

        # Create a lookup table of vec-->nobs.
        v2nobs = {}
        for 0 <= i < len(self._parts):
            part = <Part>self._parts[i]
            v2nobs[part.vec] = part.getNobs()
        # Leaf edges are in every tree.
        for vec in v2e:
            if vec not in v2nobs:
                v2nobs[vec] = N

        # Set branch lengths.
        nruns = len(self._treeLists)
        for 0 <= i < nruns:
            trees = self._treeLists[i]
            for 0 <= k < len(trees):
                bipart = Bipart(<Tree>trees[k], True)
                for 0 <= j < len(bipart.edgeVecs):
                    vec = <Vec>bipart.edgeVecs[j]
                    if vec in v2e:
                        edge = <Edge>v2e[vec]
                        edge.length += vec.edge.length / <double>v2nobs[vec]

        return tree
