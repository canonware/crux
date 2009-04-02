from libc cimport *
from Crux.Tree cimport Tree, Node, Edge
from Crux.Tree.Bipart cimport Vec

cdef class Trprob:
    cdef Tree _tree
    cdef list _nobs       # Number of times observed in each run.
    cdef list _trees      # List of trees, used to compute mean brlens.

    cdef void _observe(self, unsigned run, Tree tree) except *

    cdef uint64_t getNobs(self)
    # property nobs
    cdef double getMse(self)
    # property mse
    cdef Tree getTree(self)
    # property tree
    cdef unsigned getNruns(self)
    # property nruns

cdef class Part:
    cdef readonly Vec vec
    cdef list _nobs       # Number of times observed in each run.
    cdef list _edges      # List of edges, used to compute mean/variance.

    cdef void _observe(self, unsigned run, Edge edge) except *

    cdef uint64_t getNobs(self)
    # property nobs
    cdef double getMse(self)
    # property mse
    cdef double getMean(self)
    # property mean
    cdef double getVar(self)
    # property var
    cdef unsigned getNruns(self)
    # property nruns

cdef class Sumt:
    cdef list _treeLists
    cdef list _trprobs
    cdef list _parts

    cdef void _summarizeTrprobs(self) except *
    cdef list getTrprobs(self)
    # property trprobs

    cdef void _summarizeParts(self) except *
    cdef list getParts(self)
    # property parts

    cdef bint _incorpPart(self, Tree tree, unsigned ntaxa, Node node, \
      dict e2v, Vec incorp, double support) except *
    cpdef Tree getConTree(self, minSupport=*)
