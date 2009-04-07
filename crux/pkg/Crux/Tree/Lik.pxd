# Forward declarations.
cdef class CL
cdef class Lik

from Crux.Character cimport Character
from Crux.Tree cimport Tree, Node, Ring
from Crux.CTMatrix cimport Alignment

from CxLik cimport *

cdef class CL:
    # Array of CL's, one for each polarity.  In the case of leaf nodes, only
    # the first element's cLMat/lnScale are used (but cache-related state is
    # separate for each polarity even for leaf nodes).
    cdef CxtLikCL cLs[2]

    cdef void prepare(self, unsigned polarity, unsigned nchars, unsigned dim, \
      unsigned ncomp) except *
    cdef void resize(self, unsigned polarity, unsigned nchars, unsigned dim, \
      unsigned ncomp) except *
    cdef void flush(self, unsigned polarity) except *

cdef class Lik:
    cdef Lik mate
    cdef readonly Character char_

    cdef readonly Tree tree
    cdef readonly Alignment alignment

    # Everything that needs to be accessible from pure C code is stored in lik,
    # so that worker threads can compute conditional likelihoods in parallel.
    cdef CxtLik *lik

    # Since the tree is actually unrooted, and truly rooting it during lnL() is
    # problematic (invalidates the tree's cache), the root CL is stored outside
    # the tree.  This happens to be why CL's parent must be CL rather than
    # Ring; there is no extant ring object associated with the root.
    cdef CL rootCL

    cdef void _init0(self, Tree tree) except *
    cdef void _init1(self, Tree tree, unsigned nchars, unsigned dim, \
      unsigned polarity) except *
    cdef void _init2(self, Alignment alignment) except *
    cdef CxtLikModel *_allocModel(self, unsigned ncat) except *
    cdef void _initModel(self, CxtLikModel *modelP, double weight, \
      bint catMedian)
    cdef void _decompModel(self, CxtLikModel *modelP) except *
    cdef void _deallocModel(self, CxtLikModel *modelP, unsigned model)

    cpdef Lik unpickle(self, str pickle)
    cdef void _dup(self, Lik lik) except *
    cpdef Lik dup(self)
    cpdef Lik clone(self)
    cpdef double getWNorm(self) except -1.0
    cpdef unsigned nmodels(self)
    cpdef unsigned addModel(self, double weight, unsigned ncat=*, \
      bint catMedian=*) except *
    cpdef delModel(self, unsigned model)
    cpdef double getWeight(self, unsigned model) except -1.0
    cpdef setWeight(self, unsigned model, double weight)
    cpdef double getRmult(self, unsigned model) except -1.0
    cpdef setRmult(self, unsigned model, double rmult)
    cpdef list getRclass(self, unsigned model)
    cpdef setRclass(self, unsigned model, list rclass, list rates=*)
    cpdef unsigned getNrates(self, unsigned model) except 0
    cpdef double getRate(self, unsigned model, unsigned i) except -1.0
    cpdef setRate(self, unsigned model, unsigned i, double rate)
    cpdef double getFreq(self, unsigned model, unsigned i) except -1.0
    cpdef setFreq(self, unsigned model, unsigned i, double freq)
    cpdef double getAlpha(self, unsigned model)
    cpdef unsigned getNcat(self, unsigned model) except *
    cpdef bint getCatMedian(self, unsigned model) except *
    cpdef setAlpha(self, unsigned model, double alpha)
    cdef void _planAppend(self, CxeLikStep variant, CL parentCL, CL childCL, \
      double edgeLen) except *
    cdef void _planRecurse(self, Ring ring, CL parent, unsigned nSibs,
      double edgeLen) except *
    cdef void _plan(self, Node root) except *
    cpdef prep(self)
    cdef void _prep(self, Node root) except *
    cpdef double lnL(self, Node root=*) except 1.0
    cpdef list siteLnLs(self, Node root=*)
