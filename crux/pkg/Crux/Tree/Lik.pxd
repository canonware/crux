from Crux.Character cimport Character
from Crux.Tree cimport Tree, Node, Ring
from Crux.CTMatrix cimport Alignment

from CxLik cimport *

cdef class CL:
    # Parent which most recently used this cache when computing its conditional
    # likelihood.  All siblings must still refer to the parent, the number of
    # siblings must still be nSibs, and the branches separating them from the
    # parent must remain the same length, in order for the parent's cache to be
    # potentially valid.  (The model of evolution must also remain unchanged.)
    cdef CL parent
    cdef unsigned nSibs
    cdef double edgeLen

    # Vector of CxtLikCL structures.  For internal nodes, there is one
    # vector element for each model in the mixture.  For leaf nodes,  there is
    # only one vector element, which contains the same encoding of the
    # character data for the associated taxon, regardless of model.
    #
    # This vector is incrementally expanded as necessary, but never shrunk,
    # under the presumption that the vector will likely re-expand in the
    # future.
    cdef CxtLikCL *vec
    cdef unsigned vecMax # Number of usable elements.

    cdef void prepare(self, unsigned nchars, unsigned dim, unsigned ncat, \
      unsigned nmodels) except *
    cdef void dupModel(self, unsigned nchars, unsigned dim, unsigned ncat, \
      unsigned to, unsigned fr) except *

cdef class Lik:
    cdef readonly Character char_

    # Every model of evolution is assigned a unique serial number, so that
    # cached conditional likelihoods can be associated with them.  This field
    # is the source for serial number assignments.  Serial numbers are in
    # [1..2^64).
    cdef uint64_t sn

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

    cdef uint64_t _assignSn(self)
    cdef void _allocModel(self, CxtLikModel *model) except *
    cdef void _initModel(self, CxtLikModel *model)
    cdef void _reassignModel(self, CxtLikModel *model) except *
    cdef void _deallocModel(self, CxtLikModel *model)

    cpdef Lik dup(self)
    cpdef unsigned getNcat(self)
    cpdef unsigned nmodels(self)
    cpdef unsigned addModel(self)
    cpdef dupModel(self, unsigned to, unsigned fr, bint dupCLs=*)
    cpdef delModel(self)
    cpdef double getWeight(self, unsigned model)
    cpdef setWeight(self, unsigned model, double weight)
    cpdef list getRclass(self, unsigned model)
    cpdef setRclass(self, unsigned model, list rclass, list rates=*)
    cpdef unsigned getNrates(self, unsigned model) except 0
    cpdef double getRate(self, unsigned model, unsigned i) except -1.0
    cpdef setRate(self, unsigned model, unsigned i, double rate)
    cpdef double getFreq(self, unsigned model, unsigned i) except -1.0
    cpdef setFreq(self, unsigned model, unsigned i, double freq)
    cpdef double getAlpha(self, unsigned model)
    cpdef setAlpha(self, unsigned model, double alpha)
    cdef void _planAppend(self, unsigned model, CxeLikStep variant, \
      CL parentCL, CL childCL, double edgeLen) except *
    cdef void _planRecurse(self, Ring ring, CL parent, unsigned nSibs,
      double edgeLen) except *
    cdef void _plan(self, Node root) except *
    cdef void _prep(self, Node root) except *
    cpdef double lnL(self, Node root=*) except 1.0
    cpdef list siteLnLs(self, Node root=*)
