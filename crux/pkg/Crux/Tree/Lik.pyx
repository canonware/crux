import Crux.Config as Config
from Crux.Character cimport Character
from Crux.Taxa cimport Taxon
from Crux.Tree cimport Tree, Node, Edge, Ring
from Crux.CTMatrix cimport Alignment

from Cx cimport CxNcpus
from libc cimport *
from libm cimport *
from CxLik cimport *

# Conditional likelihoods are allocated such that they are aligned to cacheline
# boundaries, and striping is configured to avoid false cacheline sharing among
# worker threads.  Over-estimating cacheline size is okay (within reason), as
# long as the estimate is a multiple of the true cacheline size.
cdef enum:
    cacheLine = 64

cdef class CL:
    """
        Conditional likelihood, associated with a ring.
    """
    def __cinit__(self):
        self.vec = NULL

    def __dealloc__(self):
        cdef unsigned i

        if self.vec != NULL:
            for 0 <= i < self.vecMax:
                if self.vec[i].mat != NULL:
                    free(self.vec[i].mat)
                    self.vec[i].mat = NULL
            free(self.vec)
            self.vec = NULL

    def __init__(self, unsigned nchars, unsigned dim, unsigned nmodels):
        self.parent = None
        self.nSibs = 0

        self.vec = <CxtLikCL *>calloc(nmodels, sizeof(CxtLikCL));
        if self.vec == NULL:
            raise MemoryError("Error allocating vector")
        self.vecMax = nmodels

        for 0 <= i < nmodels:
            if posix_memalign(<void **>&self.vec[i].mat, cacheLine, nchars * \
              dim * sizeof(double)):
                raise MemoryError("Error allocating matrix")

    cdef void prepare(self, unsigned nchars, unsigned dim, unsigned nmodels) \
      except *:
        cdef CxtLikCL *vec
        cdef unsigned i

        if self.vecMax < nmodels:
            vec = <CxtLikCL *>realloc(self.vec, nmodels * \
              sizeof(CxtLikCL))
            if vec == NULL:
                raise MemoryError("Error reallocating vector")
            self.vec = vec
            memset(&self.vec[self.vecMax], 0, (nmodels - self.vecMax) * \
              sizeof(CxtLikCL))
            for self.vecMax <= i < nmodels:
                if posix_memalign(<void **>&self.vec[i].mat, cacheLine, \
                  nchars * dim * sizeof(double)):
                    raise MemoryError("Error allocating matrix")
            self.vecMax = nmodels

    cdef void dupModel(self, unsigned nchars, unsigned dim, unsigned to,
      unsigned fr) except *:
        if to >= self.vecMax:
            self.prepare(nchars, dim, to + 1)
        if fr >= self.vecMax:
            self.prepare(nchars, dim, fr + 1)

        self.vec[to].sn = self.vec[fr].sn
        memcpy(&self.vec[to].mat, &self.vec[fr].mat, self.nchars * self.dim * \
          sizeof(double))

cdef class Lik:
    """
        Lik manages log-likelihood computation using a tree, an alignment, and
        a weighted mixture of one or more models of evolution.  Initially, the
        models are equally weighted.

        Each model of evolution consists of a Q matrix and a Gamma-distributed
        mutation rates parameter (alpha).  Q is conceptually a matrix product
        of the following form:

            Q = R * Pi
                 _            _     _            _
                |  -  a  b  c  |   | pA  0  0  0  |
                |              |   |              |
                |  a  -  d  e  |   |  0 pC  0  0  |
              = |              | * |              |,
                |  b  d  -  f  |   |  0  0 pG  0  |
                |              |   |              |
                |_ c  e  f  - _|   |_ 0  0  0 pT _|

        where {a,b,c,d,e,f} are relative substitution rates, and {pA,pC,pG,pT}
        are frequencies.

        During Q matrix computation:
        * The relative rates are rescaled to fix one rate to 1.0 (f in the
          example above).
        * The elements on R's diagonal are set such that each row sums to 0.
        * The frequencies are rescaled to sum to 1.

        By default:
        * The relative rates are all equal (1.0).
        * The frequencies are all equal (0.25 in the example above).
        * Gamma-distributed rates are disabled (alpha is infinite).

        When computing substitution probabilities for a branch of a particular
        length, numerical methods are used to compute Q's
        eigenvector/eigenvalue decomposition, which can then be used to
        compute:

                     Q*mu*t
        P(Q,mu,t) = e
    """

    def __cinit__(self):
        self.lik = NULL

    def __dealloc__(self):
        cdef CxtLik *lik
        cdef CxtLikModel *model
        cdef unsigned i

        if self.lik != NULL:
            lik = self.lik
            for 0 <= i < lik.modelsMax:
                model = &lik.models[i]
                self._deallocModel(model)
            free(lik)
            self.lik = NULL

    def __init__(self, Tree tree, Alignment alignment, unsigned nmodels=1,
      unsigned ncat=4):
        cdef unsigned stripeWidth, stepsMax, i

        if tree.rooted:
            raise ValueError("Tree must be unrooted")

        self.char_ = alignment.charType.get()

        # Configure striping.
        if Config.threaded:
            # Stripe width should always be a multiple of 2 (assuming 64-byte
            # cachelines) in order to avoid false cacheline sharing between
            # threads.
            stripeWidth = cacheLine / (4 * sizeof(double))

            # Limit the number of stripes to the number of CPUs.
            while alignment.nchars / stripeWidth > CxNcpus:
                stripeWidth *= 2
        else:
            stripeWidth = alignment.nchars

        # Pad the alignment, if necessary.
        if alignment.nchars % stripeWidth != 0:
            alignment.pad(self.char_.val2code(self.char_.any), stripeWidth - \
              (alignment.nchars % stripeWidth))

        self.sn = 0

        self.tree = tree
        self.alignment = alignment

        self.lik = <CxtLik *>calloc(1, sizeof(CxtLik))
        if self.lik == NULL:
            raise MemoryError("Error allocating lik")

        self.lik.dim = self.char_.nstates()
        self.lik.ncat = ncat
        self.lik.nchars = self.alignment.nchars
        self.lik.charFreqs = self.alignment.freqs
        self.lik.stripeWidth = stripeWidth
        self.lik.nstripes = self.lik.nchars / self.lik.stripeWidth
        assert self.lik.nstripes <= CxNcpus

        self.lik.models = <CxtLikModel *>calloc(nmodels, sizeof(CxtLikModel))
        if self.lik.models == NULL:
            raise MemoryError("Error allocating models")
        self.lik.modelsLen = nmodels
        self.lik.modelsMax = nmodels
        for 0 <= i < self.lik.modelsMax:
            self._allocModel(&self.lik.models[i])
            self._initModel(&self.lik.models[i])

        stepsMax = ((2 * self.alignment.ntaxa) - 2) * self.lik.modelsMax
        self.lik.steps = <CxtLikStep *>malloc(stepsMax * sizeof(CxtLikStep))
        if self.lik.steps == NULL:
            raise MemoryError("Error allocating steps")
        self.lik.stepsLen = 0
        self.lik.stepsMax = stepsMax

        self.rootCL = CL(self.lik.nchars, self.lik.dim, self.lik.modelsMax)
        # Set each model's pointer to the CxtLikCL that is actually stored in
        # rootCL, so that lnL computation can be completed in the C code.  This
        # is also done in _plan() for models that are added later on.
        for 0 <= i < self.lik.modelsMax:
            self.lik.models[i].cL = &self.rootCL.vec[i]

    cdef uint64_t _assignSn(self):
        self.sn += 1
        return self.sn

    cdef void _allocModel(self, CxtLikModel *model) except *:
        cdef unsigned i

        model.rMat = <double *>calloc(self.lik.dim * self.lik.dim, \
          sizeof(double))
        if model.rMat == NULL:
            raise MemoryError("Error allocating rMat")

        model.piDiag = <double *>calloc(self.lik.dim, sizeof(double))
        if model.piDiag == NULL:
            raise MemoryError("Error allocating piDiag")

        model.qEigVecCube = <double *>malloc(self.lik.dim * self.lik.dim * \
          self.lik.dim * sizeof(double))
        if model.qEigVecCube == NULL:
            raise MemoryError("Error allocating qEigVecCube")

        model.qEigVals = <double *>malloc(self.lik.dim * sizeof(double))
        if model.qEigVals == NULL:
            raise MemoryError("Error allocating qEigVals")

        if self.lik.ncat != 0:
            model.gammas = <double *>malloc(self.lik.ncat * sizeof(double))
            if model.gammas == NULL:
                raise MemoryError("Error allocating gammas")

        model.stripeLnL = <double *>malloc(self.lik.nstripes * sizeof(double))
        if model.stripeLnL == NULL:
            raise MemoryError("Error allocating stripeLnL")

    cdef void _initModel(self, CxtLikModel *model):
        cdef unsigned i

        model.sn = 0
        model.reassign = True
        model.weight = 0.0

        # Initialize to default equal relative rates.
        for 0 <= i < self.lik.dim:
            for i < j < self.lik.dim:
                model.rMat[i*self.lik.dim + j] = 1.0

        # Initialize to default equal frequencies.
        for 0 <= i < self.lik.dim:
            model.piDiag[i] = 1.0 / <double>self.lik.dim

        model.alpha = INFINITY

    cdef void _reassignModel(self, CxtLikModel *model) except *:
        model.sn = self._assignSn()

        # Rescale R and Pi as needed, and decompose Q.
        if CxLikQDecomp(self.lik.dim, model.rMat, model.piDiag, \
          model.qEigVecCube, model.qEigVals):
            raise ValueError("Error decomposing Q")

        model.reassign = False

    cdef void _deallocModel(self, CxtLikModel *model):
        if model.rMat != NULL:
            free(model.rMat)
            model.rMat = NULL
        if model.piDiag != NULL:
            free(model.piDiag)
            model.piDiag = NULL
        if model.qEigVecCube != NULL:
            free(model.qEigVecCube)
            model.qEigVecCube = NULL
        if model.qEigVals != NULL:
            free(model.qEigVals)
            model.qEigVals = NULL
        if model.gammas != NULL:
            free(model.gammas)
            model.gammas = NULL
        if model.stripeLnL != NULL:
            free(model.stripeLnL)
            model.stripeLnL = NULL

    cpdef unsigned getNcat(self):
        """
            Get the number of discrete rates for Gamma-distributed mutation
            rate variation.
        """
        return self.lik.ncat

    cpdef unsigned nmodels(self):
        """
            Return the current number of models being mixed.  Models can be
            referenced by their offset within the mixture model vector.
        """
        return self.lik.modelsLen

    cpdef unsigned addModel(self):
        """
            Append a model to the mixture, and return its index in the model
            vector.  The new model is created with a weight of 0, so it will
            have no impact on likelihood computations unless the weight is
            changed.
        """
        cdef CxtLikModel *models

        if self.lik.modelsLen == self.lik.modelsMax:
            models = <CxtLikModel *>realloc(self.lik.models, \
              (self.lik.modelsMax + 1) * sizeof(CxtLikModel))
            if models == NULL:
                raise MemoryError("Error reallocating models")
            self.lik.models = models
            self._allocModel(&self.lik.models[self.lik.modelsMax])
            self.lik.modelsMax += 1
        self._initModel(&self.lik.models[self.lik.modelsLen])
        self.lik.modelsLen += 1

    cpdef dupModel(self, unsigned to, unsigned fr, bint dupCLs=False):
        """
            Copy all model parameters from the model at offset 'fr' within the
            mixture model vector to the model at offset 'to' (but leave their
            weight unmodified).  If dupCLs is True, copy all valid cached
            conditional likelihoods, with the expectation that those data will
            be useful for future likelihood computations.

            Among other uses, this method makes it possible to re-order the
            mixture vector, in preparation for removing the last in the vector,
            via delModel().
        """
        cdef CxtLikModel *toP, *frP
        cdef Edge edge
        cdef Ring ring
        cdef CL cL

        assert to < self.lik.modelsLen
        assert fr < self.lik.modelsLen
        assert to != fr

        toP = &self.lik.models[to]
        frP = &self.lik.models[fr]

        toP.sn = frP.sn
        toP.reassign = frP.reassign
        memcpy(toP.rMat, frP.rMat, self.lik.dim * self.lik.dim * sizeof(double))
        memcpy(toP.piDiag, frP.piDiag, self.lik.dim * sizeof(double))
        memcpy(toP.qEigVecCube, frP.qEigVecCube, self.lik.dim * self.lik.dim * \
          self.lik.dim * sizeof(double))
        memcpy(toP.qEigVals, frP.qEigVals, self.lik.dim * self.lik.dim * \
          sizeof(double))
        if self.lik.ncat != 0:
            memcpy(toP.gammas, frP.gammas, self.lik.ncat * sizeof(double))
        toP.cL.sn = frP.cL.sn
        memcpy(toP.cL.mat, frP.cL.mat, self.lik.dim * self.lik.nchars * \
          sizeof(double))
        memcpy(toP.stripeLnL, frP.stripeLnL, self.lik.nstripes * sizeof(double))
        toP.lnL = frP.lnL

        if dupCLs:
            for edge in self.tree.getEdges():
                for ring in (edge.ring, edge.ring.other):
                    cL = <CL>ring.aux
                    if cL is None:
                        cL = CL(self.lik.nchars, self.lik.dim, \
                          self.lik.modelsMax)
                        ring.aux = cL
                    cL.dupModel(self.lik.nchars, self.lik.dim, to, fr)

    cpdef delModel(self):
        """
            Remove the last model in the mixture vector.  Weights for the other
            models are not modified, which if nothing else is changed will
            result in automatic proportional rescaling during the next lnL()
            computation.
        """
        cdef CxtLikModel *models

        if self.lik.modelsLen == 1:
            raise ValueError("At least one model must remain")
        self.lik.modelsLen -= 1

    cpdef double getWeight(self, unsigned model):
        """
            Get the weight for the specified model in the mixture vector.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]

        return modelP.weight

    cpdef setWeight(self, unsigned model, double weight):
        """
            Set the weight for the specified model in the mixture vector.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]
        assert weight >= 0.0

        modelP.weight = weight

    cpdef double getRate(self, unsigned model, unsigned i, unsigned j) \
      except -1.0:
        """
            Get the relative mutation rate at row i, column j, for the
            specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]
        assert i < self.lik.dim
        assert j < self.lik.dim
        assert i != j

        # Only the upper triangle of rMat is maintained, so swap i and j if
        # necessary.
        if i > j:
            i, j = j, i

        return modelP.rMat[i * self.lik.dim + j]

    cpdef setRate(self, unsigned model, unsigned i, unsigned j, double rate):
        """
            Set the relative mutation rate at row i, column j, for the
            specified model.  Relative rates are automatically rescaled by
            lnL() such that the bottommost rate in the upper triangle of the
            rate matrix is fixed at 1.0.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]
        assert i < self.lik.dim
        assert j < self.lik.dim
        assert i != j
        assert rate > 0.0

        # Only the upper triangle of rMat is maintained, so swap i and j if
        # necessary.
        if i > j:
            i, j = j, i

        modelP.rMat[i * self.lik.dim + j] = rate
        modelP.reassign = True

    cpdef double getFreq(self, unsigned model, unsigned i) except -1.0:
        """
            Get the frequency parameter for state i, for the specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]
        assert i < self.lik.dim

        return modelP.piDiag[i]

    cpdef setFreq(self, unsigned model, unsigned i, double freq):
        """
            Set the frequency parameter for state i, for the specified model.
            Unless compensatory changes are made to other frequency parameters,
            lnL() will automatically rescale the frequencies to sum to 1.0.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]
        assert i < self.lik.dim
        assert freq > 0.0

        modelP.piDiag[i] = freq
        modelP.reassign = True

    cpdef double getAlpha(self, unsigned model):
        """
            Get the alpha shape parameter controlling Gamma-distributed
            mutation rate variation, for the specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]

        return modelP.alpha

    cpdef setAlpha(self, unsigned model, double alpha):
        """
            Set the alpha shape parameter controlling Gamma-distributed
            mutation rate variation, for the specified model.
        """
        cdef CxtLikModel *modelP
        cdef unsigned i

        assert model < self.lik.modelsLen
        modelP = &self.lik.models[model]

        modelP.alpha = alpha
        # Discretize.
        for 0 <= i < self.lik.ncat:
            assert False # XXX Not implemented.
        modelP.reassign = True

    cdef void _planAppend(self, unsigned model, CxeLikStep variant, \
      CL parentCL, CL childCL, double edgeLen) except *:
        cdef CxtLikStep *step
        cdef CxtLikModel *modelP

        step = &self.lik.steps[self.lik.stepsLen]
        self.lik.stepsLen += 1

        modelP = &self.lik.models[model]
        # Clear the 'entire' attribute, in order to force re-aggregation of
        # conditional likelihoods for the model.
        modelP.entire = False

        step.variant = variant
        step.model = modelP
        step.parentMat = parentCL.vec[model].mat
        if childCL.vecMax == 1:
            # Be careful with leaf nodes to always use the first (and only)
            # vector element.
            step.childMat = childCL.vec[0].mat
        else:
            step.childMat = childCL.vec[model].mat
        if edgeLen < 0.0:
            raise ValueError("Negative branch length")
        step.edgeLen = edgeLen

    cdef void _planRecurse(self, Ring ring, CL parent, unsigned nSibs,
      double edgeLen) except *:
        cdef CL cL
        cdef Taxon taxon
        cdef unsigned degree, i, j
        cdef char *chars
        cdef double *mat
        cdef int ind, val
        cdef Ring r

        cL = <CL>ring.aux
        degree = ring.node.getDegree()

        # Initialize/expand cL as necessary.  Leaf nodes are initialized here
        # rather than in the Lik constructor so that sequential addition
        # methods will work correctly.
        if degree == 1:
            if cL is None:
                taxon = ring.node.taxon
                if taxon is None:
                    raise ValueError("Leaf node missing taxon")
                # Leaf nodes only need one vector element, since character data
                # can be shared by all mixture models.
                cL = CL(self.lik.nchars, self.lik.dim, 1)
                ring.aux = cL

                ind = self.alignment.taxaMap.indGet(taxon)
                if ind == -1:
                    raise ValueError( \
                      "Taxon %r missing from alignment's taxa map" % \
                      taxon.label)
                chars = self.alignment.getRow(ind)
                mat = cL.vec[0].mat
                for 0 <= i < self.lik.nchars:
                    val = self.char_.code2val(chr(chars[i]))
                    for 0 <= j < self.lik.dim:
                        if val & (1 << j):
                            mat[i*self.lik.dim + j] = 1.0
                        else:
                            mat[i*self.lik.dim + j] = 0.0
        else:
            if cL is None:
                cL = CL(self.lik.nchars, self.lik.dim, self.lik.modelsMax)
                ring.aux = cL
            else:
                cL.prepare(self.lik.nchars, self.lik.dim, self.lik.modelsMax)

        # Recurse.
        for r in ring.siblings():
            self._planRecurse(r.other, cL, degree, r.edge.length)

        # Check whether the current tree topology is compatible with the
        # parent's cache.  If not, invalidate all of the parent's caches.
        if cL.parent is not parent or cL.nSibs != nSibs or \
          cL.edgeLen != edgeLen:
            for 0 <= i < self.lik.modelsLen:
                if self.lik.models[i].weight != 0.0:
                    parent.vec[i].sn = 0
            cL.parent = parent
            cL.nSibs = nSibs
            cL.edgeLen = edgeLen

        # Check cache validity for each model.
        for 0 <= i < self.lik.modelsLen:
            if self.lik.models[i].weight != 0.0 and \
              cL.vec[i].sn != self.lik.models[i].sn:
                r = ring.next
                for 1 <= j < degree:
                    assert r != ring
                    self._planAppend(i, (CxeLikStepComputeCL if j == 1 else \
                      CxeLikStepMergeCL), cL, <CL>r.other.aux, r.edge.length)
                    r = r.next
                cL.vec[i].sn = self.lik.models[i].sn
                # Propagate invalidation to the parent.
                parent.vec[i].sn = 0

    cdef void _plan(self, Node root) except *:
        cdef Ring ring, r
        cdef unsigned degree, i

        # Make sure rootCL has space for all the mixture models.
        if self.rootCL.vecMax < self.lik.modelsMax:
            i = self.rootCL.vecMax
            self.rootCL.prepare(self.lik.nchars, self.lik.dim, \
              self.lik.modelsMax)
            # Set each model's pointer to the CxtLikCL that is actually stored
            # in rootCL, so that lnL computation can be completed in the C
            # code.  This done for the initial models in __init__().
            for i <= i < self.lik.modelsMax:
                self.lik.models[i].cL = &self.rootCL.vec[i]

        if root is None:
            root = self.tree.base
        degree = root.getDegree()
        ring = root.ring
        if degree == 1:
            # The root is a leaf node, which requires special handling.  Treat
            # the root as if it were a child separated from the root by a
            # 0-length branch.
            self._planRecurse(ring.other, self.rootCL, 2, ring.edge.length)
            self._planRecurse(ring, self.rootCL, 2, 0.0)

            for 0 <= i < self.lik.modelsLen:
                if self.lik.models[i].weight != 0.0 and \
                  self.rootCL.vec[i].sn != self.lik.models[i].sn:
                    self._planAppend(i, CxeLikStepComputeCL, self.rootCL, \
                      <CL>ring.other.aux, ring.edge.length)
                    self._planAppend(i, CxeLikStepMergeCL, self.rootCL, \
                      <CL>ring.aux, 0.0)
                    self.rootCL.vec[i].sn = self.lik.models[i].sn
        else:
            for r in ring:
                self._planRecurse(r.other, self.rootCL, degree, r.edge.length)

            for 0 <= i < self.lik.modelsLen:
                if self.lik.models[i].weight != 0.0 and \
                  self.rootCL.vec[i].sn != self.lik.models[i].sn:
                    self._planAppend(i, CxeLikStepComputeCL, self.rootCL, \
                      <CL>ring.other.aux, ring.edge.length)
                    for r in ring.siblings():
                        self._planAppend(i, CxeLikStepMergeCL, self.rootCL, \
                          <CL>r.other.aux, r.edge.length)
                    self.rootCL.vec[i].sn = self.lik.models[i].sn

    cpdef double lnL(self, Node root=None) except 1.0:
        """
            Compute the log-likelihood.  Use the tree's base node as the root
            for computation, unless a root is specified.
        """
        cdef double ret, wsum
        cdef unsigned i, stepsMax
        cdef CxtLikModel *modelP
        cdef CxtLikStep *steps

        # Rescale weights, if necessary.  If all weights are 0, weight all
        # models evenly.
        wsum = 0.0
        for 0 <= i < self.lik.modelsLen:
            wsum += self.lik.models[i].weight
        if wsum != 1.0:
            if wsum != 0.0:
                for 0 <= i < self.lik.modelsLen:
                    self.lik.models[i].weight /= wsum
            else:
                # Evenly weight models.
                for 0 <= i < self.lik.modelsLen:
                    self.lik.models[i].weight = 1.0

        # Make sure that the mixture models are up to date.
        for 0 <= i < self.lik.modelsLen:
            modelP = &self.lik.models[i]
            if modelP.weight != 0.0:
                # Optimistically set the model's 'entire' attribute.
                # _planAppend() clears the attribute if any computation at all
                # is necessary to determine the tree's lnL under the model.
                modelP.entire = True
                # Update model parameters, if necessary.
                if modelP.reassign:
                    self._reassignModel(modelP)

        # Expand steps, if necessary.
        stepsMax = ((2 * self.alignment.ntaxa) - 2) * self.lik.modelsMax
        if self.lik.stepsMax < stepsMax:
            steps = <CxtLikStep *>realloc(self.lik.steps, stepsMax * \
              sizeof(CxtLikStep))
            if self.lik.steps == NULL:
                raise MemoryError("Error reallocating steps")
            self.lik.stepsMax = stepsMax
        self.lik.stepsLen = 0

        # Generate the execution plan via post-order tree traversal.
        self._plan(root)

        # Execute the plan.
        CxLikExecute(self.lik, &ret)

        return ret
