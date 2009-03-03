import math
import os
import random
import sys
import time

import Crux.Config

from libc cimport *
from libm cimport *
from SFMT cimport *
from Crux.Tree cimport Tree, Node, Edge, Ring
from Crux.Tree.Lik cimport Lik

# Create a formatted string that represents a model's relative mutation rate
# class.
cdef str formatRclass(Lik lik, unsigned model):
    cdef list strs, rclass
    cdef unsigned rlen, i

    strs = []
    strs.append("( ")
    rclass = lik.getRclass(model)
    rlen = len(rclass)
    for 0 <= i < rlen:
        if i > 0 and rlen > 10:
            strs.append(",")
        strs.append("%d" % rclass[i])
    strs.append(" )")
    return "".join(strs)

# Create a formatted string that represents a model's relative mutation rates.
cdef str formatRates(Lik lik, unsigned model, str fmt):
    cdef list strs, rclass
    cdef unsigned nrates, rlen, i
    cdef double r

    strs = []
    strs.append("[ ")
    nrates = lik.getNrates(model)
    rclass = lik.getRclass(model)
    rlen = len(rclass)
    for 0 <= i < rlen:
        if i > 0:
            strs.append(" ")
        r = lik.getRate(model, <unsigned>rclass[i])
        strs.append(fmt % r)
    strs.append(" ]")
    return "".join(strs)

# Create a formatted string that represents a model's state frequencies.
cdef str formatFreqs(Lik lik, unsigned model, str fmt):
    cdef list strs
    cdef unsigned nstates, i
    cdef double fSum, f

    strs = []
    strs.append("[ ")

    nstates = lik.char_.nstates()
    fSum = 0.0
    for 0 <= i < nstates:
        f = lik.getFreq(model, i)
        fSum += f
    for 0 <= i < nstates:
        if i > 0:
            strs.append(" ")
        f = lik.getFreq(model, i) / fSum
        strs.append(fmt % f)

    strs.append(" ]")
    return "".join(strs)

# Used to search for the appropriate proposal, based on k.
cdef int _propsCdfCmp(void *vk, void *vx):
    cdef double k = (<double *>vk)[0]
    cdef double x = (<double *>vx)[0]
    cdef double x1

    if k >= x:
        x1 = (<double *>vx)[1]
        if k < x1:
            return 0
        else:
            return 1
    else:
        return -1

cdef int _lnLCmp(void *va, void *vb):
    cdef double a = (<double *>va)[0]
    cdef double b = (<double *>vb)[0]

    return (a > b) - (a < b)

cdef class Mc3Chain:
    def __cinit__(self):
        self.swapPrng = NULL
        self.prng = NULL

    def __dealloc__(self):
        if self.swapPrng != NULL:
            fini_gen_rand(self.swapPrng)
            self.swapPrng = NULL
        if self.prng != NULL:
            fini_gen_rand(self.prng)
            self.prng = NULL

    def __init__(self, Mc3 master, unsigned run, unsigned ind, \
      uint32_t swapSeed, uint32_t seed):
        cdef Edge edge
        cdef unsigned i, nstates, rlen, m

        self.master = master
        self.run = run
        self.ind = ind
        self.nswap = 0
        for 0 <= i < Mc3Prop:
            self.accepts[i] = 0
            self.rejects[i] = 0
        self.heat = 1.0 / (1.0 + (ind * self.master._heatDelta))
        self.swapInd = ind
        # Use a separate PRNG for Metropolis-coupled chain swaps, so that all
        # chains can independently compute the same sequence of swaps.
        self.swapPrng = init_gen_rand(swapSeed)
        if self.swapPrng == NULL:
            raise MemoryError("Error allocating swapPrng")
        self.prng = init_gen_rand(seed)
        if self.prng == NULL:
            raise MemoryError("Error allocating prng")

        # Generate a random fully resolved tree.  Polytomous starting trees are
        # never generated, but that is okay since there is no requirement to
        # draw the starting tree from the prior.  The main goal here is to
        # avoid systematic starting point dependence.
        self.tree = Tree(self.master.alignment.taxaMap.ntaxa, \
          self.master.alignment.taxaMap)
        self.tree.deroot()

        self.lik = Lik(self.tree, self.master.alignment, self.master._nmodels, \
          self.master._ncat, self.master._catMedian)

        # Randomly draw model parameters from their prior distributions.
        nstates = self.lik.char_.nstates()
        rlen = nstates * (nstates-1) / 2
        for 0 <= m < self.master._nmodels:
            if self.master._nmodels > 1:
                # Relative weight.
                if self.master.props[Mc3WeightProp] > 0.0:
                    self.lik.setWeight(m, -log(1.0 - genrand_res53(self.prng)))

                # Rate multiplier.
                if self.master.props[Mc3RmultProp] > 0.0:
                    self.lik.setRmult(m, -log(1.0 - genrand_res53(self.prng)))

            # State frequencies.
            if self.master.props[Mc3FreqProp] > 0.0:
                for 0 <= i < self.lik.char_.nstates():
                    self.lik.setFreq(m, i, -log(1.0 - genrand_res53(self.prng)))

            # Relative mutation rates.
            if self.master.props[Mc3RateProp] > 0.0:
                self.lik.setRclass(m, range(rlen))
                for 0 <= i < rlen:
                    self.lik.setRate(m, i, -log(1.0 - genrand_res53(self.prng)))

            if self.master.props[Mc3RateShapeInvProp] > 0.0:
                # Gamma-distributed rates shape parameter.
                if self.master._ncat > 1:
                    self.lik.setAlpha(m, -log(1.0 - genrand_res53(self.prng)) \
                      * self.master._rateShapeInvPrior)

        # Branch lengths.
        for edge in self.tree.getEdges():
            edge.length = -log(1.0 - genrand_res53(self.prng)) / \
              self.master._brlenPrior

        self.lnL = self.lik.lnL()

        self.step = 0
        self.master.sendSample(self.run, self.step, self.heat, self.nswap, \
          self.accepts, self.rejects, self.lik, self.lnL)

    cdef bint weightPropose(self) except *:
        cdef unsigned mInd
        cdef double u, lnM, m, w0, w1, lnL1, p

        assert self.master._nmodels > 1

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, self.lik.nmodels())

        # Generate the weight multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._weightLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified weight.
        w0 = self.lik.getWeight(mInd)
        w1 = w0 * m
        self.lik.setWeight(mInd, w1)
        lnL1 = self.lik.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # The relative weight prior paremeter is arbitrarily fixed at 1.  Due
        # to the Lik code automatically normalizing weights to sum to 1,
        # changing the prior has no effect; in all cases, the normalization
        # causes the set of weight parameters to effectively be drawn from a
        # flat Dirichlet distribution.
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + (-(w1-w0))) * self.heat + lnM)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.accepts[Mc3WeightProp] += 1
        else:
            # Reject.
            self.lik.setWeight(mInd, w0)
            self.rejects[Mc3WeightProp] += 1

        return False

    cdef bint freqPropose(self) except *:
        cdef unsigned m0Ind, m1Ind, fInd
        cdef double w, f0, f1, u, lnM, m, lnL1

        # Uniformly choose a random model within the mixture.
        m0Ind = gen_rand64_range(self.prng, self.lik.nmodels())
        # Create a duplicate scratch model.
        w = self.lik.getWeight(m0Ind)
        m1Ind = self.lik.addModel(w)
        self.lik.dupModel(m1Ind, m0Ind, False)
        self.lik.setWeight(m0Ind, 0.0)

        # Uniformly choose a random frequency parameter.
        fInd = gen_rand64_range(self.prng, self.lik.char_.nstates())

        # Generate the frequency multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._freqLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified frequency.
        f0 = self.lik.getFreq(m1Ind, fInd)
        f1 = f0 * m
        self.lik.setFreq(m1Ind, fInd, f1)
        lnL1 = self.lik.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # The frequency prior paramater is arbitrarily fixed at 1.  Due to the
        # Lik code automatically normalizing frequencies to sum to 1, changing
        # the prior has no effect; in all cases, the normalization causes the
        # set of frequency parameters to effectively be drawn from a flat
        # Dirichlet distribution.
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + (-(f1-f0))) * self.heat + lnM)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik.dupModel(m0Ind, m1Ind, True)
            self.accepts[Mc3FreqProp] += 1
        else:
            self.rejects[Mc3FreqProp] += 1
        self.lik.setWeight(m0Ind, w)
        self.lik.delModel()

        return False

    cdef bint rmultPropose(self) except *:
        cdef unsigned mInd
        cdef double u, lnM, m, rm0, rm1, lnL1, p

        assert self.master._nmodels > 1

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, self.lik.nmodels())

        # Generate the rmult multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._rmultLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified rmult.
        rm0 = self.lik.getRmult(mInd)
        rm1 = rm0 * m
        self.lik.setRmult(mInd, rm1)
        lnL1 = self.lik.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # The rate multiplier prior paremeter is arbitrarily fixed at 1.  Due
        # to the Lik code automatically normalizing relative rates to a mean
        # rate of 1, changing the prior has no effect; in all cases, the
        # normalization causes the set of rate multiplier parameters to
        # effectively be drawn from a flat Dirichlet distribution.
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + (-(rm1-rm0))) * self.heat + lnM)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.accepts[Mc3RmultProp] += 1
        else:
            # Reject.
            self.lik.setRmult(mInd, rm0)
            self.rejects[Mc3RmultProp] += 1

        return False

    cdef bint ratePropose(self) except *:
        cdef unsigned m0Ind, m1Ind, nrates, r, rInd
        cdef double w, r0, r1, u, lnM, m, lnL1
        cdef list rclass

        # Uniformly choose a random model within the mixture.
        m0Ind = gen_rand64_range(self.prng, self.lik.nmodels())
        nrates = self.lik.getNrates(m0Ind)
        if nrates == 1:
            # Only one rate class, so this proposal is bogus.
            return True
        # Create a duplicate scratch model.
        w = self.lik.getWeight(m0Ind)
        m1Ind = self.lik.addModel(w)
        self.lik.dupModel(m1Ind, m0Ind, False)
        self.lik.setWeight(m0Ind, 0.0)

        # Uniformly choose a random rate parameter.
        rclass = self.lik.getRclass(m1Ind)
        rInd = gen_rand64_range(self.prng, nrates)

        # Generate the rate multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._rateLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified rate.
        r0 = self.lik.getRate(m1Ind, rInd)
        r1 = r0 * m
        self.lik.setRate(m1Ind, rInd, r1)
        lnL1 = self.lik.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # The frequency prior paramater is arbitrarily fixed at 1.  Due to the
        # Lik code automatically normalizing frequencies to sum to 1, changing
        # the prior has no effect; in all cases, the normalization causes the
        # set of frequency parameters to effectively be drawn from a flat
        # Dirichlet distribution.
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + (-(r1-r0))) * self.heat + lnM)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik.dupModel(m0Ind, m1Ind, True)
            self.accepts[Mc3RateProp] += 1
        else:
            self.rejects[Mc3RateProp] += 1
        self.lik.setWeight(m0Ind, w)
        self.lik.delModel()

        return False

    cdef bint rateShapeInvPropose(self) except *:
        cdef unsigned m0Ind, m1Ind
        cdef double w, u, lnM, m, a0, a1, a0Inv, a1Inv, lnPrior, lnL1, p

        # Uniformly choose a random model within the mixture.
        m0Ind = gen_rand64_range(self.prng, self.lik.nmodels())
        # Create a duplicate scratch model.
        w = self.lik.getWeight(m0Ind)
        m1Ind = self.lik.addModel(w)
        self.lik.dupModel(m1Ind, m0Ind, False)
        self.lik.setWeight(m0Ind, 0.0)

        # Generate the inverse rate shape multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._rateShapeInvLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified rate shape.
        a0 = self.lik.getAlpha(m1Ind)
        a0Inv = 1.0 / a0
        a1Inv = a0Inv * m
        a1 = 1.0 / a1Inv
        self.lik.setAlpha(m1Ind, a1)
        lnL1 = self.lik.lnL()

        lnPrior = -self.master._rateShapeInvPrior * (a1Inv-a0Inv)

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnM)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik.dupModel(m0Ind, m1Ind, True)
            self.accepts[Mc3RateShapeInvProp] += 1
        else:
            self.rejects[Mc3RateShapeInvProp] += 1
        self.lik.setWeight(m0Ind, w)
        self.lik.delModel()

        return False

    cdef bint brlenPropose(self) except *:
        cdef Edge edge
        cdef double v0, v1, u, lnM, m, lnL1, lnPrior, lnProp, p

        # Uniformly choose a random edge.
        edge = <Edge>self.tree.getEdges()[gen_rand64_range(self.prng, \
          self.tree.getNedges())]

        # Generate the branch length multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._brlenLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified branch length.
        v0 = edge.length
        v1 = v0 * m
        edge.length = v1
        lnL1 = self.lik.lnL(edge.ring.node)

        lnPrior = -self.master._brlenPrior * (v1-v0)
        lnProp = lnM

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            # Re-base to make the cache fully usable if the next lnL() call
            # does not specify a rooting.
            self.lik.tree.setBase(edge.ring.node)
            self.accepts[Mc3BrlenProp] += 1
        else:
            # Reject.
            edge.length = v0
            self.rejects[Mc3BrlenProp] += 1

        return False

    cdef bint etbrPropose(self) except *:
        cdef Edge e, eA, eX, eR, eY, eS
        cdef Node n, nX0, nX1, nR0, nR1, nY0, nY1, nS0, nS1
        cdef Ring r, rnext, rFrom
        cdef list nX0Sibs, nY0Sibs
        cdef bint nX0Uncon, nR0Uncon, nY0Uncon, nS0Uncon
        cdef unsigned uI, i
        cdef double vA0, vA1, lnMA, mA
        cdef double vX0, vX1, lnMX, mX
        cdef double vY0, vY1, lnMY, mY
        cdef double u, lnL1, lnPrior, lnProp, p

        assert self.lik.tree.getNtaxa() > 3

        # Uniformly choose a random edge.
        eA = <Edge>self.tree.getEdges()[gen_rand64_range(self.prng, \
          self.tree.getNedges())]

        # For each end of eA, choose a random starting direction away from eA.
        # With probability etbrPExt, "extend" (move the end point past the next
        # node along that path.  Randomly choose a direction away from the
        # previous edge.  Continue iterative extension until extension failure
        # due to randomness (the "unconstrained" case), or a leaf edge is
        # reached (the "constrained" case).  Re-arrange the tree by extracting
        # *==eA==nX0==eX== and moving it to the extension point.  The following
        # figure depicts the transformation due to extending the left end of eA
        # to eR:
        #
        #         cH               cI                            cI            #
        #          \              /    ===>                        \           #
        #           \            /                                  \          #
        #           nR1--------nR0                                  nR0--cJ    #
        #    cG     /     eR     \                cH                /          #
        #      \   /              \    ===>         \            eR/           #
        #       \ /                cJ                \            /            #
        #   cF--nX1                                  nR1========nX0            #
        #        \\               cA          cG     /     eX     \\           #
        #         \\eX            /    ===>     \   /            eA\\          #
        #          \\     eA     /               \ /                \\         #
        #       cE--nX0=========*            cF--nX1--cC             *--cA     #
        #           / \          \               / \                /          #
        #          /   \          \    ===>     /   \              /           #
        #         cD    cC         cB          cE    cD          cB            #

        #=======================================================================
        # Extend X.
        rFrom = eA.ring
        nX0 = rFrom.node
        nR1 = nX0
        nX0Uncon = (nX0.getDegree() > 1)
        if nX0Uncon:
            # Randomly choose direction.
            uI = gen_rand64_range(self.prng, nX0.getDegree()-1)
            r = rFrom.next
            for 0 <= i < uI:
                r = r.next
            # Finish initializing so that nodes are aliased as such:
            #
            #   nR1      nR0
            #   nX0==eX==nX1
            eX = r.edge
            rFrom = r.other
            nR0 = nX1 = rFrom.node
            # Iteratively extend.
            nR0Uncon = (nR0.getDegree() > 1)
            while nR0Uncon:
                # Randomly determine whether to extend.
                u = genrand_res53(self.prng)
                if self.master._etbrPExt < u:
                    break
                # Randomly choose direction.
                uI = gen_rand64_range(self.prng, nR0.getDegree()-1)
                r = rFrom.next
                for 0 <= i < uI:
                    r = r.next
                # Extend.
                rFrom = r.other
                nR1 = nR0
                nR0 = rFrom.node
                nR0Uncon = (nR0.getDegree() > 1)
            # Perform rearrangement unless it would be a no-op.
            if nR0 is not nX1:
                # Move nX0's siblings to nX1.  Keep track of them in case the
                # proposal is rejected.
                nX0Sibs = []
                r = eA.ring.next
                while r is not eA.ring:
                    rnext = r.next
                    e = r.edge
                    if e is not eX:
                        nX0Sibs.append(e)
                        n = r.other.node
                        # Preserve e's polarity when detaching/reattaching in
                        # order to keep from disturbing CL caches.
                        if e.ring is r:
                            e.detach()
                            e.attach(nX1, n)
                        else:
                            e.detach()
                            e.attach(n, nX1)
                    r = rnext
                # Move eX.  Preserve polarity.
                if eX.ring.node is nX1:
                    eX.detach()
                    eX.attach(nR1, nX0)
                else:
                    eX.detach()
                    eX.attach(nX0, nR1)
                # Move eR.  Preserve polarity.
                eR = rFrom.edge
                if eR.ring.node is nR1:
                    eR.detach()
                    eR.attach(nX0, nR0)
                else:
                    eR.detach()
                    eR.attach(nR0, nX0)

        #=======================================================================
        # Extend Y.
        rFrom = eA.ring.other
        nY0 = rFrom.node
        nS1 = nY0
        nY0Uncon = (nY0.getDegree() > 1)
        if nY0Uncon:
            # Randomly choose direction.
            uI = gen_rand64_range(self.prng, nY0.getDegree()-1)
            r = rFrom.next
            for 0 <= i < uI:
                r = r.next
            # Finish initializing so that nodes are aliased as such:
            #
            #   nS1      nS0
            #   nY0==eY==nY1
            eY = r.edge
            rFrom = r.other
            nS0 = nY1 = rFrom.node
            # Iteratively extend.
            nS0Uncon = (nS0.getDegree() > 1)
            while nS0Uncon:
                # Randomly determine whether to extend.
                u = genrand_res53(self.prng)
                if self.master._etbrPExt < u:
                    break
                # Randomly choose direction.
                uI = gen_rand64_range(self.prng, nS0.getDegree()-1)
                r = rFrom.next
                for 0 <= i < uI:
                    r = r.next
                # Extend.
                rFrom = r.other
                nS1 = nS0
                nS0 = rFrom.node
                nS0Uncon = (nS0.getDegree() > 1)
            # Perform rearrangement unless it would be a no-op.
            if nS0 is not nY1:
                # Move nY0's siblings to nY1.  Keep track of them in case the
                # proposal is rejected.
                nY0Sibs = []
                r = eA.ring.other.next
                while r is not eA.ring.other:
                    rnext = r.next
                    e = r.edge
                    if e is not eY:
                        nY0Sibs.append(e)
                        n = r.other.node
                        # Preserve e's polarity when detaching/reattaching in
                        # order to keep from disturbing CL caches.
                        if e.ring is r:
                            e.detach()
                            e.attach(nY1, n)
                        else:
                            e.detach()
                            e.attach(n, nY1)
                    r = rnext
                # Move eY.  Preserve polarity.
                if eY.ring.node is nY1:
                    eY.detach()
                    eY.attach(nS1, nY0)
                else:
                    eY.detach()
                    eY.attach(nY0, nS1)
                # Move eS.  Preserve polarity.
                eS = rFrom.edge
                if eS.ring.node is nS1:
                    eS.detach()
                    eS.attach(nY0, nS0)
                else:
                    eS.detach()
                    eS.attach(nS0, nY0)

        #=======================================================================

        lnProp = 0.0

        # Generate branch length multipliers and set new branch lengths.
        # eA.
        u = genrand_res53(self.prng)
        lnMA = self.master._brlenLambda * (u - 0.5)
        lnProp += lnMA
        mA = exp(lnMA)
        vA0 = eA.length
        vA1 = vA0 * mA
        eA.length = vA1
        # eX.
        if nX0Uncon:
            u = genrand_res53(self.prng)
            lnMX = self.master._brlenLambda * (u - 0.5)
            lnProp += lnMX
            mX = exp(lnMX)
            vX0 = eX.length
            vX1 = vX0 * mX
            eX.length = vX1
        # eY.
        if nY0Uncon:
            u = genrand_res53(self.prng)
            lnMY = self.master._brlenLambda * (u - 0.5)
            lnProp += lnMY
            mY = exp(lnMY)
            vY0 = eY.length
            vY1 = vY0 * mY
            eY.length = vY1

        #=======================================================================

        # Compute lnL with modified branch lengths and (possibly) modified
        # topology.
        lnL1 = self.lik.lnL(nX0)

        # The prior ratio is the product of the prior ratios for each modified
        # branch length.  The number of internal branches does not change
        # (though the number of polytomies may change), so the topology prior
        # ratio is always 1.
        lnPrior = -self.master._brlenPrior * (vA1-vA0)
        if nX0Uncon:
            lnPrior += -self.master._brlenPrior * (vX1-vX0)
        if nY0Uncon:
            lnPrior += -self.master._brlenPrior * (vY1-vY0)

        # The proposal ratio is the product of the proposal ratios for
        # extension of each end of eA, as well as the branch multipliers.  The
        # ratio is 1 for the constrained/constrained and
        # unconstrained/unconstrained extension cases.
        #
        # nY0/nR0.
        if nX0Uncon and nR0 is not nX1:
            if nY0Uncon:
                if not nR0Uncon:
                    lnProp += log(1.0 - self.master._etbrPExt)
            elif nR0Uncon:
                lnProp += log(1.0 / (1.0 - self.master._etbrPExt))
        # nX0/nS0.
        if nY0Uncon and nS0 is not nY1:
            if nX0Uncon:
                if not nS0Uncon:
                    lnProp += log(1.0 - self.master._etbrPExt)
            elif nS0Uncon:
                lnProp += log(1.0 / (1.0 - self.master._etbrPExt))

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            # Re-base to make the cache fully usable if the next lnL() call
            # does not specify a rooting.
            self.lik.tree.setBase(nX0)
            self.accepts[Mc3EtbrProp] += 1
        else:
            # Reject.
            eA.length = vA0
            #===================================================================
            # Revert X extension and eX branch length change.
            if nX0Uncon:
                eX.length = vX0
                if nR0 is not nX1:
                    # Move eR.  Preserve polarity.
                    if eR.ring.node is nX0:
                        eR.detach()
                        eR.attach(nR1, nR0)
                    else:
                        eR.detach()
                        eR.attach(nR0, nR1)
                    # Move eX.  Preserve polarity.
                    if eX.ring.node is nR1:
                        eX.detach()
                        eX.attach(nX1, nX0)
                    else:
                        eX.detach()
                        eX.attach(nX0, nX1)
                    # Restore nX0's siblings (move them from nX1 to nX0).
                    for 0 <= i < len(nX0Sibs):
                        e = <Edge>nX0Sibs[i]
                        r = e.ring
                        # Preserve polarity.
                        if r.node is nX1:
                            n = r.other.node
                            e.detach()
                            e.attach(nX0, n)
                        else:
                            n = r.node
                            e.detach()
                            e.attach(n, nX0)
            #===================================================================
            # Revert Y extension and eY branch length change.
            if nY0Uncon:
                eY.length = vY0
                if nS0 is not nY1:
                    # Move eS.  Preserve polarity.
                    if eS.ring.node is nY0:
                        eS.detach()
                        eS.attach(nS1, nS0)
                    else:
                        eS.detach()
                        eS.attach(nS0, nS1)
                    # Move eY.  Preserve polarity.
                    if eY.ring.node is nS1:
                        eY.detach()
                        eY.attach(nY1, nY0)
                    else:
                        eY.detach()
                        eY.attach(nY0, nY1)
                    # Restore nY0's siblings (move them from nY1 to nY0).
                    for 0 <= i < len(nY0Sibs):
                        e = <Edge>nY0Sibs[i]
                        r = e.ring
                        # Preserve polarity.
                        if r.node is nY1:
                            n = r.other.node
                            e.detach()
                            e.attach(nY0, n)
                        else:
                            n = r.node
                            e.detach()
                            e.attach(n, nY0)
            #===================================================================
            self.rejects[Mc3EtbrProp] += 1

        return False

    cdef void rateMergePropose(self, unsigned m0Ind, unsigned m1Ind, \
      list rclass, unsigned nrates) except *:
        cdef unsigned a, b, na, nb, r, i, revSplit, n
        cdef double rateA, rateB, rateAB, lnL1, pSplit, pMerge, lnPrior
        cdef double lnProp, u, p
        cdef list rates, ns

        assert nrates > 1

        # Randomly choose two rates to merge.
        a = gen_rand64_range(self.prng, nrates)
        b = gen_rand64_range(self.prng, nrates-1)
        if b >= a:
            b += 1
        if a > b:
            a, b, = b, a

        # Compute merged relative mutation rate.
        na = 0
        nb = 0
        for 0 <= i < len(rclass):
            r = <unsigned>rclass[i]
            if r == a:
                na += 1
            elif r == b:
                nb += 1
        rateA = self.lik.getRate(m1Ind, a)
        rateB = self.lik.getRate(m1Ind, b)
        rateAB = (<double>na*rateA + <double>nb*rateB) / <double>(na+nb)

        # Update rate class list.
        for 0 <= i < len(rclass):
            r = <unsigned>rclass[i]
            if r == b:
                rclass[i] = a
            elif r > b:
                rclass[i] -= 1

        # Create new list of rates.
        rates = [self.lik.getRate(m1Ind, r) for r in xrange(nrates)]
        rates[a] = rateAB
        rates.pop(b)

        # Compute lnL with modified rclass.
        self.lik.setRclass(m1Ind, rclass, rates)
        lnL1 = self.lik.lnL()

        if nrates == 2:
            # After this merge, rateJumpPropose() can only split.  Note also
            # that rate change proposals become invalid, since there will be no
            # free rate parameters.
            pSplit = self.master.props[Mc3RateJumpProp] / \
              (1.0 - self.master.props[Mc3RateProp])
            pMerge = self.master.props[Mc3RateJumpProp] / 2.0
        elif nrates == len(rclass):
            # Splitting was not an option in rateJumpPropose() for this step.
            # (pSplit/pMerge is always 0.5 .)
            pSplit = self.master.props[Mc3RateJumpProp] / 2.0
            pMerge = self.master.props[Mc3RateJumpProp]
        else:
            # (pSplit/pMerge is always 1.0 .)
            pSplit = self.master.props[Mc3RateJumpProp] / 2.0
            pMerge = self.master.props[Mc3RateJumpProp] / 2.0
        lnPrior = log(self.master._rateJumpPrior) + log(pSplit / pMerge)

        # Count how many rate classes could be split during a reverse jump.
        revSplit = 0
        ns = [0 for i in xrange(nrates)]
        for 0 <= i < len(rclass):
            r = <unsigned>rclass[i]
            n = <unsigned>ns[r]
            if n == 1:
                revSplit += 1
            n += 1
            ns[r] = n

        lnProp = log((<double>len(rclass) * <double>nrates) / \
          (<double>revSplit * (pow(2.0, na+nb)-2.0) * rateAB*<double>(na+nb)))

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik.dupModel(m0Ind, m1Ind, True)
            self.accepts[Mc3RateJumpProp] += 1
        else:
            self.rejects[Mc3RateJumpProp] += 1

    cdef void rateSplitPropose(self, unsigned m0Ind, unsigned m1Ind, \
      list rclass, unsigned nrates) except *:
        cdef list ns, splittable, reorder, rates
        cdef unsigned i, r0, n0, r, n, rMax, ra, rb, na, nb, a0, b0
        cdef double rate0, u, rateA, rateB, pMerge, pSplit, lnPrior, lnProp
        cdef double lnL1, p
        cdef bint inA

        # Count how many rates are in each class, then create a list of
        # splittable rate classes.
        ns = [0 for i in xrange(nrates)]
        for 0 <= i < len(rclass):
            r = <unsigned>rclass[i]
            n = <unsigned>ns[r]
            n += 1
            ns[r] = n

        splittable = []
        for 0 <= i < nrates:
            assert ns[i] >= 1
            if ns[i] > 1:
                splittable.append(i)

        # Uniformly choose a random splittable rate class.
        r0 = splittable[gen_rand64_range(self.prng, len(splittable))]
        n0 = ns[r0]
        rate0 = self.lik.getRate(m1Ind, r0)

        # Randomly partition the rate class.  Rate class re-numbering is
        # performed later.
        ra = r0
        rb = nrates
        na = nb = 0
        # Randomly assign two rates to the partitions, in order to
        # guarantee that they are not empty.
        a0 = gen_rand64_range(self.prng, n0)
        b0 = gen_rand64_range(self.prng, n0-1)
        if b0 >= a0:
            b0 += 1
        for 0 <= i < len(rclass):
            r = <unsigned>rclass[i]
            if r == r0:
                if na + nb == a0:
                    # This rate was pre-assigned to a.
                    inA = True
                elif na + nb == b0:
                    inA = False
                else:
                    inA = (gen_rand64_range(self.prng, 2) == 0)

                if inA:
                    rclass[i] = ra
                    na += 1
                else:
                    rclass[i] = rb
                    nb += 1
        assert na != 0
        assert nb != 0
        assert na + nb == n0

        # Determine rates for new rate classes.
        u = (<double>n0 * genrand_res53(self.prng) - <double>na) * rate0
        rateA = rate0 + u / <double>na
        rateB = rate0 - u / <double>nb

        rates0 = [self.lik.getRate(m1Ind, i) for i in xrange(nrates)]
        rates0.append(-1.0) # Place holder for new rate.

        rates0[ra] = rateA
        rates0[rb] = rateB

        # Create rate re-ordering lookup.
        reorder = [None for i in xrange(nrates+1)]
        rMax = 0
        for 0 <= i < len(rclass):
            r = rclass[i]
            if reorder[r] is None:
                reorder[r] = rMax
                rMax += 1
        # Re-order.
        for 0 <= i < len(rclass):
            r = <unsigned>rclass[i]
            rclass[i] = reorder[rclass[i]]
        rates1 = [None for i in xrange(nrates+1)]
        for 0 <= i < nrates+1:
            rates1[reorder[i]] = rates0[i]

        # Compute lnL with modified rclass.
        self.lik.setRclass(m1Ind, rclass, rates1)
        lnL1 = self.lik.lnL()

        if nrates == len(rclass) - 1:
            # After this split, rateJumpPropose() can only merge.
            # (pMerge/pSplit is always 2.0 .)
            pMerge = self.master.props[Mc3RateJumpProp]
            pSplit = self.master.props[Mc3RateJumpProp] / 2.0
        elif nrates == 1:
            # Merging was not an option in rateJumpPropose() for this step.
            # Note also that rate change proposals were invalid, since there
            # were no free rate parameters.
            pMerge = self.master.props[Mc3RateJumpProp] / 2.0
            pSplit = self.master.props[Mc3RateJumpProp] / \
              (1.0 - self.master.props[Mc3RateProp])
        else:
            # (pMerge/pSplit is always 1.0 .)
            pMerge = self.master.props[Mc3RateJumpProp] / 2.0
            pSplit = self.master.props[Mc3RateJumpProp] / 2.0
        lnPrior = -log(self.master._rateJumpPrior) + log(pMerge / pSplit)

        lnProp = log((<double>len(splittable) * (pow(2.0, n0)-2.0) * \
          rate0*<double>n0) / (<double>len(rclass) * <double>(nrates+1)))

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik.dupModel(m0Ind, m1Ind, True)
            self.accepts[Mc3RateJumpProp] += 1
        else:
            self.rejects[Mc3RateJumpProp] += 1

    cdef bint rateJumpPropose(self) except *:
        cdef unsigned m0Ind, m1Ind, nrates, r
        cdef double w
        cdef list rclass
        cdef bint merge

        assert self.lik.char_.nstates() > 2

        # Uniformly choose a random model within the mixture.
        m0Ind = gen_rand64_range(self.prng, self.lik.nmodels())
        # Create a duplicate scratch model.
        w = self.lik.getWeight(m0Ind)
        m1Ind = self.lik.addModel(w)
        self.lik.dupModel(m1Ind, m0Ind, False)
        self.lik.setWeight(m0Ind, 0.0)

        # Determine the number of rate parameters.
        rclass = self.lik.getRclass(m1Ind)
        nrates = self.lik.getNrates(m1Ind)
        if nrates == 1:
            # Can't merge.
            merge = False
        elif nrates == len(rclass):
            # Can't split.
            merge = True
        else:
            # Randomly decide whether to try a merge or a split.
            merge = <bint>gen_rand64_range(self.prng, 2)

        if merge:
            self.rateMergePropose(m0Ind, m1Ind, rclass, nrates)
        else:
            self.rateSplitPropose(m0Ind, m1Ind, rclass, nrates)
        self.lik.setWeight(m0Ind, w)
        self.lik.delModel()

        return False

    cdef void polytomyMergePropose(self, Tree tree, unsigned nedges, \
      unsigned ntaxa) except *:
        cdef unsigned uI, i, j, nodeADeg, nodeBDeg, np1
        cdef list edges, sibsB, nodes
        cdef Edge edge
        cdef Node n, nodeA, nodeB, base
        cdef Ring r, rnext
        cdef double lnL1, lnPrior, lnGam, lnJacob, lnProp, u, p

        # Uniformly choose a random internal edge to remove.
        edges = tree.getEdges()
        uI = gen_rand64_range(self.prng, nedges-ntaxa)
        j = 0
        assert nedges == len(edges)
        for 0 <= i < nedges:
            edge = <Edge>edges[i]
            nodeA = edge.ring.node
            nodeB = edge.ring.other.node
            nodeADeg = nodeA.getDegree()
            nodeBDeg = nodeB.getDegree()
            if nodeADeg > 1 and nodeBDeg > 1:
                if j == uI:
                    break
                j += 1
        assert nodeADeg > 1 and nodeBDeg > 1

        # Detach edge that is to be collapsed.
        edge.detach()

        # Move nodeB's siblings to nodeA.
        sibsB = []
        r = nodeB.ring
        assert nodeBDeg-1 == nodeB._degreeGet(True)
        for 0 <= i < nodeBDeg-1:
            rnext = r.next
            e = r.edge
            # Record the edges that connect the siblings of nodeB so that if
            # the proposal is rejected, it is possible to reverse the tree
            # modification.
            sibsB.append(e)
            # Preserve e's polarity when detaching/reattaching in order to keep
            # from disturbing CL caches.
            n = r.other.node
            e.detach()
            if e.ring is r:
                e.attach(nodeA, n)
            else:
                e.attach(n, nodeA)
            r = rnext

        base = tree.getBase()
        if base is nodeB:
            # The excised nodeB happened to be the tree base, so fix that
            # before doing any tree operations that depend on tree traversal.
            tree.setBase(nodeA)

        # Compute lnL with new polytomy.
        lnL1 = self.lik.lnL(nodeA)

        lnPrior = log(self.master._polytomyJumpPrior) \
          - log(self.master._brlenPrior) - \
          (-self.master._brlenPrior*edge.length)

        if nedges == (2*ntaxa)-3 and nedges-1 != ntaxa:
            # The current tree is fully resolved and the proposed tree is not
            # the star tree.
            lnGam = log(0.5)
        elif nedges != (2*ntaxa)-3 and nedges-1 == ntaxa:
            # The current tree is not fully resolved and the proposed tree is
            # the star tree.
            lnGam = log(2.0)
        else:
            lnGam = 0.0 # log(1.0)

        nodes = tree.getNodes()
        np1 = 0
        for 0 <= i < len(nodes):
            n = <Node>nodes[i]
            if n.getDegree() > 3:
                np1 += 1

        lnJacob = log(self.master._brlenPrior) + \
          (-self.master._brlenPrior*edge.length)
        lnProp = lnGam + log(<double>(nedges-ntaxa)) - log(<double>np1) \
          - log(pow(2.0, <double>(nodeADeg+nodeBDeg-2)-1.0) - \
          <double>(nodeADeg+nodeBDeg-2) - 1.0) + lnJacob

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.accepts[Mc3PolytomyJumpProp] += 1
        else:
            # Reject.
            for 0 <= i < len(sibsB):
                e = <Edge>sibsB[i]
                # Preserve e's polarity when detaching/reattaching in order to
                # keep from disturbing CL caches.
                n = e.ring.node
                if n is nodeA:
                    n = e.ring.other.node
                    e.detach()
                    e.attach(nodeB, n)
                else:
                    e.detach()
                    e.attach(n, nodeB)
            edge.attach(nodeA, nodeB)
            if base is not tree.getBase():
                tree.setBase(base)
            self.rejects[Mc3PolytomyJumpProp] += 1

    cdef void polytomySplitPropose(self, Tree tree, unsigned nedges, \
      unsigned ntaxa) except *:
        cdef list nodes, polys
        cdef Node n, nodeA, nodeB
        cdef Ring r, rnext
        cdef Edge e
        cdef unsigned i, deg0, a0, aa0, b0, bb0, degA, degB
        cdef bint inA
        cdef double lnL1, lnPrior, lnGam, lnJacob, lnProp, u, p

        # Uniformly choose a random polytomy to split.
        nodes = tree.getNodes()
        polys = []
        for 0 <= i < len(nodes):
            n = <Node>nodes[i]
            if n.getDegree() > 3:
                polys.append(n)
        assert len(polys) > 0
        nodeA = <Node>polys[gen_rand64_range(self.prng, len(polys))]
        nodeB = Node(tree)

        # Create the edge that will induce the new bipartition, and draw its
        # length from the prior distribution.
        edge = Edge(tree)
        edge.length = -log(1.0 - genrand_res53(self.prng)) / \
          self.master._brlenPrior

        # Randomly assign four siblings to the partitions, in order to guarantee
        # that each partition contains at least two edges.
        deg0 = nodeA.getDegree()
        (a0, aa0, b0, bb0) = random.sample(xrange(deg0), 4)

        # Partition.
        degA = degB = 0
        r = nodeA.ring
        for 0 <= i < deg0:
            rnext = r.next
            if i == a0 or i == aa0:
                # This sibling was pre-assigned to nodeA.
                inA = True
            elif i == b0 or i == bb0:
                # This sibling was pre-assigned to nodeB.
                inA = False
            else:
                inA = (gen_rand64_range(self.prng, 2) == 0)

            if inA:
                # Leave edge/node associated with r attached to nodeA.
                degA += 1
            else:
                # Move edge/node associated with r to nodeB.
                e = r.edge
                n = r.other.node
                # Preserve e's polarity when detaching/reattaching in order to
                # keep from disturbing CL caches.
                e.detach()
                if e.ring is r:
                    e.attach(nodeB, n)
                else:
                    e.attach(n, nodeB)
                degB += 1
            r = rnext
        assert degA != 0
        assert degB != 0
        assert degA + degB == deg0
        edge.attach(nodeA, nodeB)

        # Compute lnL with new polytomy.
        lnL1 = self.lik.lnL(nodeA)

        lnPrior = -log(self.master._polytomyJumpPrior) \
          + log(self.master._brlenPrior) + \
          (-self.master._brlenPrior*edge.length)

        if nedges == ntaxa and nedges+1 != (2*ntaxa)-3:
            # The current tree is the star tree, and the proposed tree is not
            # fully resolved.
            lnGam = log(0.5)
        elif nedges != ntaxa and nedges+1 == (2*ntaxa)-3:
            # The current tree is not the star tree, and the proposed tree is
            # fully resolved.
            lnGam = log(2.0)
        else:
            lnGam = 0.0 # log(1.0)

        lnJacob = -log(self.master._brlenPrior) \
          - (-self.master._brlenPrior*edge.length)
        lnProp = lnGam + log(<double>len(polys)) \
          + log(pow(2.0, <double>(deg0-1)) - <double>deg0 - 1.0) \
          - log(<double>(nedges-ntaxa+1))

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.accepts[Mc3PolytomyJumpProp] += 1
        else:
            # Reject.
            edge.detach()
            r = nodeB.ring
            for 0 <= i < degB:
                rnext = r.next
                e = r.edge
                n = r.other.node
                # Preserve e's polarity when detaching/reattaching in order to
                # keep from disturbing CL caches.
                e.detach()
                if e.ring is r:
                    e.attach(nodeA, n)
                else:
                    e.attach(n, nodeA)
                r = rnext
            self.rejects[Mc3PolytomyJumpProp] += 1

    cdef bint polytomyJumpPropose(self) except *:
        cdef Tree tree
        cdef unsigned nedges, ntaxa
        cdef bint merge

        tree = self.lik.tree
        assert tree.getNtaxa() > 3
        nedges = tree.getNedges()
        ntaxa = tree.getNtaxa()
        if nedges == ntaxa:
            # Can't merge.
            merge = False
        elif nedges == (2*ntaxa)-3:
            # Can't split.
            merge = True
        else:
            # Randomly decide whether to try a merge or a split.
            merge = <bint>gen_rand64_range(self.prng, 2)

        if merge:
            self.polytomyMergePropose(tree, nedges, ntaxa)
        else:
            self.polytomySplitPropose(tree, nedges, ntaxa)

        return False

    cdef void rateShapeInvRemovePropose(self, unsigned m0Ind, unsigned m1Ind, \
      double alpha0) except *:
        cdef double rateShapeInv0
        cdef double lnL1, lnPrior, lnJacob, lnProp, u, p

        rateShapeInv0 = 1.0 / alpha0

        # Disable +G.
        self.lik.setAlpha(m1Ind, INFINITY)

        # Compute lnL without +G.
        lnL1 = self.lik.lnL()

        lnPrior = log(self.master._rateShapeInvJumpPrior) \
          - log(self.master._rateShapeInvPrior) - \
          (-self.master._rateShapeInvPrior*rateShapeInv0)

        lnJacob = log(self.master._rateShapeInvPrior) + \
          (-self.master._rateShapeInvPrior*rateShapeInv0)
        lnProp = lnJacob

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik.dupModel(m0Ind, m1Ind, True)
            self.accepts[Mc3RateShapeInvJumpProp] += 1
        else:
            self.rejects[Mc3RateShapeInvJumpProp] += 1

    cdef void rateShapeInvAddPropose(self, unsigned m0Ind, unsigned m1Ind) \
      except *:
        cdef double rateShapeInv1, alpha1
        cdef double lnL1, lnPrior, lnJacob, lnProp, u, p

        # Draw the inverse shape parameter from the prior distribution.
        rateShapeInv1 = -log(1.0 - genrand_res53(self.prng)) / \
          self.master._rateShapeInvPrior
        alpha1 = 1.0 / rateShapeInv1

        # Enable +G.
        self.lik.setAlpha(m1Ind, alpha1)

        # Compute lnL with +G.
        lnL1 = self.lik.lnL()

        lnPrior = -log(self.master._rateShapeInvJumpPrior) \
          + log(self.master._rateShapeInvPrior) + \
          (-self.master._rateShapeInvPrior*rateShapeInv1)

        lnJacob = -log(self.master._rateShapeInvPrior) \
          - (-self.master._rateShapeInvPrior*rateShapeInv1)
        lnProp = lnJacob

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik.dupModel(m0Ind, m1Ind, True)
            self.accepts[Mc3RateShapeInvJumpProp] += 1
        else:
            self.rejects[Mc3RateShapeInvJumpProp] += 1

    cdef bint rateShapeInvJumpPropose(self) except *:
        cdef unsigned m0Ind, m1Ind
        cdef double w, alpha

        assert self.master._ncat > 1

        # Uniformly choose a random model within the mixture.
        m0Ind = gen_rand64_range(self.prng, self.lik.nmodels())
        # Create a duplicate scratch model.
        w = self.lik.getWeight(m0Ind)
        m1Ind = self.lik.addModel(w)
        self.lik.dupModel(m1Ind, m0Ind, False)
        self.lik.setWeight(m0Ind, 0.0)

        # Determine whether the model is currently +G.
        alpha = self.lik.getAlpha(m1Ind)

        if alpha != INFINITY:
            self.rateShapeInvRemovePropose(m0Ind, m1Ind, alpha)
        else:
            self.rateShapeInvAddPropose(m0Ind, m1Ind)
        self.lik.setWeight(m0Ind, w)
        self.lik.delModel()

        return False

    cdef void advance(self) except *:
        cdef bint again
        cdef unsigned propInd, a, b, other
        cdef double rcvHeat, rcvLnL, u, p

        self.step += 1

        # Finish handling a pending potential heat swap.
        if self.swapInd != self.ind:
            self.master.recvSwapInfo(self.run, self.ind, self.swapInd, \
              self.step-1, &rcvHeat, &rcvLnL)
            assert rcvHeat != self.heat
            p = exp((rcvLnL - self.lnL) * self.heat \
              + (self.lnL - rcvLnL) * rcvHeat)
            if p >= self.swapProb:
                self.heat = rcvHeat
                self.nswap += 1
            self.swapInd = self.ind

        again = True
        while again:
            # Choose proposal type.
            u = genrand_res53(self.prng)
            propInd = (<uintptr_t>bsearch(&u, self.master.propsCdf, \
              (sizeof(self.master.propsCdf) / sizeof(double)) - 1, \
              sizeof(double), _propsCdfCmp) - <uintptr_t>self.master.propsCdf) \
              / sizeof(double)

            # Propose new state.
            if propInd == Mc3WeightProp:
                again = self.weightPropose()
            elif propInd == Mc3FreqProp:
                again = self.freqPropose()
            elif propInd == Mc3RmultProp:
                again = self.rmultPropose()
            elif propInd == Mc3RateProp:
                again = self.ratePropose()
            elif propInd == Mc3RateShapeInvProp:
                again = self.rateShapeInvPropose()
            elif propInd == Mc3BrlenProp:
                again = self.brlenPropose()
            elif propInd == Mc3EtbrProp:
                again = self.etbrPropose()
            elif propInd == Mc3RateJumpProp:
                again = self.rateJumpPropose()
            elif propInd == Mc3PolytomyJumpProp:
                again = self.polytomyJumpPropose()
            elif propInd == Mc3RateShapeInvJumpProp:
                again = self.rateShapeInvJumpPropose()
            else:
                assert False

        # Sample if this step is a multiple of the sample stride.
        if self.step % self.master._stride == 0:
            self.master.sendSample(self.run, self.step, self.heat, self.nswap, \
              self.accepts, self.rejects, self.lik, self.lnL)
            # Clear swap stats.
            self.nswap = 0
            # Clear accept/reject stats.
            for 0 <= i < Mc3Prop:
                self.accepts[i] = 0
                self.rejects[i] = 0

        # Potentially try a heat swap with another chain.
        if self.master._ncoupled > 1 and \
          self.step % self.master._swapStride == 0:
            # Randomly choose two chains that should attempt a swap.  All
            # chains use the same PRNG so that the swap sequence can be
            # independently determined on the fly.
            a = gen_rand64_range(self.swapPrng, self.master._ncoupled)
            b = gen_rand64_range(self.swapPrng, self.master._ncoupled - 1)
            if b >= a:
                b += 1
            if a == self.ind:
                self.swapInd = b
            elif b == self.ind:
                self.swapInd = a

            if self.swapInd != self.ind:
                self.master.sendSwapInfo(self.run, self.ind, self.swapInd, \
                  self.step, self.heat, self.lnL)

            # All chains must advance the PRNG identically, so swapProb must be
            # drawn here regardless of whether it is used by this chain later.
            self.swapProb = genrand_res53(self.swapPrng)

cdef class Mc3:
    """
        Mc3 uses reversible jump Metropolis-coupled Markov chain Monte Carlo to
        create a representative sample (assuming convergence) of the stationary
        posterior distribution, based on the GTR+Gamma family of sequence
        evolution models.

        Polytomous tree topologies are sampled using the methods described by
        Lewis et al. (2005), and restrictions of the fully parameterized GTR
        model are sampled using the methods described by Huelsenbeck et al.
        (2004).

        Extending tree bisection and reconnection (eTBR) is used to change
        topology and to modify branch lengths (Lakner et al. 2008).  This
        implementation is generalized to support polytomies, as well as to
        apply to all edges (including leaf edges).  The latter generalization
        allows topology tranformations even on trees that have a single
        internal edge.

        If MPI (http://www.mpi-forum.org/) support is enabled, the methods
        described by Altekar et al. (2004) are used to run chains in parallel.
        Ideally, the total number of chains (number of independent runs times
        number of Metropolis-coupled chains per run) should be an even multiple
        of the number of MPI nodes.

        Convergence (based on log-likelihoods) is monitored using the
        interval-based coverage ratio diagnostic described at the bottom of
        page 441 in Brooks and Gelman (1998).  The diagnostic takes samples
        from the latter halves of multiple independent MCMC chains, then for
        each chain erects a credibility interval and computes the proportion of
        the concatenation of samples from all chains that are covered by the
        credibility interval.  The chains are considered to have converged when
        the mean coverage is within some epsilon of the credibility interval.
        This convergence condition indicates that the chains are sampling from
        approximately the same distribution.

        Constructor parameters:

          alignment
            Aligned character-by-taxon matrix input data.

          outPrefix
            Filename prefix, used for naming output files.

        References:

          Altekar, G., S. Dwarkadas, J.P. Huelsenbeck, F. Ronquist (2004)
          Parallel Metropolis coupled Markov chain Monte Carlo for Bayesian
          phylogentic inference.  Bioinformatics 20(3):407-415.

          Brooks, S.P., A. Gelman (1998) General Methods for Monitoring
          Convergence of Iterative Simulations.  J. Comput. Graph. Stat.
          7(4):434-455.

          Huelsenbeck, J.P., B. Larget, M.E. Alfaro (2004) Bayesian
          phylogenetic model selection using reversible jump Markov chain Monte
          Carlo.  Mol. Biol. Evol. 21(6):1123-1133.

          Lakner, C., P.v.d Mark, J.P. Huelsenbeck, B. Larget, F. Ronquist
          (2008) Efficiency of Markov chain Monte Carlo tree proposals in
          Bayesian phylogenetics.  Syst. Biol. 57(1):86-103.

          Lewis, P.O., M.T. Holder, K.E. Holsinger (2005) Polytomies and
          Bayesian phylogenetic inference.  Sys. Biol. 54(2):241-253.
    """
    def __cinit__(self):
        cdef unsigned i

        self.swapInfo = NULL
        self.swapStats = NULL
        for 0 <= i < Mc3Prop:
            self.propStats[i] = NULL
        self.lnLs = NULL

    def __dealloc__(self):
        cdef unsigned i

        if self.swapInfo != NULL:
            free(self.swapInfo)
            self.swapInfo = NULL
        if self.swapStats != NULL:
            free(self.swapStats)
            self.swapStats = NULL
        for 0 <= i < Mc3Prop:
            if self.propStats[i] != NULL:
                free(self.propStats[i])
                self.propStats[i] = NULL
        if self.lnLs != NULL:
            for 0 <= i < self._nruns+1:
                if self.lnLs[i] != NULL:
                    free(self.lnLs[i])
                    self.lnLs[i] = NULL
            free(self.lnLs)
            self.lnLs = NULL

    def __init__(self, Alignment alignment, str outPrefix):
        self.alignment = alignment
        self.outPrefix = outPrefix

        # Set defaults.
        self._graphDelay = -1.0
        self._emaAlpha = 1.0
        self._cvgSampStride = 1
        self._cvgAlpha = 0.05
        self._cvgEpsilon = 0.025
        self._minStep = 100000
        self._maxStep = ULLONG_MAX
        self._stride = 100
        self._nruns = 2
        self._ncoupled = 1
        self._heatDelta = 0.05
        self._swapStride = 1
        self._nmodels = 1
        self._ncat = 1
        self._catMedian = False
        self._weightLambda = 2.0 * log(1.6)
        self._freqLambda = 2.0 * log(1.6)
        self._rmultLambda = 2.0 * log(1.6)
        self._rateLambda = 2.0 * log(1.6)
        self._rateShapeInvLambda = 2.0 * log(1.6)
        self._brlenLambda = 2.0 * log(1.6)
        self._etbrPExt = 0.8
        self._etbrLambda = 2.0 * log(1.6)
        self._rateShapeInvPrior = 1.0
        self._brlenPrior = 10.0
        self._rateJumpPrior = 1.0
        self._polytomyJumpPrior = 1.0
        self._rateShapeInvJumpPrior = 1.0
        self.props[Mc3WeightProp] = 3.0
        self.props[Mc3FreqProp] = 1.0
        self.props[Mc3RmultProp] = 1.0
        self.props[Mc3RateProp] = 3.0
        self.props[Mc3RateShapeInvProp] = 1.0
        self.props[Mc3BrlenProp] = 2.0
        self.props[Mc3EtbrProp] = 8.0
        self.props[Mc3RateJumpProp] = 1.0
        self.props[Mc3PolytomyJumpProp] = 4.0
        self.props[Mc3RateShapeInvJumpProp] = 1.0

    cdef void sendSample(self, unsigned runInd, uint64_t step, double heat, \
      uint64_t nswap, uint64_t *accepts, uint64_t *rejects, Lik lik, \
      double lnL) except *:
        cdef uint64_t sample, lnLsMax
        cdef double *lnLs
        cdef list run
        cdef Mc3Chain chain
        cdef double swaprate, wsum, w
        cdef unsigned i, nmodels, m

        sample = step / self._stride

        # Record swap statistics.
        self.swapStats[runInd].n += nswap

        if heat != 1.0:
            # The remainder of this method applies only to cold chains.
            return

        # Record accept/reject rate statistics.
        for 0 <= i < Mc3Prop:
            self.propStats[i][runInd].n = accepts[i]
            self.propStats[i][runInd].d = accepts[i] + rejects[i]

        # Allocate/expand lnLs as necessary.
        if self.lnLs == NULL:
            self.lnLs = <double **>calloc(self._nruns+1, sizeof(double *))
            if self.lnLs == NULL:
                raise MemoryError("Error allocating lnLs")

            if self._minStep == 0:
                self.lnLsMax = 8 # Arbitrary starting size.
            else:
                self.lnLsMax = self._minStep / self._stride
            for 0 <= i < self._nruns+1:
                self.lnLs[i] = <double *>malloc(self.lnLsMax * sizeof(double))
                if self.lnLs[i] == NULL:
                    for 0 <= j < i:
                        free(self.lnLs[j])
                    free(self.lnLs)
                    self.lnLs = NULL
                    raise MemoryError("Error allocating lnLs[%d]" % i)
        elif self.lnLsMax <= sample:
            # Grow matrix using a slowly growing exponential series.
            lnLsMax = self.lnLsMax + 4 + (self.lnLsMax >> 3) + \
              (self.lnLsMax >> 4)
            for 0 <= i < self._nruns+1:
                lnLs = <double *>realloc(self.lnLs[i], lnLsMax * sizeof(double))
                if lnLs == NULL:
                    raise MemoryError("Error reallocating lnLs[%d]" % i)
                self.lnLs[i] = lnLs
            self.lnLsMax = lnLsMax

        # Store lnL.
        self.lnLs[runInd][sample] = lnL

        # Write to log files.
        # .t
        self.tFile.write("[%d %d] " % (runInd, step))
        lik.tree.render(True, "%.12e", None, self.tFile)
        self.tFile.flush()
        # .p
        wsum = 0.0
        nmodels = lik.nmodels()
        for 0 <= m < nmodels:
            wsum += lik.getWeight(m)
        for 0 <= m < nmodels:
            w = (lik.getWeight(m)/wsum if wsum > 0.0 else 1.0/nmodels)
            self.pFile.write("%d\t%d\t%d\t%.5f\t%.5f\t%s%s %.5e %.5e %s\n" % \
              (runInd, step, m, w, lik.getRmult(m), formatRclass(lik, m), \
              formatRates(lik, m, "%.5e"), lik.getWNorm(m), lik.getAlpha(m), \
              formatFreqs(lik, m, "%.5e")))
            if self.verbose:
                sys.stdout.write( \
                  "p\t%d\t%d\t%d\t%.5f\t%.5f\t%s%s %.5f %.5e %s\n" % \
                  (runInd, step, m, w, lik.getRmult(m), formatRclass(lik, m), \
                  formatRates(lik, m, "%.4e"), lik.getWNorm(m), \
                  lik.getAlpha(m), formatFreqs(lik, m, "%.5f")))
        self.pFile.flush()

    cdef void sendSwapInfo(self, unsigned runInd, unsigned srcChainInd, \
      unsigned dstChainInd, uint64_t step, double heat, double lnL) except *:
        cdef Mc3SwapInfo *swapInfo

        assert step % self._swapStride == 0
        swapInfo = &self.swapInfo[((runInd*self._ncoupled + \
          srcChainInd)*self._ncoupled + dstChainInd)*2 + \
          ((step/self._swapStride) % 2)]
        assert swapInfo.step == 0

        swapInfo.step = step
        swapInfo.heat = heat
        swapInfo.lnL = lnL

    cdef void recvSwapInfo(self, unsigned runInd, unsigned dstChainInd, \
      unsigned srcChainInd, uint64_t step, double *heat, double *lnL) except *:
        cdef Mc3SwapInfo *swapInfo

        assert step % self._swapStride == 0
        swapInfo = &self.swapInfo[((runInd*self._ncoupled + \
          srcChainInd)*self._ncoupled + dstChainInd)*2 + \
          ((step/self._swapStride) % 2)]
        assert swapInfo.step == step

        heat[0] = swapInfo.heat
        lnL[0] = swapInfo.lnL

        swapInfo.step = 0

    cdef void initLogs(self) except *:
        cdef file f

        self.lFile = open("%s.l" % self.outPrefix, "w")
        self.lFile.write("Begin run: %s\n" % \
          time.strftime("%Y/%m/%d %H:%M:%S (%Z)", time.localtime(time.time())))
        self.lFile.write("Host machine: %r\n" % (os.uname(),))
        for f in ((self.lFile, sys.stdout) if self.verbose else (self.lFile,)):
            f.write("PRNG seed: %d\n" % Crux.Config.seed)
            f.write("Compact alignment:\n")
            self.alignment.render(50, f)

            f.write("Configuration parameters:\n")
            f.write("  outPrefix: %r\n" % self.outPrefix)
            f.write("  graphDelay: %r\n" % self._graphDelay)
            f.write("  emaAlpha: %r\n" % self._emaAlpha)
            f.write("  cvgSampStride: %r\n" % self._cvgSampStride)
            f.write("  cvgAlpha: %r\n" % self._cvgAlpha)
            f.write("  cvgEpsilon: %r\n" % self._cvgEpsilon)
            f.write("  minStep: %r\n" % self._minStep)
            f.write("  maxStep: %r\n" % self._maxStep)
            f.write("  stride: %r\n" % self._stride)
            f.write("  nruns: %r\n" % self._nruns)
            f.write("  ncoupled: %r\n" % self._ncoupled)
            f.write("  heatDelta: %r\n" % self._heatDelta)
            f.write("  swapStride: %r\n" % self._swapStride)
            f.write("  nmodels: %r\n" % self._nmodels)
            f.write("  ncat: %r\n" % self._ncat)
            f.write("  catMedian: %r\n" % self._catMedian)
            f.write("  weightLambda: %r\n" % self._weightLambda)
            f.write("  freqLambda: %r\n" % self._freqLambda)
            f.write("  rmultLambda: %r\n" % self._rmultLambda)
            f.write("  rateLambda: %r\n" % self._rateLambda)
            f.write("  rateShapeInvLambda: %r\n" % self._rateShapeInvLambda)
            f.write("  brlenLambda: %r\n" % self._brlenLambda)
            f.write("  etbrPExt: %r\n" % self._etbrPExt)
            f.write("  etbrLambda: %r\n" % self._etbrLambda)
            f.write("  rateShapeInvPrior: %r\n" % self._rateShapeInvPrior)
            f.write("  brlenPrior: %r\n" % self._brlenPrior)
            f.write("  rateJumpPrior: %r\n" % self._rateJumpPrior)
            f.write("  polytomyJumpPrior: %r\n" % self._polytomyJumpPrior)
            f.write("  rateShapeInvJumpPrior: %r\n" % \
              self._rateShapeInvJumpPrior)
            f.write("  weightProp: %r\n" % self.props[Mc3WeightProp])
            f.write("  freqProp: %r\n" % self.props[Mc3FreqProp])
            f.write("  rmultProp: %r\n" % self.props[Mc3RmultProp])
            f.write("  rateProp: %r\n" % self.props[Mc3RateProp])
            f.write("  rateShapeInvProp: %r\n" % \
              self.props[Mc3RateShapeInvProp])
            f.write("  brlenProp: %r\n" % self.props[Mc3BrlenProp])
            f.write("  etbrProp: %r\n" % self.props[Mc3EtbrProp])
            f.write("  rateJumpProp: %r\n" % self.props[Mc3RateJumpProp])
            f.write("  polytomyJumpProp: %r\n" % \
              self.props[Mc3PolytomyJumpProp])
            f.write("  rateShapeInvJumpProp: %r\n" % \
              self.props[Mc3RateShapeInvJumpProp])
        self.lFile.flush()

        self.tFile = open("%s.t" % self.outPrefix, "w")
        self.tFile.write("[[run step] tree]\n")
        self.tFile.flush()

        self.pFile = open("%s.p" % self.outPrefix, "w")
        self.pFile.write( \
          "run\tstep\tmodel\tweight\trmult\t( rclass )[ R ] wNorm alpha" \
          " [ Pi ]\n")
        self.pFile.flush()
        if self.verbose:
            sys.stdout.write( \
              "p\trun\tstep\tmodel\tweight\trmult\t( rclass )[ R ] wNorm" \
              " alpha [ Pi ]\n")

        self.sFile = open("%s.s" % self.outPrefix, "w")
        self.sFile.write("step\t[ lnLs ] Rcov [ swapRates ] {" \
          " [ weightPropRates ] [ freqPropRates ] [ rmultPropRates ]" \
          " [ ratePropRates ] [ rateShapeInvPropRates ] [ brlenPropRates ]" \
          " [ etbrPropRates ] [ rateJumpPropRates ] [ polytomyJumpPropRates ]" \
          " [ rateShapeInvJumpPropRates ] }\n")
        self.sFile.flush()
        if self.verbose:
            sys.stdout.write("s\tstep\t[ lnLs ] Rcov\n")

    # Compute Rcov convergence diagnostic.
    cdef double computeRcov(self, uint64_t last) except *:
        cdef ret
        cdef uint64_t first, past, lower, upper, i, j, n
        cdef unsigned runInd
        cdef double alphaEmp
        cdef double *rcovScratch = self.lnLs[self._nruns]

        # Utilize (at most) the last half of each chain.
        past = last + 1
        first = (last/2)+1
        if first == past:
            return 0.0

        # Compute the lower and upper bounds for the credibility intervals.  In
        # cases of rounding, expand rather than contract the credibility
        # intervals, in order to be conservative.  Expansion is conservative
        # because it is (for example) like computing the 97% credibility
        # interval instead of the 95% credibility interval, which is a more
        # stringent measure of agreement between sample sets.
        lower = <uint64_t>floor((self._cvgAlpha/2.0) * (past-first))
        upper = (past-first-1) - lower

        ret = 0.0
        for 0 <= runInd < self._nruns:
            # Copy one run's samples into rcovScratch, and sort to make the
            # bounds of the credibility interval available.
            memcpy(rcovScratch, &self.lnLs[runInd][first], \
              (past-first)*sizeof(double))
            qsort(rcovScratch, past - first, sizeof(double), _lnLCmp)

            # Compute the proportion of the relevant lnLs that fall within the
            # credibility interval.
            n = 0
            for 0 <= i < self._nruns:
                for first <= j < past:
                    if self.lnLs[i][j] >= rcovScratch[lower] and \
                      self.lnLs[i][j] <= rcovScratch[upper]:
                        n += 1
            ret += <double>n / <double>(self._nruns * (past - first))
        ret /= <double>self._nruns
        # Compute the empirical alpha, based on the credibility interval
        # proportion, and use this to rescale Rcov, so that results can be
        # interpreted in the context of the intended alpha.
        alphaEmp = <double>(lower*2) / <double>(past-first)
        ret *= (1.0 - self._cvgAlpha) / (1.0 - alphaEmp)
        return ret

    cdef bint writeGraph(self, uint64_t sample) except *:
        cdef file gfile
        cdef list cols0, cols1, lnLs0, lnLs1, colors
        cdef uint64_t past, xsp, i
        cdef unsigned runInd
        cdef double lnL, rcov

        past = sample+1
        if past < 4:
            # The code below does not generate an informative graph until there
            # are at least four samples to work with.
            return True

        # Compute which sample will be at the edge of the right graph.
        xsp = (sample/2)+1
        assert xsp > 0

        # Write to a temporary file, in order to be able to make the complete
        # file appear atomically in its final location.
        gfile = open("%s.lnL.R.tmp" % self.outPrefix, "w")

        cols0 = []
        cols1 = []
        for 0 <= runInd < self._nruns:
            lnLs0 = []
            lnLs1 = []
            for 0 <= i < xsp:
                lnLs0.append(self.lnLs[runInd][i])
            for xsp <= i < past:
                lnLs1.append(self.lnLs[runInd][i])
            cols0.append(", ".join(["%.6e" % lnL for lnL in lnLs0]))
            cols1.append(", ".join(["%.6e" % lnL for lnL in lnLs1]))
        gfile.write("lnLs0 = matrix(c(%s), ncol=%d)\n" % \
          (",\n".join(cols0), self._nruns))
        gfile.write("lnLs1 = matrix(c(%s), ncol=%d)\n" % \
          (",\n".join(cols1), self._nruns))

        gfile.write("x0 = seq(0, %d, by=%d)\n" % \
          (self._stride * (xsp-1), self._stride))
        gfile.write("x1 = seq(%d, %d, by=%d)\n" % \
          (self._stride * xsp, self._stride * sample, self._stride))

        gfile.write("par(mfrow=c(1, 2))\n")

        colors = ["black", "red", "green", "blue", "orange", "brown", "cyan", \
          "purple"]

        gfile.write('plot(c(min(x0), max(x0)), c(min(lnLs0), max(lnLs0)), pch="", xlab="step", ylab="lnL")\n')
        for 0 <= i < self._nruns:
            gfile.write('lines(x0, lnLs0[,%d], col="%s")\n' % \
              (i+1, colors[i % len(colors)]))

        gfile.write('plot(c(min(x1), max(x1)), c(min(lnLs1), max(lnLs1)), pch="", xlab="step", ylab="")\n')
        for 0 <= i < self._nruns:
            gfile.write('lines(x1, lnLs1[,%d], col="%s")\n' % \
              (i+1, colors[i % len(colors)]))

        gfile.close()
        os.rename("%s.lnL.R.tmp" % self.outPrefix, "%s.lnL.R" % self.outPrefix)

        return False

    cdef str formatLnLs(self, uint64_t sample, str fmt):
        cdef list strs
        cdef unsigned i

        strs = []
        strs.append("[ ")
        for 0 <= i < self._nruns:
            if i > 0:
                strs.append(" ")
            strs.append(fmt % self.lnLs[i][sample])
        strs.append(" ]")
        return "".join(strs)

    cdef str formatRateStats(self, Mc3RateStats *rateStats):
        cdef list strs
        cdef unsigned i

        strs = []
        strs.append("[ ")
        for 0 <= i < self._nruns:
            if i > 0:
                strs.append(" ")
            strs.append("%.6f" % rateStats[i].ema)
        strs.append(" ]")
        return "".join(strs)

    cdef str formatPropStats(self):
        cdef list strs
        cdef unsigned i, j

        strs = []
        strs.append("{ ")
        for 0 <= i < Mc3Prop:
            if i > 0:
                strs.append(" ")
            if self.props[i] != 0.0:
                strs.append(self.formatRateStats(self.propStats[i]))
            else:
                strs.append("[ ")
                for 0 <= j < self._nruns:
                    if j > 0:
                        strs.append(" ")
                    strs.append("---")
                strs.append(" ]")
        strs.append(" }")

        return "".join(strs)

    cdef void updateDiags(self, uint64_t step):
        cdef unsigned runInd, i
        cdef Mc3RateStats *rateStats
        cdef double xt

        for 0 <= runInd < self._nruns:
            # Compute the exponential moving averages for heat swap rates.
            # Swaps are double-counted, since both participants count them.
            rateStats = &self.swapStats[runInd]
            rateStats.ema += self._emaAlpha * \
              (<double>rateStats.n/(2.0*<double>self._stride) - \
              rateStats.ema)
            rateStats.n = 0

            # Compute the exponential moving averages for proposal acceptance
            # rates.
            for 0 <= i < Mc3Prop:
                rateStats = &self.propStats[i][runInd]
                if rateStats.d == 0:
                    xt = 0.0
                else:
                    xt = (<double>rateStats.n/<double>rateStats.d)
                rateStats.ema += self._emaAlpha * (xt - rateStats.ema)

    cpdef bint run(self, bint verbose=False) except *:
        """
            Run until convergence is reached, or the maximum number of steps
            is reached, whichever comes first.  Collate the results from the
            unheated chain(s) and write the raw results to disk.
        """
        cdef double propsSum, rcov, graphT0, graphT1
        cdef Mc3Chain chain
        cdef uint32_t seed, swapSeed
        cdef unsigned i, j
        cdef uint64_t step, sample
        cdef list run
        cdef bint converged
        cdef str swapStats, propStats, rcovStr

        self.verbose = verbose

        self.initLogs()

        # Allocate swapInfo matrix.
        if self._ncoupled > 1:
            if self.swapInfo != NULL:
                free(self.swapInfo)
            self.swapInfo = <Mc3SwapInfo *>calloc(self._nruns * \
              self._ncoupled * self._ncoupled * 2, sizeof(Mc3SwapInfo))
            if self.swapInfo == NULL:
                raise MemoryError("Error allocating swapInfo")

        # Allocate swapStats vector.
        if self.swapStats != NULL:
            free(self.swapStats)
        self.swapStats = <Mc3RateStats *>calloc(self._nruns, \
          sizeof(Mc3RateStats))
        if self.swapStats == NULL:
            raise MemoryError("Error allocating swapStats")

        # Allocate propStats vectors.
        for 0 <= i < Mc3Prop:
            if self.propStats[i] != NULL:
                free(self.propStats[i])
            self.propStats[i] = <Mc3RateStats *>calloc(self._nruns, \
              sizeof(Mc3RateStats))
            if self.propStats[i] == NULL:
                raise MemoryError("Error allocating propStats[%d]" % i)

        # Initialize propsCdf, which is used to choose proposal types.
        assert sizeof(self.props)/sizeof(double) == Mc3Prop
        assert sizeof(self.propsCdf)/sizeof(double) == Mc3Prop + 1
        if self._nmodels < 2:
            # There are not enough models in the mixture for weights or rate
            # multipliers to be relevant.
            self.props[Mc3WeightProp] = 0.0
            self.props[Mc3RmultProp] = 0.0
        if self.alignment.charType.get().nstates() <= 2:
            # There are not enough states to allow rate class grouping.
            self.props[Mc3RateJumpProp] = 0.0
        if self.alignment.ntaxa <= 3:
            # There are not enough taxa to allow topology changes.
            self.props[Mc3EtbrProp] = 0.0
            self.props[Mc3PolytomyJumpProp] = 0.0
        if self._ncat < 2:
            # There aren't multiple rate categories, so Gamma-distributed rates
            # are irrelevant.
            self.props[Mc3RateShapeInvProp] = 0.0
            self.props[Mc3RateShapeInvJumpProp] = 0.0
        propsSum = 0.0
        for 0 <= i < sizeof(self.props) / sizeof(double):
            propsSum += self.props[i]
        if propsSum == 0.0:
            raise ValueError("No proposals are enabled")
        self.propsCdf[0] = 0.0
        for 1 <= i < sizeof(self.propsCdf) / sizeof(double):
            self.propsCdf[i] = self.propsCdf[i-1] + (self.props[i-1] / propsSum)

        # Create/initialize chains.
        seed = random.randint(0, 0xffffffffU)
        self.runs = []
        for 0 <= i < self._nruns:
            swapSeed = seed
            seed += 1
            run = []
            self.runs.append(run)
            for 0 <= j < self._ncoupled:
                # Seed every chain differently.
                chain = Mc3Chain(self, i, j, swapSeed, seed)
                seed += 1
                run.append(chain)

        if self._graphDelay >= 0.0:
            graphT0 = 0.0

        try:
            # Run the chains to at least _minStep.  Step 0 was created during
            # chain initialization.
            for 1 <= step < self._minStep:
                sample = step / self._stride
                for 0 <= i < self._nruns:
                    run = <list>self.runs[i]
                    for 0 <= j < self._ncoupled:
                        chain = <Mc3Chain>run[j]
                        chain.advance()
                if step % self._stride == 0:
                    self.updateDiags(step)
                    # Write to .s log file.
                    swapStats = self.formatRateStats(self.swapStats)
                    propStats = self.formatPropStats()
                    self.sFile.write("%d\t%s -------- %s %s\n" % (step, \
                      self.formatLnLs(sample, "%.11e"), swapStats, propStats))
                    self.sFile.flush()
                    if self.verbose:
                        sys.stdout.write("s\t%d\t%s --------\n" % \
                          (step, self.formatLnLs(sample, "%.6f")))

                    # Write graph.
                    if self._graphDelay >= 0.0:
                        graphT1 = time.time()
                        if graphT1 - graphT0 >= self._graphDelay:
                            if not self.writeGraph(sample):
                                graphT0 = time.time()

            # Run the chains no further than than _maxStep.
            for step <= step <= self._maxStep:
                sample = step / self._stride
                for 0 <= i < self._nruns:
                    run = <list>self.runs[i]
                    for 0 <= j < self._ncoupled:
                        chain = <Mc3Chain>run[j]
                        chain.advance()
                if step % self._stride == 0:
                    self.updateDiags(step)
                    if self._nruns > 1 and \
                      step % (self._stride * self._cvgSampStride) == 0:
                        # Check for convergence.
                        rcov = self.computeRcov(sample)
                        converged = (rcov + self._cvgAlpha + self._cvgEpsilon \
                          >= 1.0)
                    else:
                        rcov = -1.0
                        converged = False

                    # Write to .s log file.
                    rcovStr = ("%.6f" % rcov if rcov != -1.0 else "--------")
                    swapStats = self.formatRateStats(self.swapStats)
                    propStats = self.formatPropStats()
                    self.sFile.write("%d\t%s %s %s %s\n" % (step, \
                      self.formatLnLs(sample, "%.11e"), rcovStr, swapStats, \
                      propStats))
                    self.sFile.flush()
                    if self.verbose:
                        sys.stdout.write("s\t%d\t%s %s\n" % (step, \
                          self.formatLnLs(sample, "%.6f"), rcovStr))

                    if converged:
                        if self._graphDelay >= 0.0:
                            # Write a graph one last time.
                            self.writeGraph(sample)
                        self.lFile.write("Runs converged\n")
                        self.lFile.flush()
                        if self.verbose:
                            sys.stdout.write("Runs converged\n")
                        break

                    if self._graphDelay >= 0.0:
                        graphT1 = time.time()
                        if graphT1 - graphT0 >= self._graphDelay:
                            if not self.writeGraph(sample):
                                graphT0 = time.time()
        except:
            error = sys.exc_info()
            self.lFile.write("Premature termination: Exception %r\n" % \
              sys.exc_info()[1])
            raise
        finally:
            self.lFile.write("Finish run: %s\n" % \
              time.strftime("%Y/%m/%d %H:%M:%S (%Z)", \
              time.localtime(time.time())))
            self.lFile.flush()

    cdef double getGraphDelay(self):
        return self._graphDelay
    cdef void setGraphDelay(self, double graphDelay):
        self._graphDelay = graphDelay
    property graphDelay:
        """
            If non-negative, periodically (no more often than every graphDelay
            seconds) output an R program to <outPrefix>.lnL.R that generates a
            graphical representation of the log-likelihoods.  To monitor
            progress, use something like the following at an interactive R
            prompt:

              while (1) {source("outPrefix.lnL.R"); Sys.sleep(10)}
        """
        def __get__(self):
            return self.getGraphDelay()
        def __set__(self, double graphDelay):
            self.setGraphDelay(graphDelay)

    cdef double getEmaAlpha(self):
        return self._emaAlpha
    cdef void setEmaAlpha(self, double emaAlpha) except *:
        if not (0.0 < emaAlpha and emaAlpha <= 1.0):
            raise ValueError("Validation failure: 0.0 < emaAlpha <= 1.0")
        self._emaAlpha = emaAlpha
    property emaAlpha:
        """
            Various diagnostic statistics (namely heat swap rates and proposal
            acceptance rates) are reported as exponential moving averages,
            using the following formulas:

              s(0) = x(0)
              s(t) = [emaAlpha * x(t)] + [s(t-1) * (1-emaAlpha)]
                   = s(t-1) + emaAlpha * [x(t) - s(t-1)]

            Discrete time (t) is measured in samples (one per stride).
            (0 < emaAlpha <= 1).  Larger values cause faster decay.
        """
        def __get__(self):
            return self.getEmaAlpha()
        def __set__(self, double emaAlpha):
            self.setEmaAlpha(emaAlpha)

    cdef uint64_t getCvgSampStride(self):
        return self._cvgSampStride
    cdef void setCvgSampStride(self, uint64_t cvgSampStride):
        if not (cvgSampStride >= 1):
            raise ValueError("Validation failure: cvgSampStride >= 1")
        self._cvgSampStride = cvgSampStride
    property cvgSampStride:
        """
            cvgSampStride limits convergence diagnostics computation to once
            every stride*cvgSampStride steps.  This is primarily useful for
            reducing convergence diagnostic overhead.
        """
        def __get__(self):
            return self.getCvgSampStride()
        def __set__(self, uint64_t cvgSampStride):
            self.setCvgSampStride(cvgSampStride)

    cdef double getCvgAlpha(self):
        return self._cvgAlpha
    cdef void setCvgAlpha(self, double cvgAlpha) except *:
        if not (0.0 <= cvgAlpha and cvgAlpha <= 1.0):
            raise ValueError("Validation failure: 0.0 <= cvgAlpha <= 1.0")
        self._cvgAlpha = cvgAlpha
    property cvgAlpha:
        """
            (1.0 - cvgAlpha) is the confidence interval used for diagnosing
            convergence.
        """
        def __get__(self):
            return self.getCvgAlpha()
        def __set__(self, double cvgAlpha):
            self.setCvgAlpha(cvgAlpha)

    cdef double getCvgEpsilon(self):
        return self._cvgEpsilon
    cdef void setCvgEpsilon(self, double cvgEpsilon) except *:
        if not (0.0 <= self._cvgEpsilon and self._cvgEpsilon <= 1.0):
            raise ValueError("Validation failure: 0.0 <= cvgEpsilon <= 1.0")
        self._cvgEpsilon = cvgEpsilon
    property cvgEpsilon:
        """
            Convergence is achieved once the mean coverage is within cvgEpsilon
            of the confidence interval.
        """
        def __get__(self):
            return self.getCvgEpsilon()
        def __set__(self, double cvgEpsilon):
            self.setCvgEpsilon(cvgEpsilon)

    cdef uint64_t getMinStep(self):
        return self._minStep
    cdef void setMinStep(self, uint64_t minStep) except *:
        if not (minStep <= self._maxStep):
            raise ValueError("Validation failure: minStep <= maxStep")
        self._minStep = minStep
    property minStep:
        """
            Minimum final step index for each chain.
        """
        def __get__(self):
            return self.getMinStep()
        def __set__(self, uint64_t minStep):
            self.setMinStep(minStep)

    cdef uint64_t getMaxStep(self):
        return self._maxStep
    cdef void setMaxStep(self, uint64_t maxStep) except *:
        if not (self._minStep <= maxStep):
            raise ValueError("Validation failure: minStep <= maxStep")
        self._maxStep = maxStep
    property maxStep:
        """
            Maximum final step index for each chain.
        """
        def __get__(self):
            return self.getMaxStep()
        def __set__(self, uint64_t maxStep):
            self.setMaxStep(maxStep)

    cdef unsigned getStride(self):
        return self._stride
    cdef void setStride(self, unsigned stride) except *:
        if not stride >= 1:
            raise ValueError("Validation failure: stride >= 1")
        self._stride = stride
    property stride:
        """
            Chain sample interval.
        """
        def __get__(self):
            return self.getStride()
        def __set__(self, unsigned stride):
            self.setStride(stride)

    cdef unsigned getNruns(self):
        return self._nruns
    cdef void setNruns(self, unsigned nruns) except *:
        if not nruns >= 1:
            raise ValueError("Validation failure: nruns >= 1")
        self._nruns = nruns
    property nruns:
        """
            Number of independent runs.  If there is only one run, convergence
            diagnostics will not be computed.
        """
        def __get__(self):
            return self.getNruns()
        def __set__(self, unsigned nruns):
            self.setNruns(nruns)

    cdef unsigned getNcoupled(self):
        return self._ncoupled
    cdef void setNcoupled(self, unsigned ncoupled) except *:
        if not ncoupled >= 1:
            raise ValueError("Validation failure: ncoupled >= 1")
        self._ncoupled = ncoupled
    property ncoupled:
        """
            Number of Metropolis-coupled chains per run.
        """
        def __get__(self):
            return self.getNcoupled()
        def __set__(self, unsigned ncoupled):
            self.setNcoupled(ncoupled)

    cdef double getHeatDelta(self):
        return self._heatDelta
    cdef void setHeatDelta(self, double heatDelta) except *:
        if not heatDelta > 0.0:
            raise ValueError("Validation failure: heatDelta > 0.0")
        self._heatDelta = heatDelta
    property heatDelta:
        """
            Temperature interval between Metropolis-coupled chains.
        """
        def __get__(self):
            return self.getHeatDelta()
        def __set__(self, double heatDelta):
            self.setHeatDelta(heatDelta)

    cdef unsigned getSwapStride(self):
        return self._swapStride
    cdef void setSwapStride(self, unsigned swapStride) except *:
        if not swapStride >= 1:
            raise ValueError("Validation failure: swapStride >= 1")
        self._swapStride = swapStride
    property swapStride:
        """
            Chain heat swap interval.
        """
        def __get__(self):
            return self.getSwapStride()
        def __set__(self, unsigned swapStride):
            self.setSwapStride(swapStride)

    cdef unsigned getNmodels(self):
        return self._nmodels
    cdef void setNmodels(self, unsigned nmodels) except *:
        if not nmodels >= 1:
            raise ValueError("Validation failure: nmodels >= 1")
        self._nmodels = nmodels
    property nmodels:
        """
            Number of mixture models.
        """
        def __get__(self):
            return self.getNmodels()
        def __set__(self, unsigned nmodels):
            self.setNmodels(nmodels)

    cdef unsigned getNcat(self):
        return self._ncat
    cdef void setNcat(self, unsigned ncat) except *:
        if not ncat >= 1:
            raise ValueError("Validation failure: ncat >= 1")
        self._ncat = ncat
    property ncat:
        """
            Number of discrete Gamma-distributed relative mutation rate
            categories.
        """
        def __get__(self):
            return self.getNcat()
        def __set__(self, unsigned ncat):
            self.setNcat(ncat)

    cdef bint getCatMedian(self):
        return self._catMedian
    cdef void setCatMedian(self, bint catMedian):
        self._catMedian = catMedian
    property catMedian:
        """
            Use category means for Gamma-distributed relative mutation rate
            discretization if false; use category medians otherwise.
        """
        def __get__(self):
            return self.getCatMedian()
        def __set__(self, bint catMedian):
            self.setCatMedian(catMedian)

    cdef double getWeightLambda(self):
        return self._weightLambda
    cdef void setWeightLambda(self, double weightLambda) except *:
        if not weightLambda >= 0.0:
            raise ValueError("Validation failure: weightLambda >= 0.0")
        self._weightLambda = weightLambda
    property weightLambda:
        """
            Relative model weight multiplier.  This controls the range of
            weight change for weight change proposals.

            Proposed weights are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = ----------------
                   weightLambda * m

                            /        1             weightLambda/2 \ 
            in the interval | ----------------- , e               | .
                            |   weightLambda/2                    |
                            \  e                                  /

            Where weightLambda = 2 * log(a), proposed weights w* are in
            the interval (w/a, aw).
        """
        def __get__(self):
            return self.getWeightLambda()
        def __set__(self, double weightLambda):
            self.setWeightLambda(weightLambda)

    cdef double getFreqLambda(self):
        return self._freqLambda
    cdef void setFreqLambda(self, double freqLambda) except *:
        if not freqLambda >= 0.0:
            raise ValueError("Validation failure: freqLambda >= 0.0")
        self._freqLambda = freqLambda
    property freqLambda:
        """
            State frequency multiplier.  This controls the range of state
            frequency change for state frequency change proposals.

            Proposed state frequencies are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = --------------
                   freqLambda * m

                            /        1           freqLambda/2 \ 
            in the interval | --------------- , e             | .
                            |   freqLambda/2                  |
                            \  e                              /

            Where freqLambda = 2 * log(a), proposed frequencies f* are in
            the interval (f/a, af).
        """
        def __get__(self):
            return self.getFreqLambda()
        def __set__(self, double freqLambda):
            self.setFreqLambda(freqLambda)

    cdef double getRmultLambda(self):
        return self._rmultLambda
    cdef void setRmultLambda(self, double rmultLambda) except *:
        if not rmultLambda >= 0.0:
            raise ValueError("Validation failure: rmultLambda >= 0.0")
        self._rmultLambda = rmultLambda
    property rmultLambda:
        """
            Rate multiplier multiplier (yes, really).  This controls the range
            of rate multiplier change for rate multiplier change proposals.

            Proposed rate multipliers are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = ---------------
                   rmultLambda * m

                            /        1           rmultLambda/2 \ 
            in the interval | --------------- , e              | .
                            |   rmultLambda/2                  |
                            \  e                               /

            Where rmultLambda = 2 * log(a), proposed rate multipliers r* are in
            the interval (r/a, ar).
        """
        def __get__(self):
            return self.getRmultLambda()
        def __set__(self, double rmultLambda):
            self.setRmultLambda(rmultLambda)

    cdef double getRateLambda(self):
        return self._rateLambda
    cdef void setRateLambda(self, double rateLambda) except *:
        if not rateLambda >= 0.0:
            raise ValueError("Validation failure: rateLambda >= 0.0")
        self._rateLambda = rateLambda
    property rateLambda:
        """
            Relative mutation rate multiplier.  This controls the range of
            rate change for rate change proposals.

            Proposed rates are drawn from a log-transformed sliding window.  A
            multiplier, m, is drawn from

                         1
            g(m) = ---------------
                   rateLambda * m

                            /        1           rateLambda/2 \ 
            in the interval | --------------- , e             | .
                            |   rateLambda/2                  |
                            \  e                              /

            Where rateLambda = 2 * log(a), proposed rates r* are in the
            interval (r/a, ar).
        """
        def __get__(self):
            return self.getRateLambda()
        def __set__(self, double rateLambda):
            self.setRateLambda(rateLambda)

    cdef double getRateShapeInvLambda(self):
        return self._rateShapeInvLambda
    cdef void setRateShapeInvLambda(self, double rateShapeInvLambda) except *:
        if not rateShapeInvLambda >= 0.0:
            raise ValueError("Validation failure: rateShapeInvLambda >= 0.0")
        self._rateShapeInvLambda = rateShapeInvLambda
    property rateShapeInvLambda:
        """
            Gamma-distributed inverse shape change multiplier.  This controls
            the range of inverse shape change for inverse shape change
            proposals.

            Proposed inverse shapes are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                             1
            g(m) = ----------------------
                   rateShapeInvLambda * m

                            /        1                  rateShapeInvLambda/2 \ 
            in the interval | ---------------------- , e                     | .
                            |   rateShapeInvLambda/2                         |
                            \  e                                             /

            Where rateShapeInvLambda = 2 * log(a), proposed inverse shapes s*
            are in the interval (s/a, as).
        """
        def __get__(self):
            return self.getRateShapeInvLambda()
        def __set__(self, double rateShapeInvLambda):
            self.setRateShapeInvLambda(rateShapeInvLambda)

    cdef double getBrlenLambda(self):
        return self._brlenLambda
    cdef void setBrlenLambda(self, double brlenLambda) except *:
        if not brlenLambda >= 0.0:
            raise ValueError("Validation failure: brlenLambda >= 0.0")
        self._brlenLambda = brlenLambda
    property brlenLambda:
        """
            Branch change length multiplier.  This controls the range of branch
            length change for basic branch length change proposals.

            Proposed branch lengths are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = ---------------
                   brlenLambda * m

                            /        1            brlenLambda/2 \ 
            in the interval | ---------------- , e              | .
                            |   brlenLambda/2                   |
                            \  e                                /

            Where brlenLambda = 2 * log(a), proposed branch lengths v* are in
            the interval (v/a, av).
        """
        def __get__(self):
            return self.getBrlenLambda()
        def __set__(self, double brlenLambda):
            self.setBrlenLambda(brlenLambda)

    cdef double getEtbrPExt(self):
        return self._etbrPExt
    cdef void setEtbrPExt(self, double etbrPExt) except *:
        if not (0.0 <= etbrPExt and etbrPExt <= 1.0):
            raise ValueError("Validation failure: 0.0 <= etbrPExt <= 1.0")
        self._etbrPExt = etbrPExt
    property etbrPExt:
        """
            eTBR extension probability.  This controls how far away from the
            bisection edge reattachment occurs on average.
        """
        def __get__(self):
            return self.getEtbrPExt()
        def __set__(self, double etbrPExt):
            self.setEtbrPExt(etbrPExt)

    cdef double getEtbrLambda(self):
        return self._etbrLambda
    cdef void setEtbrLambda(self, double etbrLambda) except *:
        if not etbrLambda >= 0.0:
            raise ValueError("Validation failure: etbrLambda >= 0.0")
        self._etbrLambda = etbrLambda
    property etbrLambda:
        """
            eTBR branch length multiplier.  This controls the range of branch
            length change.
        """
        def __get__(self):
            return self.getEtbrLambda()
        def __set__(self, double etbrLambda):
            self.setEtbrLambda(etbrLambda)

    cdef double getRateShapeInvPrior(self):
        return self._rateShapeInvPrior
    cdef void setRateShapeInvPrior(self, double rateShapeInvPrior) except *:
        if not rateShapeInvPrior > 0.0:
            raise ValueError("Validation failure: rateShapeInvPrior > 0.0")
        self._rateShapeInvPrior = rateShapeInvPrior
    property rateShapeInvPrior:
        """
            Prior for exponentially distributed inverse (variance) of shape
            parameter for Gamma-distributed relative mutation rates.
        """
        def __get__(self):
            return self.getRateShapeInvPrior()
        def __set__(self, double rateShapeInvPrior):
            self.setRateShapeInvPrior(rateShapeInvPrior)

    cdef double getBrlenPrior(self):
        return self._brlenPrior
    cdef void setBrlenPrior(self, double brlenPrior) except *:
        if not brlenPrior > 0.0:
            raise ValueError("Validation failure: brlenPrior > 0.0")
        self._brlenPrior = brlenPrior
    property brlenPrior:
        """
            Prior for exponentially distributed branch lengths.
        """
        def __get__(self):
            return self.getBrlenPrior()
        def __set__(self, double brlenPrior):
            self.setBrlenPrior(brlenPrior)

    cdef double getRateJumpPrior(self):
        return self._rateJumpPrior
    cdef void setRateJumpPrior(self, double rateJumpPrior) except *:
        if not rateJumpPrior > 0.0:
            raise ValueError("Validation failure: rateJumpPrior > 0.0")
        self._rateJumpPrior = rateJumpPrior
    property rateJumpPrior:
        """
            Prior for number of relative mutation rate parameters.  1.0
            indicates a flat prior; <1.0 favors more complex models, and >1.0
            favors less complex models.
        """
        def __get__(self):
            return self.getRateJumpPrior()
        def __set__(self, double rateJumpPrior):
            self.setRateJumpPrior(rateJumpPrior)

    cdef double getPolytomyJumpPrior(self):
        return self._polytomyJumpPrior
    cdef void setPolytomyJumpPrior(self, double polytomyJumpPrior) except *:
        if not polytomyJumpPrior > 0.0:
            raise ValueError("Validation failure: polytomyJumpPrior > 0.0")
        self._polytomyJumpPrior = polytomyJumpPrior
    property polytomyJumpPrior:
        """
            Prior for number of internal branches.  1.0 indicates a flat prior
            (all tree topologies, polytomous and fully resolved are given equal
            weight), <1.0 favors more resolved trees, and >1.0 favors less
            resolved trees.
        """
        def __get__(self):
            return self.getPolytomyJumpPrior()
        def __set__(self, double polytomyJumpPrior):
            self.setPolytomyJumpPrior(polytomyJumpPrior)

    cdef double getRateShapeInvJumpPrior(self):
        return self._rateShapeInvJumpPrior
    cdef void setRateShapeInvJumpPrior(self, double rateShapeInvJumpPrior) \
      except *:
        if not rateShapeInvJumpPrior > 0.0:
            raise ValueError("Validation failure: rateShapeInvJumpPrior > 0.0")
        self._rateShapeInvJumpPrior = rateShapeInvJumpPrior
    property rateShapeInvJumpPrior:
        """
            Prior for presence/absence of the Gamma-distributed relative
            mutation rates shape parameter.  Such models are commonly referred
            to as +G (plus Gamma) models (ex: GTR+G).  1.0 indicates a flat
            prior (no preference for/against +G models), <1.0 favors +G models,
            and >1.0 favors +G-less models.

            The rateShapeInvJumpPrior parameterization is chosen to be
            consistent with rateJumpPrior and polytomyJumpPrior.  Conceptually,
            all three are structured the same way, but rateShapeInvJumpPrior
            only applies to two levels of nested models, whereas the other
            priors may apply to many levels.
        """
        def __get__(self):
            return self.getRateShapeInvJumpPrior()
        def __set__(self, double rateShapeInvJumpPrior):
            self.setRateShapeInvJumpPrior(rateShapeInvJumpPrior)

    cdef double getWeightProp(self):
        return self.props[Mc3WeightProp]
    cdef void setWeightProp(self, double weightProp) except *:
        if not weightProp >= 0.0:
            raise ValueError("Validation failure: weightProp >= 0.0")
        self.props[Mc3WeightProp] = weightProp
    property weightProp:
        """
            Relative proportion of proposals that modify relative model weights.
        """
        def __get__(self):
            return self.getWeightProp()
        def __set__(self, double weightProp):
            self.setWeightProp(weightProp)

    cdef double getFreqProp(self):
        return self.props[Mc3FreqProp]
    cdef void setFreqProp(self, double freqProp) except *:
        if not freqProp >= 0.0:
            raise ValueError("Validation failure: freqProp >= 0.0")
        self.props[Mc3FreqProp] = freqProp
    property freqProp:
        """
            Relative proportion of proposals that modify state frequencies.
        """
        def __get__(self):
            return self.getFreqProp()
        def __set__(self, double freqProp):
            self.setFreqProp(freqProp)

    cdef double getRmultProp(self):
        return self.props[Mc3RmultProp]
    cdef void setRmultProp(self, double rmultProp) except *:
        if not rmultProp >= 0.0:
            raise ValueError("Validation failure: rmultProp >= 0.0")
        self.props[Mc3RmultProp] = rmultProp
    property rmultProp:
        """
            Relative proportion of proposals that modify rate multipliers.
        """
        def __get__(self):
            return self.getRmultProp()
        def __set__(self, double rmultProp):
            self.setRmultProp(rmultProp)

    cdef double getRateProp(self):
        return self.props[Mc3RateProp]
    cdef void setRateProp(self, double rateProp) except *:
        if not rateProp >= 0.0:
            raise ValueError("Validation failure: rateProp >= 0.0")
        self.props[Mc3RateProp] = rateProp
    property rateProp:
        """
            Relative proportion of proposals that modify relative mutation
            rates.
        """
        def __get__(self):
            return self.getRateProp()
        def __set__(self, double rateProp):
            self.setRateProp(rateProp)

    cdef double getRateShapeInvProp(self):
        return self.props[Mc3RateShapeInvProp]
    cdef void setRateShapeInvProp(self, double rateShapeInvProp) except *:
        if not rateShapeInvProp >= 0.0:
            raise ValueError("Validation failure: rateShapeInvProp >= 0.0")
        self.props[Mc3RateShapeInvProp] = rateShapeInvProp
    property rateShapeInvProp:
        """
            Relative proportion of proposals that modify inverse (variance) of
            shape parameters for Gamma-distributed relative mutation rates.
        """
        def __get__(self):
            return self.getRateShapeInvProp()
        def __set__(self, double rateShapeInvProp):
            self.setRateShapeInvProp(rateShapeInvProp)

    cdef double getBrlenProp(self):
        return self.props[Mc3BrlenProp]
    cdef void setBrlenProp(self, double brlenProp) except *:
        if not brlenProp >= 0.0:
            raise ValueError("Validation failure: brlenProp >= 0.0")
        self.props[Mc3BrlenProp] = brlenProp
    property brlenProp:
        """
            Relative proportion of proposals that modify branch lengths.
        """
        def __get__(self):
            return self.getBrlenProp()
        def __set__(self, double brlenProp):
            self.setBrlenProp(brlenProp)

    cdef double getEtbrProp(self):
        return self.props[Mc3EtbrProp]
    cdef void setEtbrProp(self, double etbrProp) except *:
        if not etbrProp >= 0.0:
            raise ValueError("Validation failure: etbrProp >= 0.0")
        self.props[Mc3EtbrProp] = etbrProp
    property etbrProp:
        """
            Relative proportion of proposals that modify tree topology via
            eTBR.
        """
        def __get__(self):
            return self.getEtbrProp()
        def __set__(self, double etbrProp):
            self.setEtbrProp(etbrProp)

    cdef double getRateJumpProp(self):
        return self.props[Mc3RateJumpProp]
    cdef void setRateJumpProp(self, double rateJumpProp) except *:
        if not rateJumpProp >= 0.0:
            raise ValueError("Validation failure: rateJumpProp >= 0.0")
        self.props[Mc3RateJumpProp] = rateJumpProp
    property rateJumpProp:
        """
            Relative proportion of proposals that split/merge relative mutation
            rate parameters.
        """
        def __get__(self):
            return self.getRateJumpProp()
        def __set__(self, double rateJumpProp):
            self.setRateJumpProp(rateJumpProp)

    cdef double getPolytomyJumpProp(self):
        return self.props[Mc3PolytomyJumpProp]
    cdef void setPolytomyJumpProp(self, double polytomyJumpProp) except *:
        if not polytomyJumpProp >= 0.0:
            raise ValueError("Validation failure: polytomyJumpProp >= 0.0")
        self.props[Mc3PolytomyJumpProp] = polytomyJumpProp
    property polytomyJumpProp:
        """
            Relative proportion of proposals that add/remove internal branches.
        """
        def __get__(self):
            return self.getPolytomyJumpProp()
        def __set__(self, double polytomyJumpProp):
            self.setPolytomyJumpProp(polytomyJumpProp)

    cdef double getRateShapeInvJumpProp(self):
        return self.props[Mc3RateShapeInvJumpProp]
    cdef void setRateShapeInvJumpProp(self, double rateShapeInvJumpProp) \
      except *:
        if not rateShapeInvJumpProp >= 0.0:
            raise ValueError("Validation failure: rateShapeInvJumpProp >= 0.0")
        self.props[Mc3RateShapeInvJumpProp] = rateShapeInvJumpProp
    property rateShapeInvJumpProp:
        """
            Relative proportion of proposals that add/remove Gamma-distributed
            ralative rates parameters (+G models).
        """
        def __get__(self):
            return self.getRateShapeInvJumpProp()
        def __set__(self, double rateShapeInvJumpProp):
            self.setRateShapeInvJumpProp(rateShapeInvJumpProp)
