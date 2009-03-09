import random

from libc cimport *
from libm cimport *
from SFMT cimport *
from Crux.Tree cimport Tree, Node, Edge, Ring
from Crux.Tree.Lik cimport Lik
from Crux.Mc3 cimport Mc3

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

cdef class Chain:
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
        for 0 <= i < PropCnt:
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
                if self.master.props[PropWeight] > 0.0:
                    self.lik.setWeight(m, -log(1.0 - genrand_res53(self.prng)))

                # Rate multiplier.
                if self.master.props[PropRmult] > 0.0:
                    self.lik.setRmult(m, -log(1.0 - genrand_res53(self.prng)))

            # State frequencies.
            if self.master.props[PropFreq] > 0.0:
                for 0 <= i < self.lik.char_.nstates():
                    self.lik.setFreq(m, i, -log(1.0 - genrand_res53(self.prng)))

            # Relative mutation rates.
            if self.master.props[PropRate] > 0.0:
                self.lik.setRclass(m, range(rlen))
                for 0 <= i < rlen:
                    self.lik.setRate(m, i, -log(1.0 - genrand_res53(self.prng)))

            if self.master.props[PropRateShapeInv] > 0.0:
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
            self.accepts[PropWeight] += 1
        else:
            # Reject.
            self.lik.setWeight(mInd, w0)
            self.rejects[PropWeight] += 1

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
            self.accepts[PropFreq] += 1
        else:
            self.rejects[PropFreq] += 1
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
            self.accepts[PropRmult] += 1
        else:
            # Reject.
            self.lik.setRmult(mInd, rm0)
            self.rejects[PropRmult] += 1

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
            self.accepts[PropRate] += 1
        else:
            self.rejects[PropRate] += 1
        self.lik.setWeight(m0Ind, w)
        self.lik.delModel()

        return False

    cdef bint rateShapeInvPropose(self) except *:
        cdef unsigned m0Ind, m1Ind
        cdef double w, u, lnM, m, a0, a1, a0Inv, a1Inv, lnPrior, lnL1, p

        # Uniformly choose a random model within the mixture.
        m0Ind = gen_rand64_range(self.prng, self.lik.nmodels())
        a0 = self.lik.getAlpha(m0Ind)
        if a0 == INFINITY:
            return True
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
            self.accepts[PropRateShapeInv] += 1
        else:
            self.rejects[PropRateShapeInv] += 1
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
            self.accepts[PropBrlen] += 1
        else:
            # Reject.
            edge.length = v0
            self.rejects[PropBrlen] += 1

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
            self.accepts[PropEtbr] += 1
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
            self.rejects[PropEtbr] += 1

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
            pSplit = self.master.props[PropRateJump] / \
              (1.0 - self.master.props[PropRate])
            pMerge = self.master.props[PropRateJump] / 2.0
        elif nrates == len(rclass):
            # Splitting was not an option in rateJumpPropose() for this step.
            # (pSplit/pMerge is always 0.5 .)
            pSplit = self.master.props[PropRateJump] / 2.0
            pMerge = self.master.props[PropRateJump]
        else:
            # (pSplit/pMerge is always 1.0 .)
            pSplit = self.master.props[PropRateJump] / 2.0
            pMerge = self.master.props[PropRateJump] / 2.0
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
            self.accepts[PropRateJump] += 1
        else:
            self.rejects[PropRateJump] += 1

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
            pMerge = self.master.props[PropRateJump]
            pSplit = self.master.props[PropRateJump] / 2.0
        elif nrates == 1:
            # Merging was not an option in rateJumpPropose() for this step.
            # Note also that rate change proposals were invalid, since there
            # were no free rate parameters.
            pMerge = self.master.props[PropRateJump] / 2.0
            pSplit = self.master.props[PropRateJump] / \
              (1.0 - self.master.props[PropRate])
        else:
            # (pMerge/pSplit is always 1.0 .)
            pMerge = self.master.props[PropRateJump] / 2.0
            pSplit = self.master.props[PropRateJump] / 2.0
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
            self.accepts[PropRateJump] += 1
        else:
            self.rejects[PropRateJump] += 1

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
            self.accepts[PropPolytomyJump] += 1
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
            self.rejects[PropPolytomyJump] += 1

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
            self.accepts[PropPolytomyJump] += 1
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
            self.rejects[PropPolytomyJump] += 1

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
            self.accepts[PropRateShapeInvJump] += 1
        else:
            self.rejects[PropRateShapeInvJump] += 1

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
            self.accepts[PropRateShapeInvJump] += 1
        else:
            self.rejects[PropRateShapeInvJump] += 1

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
            if propInd == PropWeight:
                again = self.weightPropose()
            elif propInd == PropFreq:
                again = self.freqPropose()
            elif propInd == PropRmult:
                again = self.rmultPropose()
            elif propInd == PropRate:
                again = self.ratePropose()
            elif propInd == PropRateShapeInv:
                again = self.rateShapeInvPropose()
            elif propInd == PropBrlen:
                again = self.brlenPropose()
            elif propInd == PropEtbr:
                again = self.etbrPropose()
            elif propInd == PropRateJump:
                again = self.rateJumpPropose()
            elif propInd == PropPolytomyJump:
                again = self.polytomyJumpPropose()
            elif propInd == PropRateShapeInvJump:
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
            for 0 <= i < PropCnt:
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

