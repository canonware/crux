import random

from libc cimport *
from libm cimport *
from SFMT cimport *
from CxRi cimport *
from Crux.Tree cimport Tree, Node, Edge, Ring
from Crux.Tree.Lik cimport Lik
from Crux.Tree.Bipart cimport Vec, Bipart
from Crux.Tree.Sumt cimport Trprob, Sumt
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
      uint32_t swapSeed, uint32_t seed, Lik lik):
        cdef unsigned i

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

        self.lik = lik
        self.tree = self.lik.tree
        self.lnL = self.lik.lnL()
        if self.master._prelim is not None:
            self.lnRisk = self.computeLnRisk()

        self.step = 0
        self.master.sendSample(self.run, self.step, self.heat, self.nswap, \
          self.accepts, self.rejects, self.lik, self.lnL)

    cdef unsigned rfUnscaled(self, Bipart a, Bipart b):
        cdef unsigned nUnique, iA, iB, lenA, lenB
        cdef double falseNegativeRate, falsePositiveRate
        cdef int rel
        cdef Vec vecA, vecB

        if a is b:
            return 0

        nUnique = iA = iB = 0
        lenA = len(a.edgeVecs)
        lenB = len(b.edgeVecs)
        while iA < lenA and iB < lenB:
            vecA = <Vec>a.edgeVecs[iA]
            vecB = <Vec>b.edgeVecs[iB]

            rel = vecA.cmp(vecB)
            if rel == -1:
                # a has a unique bipartition.
                nUnique += 1
                iA += 1
            elif rel == 0:
                iA += 1
                iB += 1
            elif rel == 1:
                # b has a unique bipartition.
                nUnique += 1
                iB += 1
            else:
                assert False

        if iA < lenA:
            nUnique += lenA - iA
        elif iB < lenB:
            nUnique += lenB - iB

        return nUnique

    cdef double computeLnRisk(self) except *:
        cdef double risk
        cdef Bipart bipart
        cdef uint64_t nsamples, i
        cdef Sumt sumt
        cdef list trprobs
        cdef Trprob trprob
        cdef Tree tree

        bipart = self.tree.getBipart()
        risk = 0.0
        nsamples = self.master._prelim.nsamples * self.master._prelim.mc3.nruns
        sumt = self.master._prelim.getSumt()
        trprobs = sumt.getTrprobs()
        for 0 <= i < len(trprobs):
            trprob = <Trprob>trprobs[i]

            tree = trprob.getTree()
            risk += (<double>self.rfUnscaled(tree.getBipart(), bipart) * \
              <double>trprob.getNobs()) / <double>nsamples

        return log(risk)

    cdef bint weightPropose(self) except *:
        cdef Lik lik1
        cdef unsigned mInd
        cdef double u, lnM, m, w0, w1, lnL1, p

        if self.lik.nmodels() == 1:
            return True

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, lik1.nmodels())

        # Generate the weight multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._weightLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified weight.
        w0 = lik1.getWeight(mInd)
        w1 = w0 * m
        lik1.setWeight(mInd, w1)
        lnL1 = lik1.lnL()

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
            self.lik = lik1
            self.accepts[PropWeight] += 1
        else:
            # Reject.
            self.rejects[PropWeight] += 1

        return False

    cdef bint freqsEqual(self, Lik lik, unsigned mInd):
        cdef unsigned i
        cdef double f

        # Determine whether a model is currently using estimated frequencies.
        # We assume that identical frequencies indicate fixed frequencies.
        # Mathematically it does not affect correctness in the (extremely
        # unlikely) event that the the frequencies happen to have been
        # estimated as identical.
        f = lik.getFreq(mInd, 0)
        for 1 <= i < lik.char_.nstates:
            if lik.getFreq(mInd, i) != f:
                return False
        return True

    cdef unsigned nModelsFreqsEstim(self, Lik lik):
        cdef unsigned ret, i

        # Count how many models contain estimated frequencies.
        ret = 0
        for 0 <= i < lik.nmodels():
            if not self.freqsEqual(lik, i):
                ret += 1
        return ret

    cdef bint freqPropose(self) except *:
        cdef Lik lik1
        cdef unsigned nMEstim, mChoice, mInd, fInd
        cdef double f0, f1, u, lnM, m, lnL1

        # Uniformly choose a random model within the mixture.
        if self.master.propsPdf[PropFreqJump] != 0.0:
            # We assume identical frequencies indicate fixed frequencies when
            # frequency jump proposals are enabled.
            nMEstim = self.nModelsFreqsEstim(self.lik)
            if nMEstim == 0:
                return True
            mChoice = gen_rand64_range(self.prng, nMEstim)
            # Map mChoice onto the set of models with estimated frequencies.
            for 0 <= mInd < self.lik.nmodels():
                if not self.freqsEqual(self.lik, mInd):
                    if mChoice == 0:
                        break
                    mChoice -= 1
        else:
            mInd = gen_rand64_range(self.prng, self.lik.nmodels())

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Uniformly choose a random frequency parameter.
        fInd = gen_rand64_range(self.prng, lik1.char_.nstates)

        # Generate the frequency multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._freqLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified frequency.
        f0 = lik1.getFreq(mInd, fInd)
        f1 = f0 * m
        lik1.setFreq(mInd, fInd, f1)
        lnL1 = lik1.lnL()

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
            self.lik = lik1
            self.accepts[PropFreq] += 1
        else:
            self.rejects[PropFreq] += 1

        return False

    cdef bint rmultPropose(self) except *:
        cdef Lik lik1
        cdef unsigned mInd
        cdef double u, lnM, m, rm0, rm1, lnL1, p

        if self.lik.nmodels() == 1:
            return True

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, lik1.nmodels())

        # Generate the rmult multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._rmultLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified rmult.
        rm0 = lik1.getRmult(mInd)
        rm1 = rm0 * m
        lik1.setRmult(mInd, rm1)
        lnL1 = lik1.lnL()

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
            self.lik = lik1
            self.accepts[PropRmult] += 1
        else:
            # Reject.
            self.rejects[PropRmult] += 1

        return False

    cdef unsigned nModelsRatesEstim(self, Lik lik):
        cdef unsigned ret, i

        # Count how many models contain estimated rates.
        ret = 0
        for 0 <= i < lik.nmodels():
            if lik.getNrates(i) != 1:
                ret += 1
        return ret

    cdef bint ratePropose(self) except *:
        cdef Lik lik1
        cdef unsigned nMEstim, mChoice, mInd, nrates, r, rInd
        cdef double r0, r1, u, lnM, m, lnL1
        cdef list rclass

        # Uniformly choose a random model within the mixture.
        nMEstim = self.nModelsRatesEstim(self.lik)
        if nMEstim == 0:
            return True
        mChoice = gen_rand64_range(self.prng, nMEstim)
        # Map mChoice onto the set of models with estimated rates.
        for 0 <= mInd < self.lik.nmodels():
            nrates = self.lik.getNrates(mInd)
            if nrates != 1:
                if mChoice == 0:
                    break
                mChoice -= 1

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        rclass = lik1.getRclass(mInd)
        rInd = gen_rand64_range(self.prng, nrates)

        # Generate the rate multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._rateLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified rate.
        r0 = lik1.getRate(mInd, rInd)
        r1 = r0 * m
        lik1.setRate(mInd, rInd, r1)
        lnL1 = lik1.lnL()

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
            self.lik = lik1
            self.accepts[PropRate] += 1
        else:
            self.rejects[PropRate] += 1

        return False

    cdef unsigned nModelsRatesGamma(self, Lik lik):
        cdef unsigned ret, i

        # Count how many models contain +G rates.
        ret = 0
        for 0 <= i < lik.nmodels():
            if lik.getAlpha(i) != INFINITY:
                ret += 1
        return ret

    cdef bint rateShapeInvPropose(self) except *:
        cdef Lik lik1
        cdef unsigned nMGamma, mChoice, mInd
        cdef double u, lnM, m, a0, a1, a0Inv, a1Inv, lnPrior, lnL1, p

        # Uniformly choose a random model within the mixture.
        nMGamma = self.nModelsRatesGamma(self.lik)
        if nMGamma == 0:
            return True
        mChoice = gen_rand64_range(self.prng, nMGamma)
        # Map mChoice onto the set of models with +G rates.
        for 0 <= mInd < self.lik.nmodels():
            a0 = self.lik.getAlpha(mInd)
            if a0 != INFINITY:
                if mChoice == 0:
                    break
                mChoice -= 1

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Generate the inverse rate shape multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._rateShapeInvLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified rate shape.
        a0Inv = 1.0 / a0
        a1Inv = a0Inv * m
        a1 = 1.0 / a1Inv
        lik1.setAlpha(mInd, a1)
        lnL1 = lik1.lnL()

        lnPrior = -self.master._rateShapeInvPrior * (a1Inv-a0Inv)

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnM)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropRateShapeInv] += 1
        else:
            self.rejects[PropRateShapeInv] += 1

        return False

    cdef unsigned nModelsInvar(self, Lik lik):
        cdef unsigned ret, i

        # Count how many models contain +I pinvar.
        ret = 0
        for 0 <= i < lik.nmodels():
            if lik.getWInvar(i) != 0:
                ret += 1
        return ret

    cdef bint invarPropose(self) except *:
        cdef Lik lik1
        cdef unsigned nMInvar, mChoice, mInd, which
        cdef double wInvar, u, lnM, m, w0, w1, lnL1, lnPrior, p

        # Uniformly choose a random model within the mixture.
        nMInvar = self.nModelsInvar(self.lik)
        if nMInvar == 0:
            return True
        mChoice = gen_rand64_range(self.prng, nMInvar)
        # Map mChoice onto the set of models with +I pinvar.
        for 0 <= mInd < self.lik.nmodels():
            wInvar = self.lik.getWInvar(mInd)
            if wInvar != INFINITY:
                if mChoice == 0:
                    break
                mChoice -= 1

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Randomly choose which weight to modify.
        which = gen_rand64_range(self.prng, 2)

        # Generate weight multiplier.
        u = genrand_res53(self.prng)
        lnM = self.master._invarLambda * (u - 0.5)
        m = exp(lnM)

        # Compute lnL with modified weight.
        if which == 0:
            w0 = lik1.getWVar(mInd)
            w1 = w0 * m
            lik1.setWVar(mInd, w1)
        else:
            w0 = lik1.getWInvar(mInd)
            w1 = w0 * m
            lik1.setWInvar(mInd, w1)
        lnL1 = lik1.lnL()

        if which == 0:
            lnPrior = -(1.0-self.master._invarPrior) * (w1-w0)
        else:
            lnPrior = -self.master._invarPrior * (w1-w0)

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnM)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropInvar] += 1
        else:
            self.rejects[PropInvar] += 1

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
            self.tree.setBase(edge.ring.node)
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
        cdef double u, lnL1, lnRisk1, lnPrior, lnProp, p

        # eTBR degenerates to a single branch length change for less than 3
        # taxa.
        assert self.tree.getNtaxa() >= 3

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

        if self.master._prelim is not None:
            # Incorporate risk.
            lnRisk1 = self.computeLnRisk()
            lnProp += (lnRisk1 - self.lnRisk)

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            if self.master._prelim is not None:
                self.lnRisk = lnRisk1
            # Re-base to make the cache fully usable if the next lnL() call
            # does not specify a rooting.
            self.tree.setBase(nX0)
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

    cdef void rateMergePropose(self, unsigned mInd, list rclass, \
      unsigned nrates) except *:
        cdef Lik lik1
        cdef unsigned a, b, na, nb, r, i, revSplit, n
        cdef double rateA, rateB, rateAB, lnL1, pSplit, pMerge, lnPrior
        cdef double lnProp, u, p
        cdef list rates, ns

        assert nrates > 1

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

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
        rateA = lik1.getRate(mInd, a)
        rateB = lik1.getRate(mInd, b)
        rateAB = (<double>na*rateA + <double>nb*rateB) / <double>(na+nb)

        # Update rate class list.
        for 0 <= i < len(rclass):
            r = <unsigned>rclass[i]
            if r == b:
                rclass[i] = a
            elif r > b:
                rclass[i] -= 1

        # Create new list of rates.
        rates = [lik1.getRate(mInd, r) for r in xrange(nrates)]
        rates[a] = rateAB
        rates.pop(b)

        # Compute lnL with modified rclass.
        lik1.setRclass(mInd, rclass, rates)
        lnL1 = lik1.lnL()

        if nrates == 2:
            # After this merge, rateJumpPropose() can only split.  Note also
            pSplit = self.master.propsPdf[PropRateJump]
            if self.nModelsRatesEstim(lik1) == 0:
                # Rate change proposals become invalid, since there will be no
                # free rate parameters.
                pSplit /= (1.0 - self.master.propsPdf[PropRate])
            pMerge = self.master.propsPdf[PropRateJump] / 2.0
        elif nrates == len(rclass):
            # Splitting was not an option in rateJumpPropose() for this step.
            # (pSplit/pMerge is always 0.5.)
            pSplit = self.master.propsPdf[PropRateJump] / 2.0
            pMerge = self.master.propsPdf[PropRateJump]
        else:
            # (pSplit/pMerge is always 1.0.)
            pSplit = self.master.propsPdf[PropRateJump] / 2.0
            pMerge = self.master.propsPdf[PropRateJump] / 2.0
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
            self.lik = lik1
            self.accepts[PropRateJump] += 1
        else:
            self.rejects[PropRateJump] += 1

    cdef void rateSplitPropose(self, unsigned mInd, list rclass, \
      unsigned nrates) except *:
        cdef Lik lik1
        cdef list ns, splittable, reorder, rates
        cdef unsigned i, r0, n0, r, n, rMax, ra, rb, na, nb, a0, b0
        cdef double rate0, u, rateA, rateB, pMerge, pSplit, lnPrior, lnProp
        cdef double lnL1, p
        cdef bint inA

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

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
        rate0 = lik1.getRate(mInd, r0)

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

        rates0 = [lik1.getRate(mInd, i) for i in xrange(nrates)]
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
        lik1.setRclass(mInd, rclass, rates1)
        lnL1 = lik1.lnL()

        if nrates == len(rclass) - 1:
            # After this split, rateJumpPropose() can only merge.
            # (pMerge/pSplit is always 2.0 .)
            pMerge = self.master.propsPdf[PropRateJump]
            pSplit = self.master.propsPdf[PropRateJump] / 2.0
        elif nrates == 1:
            # Merging was not an option in rateJumpPropose() for this step.
            pMerge = self.master.propsPdf[PropRateJump] / 2.0
            pSplit = self.master.propsPdf[PropRateJump]
            if self.nModelsRatesEstim(lik1) == 1:
                # Rate change proposals were invalid, since there were no free
                # rate parameters.
                pSplit /= (1.0 - self.master.propsPdf[PropRate])
        else:
            # (pMerge/pSplit is always 1.0 .)
            pMerge = self.master.propsPdf[PropRateJump] / 2.0
            pSplit = self.master.propsPdf[PropRateJump] / 2.0
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
            self.lik = lik1
            self.accepts[PropRateJump] += 1
        else:
            self.rejects[PropRateJump] += 1

    cdef bint rateJumpPropose(self) except *:
        cdef unsigned mInd, nrates, r
        cdef list rclass
        cdef bint merge

        assert self.lik.char_.nstates > 2

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, self.lik.nmodels())

        # Determine the number of rate parameters.
        rclass = self.lik.getRclass(mInd)
        nrates = self.lik.getNrates(mInd)
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
            self.rateMergePropose(mInd, rclass, nrates)
        else:
            self.rateSplitPropose(mInd, rclass, nrates)

        return False

    cdef void polytomyMergePropose(self, Tree tree, unsigned nedges, \
      unsigned ntaxa) except *:
        cdef unsigned uI, i, j, nodeADeg, nodeBDeg, np1
        cdef list edges, sibsB, nodes
        cdef Edge edge
        cdef Node n, nodeA, nodeB, base
        cdef Ring r, rnext
        cdef double lnL1, lnRisk1, lnPrior, lnGam, lnJacob, lnProp, u, p

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

        if self.master._prelim is not None:
            # Incorporate risk.
            lnRisk1 = self.computeLnRisk()
            lnProp += (lnRisk1 - self.lnRisk)

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            if self.master._prelim is not None:
                self.lnRisk = lnRisk1
            self.accepts[PropPolytomyJump] += 1
            # Dissociate edge's cached CL's in order to allow them to be
            # reclaimed by GC sooner.
            edge.ring.aux = None
            edge.ring.other.aux = None
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
        cdef double lnL1, lnRisk1, lnPrior, lnGam, lnJacob, lnProp, u, p

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
        cdef CxtRi ri
        CxRiNew(&ri, self.prng)
        if CxRiInit(&ri, deg0):
            CxRiDelete(&ri)
            raise MemoryError("Error initializing ri")
        a0 = CxRiRandomGet(&ri)
        aa0 = CxRiRandomGet(&ri)
        b0 = CxRiRandomGet(&ri)
        bb0 = CxRiRandomGet(&ri)
        CxRiDelete(&ri)

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

        if self.master._prelim is not None:
            # Incorporate risk.
            lnRisk1 = self.computeLnRisk()
            lnProp += (lnRisk1 - self.lnRisk)

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            if self.master._prelim is not None:
                self.lnRisk = lnRisk1
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

        tree = self.tree
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

    cdef void rateShapeInvRemovePropose(self, unsigned mInd, double alpha0) \
      except *:
        cdef Lik lik1
        cdef double rateShapeInv0
        cdef double lnL1, lnPrior, lnJacob, lnProp, u, p

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        rateShapeInv0 = 1.0 / alpha0

        # Disable +G.
        lik1.setAlpha(mInd, INFINITY)

        # Compute lnL without +G.
        lnL1 = lik1.lnL()

        lnPrior = log(self.master._rateShapeInvJumpPrior) \
          - log(self.master._rateShapeInvPrior) - \
          (-self.master._rateShapeInvPrior*rateShapeInv0)

        lnJacob = log(self.master._rateShapeInvPrior) + \
          (-self.master._rateShapeInvPrior*rateShapeInv0)
        lnProp = lnJacob

        if self.nModelsRatesGamma(lik1) == 0:
            # This proposal implicitly disables the PropRateShapeInv proposal.
            lnProp += -log(1.0 - self.master.propsPdf[PropRateShapeInv])

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropRateShapeInvJump] += 1
        else:
            self.rejects[PropRateShapeInvJump] += 1

    cdef void rateShapeInvAddPropose(self, unsigned mInd) except *:
        cdef Lik lik1
        cdef double rateShapeInv1, alpha1
        cdef double lnL1, lnPrior, lnJacob, lnProp, u, p

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Draw the inverse shape parameter from the prior distribution.
        rateShapeInv1 = -log(1.0 - genrand_res53(self.prng)) / \
          self.master._rateShapeInvPrior
        alpha1 = 1.0 / rateShapeInv1

        # Enable +G.
        lik1.setAlpha(mInd, alpha1)

        # Compute lnL with +G.
        lnL1 = lik1.lnL()

        lnPrior = -log(self.master._rateShapeInvJumpPrior) \
          + log(self.master._rateShapeInvPrior) + \
          (-self.master._rateShapeInvPrior*rateShapeInv1)

        lnJacob = -log(self.master._rateShapeInvPrior) \
          - (-self.master._rateShapeInvPrior*rateShapeInv1)
        lnProp = lnJacob

        if self.nModelsRatesGamma(lik1) == 1:
            # This proposal implicitly enables the PropRateShapeInv proposal.
            lnProp += log(1.0 - self.master.propsPdf[PropRateShapeInv])

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnProp)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropRateShapeInvJump] += 1
        else:
            self.rejects[PropRateShapeInvJump] += 1

    cdef bint rateShapeInvJumpPropose(self) except *:
        cdef unsigned mInd
        cdef double alpha

        assert self.master._ncat > 1

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, self.lik.nmodels())

        # Determine whether the model is currently +G.
        alpha = self.lik.getAlpha(mInd)

        if alpha != INFINITY:
            self.rateShapeInvRemovePropose(mInd, alpha)
        else:
            self.rateShapeInvAddPropose(mInd)

        return False

    cdef void invarRemovePropose(self, unsigned mInd, double wVar, \
      double wInvar) except *:
        cdef Lik lik1
        cdef double lnPrior, lnHast, lnL1, u, p

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        lnPrior = log(self.master._invarJumpPrior)
        lnHast = 0.0

        lnPrior += -log(1.0 - self.master._invarPrior) + \
          (1.0-self.master._invarPrior)*wVar
        lnHast += log(1.0 - self.master._invarPrior) - \
          self.master._invarPrior*wVar
        lik1.setWVar(mInd, 1.0)

        lnPrior += log(self.master._invarPrior) + wInvar
        lnHast += -log(self.master._invarPrior) - wInvar
        lik1.setWInvar(mInd, 0.0)

        if self.nModelsInvar(lik1) == 0:
            # This proposal implicitly disables the PropInvar proposal.
            lnHast += -log(1.0 - self.master.propsPdf[PropInvar])

        # Compute lnL with no invariable sites.
        lnL1 = lik1.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnHast)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropInvarJump] += 1
        else:
            self.rejects[PropInvarJump] += 1

    cdef void invarAddPropose(self, unsigned mInd) except *:
        cdef Lik lik1
        cdef double lnPrior, lnHast, wVar, wInvar, lnL1, u, p

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Draw weights from the prior and compute proposal ratio factors.
        lnPrior = -log(self.master._invarJumpPrior)
        lnHast = 0.0

        wVar = -log(1.0 - genrand_res53(self.prng)) * \
          (1.0 - self.master._invarPrior)
        lnPrior += -log(1.0-self.master._invarPrior) - wVar
        lnHast += log(1.0-self.master._invarPrior) + wVar
        lik1.setWVar(mInd, wVar)

        wInvar = -log(1.0 - genrand_res53(self.prng)) * self.master._invarPrior
        lnPrior += log(self.master._invarPrior) - wInvar
        lnHast += -log(self.master._invarPrior) + wInvar
        lik1.setWInvar(mInd, wInvar)

        if self.nModelsInvar(lik1) == 1:
            # This proposal implicitly enables the PropInvar proposal.
            lnHast += log(1.0 - self.master.propsPdf[PropInvar])

        # Compute lnL with invariable sites.
        lnL1 = lik1.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnHast)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropInvarJump] += 1
        else:
            self.rejects[PropInvarJump] += 1

    cdef bint invarJumpPropose(self) except *:
        cdef unsigned mInd
        cdef double wInvar

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, self.lik.nmodels())

        # Determine whether the model is currently +I.
        wInvar = self.lik.getWInvar(mInd)

        if wInvar != 0.0:
            self.invarRemovePropose(mInd, self.lik.getWVar(mInd), wInvar)
        else:
            self.invarAddPropose(mInd)

        return False

    cdef void freqEqualPropose(self, unsigned mInd) except *:
        cdef Lik lik1
        cdef unsigned i, nstates
        cdef double f, lnPrior, lnHast, lnL1, u, p

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Fix frequencies and compute proposal ratio factors.
        lnPrior = log(self.master._freqJumpPrior)
        lnHast = 0.0
        nstates = lik1.char_.nstates
        for 0 <= i < nstates:
            f = lik1.getFreq(mInd, i)
            lnPrior += f
            lnHast += -f
            lik1.setFreq(mInd, i, 1.0 / <double>nstates)

        if self.nModelsFreqsEstim(lik1) == 0:
            # This proposal implicitly disables the PropFreq proposal.
            lnHast += -log(1.0 - self.master.propsPdf[PropFreq])

        # Compute lnL with equal frequencies.
        lnL1 = lik1.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnHast)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropFreqJump] += 1
        else:
            self.rejects[PropFreqJump] += 1

    cdef void freqEstimPropose(self, unsigned mInd) except *:
        cdef Lik lik1
        cdef unsigned i
        cdef double f, lnPrior, lnHast, lnL1, u, p

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Draw frequencies from the prior and compute proposal ratio factors.
        lnPrior = -log(self.master._freqJumpPrior)
        lnHast = 0.0
        for 0 <= i < lik1.char_.nstates:
            f = -log(1.0 - genrand_res53(self.prng))
            lnPrior += -f
            lnHast += f
            lik1.setFreq(mInd, i, f)

        if self.nModelsFreqsEstim(lik1) == 1:
            # This proposal implicitly enables the PropFreq proposal.
            lnHast += log(1.0 - self.master.propsPdf[PropFreq])

        # Compute lnL with estimated frequencies.
        lnL1 = lik1.lnL()

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnHast)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropFreqJump] += 1
        else:
            self.rejects[PropFreqJump] += 1

    cdef bint freqJumpPropose(self) except *:
        cdef unsigned mInd

        # Uniformly choose a random model within the mixture.
        mInd = gen_rand64_range(self.prng, self.lik.nmodels())

        if self.freqsEqual(self.lik, mInd):
            self.freqEstimPropose(mInd)
        else:
            self.freqEqualPropose(mInd)

        return False

    cdef void mixtureRemovePropose(self, unsigned nmodels) except *:
        cdef Lik lik1
        cdef double w, rmult, lnL1, freq, rate, alpha, omega
        cdef double lnPrior, lnHast, u, p
        cdef unsigned mInd, nfreqs, nrates, i
        cdef list freqs, rates

        assert nmodels > 1

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Remove a random Q.
        mInd = gen_rand64_range(self.prng, nmodels)
        w = lik1.getWeight(mInd)
        rmult = lik1.getRmult(mInd)
        nfreqs = (1 if self.freqsEqual(lik1, mInd) else self.lik.char_.nstates)
        freqs = [lik1.getFreq(mInd, i) for i in xrange(nfreqs)]
        nrates = lik1.getNrates(mInd)
        rates = [lik1.getRate(mInd, i) for i in xrange(nrates)]
        alpha = lik1.getAlpha(mInd)
        lik1.delModel(mInd)

        # Compute lnL with removed Q.
        lnL1 = lik1.lnL()

        lnPrior = w
        lnHast = -w
        lnPrior += rmult
        lnHast += -rmult
        if nfreqs > 1:
            for 0 <= i < nfreqs:
                freq = <double>freqs[i]
                lnPrior += freq
                lnHast += -freq
        if nrates > 1:
            for 0 <= i < nrates:
                rate = <double>rates[i]
                lnPrior += rate
                lnHast += -rate
        if alpha != INFINITY:
            omega = 1.0 / alpha
            lnPrior += log(self.master._rateShapeInvPrior) + \
              (omega * self.master._rateShapeInvPrior)
            lnHast += -log(self.master._rateShapeInvPrior) - \
              (omega * self.master._rateShapeInvPrior)

        lnPrior += -log(1.0 - self.master._mixtureJumpPrior)

        if nmodels == 2:
            lnHast += log(2.0)
            # This proposal implicitly disables the PropWeight and PropRmult
            # proposals.
            lnHast += log(1.0 - self.master.propsPdf[PropWeight] - \
              self.master.propsPdf[PropRmult])
        #else:
        #    lnHast += log(1.0)

        if self.nModelsFreqsEstim(self.lik) == 1 and \
          self.nModelsFreqsEstim(lik1) == 0:
            # This proposal implicitly disables the PropFreq proposal.
            lnHast += -log(1.0 - self.master.propsPdf[PropFreq])
        if self.nModelsRatesEstim(self.lik) == 1 and \
          self.nModelsRatesEstim(lik1) == 0:
            # This proposal implicitly disables the PropRate proposal.
            lnHast += -log(1.0 - self.master.propsPdf[PropRate])
        if self.nModelsRatesGamma(self.lik) == 1 and \
          self.nModelsRatesGamma(lik1) == 0:
            # This proposal implicitly disables the PropRateShapeInv proposal.
            lnHast += -log(1.0 - self.master.propsPdf[PropRateShapeInv])

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnHast)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropMixtureJump] += 1
        else:
            self.rejects[PropMixtureJump] += 1

    cdef void mixtureAddPropose(self, unsigned nmodels) except *:
        cdef Lik lik1
        cdef double w, rmult, lnL1, freq, rate, alpha, omega
        cdef double lnPrior, lnHast, u, p
        cdef unsigned mInd, nrates, i
        cdef list freqs, rates

        # Create a cloned scratch Lik.
        lik1 = self.lik.clone()

        # Draw a new Q from the prior.
        mInd = lik1.addModel(1.0, self.master._ncat, self.master._catMedian, \
          self.master._invar)
        self.master.randomDnaQ(lik1, mInd, self.prng)
        w = lik1.getWeight(mInd)
        rmult = lik1.getRmult(mInd)
        nfreqs = (1 if self.freqsEqual(lik1, mInd) else self.lik.char_.nstates)
        freqs = [lik1.getFreq(mInd, i) for i in xrange(nfreqs)]
        nrates = lik1.getNrates(mInd)
        rates = [lik1.getRate(mInd, i) for i in xrange(nrates)]
        alpha = lik1.getAlpha(mInd)

        # Compute lnL with added Q.
        lnL1 = lik1.lnL()

        lnPrior = -w
        lnHast = w
        lnPrior += -rmult
        lnHast += rmult
        if nfreqs > 1:
            for 0 <= i < nfreqs:
                freq = <double>freqs[i]
                lnPrior += -freq
                lnHast += freq
        if nrates > 1:
            for 0 <= i < nrates:
                rate = <double>rates[i]
                lnPrior += -rate
                lnHast += rate
        if alpha != INFINITY:
            omega = 1.0 / alpha
            lnPrior += -log(self.master._rateShapeInvPrior) - \
              (omega * self.master._rateShapeInvPrior)
            lnHast += log(self.master._rateShapeInvPrior) + \
              (omega * self.master._rateShapeInvPrior)

        lnPrior += log(1.0 - self.master._mixtureJumpPrior)

        if nmodels == 1:
            lnHast += log(0.5)
            # This proposal implicitly enables the PropWeight and PropRmult
            # proposals.
            lnHast += log(1.0 - self.master.propsPdf[PropWeight] - \
              self.master.propsPdf[PropRmult])
        #else:
        #    lnHast += log(1.0)

        if self.nModelsFreqsEstim(self.lik) == 0 and \
          self.nModelsFreqsEstim(lik1) == 1:
            # This proposal implicitly enables the PropFreq proposal.
            lnHast += log(1.0 - self.master.propsPdf[PropFreq])
        if self.nModelsRatesEstim(self.lik) == 0 and \
          self.nModelsRatesEstim(lik1) == 1:
            # This proposal implicitly enables the PropRate proposal.
            lnHast += log(1.0 - self.master.propsPdf[PropRate])
        if self.nModelsRatesGamma(self.lik) == 0 and \
          self.nModelsRatesGamma(lik1) == 1:
            # This proposal implicitly enables the PropRateShapeInv proposal.
            lnHast += log(1.0 - self.master.propsPdf[PropRateShapeInv])

        # Determine whether to accept proposal.
        u = genrand_res53(self.prng)
        # p = [(likelihood ratio) * (prior ratio)]^heat * (Hastings ratio)
        p = exp(((lnL1 - self.lnL) + lnPrior) * self.heat + lnHast)
        if p >= u:
            # Accept.
            self.lnL = lnL1
            self.lik = lik1
            self.accepts[PropMixtureJump] += 1
        else:
            self.rejects[PropMixtureJump] += 1

    cdef bint mixtureJumpPropose(self) except *:
        cdef unsigned nmodels
        cdef bint remove

        nmodels = self.lik.nmodels()
        if nmodels == 1:
            # Can't remove a model.
            remove = False
        else:
            remove = <bint>gen_rand64_range(self.prng, 2)

        if remove:
            self.mixtureRemovePropose(nmodels)
        else:
            self.mixtureAddPropose(nmodels)

        return False

    cdef void advance0(self) except *:
        cdef bint again
        cdef unsigned propInd, a, b, other
        cdef double u

        self.step += 1

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
            elif propInd == PropInvar:
                again = self.invarPropose()
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
            elif propInd == PropInvarJump:
                again = self.invarJumpPropose()
            elif propInd == PropFreqJump:
                again = self.freqJumpPropose()
            elif propInd == PropMixtureJump:
                again = self.mixtureJumpPropose()
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

    cdef void advance1(self) except *:
        cdef double rcvHeat, rcvLnL, p

        # Finish handling a pending potential heat swap.
        if self.swapInd != self.ind:
            self.master.recvSwapInfo(self.run, self.ind, self.swapInd, \
              self.step, &rcvHeat, &rcvLnL)
            assert rcvHeat != self.heat
            p = exp((rcvLnL - self.lnL) * self.heat \
              + (self.lnL - rcvLnL) * rcvHeat)
            if p >= self.swapProb:
                self.heat = rcvHeat
                self.nswap += 1
            self.swapInd = self.ind
