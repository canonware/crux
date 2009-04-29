"""
    Posterior distribution statistics for Mc3.

    Post reads the log files generated by Mc3 and computes various statistics
    on their contents.
"""

import re
import sys

from libc cimport *
from libm cimport *
from Crux.Mc3 cimport Mc3
from Crux.CTMatrix cimport Alignment
from Crux.Tree cimport Tree

cdef class Msamp:
    def __init__(self, double weight, double rmult, str rclass, list rates, \
      double alpha, list freqs):
        self.weight = weight
        self.rmult = rmult
        self.rclass = rclass
        self.rates = rates
        self.alpha = alpha
        self.freqs = freqs

cdef class Samp:
    def __init__(self, double lnL):
        self.lnL = lnL
        self.msamps = []
        self.wNorm = -1.0
        self.tree = None

cdef class Post:
    """
        alignment
          Aligned character-by-taxon matrix input data.

        outPrefix
          Filename prefix, used for naming output files.

        burnin
          Number of samples to discard, including sample 0.  By default, the
          first half of each chain is discarded as burn-in.

        verbose
          If True, print diagnostic output.
    """
    def __init__(self, Alignment alignment, str outPrefix, \
      uint64_t burnin=ULLONG_MAX, bint verbose=False):
        self.mc3 = Mc3(alignment, outPrefix)
        self.verbose = verbose

        self.nsamples = 0
        self.burnin = burnin
        self.stepFirst = 0
        self.stepLast = 0
        self._sDone = False
        self._pDone = False
        self._tDone = False
        self.maxModels = 0

        self._rcovLnL = -1.0
        self._rcovRclass = -1.0

        self._lnLMinCred = -1.0
        self._lnLMaxCred = -1.0
        self._lnLMean = -1.0
        self._lnLVar = -1.0
        self._lnLMed = -1.0

        self._nmodelsMinCred = -1.0
        self._nmodelsMaxCred = -1.0
        self._nmodelsMean = -1.0
        self._nmodelsVar = -1.0
        self._nmodelsMed = -1.0

        self._ratesMinCred = None
        self._ratesMaxCred = None
        self._ratesMean = None
        self._ratesVar = None
        self._ratesMed = None

        self._freqsMinCred = None
        self._freqsMaxCred = None
        self._freqsMean = None
        self._freqsVar = None
        self._freqsMed = None

        self._alphaMinCred = -1.0
        self._alphaMaxCred = -1.0
        self._alphaMean = -1.0
        self._alphaVar = -1.0
        self._alphaMed = -1.0

        self._sumt = None

        self.runs = None

        self._parseL()

    cdef void _parseL(self) except *:
        cdef file lFile
        cdef bint inParms
        cdef str line, k, v
        cdef object m

        lFile = open("%s.l" % self.mc3.outPrefix, "r")
        inParms = False
        for line in lFile:
            if not inParms:
                m = re.compile(r'^PRNG seed: ([0-9]+)').match(line)
                if m:
                    self.seed = int(m.group(1))
                    if self.verbose:
                        sys.stdout.write("PRNG seed: %d\n" % self.seed)
                    continue
                m = re.compile(r'^Configuration parameters:$').match(line)
                if m:
                    inParms = True
                    if self.verbose:
                        sys.stdout.write("Configuration parameters:\n")
            else:
                m = re.compile(r'^  ([A-Za-z]+):\s+(.*)$').match(line)
                if not m:
                    return
                k = m.group(1)
                v = m.group(2)
                if k == 'outPrefix':
                    pass
                elif k == 'graphDelay':
                    self.mc3.setGraphDelay(float(v))
                elif k == 'emaAlpha':
                    self.mc3.setEmaAlpha(float(v))
                elif k == 'cvgSampStride':
                    self.mc3.setCvgSampStride(long(v))
                elif k == 'cvgAlpha':
                    self.mc3.setCvgAlpha(float(v))
                elif k == 'cvgEpsilon':
                    self.mc3.setCvgEpsilon(float(v))
                elif k == 'minStep':
                    self.mc3.setMinStep(long(v))
                elif k == 'maxStep':
                    self.mc3.setMaxStep(long(v))
                elif k == 'stride':
                    self.mc3.setStride(long(v))
                elif k == 'nruns':
                    self.mc3.setNruns(long(v))
                elif k == 'ncoupled':
                    self.mc3.setNcoupled(long(v))
                elif k == 'heatDelta':
                    self.mc3.setHeatDelta(float(v))
                elif k == 'swapStride':
                    self.mc3.setSwapStride(long(v))
                elif k == 'nmodels':
                    self.mc3.setNmodels(long(v))
                elif k == 'ncat':
                    self.mc3.setNcat(long(v))
                elif k == 'catMedian':
                    self.mc3.setCatMedian(v == 'True')
                elif k == 'weightLambda':
                    self.mc3.setWeightLambda(float(v))
                elif k == 'freqLambda':
                    self.mc3.setFreqLambda(float(v))
                elif k == 'rmultLambda':
                    self.mc3.setRmultLambda(float(v))
                elif k == 'rateLambda':
                    self.mc3.setRateLambda(float(v))
                elif k == 'rateShapeInvLambda':
                    self.mc3.setRateShapeInvLambda(float(v))
                elif k == 'brlenLambda':
                    self.mc3.setBrlenLambda(float(v))
                elif k == 'etbrPExt':
                    self.mc3.setEtbrPExt(float(v))
                elif k == 'etbrLambda':
                    self.mc3.setEtbrLambda(float(v))
                elif k == 'rateShapeInvPrior':
                    self.mc3.setRateShapeInvPrior(float(v))
                elif k == 'brlenPrior':
                    self.mc3.setBrlenPrior(float(v))
                elif k == 'rateJumpPrior':
                    self.mc3.setRateJumpPrior(float(v))
                elif k == 'polytomyJumpPrior':
                    self.mc3.setPolytomyJumpPrior(float(v))
                elif k == 'rateShapeInvJumpPrior':
                    self.mc3.setRateShapeInvJumpPrior(float(v))
                elif k == 'freqJumpPrior':
                    self.mc3.setFreqJumpPrior(float(v))
                elif k == 'mixtureJumpPrior':
                    self.mc3.setMixtureJumpPrior(float(v))
                elif k == 'weightProp':
                    self.mc3.setWeightProp(float(v))
                elif k == 'freqProp':
                    self.mc3.setFreqProp(float(v))
                elif k == 'rmultProp':
                    self.mc3.setRmultProp(float(v))
                elif k == 'rateProp':
                    self.mc3.setRateProp(float(v))
                elif k == 'rateShapeInvProp':
                    self.mc3.setRateShapeInvProp(float(v))
                elif k == 'brlenProp':
                    self.mc3.setBrlenProp(float(v))
                elif k == 'etbrProp':
                    self.mc3.setEtbrProp(float(v))
                elif k == 'rateJumpProp':
                    self.mc3.setRateJumpProp(float(v))
                elif k == 'polytomyJumpProp':
                    self.mc3.setPolytomyJumpProp(float(v))
                elif k == 'rateShapeInvJumpProp':
                    self.mc3.setRateShapeInvJumpProp(float(v))
                elif k == 'freqJumpProp':
                    self.mc3.setFreqJumpProp(float(v))
                elif k == 'mixtureJumpProp':
                    self.mc3.setMixtureJumpProp(float(v))
                else:
                    assert False

                if self.verbose:
                    sys.stdout.write("  %s: %s\n" % (k, v))

    cpdef parseS(self):
        """
            Parse <self.outPrefix>.s in order to acquire lnL sample data.

            This method is called implicitly as necessary when computing
            various statistics, but it must be called manually prior to
            directly accessing Samp.lnL (via Cython code), in order to assure
            the field's validity.
        """
        cdef file sFile
        cdef int nruns
        cdef uint64_t stride, nlines, i, j
        cdef str line
        cdef list toks
        cdef Samp samp

        if self._sDone:
            return

        sFile = open("%s.s" % self.mc3.outPrefix, "r")

        nruns = self.mc3.getNruns()
        self.runs = [[] for i in xrange(nruns)]
        stride = self.mc3.getStride()

        # Get the number of lines and the step number for the last sample.
        nlines = 0
        for line in sFile:
            nlines += 1
        self.stepLast = long(line.split(None, 1)[0])

        if self.burnin == ULLONG_MAX:
            # Use the last half of each chain by default.
            self.stepFirst = ((self.stepLast/stride/2)+1)*stride
        else:
            if self.burnin > self.stepLast / stride:
                raise ValueError("Burnin (%d) consumes all samples (%d)" % \
                  (self.burnin, self.stepLast/stride + 1))
            self.stepFirst = self.burnin * stride
        self.nsamples = (self.stepLast - self.stepFirst + stride) / stride

        # Rewind to beginning of file, then skip the header and burnin.
        sFile.seek(0)
        for 0 <= i < nlines - self.nsamples:
            sFile.readline()

        # Parse posterior samples.
        for i <= i < nlines:
            line = sFile.readline()
            toks = line.split(None, nruns+2)
            for 0 <= j < nruns:
                samp = Samp(float(toks[j+2]))
                self.runs[j].append(samp)

        self._sDone = True

    cpdef parseP(self):
        """
            Parse <self.outPrefix>.p in order to acquire model parameter sample
            data.

            This method is called implicitly as necessary when computing
            various statistics, but it must be called manually prior to
            directly accessing Samp.msamps (via Cython code), in order to
            assure the field's validity.
        """
        cdef file pFile
        cdef unsigned nfreqs, nrates, run, model
        cdef str line, rclass
        cdef list toks, rates, freqs
        cdef uint64_t stride, step
        cdef double weight, rmult, wNorm, alpha
        cdef Samp samp
        cdef Msamp msamp

        if self._pDone:
            return
        self.parseS()

        pFile = open("%s.p" % self.mc3.outPrefix, "r")

        nfreqs = self.mc3.alignment.charType.get().nstates
        nrates = (nfreqs*(nfreqs-1))/2
        stride = self.mc3.getStride()

        # Skip header.
        pFile.readline()

        for line in pFile:
            toks = line.split()
            step = long(toks[1])
            if step >= self.stepFirst:
                assert step <= self.stepLast
                run = int(toks[0])
                model = int(toks[2])
                weight = float(toks[3])
                rmult = float(toks[4])
                rclass = toks[6]
                rates = []
                for 0 <= i < nrates:
                    rates.append(float(toks[i+8]))
                wNorm = float(toks[nrates+9])
                alpha = float(toks[nrates+10])
                freqs = []
                for 0 <= i < nfreqs:
                    freqs.append(float(toks[i+nrates+12]))
                msamp = Msamp(weight, rmult, rclass, rates, alpha, freqs)
                samp = \
                  <Samp>(<list>self.runs[run])[(step-self.stepFirst)/stride]
                samp.msamps.append(msamp)
                assert samp.wNorm == -1.0 or samp.wNorm == wNorm
                samp.wNorm = wNorm
                if len(samp.msamps) > self.maxModels:
                    self.maxModels = len(samp.msamps)

        self._pDone = True

    cpdef parseT(self):
        """
            Parse <self.outPrefix>.t in order to acquire tree sample data.

            This method is called implicitly as necessary when computing
            various statistics, but it must be called manually prior to
            directly accessing Samp.tree (via Cython code), in order to assure
            the field's validity.
        """
        cdef file tFile
        cdef uint64_t stride, step
        cdef str line
        cdef list toks
        cdef unsigned run
        cdef Tree tree
        cdef Samp samp

        if self._tDone:
            return
        self.parseS()

        tFile = open("%s.t" % self.mc3.outPrefix, "r")

        stride = self.mc3.getStride()

        # Skip header.
        tFile.readline()

        for line in tFile:
            toks = line.split(None, 2)
            step = long(toks[1][:-1])
            if step >= self.stepFirst:
                assert step <= self.stepLast
                run = int(toks[0][1:])
                tree = Tree(toks[2], None, False)
                tree.deroot()
                samp = \
                  <Samp>(<list>self.runs[run])[(step-self.stepFirst)/stride]
                samp.tree = tree

        self._pDone = True

    # This method is conceptually identical to Crux.Mc3.Mc3.computeRcov(), but
    # it does not have to deal with discarding burn-in, nor does speed matter
    # as much.
    cdef void _computeRcovLnL(self) except *:
        cdef double rcovLnL, cvgAlpha, minLnL, maxLnL, cvgAlphaEmp
        cdef uint64_t lower, upper, j, k, n
        cdef unsigned nruns, i
        cdef list lnLs
        cdef Samp samp

        self.parseS()

        nruns = self.mc3.getNruns()
        if nruns < 2:
            raise ValueError("Rcov computation requires at least 2 runs")

        if self.nsamples == 0:
            raise ValueError("No samples")

        rcovLnL = 0.0

        # Compute the lower and upper bounds for the credibility intervals.
        cvgAlpha = self.mc3.getCvgAlpha()
        lower = <uint64_t>floor(cvgAlpha/2.0 * <double>self.nsamples)
        upper = self.nsamples - 1 - lower
        print "[%d..%d] (%d)" % (lower, upper, self.nsamples)

        for 0 <= i < nruns:
            # Accumulate sorted lnLs for run i, in order to make the bounds of
            # the credibility interval available.
            lnLs = []
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                lnLs.append(samp.lnL)
            lnLs.sort()
            minLnL = <double>lnLs[lower]
            maxLnL = <double>lnLs[upper]

            # Compute the proportion of the lnLs that fall within the
            # credibility interval.
            n = 0
            for 0 <= j < nruns:
                for 0 <= k < self.nsamples:
                    samp = <Samp>(<list>self.runs[j])[k]
                    if samp.lnL >= minLnL and samp.lnL <= maxLnL:
                        n += 1
            rcovLnL += <double>n / <double>(nruns * self.nsamples)
        rcovLnL /= <double>nruns
        # Compute the empirical alpha, based on the credibility interval
        # proportion, and use this to rescale Rcov, so that results can be
        # interpreted in the context of the intended alpha.
        cvgAlphaEmp = <double>(lower*2) / <double>self.nsamples
        rcovLnL *= (1.0-cvgAlpha) / (1.0-cvgAlphaEmp)
        self._rcovLnL = rcovLnL

    cdef double getRcovLnL(self) except -1.0:
        if self._rcovLnL == -1.0:
            self._computeRcovLnL()
        return self._rcovLnL
    property rcovLnL:
        """
            Rcov(lnL) coverage ratio, based on the non-burn-in samples from
            multiple independent runs.  For more details on Rcov, see the
            Crux.Mc3 documentation.
        """
        def __get__(self):
            return self.getRcovLnL()

    cdef void _computeRcovRclass(self) except *:
        cdef double rcovRclass, cvgAlpha, cvgAlphaEmp, r
        cdef unsigned nruns, i
        cdef uint64_t tail, j, k, n, nCum
        cdef dict rclasses
        cdef Samp samp
        cdef Msamp msamp
        cdef str rclass
        cdef list deco

        self.parseP()
        if self.maxModels != 1:
            raise ValueError("Rcov(rclass) only applies to nmodels=1")

        nruns = self.mc3.getNruns()
        if nruns < 2:
            raise ValueError("Rcov computation requires at least 2 runs")

        if self.nsamples == 0:
            raise ValueError("No samples")

        rcovRclass = 0.0

        # Compute length of tail.
        cvgAlpha = self.mc3.getCvgAlpha()
        tail = <uint64_t>floor(cvgAlpha * <double>self.nsamples)

        for 0 <= i < nruns:
            # Create a rclass-->count mapping of rclasses.
            rclasses = {}
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                msamp = <Msamp>samp.msamps[0]
                rclass = msamp.rclass
                if rclass in rclasses:
                    rclasses[rclass] += 1
                else:
                    rclasses[rclass] = 1

            deco = [(rclasses[rclass], rclass) for rclass in rclasses]
            deco.sort()

            # Find the boundary between the tail and the credible set of
            # rclasses.  If an rclass falls partly within the credible set,
            # include it (which increases cvgAlphaEmp).  If multiple rclasses
            # tie for inclusion in the credible set, include all of them (which
            # increases cvgAlphaEmp).
            nCum = 0
            for 0 <= j < len(deco):
                (n, rclass) = deco[j]
                if nCum + n > tail:
                    break
                nCum += n
            for j >= j > 0:
                if deco[j-1][0] != deco[j][0]:
                    # Not a tie.
                    break
            # Remove tail from rclasses.
            nCum = 0
            for 0 <= k < j:
                (n, rclass) = deco[k]
                nCum += n
                del rclasses[rclass]
            cvgAlphaEmp = <double>nCum / <double>self.nsamples

            # Compute the proportion of the rclasses that fall within the
            # credibility interval.
            n = 0
            for 0 <= j < nruns:
                for 0 <= k < self.nsamples:
                    samp = <Samp>(<list>self.runs[j])[k]
                    msamp = <Msamp>samp.msamps[0]
                    rclass = msamp.rclass
                    if rclass in rclasses:
                        n += 1

            r = <double>n / <double>(nruns * self.nsamples)
            # Rescale so that results can be interpreted in the context of the
            # intended alpha.
            r *= (1.0-cvgAlpha) / (1.0-cvgAlphaEmp)
            rcovRclass += r
        rcovRclass /= <double>nruns
        self._rcovRclass = rcovRclass

    cdef double getRcovRclass(self) except -1.0:
        if self._rcovRclass == -1.0:
            self._computeRcovRclass()
        return self._rcovRclass
    property rcovRclass:
        """
            Rcov(rclass) coverage ratio, based on non-burn-in frequency-ordered
            samples from multiple independent runs.  This statistic can only be
            computed if nmodels=1.  For more details on Rcov, see the Crux.Mc3
            documentation.
        """
        def __get__(self):
            return self.getRcovRclass()

    cdef void _summarizeLnL(self) except *:
        cdef list lnLs
        cdef unsigned nruns
        cdef uint64_t i, j, N, lower, upper
        cdef Samp samp
        cdef double cvgAlpha, lnL, lnLMean, diff, lnLVar

        self.parseS()

        if self.nsamples == 0:
            raise ValueError("No samples")

        lnLs = []
        nruns = self.mc3.getNruns()
        for 0 <= i < nruns:
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                lnLs.append(samp.lnL)
        lnLs.sort()
        N = len(lnLs)
        # Compute the lower and upper bounds for the credibility interval.
        cvgAlpha = self.mc3.getCvgAlpha()
        lower = <uint64_t>floor(cvgAlpha/2.0 * N)
        upper = N - 1 - lower

        self._lnLMinCred = <double>lnLs[lower]
        self._lnLMaxCred = <double>lnLs[upper]

        lnLMean = 0.0
        for 0 <= i < N:
            lnL = <double>lnLs[i]
            lnLMean += lnL / <double>N
        self._lnLMean = lnLMean

        lnLVar = 0.0
        for 0 <= i < N:
            lnL = <double>lnLs[i]
            diff = lnL - lnLMean
            lnLVar += (diff * diff) / <double>N
        self._lnLVar = lnLVar

        if N % 2 == 1:
            self._lnLMed = <double>lnLs[N/2]
        else:
            # Average the middle two samples.
            self._lnLMed = (<double>lnLs[N/2 - 1] + <double>lnLs[N/2]) / 2.0

    cdef double getLnLMinCred(self) except -1.0:
        if self._lnLMinCred == -1.0:
            self._summarizeLnL()
        return self._lnLMinCred
    property lnLMinCred:
        """
            Minimum lnL within the credibility interval, as defined by the
            cvgAlpha parameter used during sampling.
        """
        def __get__(self):
            return self.getLnLMinCred()

    cdef double getLnLMaxCred(self) except -1.0:
        if self._lnLMaxCred == -1.0:
            self._summarizeLnL()
        return self._lnLMaxCred
    property lnLMaxCred:
        """
            Maximum lnL within the credibility interval, as defined by the
            cvgAlpha parameter used during sampling.
        """
        def __get__(self):
            return self.getLnLMaxCred()

    cdef double getLnLMean(self) except -1.0:
        if self._lnLMean == -1.0:
            self._summarizeLnL()
        return self._lnLMean
    property lnLMean:
        """
            lnL sample mean.
        """
        def __get__(self):
            return self.getLnLMean()

    cdef double getLnLVar(self) except -1.0:
        if self._lnLVar == -1.0:
            self._summarizeLnL()
        return self._lnLVar
    property lnLVar:
        """
            lnL sample variance.
        """
        def __get__(self):
            return self.getLnLVar()

    cdef double getLnLMed(self) except -1.0:
        if self._lnLMed == -1.0:
            self._summarizeLnL()
        return self._lnLMed
    property lnLMed:
        """
            lnL sample median.
        """
        def __get__(self):
            return self.getLnLMed()

    cdef void _summarizeNmodels(self) except *:
        cdef list modelCnts
        cdef unsigned nruns
        cdef uint64_t i, j, N, lower, upper
        cdef Samp samp
        cdef double cvgAlpha, nmodels, nmodelsMean, diff, nmodelsVar

        self.parseP()

        if self.nsamples == 0:
            raise ValueError("No samples")

        modelCnts = []
        nruns = self.mc3.getNruns()
        for 0 <= i < nruns:
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                modelCnts.append(<double>len(samp.msamps))
        modelCnts.sort()
        N = len(modelCnts)
        # Compute the lower and upper bounds for the credibility interval.
        cvgAlpha = self.mc3.getCvgAlpha()
        lower = <uint64_t>floor(cvgAlpha/2.0 * N)
        upper = N - 1 - lower

        self._nmodelsMinCred = <double>modelCnts[lower]
        self._nmodelsMaxCred = <double>modelCnts[upper]

        nmodelsMean = 0.0
        for 0 <= i < N:
            nmodels = <double>modelCnts[i]
            nmodelsMean += nmodels / <double>N
        self._nmodelsMean = nmodelsMean

        nmodelsVar = 0.0
        for 0 <= i < N:
            nmodels = <double>modelCnts[i]
            diff = nmodels - nmodelsMean
            nmodelsVar += (diff * diff) / <double>N
        self._nmodelsVar = nmodelsVar

        if N % 2 == 1:
            self._nmodelsMed = <double>modelCnts[N/2]
        else:
            # Average the middle two samples.
            self._nmodelsMed = (<double>modelCnts[N/2 - 1] + \
              <double>modelCnts[N/2]) / 2.0

    cdef double getNmodelsMinCred(self) except -1.0:
        if self._nmodelsMinCred == -1.0:
            self._summarizeNmodels()
        return self._nmodelsMinCred
    property nmodelsMinCred:
        """
            Minimum nmodels within the credibility interval, as defined by the
            cvgAlpha parameter used during sampling.
        """
        def __get__(self):
            return self.getNmodelsMinCred()

    cdef double getNmodelsMaxCred(self) except -1.0:
        if self._nmodelsMaxCred == -1.0:
            self._summarizeNmodels()
        return self._nmodelsMaxCred
    property nmodelsMaxCred:
        """
            Maximum nmodels within the credibility interval, as defined by the
            cvgAlpha parameter used during sampling.
        """
        def __get__(self):
            return self.getNmodelsMaxCred()

    cdef double getNmodelsMean(self) except -1.0:
        if self._nmodelsMean == -1.0:
            self._summarizeNmodels()
        return self._nmodelsMean
    property nmodelsMean:
        """
            nmodels sample mean.
        """
        def __get__(self):
            return self.getNmodelsMean()

    cdef double getNmodelsVar(self) except -1.0:
        if self._nmodelsVar == -1.0:
            self._summarizeNmodels()
        return self._nmodelsVar
    property nmodelsVar:
        """
            nmodels sample variance.
        """
        def __get__(self):
            return self.getNmodelsVar()

    cdef double getNmodelsMed(self) except -1.0:
        if self._nmodelsMed == -1.0:
            self._summarizeNmodels()
        return self._nmodelsMed
    property nmodelsMed:
        """
            nmodels sample median.
        """
        def __get__(self):
            return self.getNmodelsMed()

    cdef void _summarizeRates(self) except *:
        cdef list rates, ratesK
        cdef Samp samp
        cdef Msamp msamp
        cdef unsigned rlen, nruns, k
        cdef uint64_t i, j, N, lower, upper
        cdef double cvgAlpha, rate, diff, rateMean, rateVar

        if self.mc3.nmodels != 1 or self.mc3.rateJumpProp != 0.0:
            raise ValueError( \
              "Rates cannot be summarized for chosen model parameters")

        self.parseP()

        if self.nsamples == 0:
            raise ValueError("No samples")

        samp = <Samp>(<list>self.runs[0])[0]
        msamp = <Msamp>samp.msamps[0]
        rlen = len(msamp.rates)
        rates = [[] for i in xrange(rlen)]

        nruns = self.mc3.getNruns()
        for 0 <= i < nruns:
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                msamp = <Msamp>samp.msamps[0]
                for 0 <= k < rlen:
                    ratesK = <list>rates[k]
                    ratesK.append(<double>msamp.rates[k] * samp.wNorm * \
                      msamp.rmult)
        for 0 <= k < rlen:
            ratesK = <list>rates[k]
            ratesK.sort()
        N = len(<list>rates[0])
        # Compute the lower and upper bounds for the credibility interval.
        cvgAlpha = self.mc3.getCvgAlpha()
        lower = <uint64_t>floor(cvgAlpha/2.0 * N)
        upper = N - 1 - lower

        self._ratesMinCred = []
        self._ratesMaxCred = []
        self._ratesMean = []
        self._ratesVar = []
        self._ratesMed = []
        for 0 <= k < rlen:
            ratesK = <list>rates[k]

            self._ratesMinCred.append(<double>ratesK[lower])
            self._ratesMaxCred.append(<double>ratesK[upper])

            rateMean = 0.0
            for 0 <= i < N:
                rate = <double>ratesK[i]
                rateMean += rate / <double>N
            self._ratesMean.append(rateMean)

            rateVar = 0.0
            for 0 <= i < N:
                rate = <double>ratesK[i]
                diff = rate - rateMean
                rateVar += (diff * diff) / <double>N
            self._ratesVar.append(rateVar)

            if N % 2 == 1:
                self._rateMed.append(<double>ratesK[N/2])
            else:
                # Average the middle two samples.
                self._ratesMed.append((<double>ratesK[N/2 - 1] + \
                  <double>ratesK[N/2]) / 2.0)
    cdef list getRatesMinCred(self):
        if self._ratesMinCred is None:
            self._summarizeRates()
        return self._ratesMinCred
    property ratesMinCred:
        """
            List of minimum rates within their respective credibility
            intervals, as defined by the cvgAlpha parameter used during
            sampling.
        """
        def __get__(self):
            return self.getRatesMinCred()
    cdef list getRatesMaxCred(self):
        if self._ratesMaxCred is None:
            self._summarizeRates()
        return self._ratesMaxCred
    property ratesMaxCred:
        """
            List of maximum rates within their respective credibility
            intervals, as defined by the cvgAlpha parameter used during
            sampling.
        """
        def __get__(self):
            return self.getRatesMaxCred()
    cdef list getRatesMean(self):
        if self._ratesMean is None:
            self._summarizeRates()
        return self._ratesMean
    property ratesMean:
        """
            List of rate sample means.
        """
        def __get__(self):
            return self.getRatesMean()
    cdef list getRatesVar(self):
        if self._ratesVar is None:
            self._summarizeRates()
        return self._ratesVar
    property ratesVar:
        """
            List of rate sample variances.
        """
        def __get__(self):
            return self.getRatesVar()
    cdef list getRatesMed(self):
        if self._ratesMed is None:
            self._summarizeRates()
        return self._ratesMed
    property ratesMed:
        """
            List of rate sample medians.
        """
        def __get__(self):
            return self.getRatesMed()

    cdef void _summarizeFreqs(self) except *:
        cdef list freqs, freqsK
        cdef Samp samp
        cdef Msamp msamp
        cdef unsigned flen, nruns, k
        cdef uint64_t i, j, N, lower, upper
        cdef double cvgAlpha, freq, diff, freqMean, freqVar

        if self.mc3.nmodels != 1:
            raise ValueError( \
              "Freqs cannot be summarized for chosen model parameters")

        self.parseP()

        if self.nsamples == 0:
            raise ValueError("No samples")

        samp = <Samp>(<list>self.runs[0])[0]
        msamp = <Msamp>samp.msamps[0]
        flen = len(msamp.freqs)
        freqs = [[] for i in xrange(flen)]

        nruns = self.mc3.getNruns()
        for 0 <= i < nruns:
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                msamp = <Msamp>samp.msamps[0]
                for 0 <= k < flen:
                    freqsK = <list>freqs[k]
                    freqsK.append(<double>msamp.freqs[k])
        for 0 <= k < flen:
            freqsK = <list>freqs[k]
            freqsK.sort()
        N = len(<list>freqs[0])
        # Compute the lower and upper bounds for the credibility interval.
        cvgAlpha = self.mc3.getCvgAlpha()
        lower = <uint64_t>floor(cvgAlpha/2.0 * N)
        upper = N - 1 - lower

        self._freqsMinCred = []
        self._freqsMaxCred = []
        self._freqsMean = []
        self._freqsVar = []
        self._freqsMed = []
        for 0 <= k < flen:
            freqsK = <list>freqs[k]

            self._freqsMinCred.append(<double>freqsK[lower])
            self._freqsMaxCred.append(<double>freqsK[upper])

            freqMean = 0.0
            for 0 <= i < N:
                freq = <double>freqsK[i]
                freqMean += freq / <double>N
            self._freqsMean.append(freqMean)

            freqVar = 0.0
            for 0 <= i < N:
                freq = <double>freqsK[i]
                diff = freq - freqMean
                freqVar += (diff * diff) / <double>N
            self._freqsVar.append(freqVar)

            if N % 2 == 1:
                self._freqMed.append(<double>freqsK[N/2])
            else:
                # Average the middle two samples.
                self._freqsMed.append((<double>freqsK[N/2 - 1] + \
                  <double>freqsK[N/2]) / 2.0)
    cdef list getFreqsMinCred(self):
        if self._freqsMinCred is None:
            self._summarizeFreqs()
        return self._freqsMinCred
    property freqsMinCred:
        """
            List of minimum frequencies within their respective credibility
            intervals, as defined by the cvgAlpha parameter used during
            sampling.
        """
        def __get__(self):
            return self.getFreqsMinCred()
    cdef list getFreqsMaxCred(self):
        if self._freqsMaxCred is None:
            self._summarizeFreqs()
        return self._freqsMaxCred
    property freqsMaxCred:
        """
            List of maximum frequencies within their respective credibility
            intervals, as defined by the cvgAlpha parameter used during
            sampling.
        """
        def __get__(self):
            return self.getFreqsMaxCred()
    cdef list getFreqsMean(self):
        if self._freqsMean is None:
            self._summarizeFreqs()
        return self._freqsMean
    property freqsMean:
        """
            List of frequency sample means.
        """
        def __get__(self):
            return self.getFreqsMean()
    cdef list getFreqsVar(self):
        if self._freqsVar is None:
            self._summarizeFreqs()
        return self._freqsVar
    property freqsVar:
        """
            List of frequency sample variances.
        """
        def __get__(self):
            return self.getFreqsVar()
    cdef list getFreqsMed(self):
        if self._freqsMed is None:
            self._summarizeFreqs()
        return self._freqsMed
    property freqsMed:
        """
            List of frequency sample medians.
        """
        def __get__(self):
            return self.getFreqsMed()

    cdef void _summarizeAlpha(self) except *:
        cdef list alphas
        cdef unsigned nruns
        cdef uint64_t i, j, N, lower, upper
        cdef Samp samp
        cdef Msamp msamp
        cdef double cvgAlpha, alpha, diff, alphaMean, alphaVar

        if self.mc3.nmodels != 1 or self.mc3.ncat <= 1:
            raise ValueError( \
              "+G alpha cannot be summarized for chosen model parameters")

        self.parseP()

        if self.nsamples == 0:
            raise ValueError("No samples")

        alphas = []
        nruns = self.mc3.getNruns()
        for 0 <= i < nruns:
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                msamp = <Msamp>samp.msamps[0]
                alphas.append(msamp.alpha)
        alphas.sort()
        N = len(alphas)
        # Compute the lower and upper bounds for the credibility interval.
        cvgAlpha = self.mc3.getCvgAlpha()
        lower = <uint64_t>floor(cvgAlpha/2.0 * N)
        upper = N - 1 - lower

        self._alphaMinCred = <double>alphas[lower]
        self._alphaMaxCred = <double>alphas[upper]

        alphaMean = 0.0
        for 0 <= i < N:
            alpha = <double>alphas[i]
            alphaMean += alpha / <double>N
        self._alphaMean = alphaMean

        alphaVar = 0.0
        for 0 <= i < N:
            alpha = <double>alphas[i]
            diff = alpha - alphaMean
            alphaVar += (diff * diff) / <double>N
        self._alphaVar = alphaVar

        if N % 2 == 1:
            self._alphaMed = <double>alphas[N/2]
        else:
            # Average the middle two samples.
            self._alphaMed = (<double>alphas[N/2 - 1] + \
              <double>alphas[N/2]) / 2.0

    cdef double getAlphaMinCred(self) except -1.0:
        if self._alphaMinCred == -1.0:
            self._summarizeAlpha()
        return self._alphaMinCred
    property alphaMinCred:
        """
            Minimum alpha within the credibility interval, as defined by the
            cvgAlpha parameter used during sampling.
        """
        def __get__(self):
            return self.getAlphaMinCred()

    cdef double getAlphaMaxCred(self) except -1.0:
        if self._alphaMaxCred == -1.0:
            self._summarizeAlpha()
        return self._alphaMaxCred
    property alphaMaxCred:
        """
            Maximum alpha within the credibility interval, as defined by the
            cvgAlpha parameter used during sampling.
        """
        def __get__(self):
            return self.getAlphaMaxCred()

    cdef double getAlphaMean(self) except -1.0:
        if self._alphaMean == -1.0:
            self._summarizeAlpha()
        return self._alphaMean
    property alphaMean:
        """
            alpha sample mean.
        """
        def __get__(self):
            return self.getAlphaMean()

    cdef double getAlphaVar(self) except -1.0:
        if self._alphaVar == -1.0:
            self._summarizeAlpha()
        return self._alphaVar
    property alphaVar:
        """
            alpha sample variance.
        """
        def __get__(self):
            return self.getAlphaVar()

    cdef double getAlphaMed(self) except -1.0:
        if self._alphaMed == -1.0:
            self._summarizeAlpha()
        return self._alphaMed
    property alphaMed:
        """
            alpha sample median.
        """
        def __get__(self):
            return self.getAlphaMed()

    cdef void _summarizeTrees(self) except *:
        cdef list treeLists, trees
        cdef unsigned nruns, i
        cdef uint64_t j
        cdef Samp samp
        cdef Tree tree

        self.parseT()

        if self.nsamples == 0:
            raise ValueError("No samples")

        treeLists = []
        nruns = self.mc3.getNruns()
        for 0 <= i < nruns:
            trees = []
            treeLists.append(trees)
            for 0 <= j < self.nsamples:
                samp = <Samp>(<list>self.runs[i])[j]
                tree = samp.tree
                trees.append(tree)
        self._sumt = Sumt(treeLists)
    cdef Sumt getSumt(self):
        if self._sumt is None:
            self._summarizeTrees()
        return self._sumt
    property sumt:
        """
            Sumt instance associated with posterior trees.
        """
        def __get__(self):
            return self.getSumt()
