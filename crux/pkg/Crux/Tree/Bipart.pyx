from Crux.Tree cimport Tree, Node, Edge, Ring

from Cx cimport CxCmp2Richcmp
from libc cimport *

cdef class Vec:
    def __cinit__(self):
        self.bits = NULL

    def __dealloc__(self):
        if self.bits != NULL:
            free(self.bits)
            self.bits = NULL

    def __init__(self, unsigned nBits):
        cdef unsigned nBytes

        nBytes = nBits >> 3
        if (nBits & 0x7) != 0:
            nBytes += 1
        nBits = nBytes << 3

        self.nBits = nBits
        self.bits = <unsigned char *>calloc(1, nBytes)
        if self.bits == NULL:
            raise MemoryError("Error allocating %d-bit vector" % nBits)

    def __richcmp__(Vec self, Vec other, int op):
        cdef int rel

        assert self.nBits == other.nBits

        rel = memcmp(self.bits, other.bits, (self.nBits >> 3))
        rel = (rel > 0) - (rel < 0)
        return CxCmp2Richcmp(rel, op)

    cdef int cmp(Vec self, Vec other) except *:
        cdef int rel

        assert self.nBits == other.nBits

        rel = memcmp(self.bits, other.bits, (self.nBits >> 3))
        rel = (rel > 0) - (rel < 0)
        return rel

    cdef reset(self):
        memset(self.bits, 0, self.nBits >> 3)

    cdef bint get(self, unsigned bit) except *:
        cdef bint ret
        cdef unsigned byteOffset, bitOffset
        cdef unsigned char byte

        assert bit < self.nBits

        byteOffset = bit >> 3
        bitOffset = bit & 0x7

        byte = self.bits[byteOffset]
        byte >>= bitOffset
        byte &= 0x1
        ret = <bint>byte

        return ret

    cdef void set(self, unsigned bit, bint val) except *:
        cdef unsigned byteOffset, bitOffset
        cdef unsigned char byte, mask

        assert bit < self.nBits

        byteOffset = bit >> 3
        bitOffset = bit & 0x7

        byte = self.bits[byteOffset]

        mask = (0xff ^ (0x1 << bitOffset))
        byte &= mask
        byte |= (val << bitOffset)

        self.bits[byteOffset] = byte

    cdef void invert(self):
        cdef unsigned i, nBytes

        nBytes = self.nBits >> 3
        for 0 <= i < nBytes:
            self.bits[i] = (~self.bits[i])

    cdef void merge(self, Vec other) except *:
        cdef unsigned i, nBytes

        assert self.nBits == other.nBits

        nBytes = self.nBits >> 3
        for 0 <= i < nBytes:
            self.bits[i] |= other.bits[i]

cdef class Bipart:
    def __init__(self, Tree tree):
        cdef list taxa

        # Create a taxon-->index translation.
        self.taxaX = {}
        taxa = tree.getTaxa()
        for 0 <= i < len(taxa):
            if taxa[i] in self.taxaX:
                raise ValueError("Duplicate taxa")
            self.taxaX[taxa[i]] = i

        self.edgeVecs = []
        self.leafVec = Vec(len(taxa))

        self._bipartitions(tree)

    cdef Vec _bipartitionsRecurse(self, Ring ring, bint calcVec):
        cdef Vec ret, vec
        cdef Node node

        node = ring.node
        if node.getDegree() <= 1:
            # Leaf node.

            # Use the special temp vector for this edge.
            ret = self.leafVec
            ret.reset()

            # Mark this taxon as being in the set.
            assert node._taxon is not None
            ret.set(self.taxaX[node._taxon], True)
        else:
            if calcVec:
                # Create a vector for this edge.
                ret = Vec(len(self.taxaX))
                self.edgeVecs.append(ret)

                # Even internal nodes may have associated taxa.
                if node._taxon is not None:
                    ret.set(self.taxaX[node._taxon], True)
            else:
                # Avoid creating a vector, because the previous node is the
                # tree base, and it's a leaf node.
                ret = None

            # Recurse.
            for r in ring.siblings():
                vec = self._bipartitionsRecurse(r.other, True)
                if calcVec:
                    ret.merge(vec)

        return ret

    cdef void _bipartitions(self, Tree tree) except *:
        cdef Node base
        cdef Ring ring, r
        cdef int baseDegree
        cdef unsigned i
        cdef Vec vec

        base = tree.getBase()
        if base is not None:
            baseDegree = base.getDegree()
            ring = base.ring
            if ring is not None:
                for r in ring:
                    self._bipartitionsRecurse(r.other, (baseDegree != 1))

        # Iterate through the bipartitions and make sure that bit 0 is never
        # set.
        for 0 <= i < len(self.edgeVecs):
            vec = <Vec>self.edgeVecs[i]
            if vec.get(0):
                vec.invert()

        # Sort the list of bipartitions, so that lists can be linearly scanned
        # for differences.
        self.edgeVecs.sort()

    cdef double rfDist(self, Bipart other) except -1.0:
        cdef double falseNegativeRate, falsePositiveRate
        cdef unsigned nUniqueA, nUniqueB, iA, iB, lenA, lenB
        cdef int rel
        cdef Vec vecA, vecB

        if self is other:
            return 0.0

        nUniqueA = nUniqueB = iA = iB = 0
        lenA = len(self.edgeVecs)
        lenB = len(other.edgeVecs)
        while iA < lenA and iB < lenB:
            vecA = <Vec>self.edgeVecs[iA]
            vecB = <Vec>other.edgeVecs[iB]

            rel = vecA.cmp(vecB)
            if rel == -1:
                # self has a unique bipartition.
                nUniqueA += 1
                iA += 1
            elif rel == 0:
                iA += 1
                iB += 1
            elif rel == 1:
                # other has a unique bipartition.
                nUniqueB += 1
                iB += 1
            else:
                assert False

        if iA < lenA:
            nUniqueA += lenA - iA
        elif iB < lenB:
            nUniqueB += lenB - iB

        # Convert counts to the Robinson-Foulds distance.
        if lenA > 0:
            falseNegativeRate = <double>nUniqueA / <double>lenA
        else:
            falseNegativeRate = 0.0

        if lenB > 0:
            falsePositiveRate = <double>nUniqueB / <double>lenB
        else:
            falsePositiveRate = 0.0

        return (falseNegativeRate + falsePositiveRate) / 2.0

cdef double rf(Tree a, Tree b) except -1.0:
    cdef list taxaA
    cdef dict taxa
    cdef unsigned i
    cdef Bipart bipartA, bipartB

    if a.getTaxa() != b.getTaxa():
        # The trees do not contain the same taxa, so treat them as being
        # infinitely distant.
        return 1.0

    bipartA = a.getBipart()
    bipartB = b.getBipart()

    return bipartA.rfDist(bipartB)

cdef list rfs(Tree a, list others):
    cdef list ret, taxaA
    cdef Tree b
    cdef Bipart bipartA, bipartB

    ret = []
    taxaA = a.getTaxa()

    bipartA = Bipart(a)
    for b in others:
        if taxaA != b.getTaxa():
            # The trees do not contain the same taxa, so treat them as being
            # infinitely distant.
            ret.append(1.0)
        else:
            bipartB = Bipart(b)
            ret.append(bipartA.rfDist(bipartB))
    return ret
