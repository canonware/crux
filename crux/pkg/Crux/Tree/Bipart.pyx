from Crux.Tree cimport Tree, Node, Edge, Ring

from Cx cimport CxCmp2Richcmp
from libc cimport *

cdef extern from "Python.h":
    cdef object PyString_FromStringAndSize(char *s, Py_ssize_t len)

cdef class Vec:
    def __cinit__(self):
        self.bits = NULL

    def __dealloc__(self):
        if self.bits != NULL:
            free(self.bits)
            self.bits = NULL

    def __init__(self, Edge edge, unsigned nBits):
        cdef unsigned nBytes

        self.edge = edge
        self.nBits = nBits

        self.nBytes = nBits >> 3
        if (nBits & 0x7) != 0:
            self.nBytes += 1

        self.bits = <unsigned char *>calloc(1, self.nBytes)
        if self.bits == NULL:
            raise MemoryError("Error allocating %d-bit vector" % nBits)

    def __repr__(self):
        cdef list strs
        cdef unsigned i

        strs = []
        for 0 <= i < self.nBits:
            strs.append("*" if self.get(i) else ".")
        return "".join(strs)

    def __hash__(self):
        cdef Py_ssize_t len
        cdef str s

        # Only hash the first VecHashMax bytes.
        DEF VecHashMax = 16
        len = (<Py_ssize_t>VecHashMax if \
          self.nBytes > <Py_ssize_t>VecHashMax else \
          self.nBytes)
        s = PyString_FromStringAndSize(<char *>self.bits, len)
        return hash(s)

    def __richcmp__(Vec self, Vec other, int op):
        cdef int rel

        assert self.nBits == other.nBits

        rel = memcmp(self.bits, other.bits, self.nBytes)
        rel = (rel > 0) - (rel < 0)
        return CxCmp2Richcmp(rel, op)

    cpdef int cmp(self, Vec other) except *:
        cdef int rel

        assert self.nBits == other.nBits

        rel = memcmp(self.bits, other.bits, self.nBytes)
        rel = (rel > 0) - (rel < 0)
        return rel

    cpdef reset(self):
        memset(self.bits, 0, self.nBytes)

    cpdef bint get(self, unsigned bit) except *:
        cdef bint ret
        cdef unsigned byteOffset, bitOffset
        cdef unsigned char byte

        assert bit < self.nBits

        byteOffset = bit >> 3
        bitOffset = 7 - (bit & 0x7)

        byte = self.bits[byteOffset]
        byte >>= bitOffset
        byte &= 0x1
        ret = <bint>byte

        return ret

    cdef _set(self, unsigned bit, bint val):
        cdef unsigned byteOffset, bitOffset
        cdef unsigned char byte, mask

        byteOffset = bit >> 3
        bitOffset = 7 - (bit & 0x7)

        byte = self.bits[byteOffset]

        mask = (0xff ^ (0x1 << bitOffset))
        byte &= mask
        byte |= (val << bitOffset)

        self.bits[byteOffset] = byte

    cpdef set(self, unsigned bit, bint val):
        assert bit < self.nBits
        self._set(bit, val)

    cpdef invert(self):
        cdef unsigned i, nBytes

        for 0 <= i < self.nBytes:
            self.bits[i] = (~self.bits[i])

        # Clear trailing bits in order to allow memcmp() to be used for
        # comparison.
        for self.nBits <= i < (self.nBytes << 3):
            self._set(i, False)

    cpdef merge(self, Vec other):
        cdef unsigned i, nBytes

        assert self.nBits == other.nBits

        for 0 <= i < self.nBytes:
            self.bits[i] |= other.bits[i]

cdef class Bipart:
    def __init__(self, Tree tree, bint leaves=False):
        cdef list taxa

        # Create a taxon-->index translation.
        self.taxaX = {}
        taxa = tree.getTaxa()
        for 0 <= i < len(taxa):
            if taxa[i] in self.taxaX:
                raise ValueError("Duplicate taxa")
            self.taxaX[taxa[i]] = i

        self.leaves = leaves
        self.edgeVecs = []
        if not leaves:
            self.leafVec = Vec(None, len(taxa))

        self._bipartitions(tree)

    def __repr__(self):
        cdef list strs
        cdef unsigned i

        strs = []
        for 0 <= i < len(self.edgeVecs):
            strs.append((<Vec>self.edgeVecs[i]).__repr__())

        return "\n".join(strs)

    def __hash__(self):
        cdef unsigned i, lim, hash

        # Only hash the first BipartHashMax vectors.
        DEF BipartHashMax = 16
        lim = (BipartHashMax if len(self.edgeVecs) > BipartHashMax \
          else len(self.edgeVecs))
        hash = 0
        for 0 <= i < lim:
            hash = (hash + (<Vec>self.edgeVecs[i]).__hash__()) & 0x7fffffff
        return <int>hash

    def __richcmp__(Bipart self, Bipart other, int op):
        cdef int rel

        rel = cmp(self.edgeVecs, other.edgeVecs)
        rel = (rel > 0) - (rel < 0)
        return CxCmp2Richcmp(rel, op)

    cpdef int cmp(Bipart self, Bipart other):
        cdef int rel

        rel = cmp(self.edgeVecs, other.edgeVecs)
        rel = (rel > 0) - (rel < 0)
        return rel

    cdef Vec _bipartitionsRecurse(self, Ring ring, bint calcVec):
        cdef Vec ret, vec
        cdef Node node

        node = ring.node
        if node.getDegree() <= 1:
            # Leaf node.

            if self.leaves:
                ret = Vec(ring.edge, len(self.taxaX))
                self.edgeVecs.append(ret)
            else:
                # Use the special temp vector for this edge.
                ret = self.leafVec
                ret.reset()

            # Mark this taxon as being in the set.
            assert node._taxon is not None
            ret.set(self.taxaX[node._taxon], True)
        else:
            if calcVec:
                # Create a vector for this edge.
                ret = Vec(ring.edge, len(self.taxaX))
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

    cpdef double rfDist(self, Bipart other) except -1.0:
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
