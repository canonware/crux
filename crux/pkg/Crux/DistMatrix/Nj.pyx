#===============================================================================
# This file implements the neigbhor joining (NJ) and relaxed neighbor joining
# (RNJ) algorithms.
#
# Matrices are stored as upper-half matrices, which means that addressing is a
# bit tricky.  The following matrix shows how addressing logically works, as
# well as the order in which matrix elements are stored in memory:
#
#   x: 0  1  2  3  4  5  6  7  8
#     --+--+--+--+--+--+--+--+--+ y:
#       | 0| 1| 2| 3| 4| 5| 6| 7| 0
#       +--+--+--+--+--+--+--+--+
#          | 8| 9|10|11|12|13|14| 1
#          +--+--+--+--+--+--+--+
#             |15|16|17|18|19|20| 2
#             +--+--+--+--+--+--+
#                |21|22|23|24|25| 3
#                +--+--+--+--+--+
#                   |26|27|28|29| 4
#                   +--+--+--+--+
#                      |30|31|32| 5
#                      +--+--+--+
#                         |33|34| 6
#                         +--+--+
#                            |35| 7
#                            +--+
#                               | 8
#
# The following formula can be used to convert from (x,y) coordinates to array
# offsets:
#
#   n : Number of nodes currently in the matrix.
#   x : Row.
#   y : Column.
#
#                        2
#                       x  + 3x
#   f(n,x,y) = nx + y - ------- - 1
#                          2
#
#===============================================================================
#
# Throughout the code, the matrix is stored as an upper-triangle symmetric
# matrix, with additional bookkeeping, as necessary for RNJ and NJ.  For
# example (n is current matrix size):
#
#                             |  r  ||
#                             | --- ||
#   | A | B | C | D | E ||  r | n-2 ||
#   +===+===+===+===+===++====+=====++===
#       | 0 | 1 | 2 | 3 ||  6 | 2.0 || A
#       +---+---+---+---++----+-----++---
#           | 4 | 5 | 6 || 15 | 5.0 || B
#           +---+---+---++----+-----++---
#               | 7 | 8 || 20 | 6.7 || C
#               +---+---++----+-----++---
#                   | 9 || 23 | 7.7 || D
#                   +---++----+-----++---
#                       || 26 | 8.7 || E
#                       ++----+-----++---
#
# is stored as:
#
#   d:
#   |---+---+---+---+---+---+---+---+---+---|
#   | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
#   |---+---+---+---+---+---+---+---+---+---|
#
#   r:
#   |------+------+------+------+------|
#   |  6   | 15   | 20   | 23   | 26   |
#   |------+------+------+------+------|
#
#   rScaled:
#   |------+------+------+------+------|
#   |  2.0 |  5.0 |  6.7 |  7.7 |  8.7 |
#   |------+------+------+------+------|
#
#   nodes:
#   |------+------+------+------+------|
#   |  A   |  B   |  C   |  D   |  E   |
#   |------+------+------+------+------|
#
#===============================================================================
#
# Since neighbor joining involves repeatedly joining two nodes and removing a
# row from the matrix, the performance of matrix collapsing is important.
# Additionally, it is important to keep the matrix compactly stored in memory,
# for cache locality reasons.  This implementation removes row x by moving row
# 0 into its place, then discarding row 0 of the array.  This has the effects
# of 1) re-ordering rows, and 2) shifting row addresses, so care is necessary
# in code that both iterates over rows and collapses the matrix.
#
#===============================================================================
#
# The key to RNJ's performance is the way in which clustering (node joining)
# decisions are made.  Rather than calculating the transformed distances for
# all possible node pairings and joining those two nodes as in NJ, RNJ
# iteratively checks to see if various possible pairings of nodes are legal,
# according to the constraints of the neighbor joining algorithm.  For example,
# consider the joining of nodes 4 and 6:
#
#   x: 0  1  2  3  4  5  6  7  8
#     --+--+--+--+--+--+--+--+--+ y:
#       |  |  |  |XX|  |YY|  |  | 0
#       +--+--+--+--+--+--+--+--+
#          |  |  |XX|  |YY|  |  | 1
#          +--+--+--+--+--+--+--+
#             |  |XX|  |YY|  |  | 2
#             +--+--+--+--+--+--+
#                |XX|  |YY|  |  | 3
#                +--+--+--+--+--+
#                   |XX|**|XX|XX| 4
#                   +--+--+--+--+
#                      |YY|  |  | 5
#                      +--+--+--+
#                         |YY|YY| 6
#                         +--+--+
#                            |  | 7
#                            +--+
#                               | 8
#
# As long as:
#   1) the transformed distance for (4,6), denoted by **, is less than or equal
#      to the transformed distances for the matrix elements marked by XX or YY,
#      and
#   2) joining rows 4 and 6 does not change the additivity of distances for the
#      nodes still represented in the distance matrix (only applicable when
#      dealing with additive distance matrices),
# then joining nodes 4 and 6 poses no correctness problems for the neighbor
# joining algorithm.  This implementation searches for such clusterings in an
# efficient manner.
#
# It is important to note that in the worst case, this implementation has
# O(n^3) performance.  However, worst case performance requires that the tree
# be very long, with a particular pattern of branch lengths, and that the taxa
# be inserted into the matrix in a particular order.  As such, worst case
# performance almost never occurs, and typical runtime is approximately
# proportional to (n^2 lg n).
#
#===============================================================================

cdef extern from "sys/types.h":
    ctypedef unsigned long size_t

cdef extern from "stdlib.h":
    cdef void *calloc(size_t nmemb, size_t size)
    cdef void *malloc(size_t size)
    cdef void free(void *ptr)

cdef extern from "math.h":
    cdef double HUGE_VAL

from CxDistMatrix cimport *
from CxDistMatrixNj cimport *

from Crux.Tree cimport Tree, Node, Edge

import random

cdef class Nj:
    def __cinit__(self):
        self.dBase = NULL
        self.d = NULL
        self.rBase = NULL
        self.r = NULL
        self.rScaledBase = NULL
        self.rScaled = NULL

    def __dealloc__(self):
        if self.dBase != NULL:
            free(self.dBase)
            self.dBase = NULL
        if self.rBase != NULL:
            free(self.rBase)
            self.rBase = NULL
        if self.rScaledBase != NULL:
            free(self.rScaledBase)
            self.rScaledBase = NULL

    cdef void _rInit(self) except *:
        cdef CxtDMDist *d, *r, dist
        cdef CxtDMSize x, y, i

        # Allocate zeroed vector.
        r = <CxtDMDist *>calloc(self.nBase, sizeof(CxtDMDist))
        if r == NULL:
            raise MemoryError("Failed to allocate %d-taxon distance vector" %
              <int>self.nBase)

        # Calculate r (sum of distances to other nodes) for each node.
        d = self.d
        i = 0
        for 0 <= x < self.nBase:
            for x + 1 <= y < self.nBase:
                dist = d[i]
                i += 1

                r[x] += dist
                r[y] += dist

        self.rBase = r
        self.r = r

    cdef void _rScaledInit(self) except *:
        cdef CxtDMDist *rScaled

        rScaled = <CxtDMDist *>malloc(sizeof(CxtDMDist) * self.nBase)
        if rScaled == NULL:
            raise MemoryError("Failed to allocate %d-taxon distance vector" %
              <int>self.nBase)

        self.rScaledBase = rScaled
        self.rScaled = rScaled

    cdef void _nodesInit(self, Taxa.Map taxaMap) except *:
        cdef int i
        cdef Node node

        assert taxaMap.ntaxa == self.nBase

        self.nodes = []
        for 0 <= i < self.nBase:
            node = Node(self.tree)
            node.taxon = taxaMap.taxonGet(i)
            self.nodes.append(node)

    cdef void _rScaledUpdate(self):
        cdef CxtDMSize x, denom

        # Calculate rScaled for each node.
        denom = self.n - 2
        for 0 <= x < self.n:
            self.rScaled[x] = self.r[x] / denom

    cdef void _njRandomMinFind(self, CxtDMSize *rXMin, CxtDMSize *rYMin):
        cdef CxtDMSize nmins, i, x, y, xMin, yMin
        cdef CxtDMDist transMin, transCur

        # Calculate the transformed distance for each pairwise distance.  Keep
        # track of the minimum transformed distance, so that the corresponding
        # nodes can be joined.
        #
        # This is by far the most time-consuming portion of NJ.
        nmins = 1
        transMin = HUGE_VAL
        i = 0
        for 0 <= x <= self.n:
            for x + 1 <= y < self.n:
                transCur = self.d[i]
                i += 1

                # Use CxDistMatrixNjDistCompare() in order to compare
                # transformed distances, so that random selection is possible.
                rel = CxDistMatrixNjDistCompare(transCur, transMin)
                if rel == -1:
                    nmins = 1
                    xMin = x
                    yMin = y
                    transMin = transCur
                elif rel == 0:
                    # Choose such that all tied distances have an equal
                    # probability of being chosen.
                    nmins += 1
                    if random.randint(0, nmins - 1) == 0:
                        xMin = x
                        yMin = y
                        transMin = transCur
                elif rel == 1:
                    pass
                else:
                    assert False

                # Since an arbitrary tie-breaking decision is being made anyway,
                # don't bother using CxDistMatrixNjDistCompare() here.
                if transCur < transMin:
                    xMin = x
                    yMin = y
                    transMin = transCur

        rXMin[0] = xMin
        rYMin[0] = yMin

    cdef void _njDeterministicMinFind(self, CxtDMSize *rXMin, CxtDMSize *rYMin):
        cdef CxtDMSize i, x, y, xMin, yMin
        cdef CxtDMDist transMin, transCur

        # Calculate the transformed distance for each pairwise distance.  Keep
        # track of the minimum transformed distance, so that the corresponding
        # nodes can be joined.
        #
        # This is by far the most time-consuming portion of NJ.
        transMin = HUGE_VAL
        i = 0
        for 0 <= x <= self.n:
            for x + 1 <= y < self.n:
                transCur = self.d[i]
                i += 1

                # Since an arbitrary tie-breaking decision is being made anyway,
                # don't bother using CxDistMatrixNjDistCompare() here.
                if transCur < transMin:
                    xMin = x
                    yMin = y
                    transMin = transCur

        rXMin[0] = xMin
        rYMin[0] = yMin

    cdef Node _njNodesJoin(self, CxtDMSize xMin, CxtDMSize yMin,
      CxtDMDist *rDistX, CxtDMDist *rDistY):
        cdef Node node
        cdef Edge edgeX, edgeY
        cdef CxtDMSize iMin
        cdef CxtDMDist distX, distY

        node = Node(self.tree)

        edgeX = Edge(self.tree)
        edgeX.attach(node, self.nodes[xMin])
        iMin = CxDistMatrixNxy2i(self.n, xMin, yMin)
        distX = (self.d[iMin] + self.rScaled[xMin] - self.rScaled[yMin]) / 2
        edgeX.length = distX

        edgeY = Edge(self.tree)
        edgeY.attach(node, self.nodes[yMin])
        distY = self.d[iMin] - distX
        edgeY.length = distY

        rDistX[0] = distX
        rDistY[0] = distY
        return node

    cdef void _njRSubtract(self, CxtDMSize xMin, CxtDMSize yMin):
        assert False # XXX

    cdef void _njCompact(self, CxtDMSize xMin, CxtDMSize yMin, Node node,
      CxtDMDist distX, CxtDMDist distY) except *:
        assert False # XXX

    cdef void _njDiscard(self):
        assert False # XXX

    cdef void prepare(self, CxtDMDist *d, CxtDMSize n, Taxa.Map taxaMap) \
      except *:
        assert self.dBase == NULL # Was self.prepare() already called?

        self.dBase = d
        self.d = d
        self.nBase = n

        self._rInit()
        self._rScaledInit()
        self.tree = Tree()
        self._nodesInit(taxaMap)

    cdef Tree nj(self, bint random):
        cdef CxtDMSize nleft, xMin, yMin
        cdef CxtDMDist distX, distY
        cdef Node node

        assert self.dBase != NULL # Was self.prepare() called?

        # Iteratively join two nodes in the matrix, until only two are left.
        for self.nBase >= nleft > 2:
            # Standard neighbor joining.
            self.n = nleft
            self._rScaledUpdate()
            if random:
                self._njRandomMinFind(&xMin, &yMin)
            else:
                self._njDeterministicMinFind(&xMin, &yMin)
            node = self._njNodesJoin(xMin, yMin, &distX, &distY)
            self._njRSubtract(xMin, yMin)
            self._njCompact(xMin, yMin, node, distX, distY)
            self._njDiscard()

        # Join last two nodes.
        self._njFinalJoin()

        return self.tree

cdef class Rnj(Nj):
    cdef Tree rnj(self, bint random, bint additive):
#        if random:
#            self._njRandomCluster()
#        else:
#            self._njDeterministicCluster()
        # XXX

        return self.tree
