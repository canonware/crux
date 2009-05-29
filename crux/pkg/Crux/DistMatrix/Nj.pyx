#===============================================================================
# This file implements the neigbhor joining (NJ) and relaxed neighbor joining
# (RNJ) algorithms.
#
#   Evans, J., L. Sheneman, J.A. Foster (2006) Relaxed neighbor joining: A
#   fast distance-based phylogenetic tree construction method.  J Mol Evol
#   62:785-792.
#
#   Saitou, N., M. Nei (1987) The neighbor-joining method: A new method for
#   reconstructing phylogenetic trees.  Mol Biol Evol 4:406-425.
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

from libc cimport *
from libm cimport *

from SFMT cimport *
from CxRi cimport *

from Crux.Tree cimport Tree, Node, Edge
from Crux.DistMatrix cimport nxy2i

DEF NjDebug = False # Note companion DEF in Nj.pxd.
IF NjDebug:
    import sys

cdef union fiUnion:
    float d
    int32_t i

cdef enum:
    MaxUlps = 0x7f

# Compare two distances, and consider them equal if they are close enough.
# MaxUlps (ULP: Units in Last Place) specifies the maximum number of ulps that
# a and b may differ by and still be considered equal.  Ideally, we would keep
# track of the maximum possible accumulated error during tree construction, but
# the error grows exponentially, so that approach ends up not being practical.
#
# Provide 5 (base 10) digits of accuracy.
cdef inline bint distEq(float a, float b):
    cdef fiUnion aU, bU

    assert sizeof(float) == sizeof(int32_t)

    # Convert a and b to lexicographically ordered ints.
    aU.d = a
    if aU.i < 0:
        aU.i = <uint32_t>0x80000000 - aU.i

    bU.d = b
    if bU.i < 0:
        bU.i = <uint32_t>0x80000000U - bU.i

    # Check if a and b are within MaxUlps of each other.
    return abs(aU.i - bU.i) <= MaxUlps

cdef class Nj:
    def __cinit__(self):
        self.prng = NULL
        self.dBase = NULL
        self.d = NULL
        self.rBase = NULL
        self.r = NULL
        self.rScaledBase = NULL
        self.rScaled = NULL

    def __dealloc__(self):
        if self.prng != NULL:
            fini_gen_rand(self.prng)
            self.prng = NULL
        if self.dBase != NULL:
            free(self.dBase)
            self.dBase = NULL
        if self.rBase != NULL:
            free(self.rBase)
            self.rBase = NULL
        if self.rScaledBase != NULL:
            free(self.rScaledBase)
            self.rScaledBase = NULL

    def __init__(self):
        import random

        self.prng = init_gen_rand(random.randint(0, 0xffffffffU))
        if self.prng == NULL:
            raise MemoryError("Error initializing prng")

    IF NjDebug:
        cdef void _njDump(self) except *:
            cdef size_t i, x, y

            sys.stdout.write( \
              "----------------------------------------" \
              "----------------------------------------\n")
            i = 0
            for 0 <= x < self.n:
                # || node
                sys.stdout.write(" " * (x * 9))
                y -= self.n - (x + 1)
                for x + 1 <= y < self.n:
                    sys.stdout.write(" " * 9)
                taxon = self.nodes[x].taxon
                if taxon is not None:
                    sys.stdout.write(" || %s\n" % taxon.label)
                else:
                    sys.stdout.write(" || %r\n" % self.nodes[x])

                # dist || r
                sys.stdout.write(" " * (x * 9))
                y -= self.n - (x + 1)
                for y <= y < self.n:
                    sys.stdout.write(" %8.4f" % self.d[i])
                    i += 1
                sys.stdout.write(" || %8.4f\n" % self.r[x])

                # tdist || rScaled
                sys.stdout.write(" " * (x * 9))
                y -= self.n - (x + 1)
                i -= self.n - (x + 1)
                for y <= y < self.n:
                    sys.stdout.write(" %8.4f" % \
                      (self.d[i] - (self.rScaled[x] + self.rScaled[y])))
                    i += 1
                sys.stdout.write(" || %8.4f\n" % self.rScaled[x])

    cdef void _rInit(self) except *:
        cdef float *d, *r, dist
        cdef size_t x, y, i

        # Allocate zeroed vector.
        r = <float *>calloc(self.nBase, sizeof(float))
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
        cdef float *rScaled

        rScaled = <float *>malloc(sizeof(float) * self.nBase)
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
        cdef size_t x, denom

        # Calculate rScaled for each node.
        denom = self.n - 2
        for 0 <= x < self.n:
            self.rScaled[x] = self.r[x] / denom

    cdef void _njRandomMinFind(self, size_t *rXMin, size_t *rYMin):
        cdef size_t nmins, i, n, x, y, xMin, yMin
        cdef float *d, *rScaled
        cdef float transMin, transCur, rScaledX

        IF @enable_cc_silence@:
            xMin = yMin = 0

        # Calculate the transformed distance for each pairwise distance.  Keep
        # track of the minimum transformed distance, so that the corresponding
        # nodes can be joined.
        #
        # This is by far the most time-consuming portion of NJ.
        nmins = 1
        transMin = HUGE_VALF
        i = 0
        d = self.d
        rScaled = self.rScaled
        n = self.n
        for 0 <= x < n:
            rScaledX = rScaled[x]
            for x + 1 <= y < n:
                transCur = d[i] - (rScaledX + rScaled[y])
                i += 1

                # Use distEq() in order to compare transformed distances,
                # so that random selection is possible.
                if distEq(transCur, transMin):
                    # Choose such that all tied distances have an equal
                    # probability of being chosen.
                    nmins += 1
                    if gen_rand64_range(self.prng, nmins) == 0:
                        xMin = x
                        yMin = y
                        transMin = transCur
                elif transCur < transMin:
                    nmins = 1
                    xMin = x
                    yMin = y
                    transMin = transCur

        rXMin[0] = xMin
        rYMin[0] = yMin

    cdef void _njDeterministicMinFind(self, size_t *rXMin, size_t *rYMin):
        cdef size_t i, n, x, y, xMin, yMin
        cdef float *d, *rScaled
        cdef float transMin, transCur, rScaledX

        IF @enable_cc_silence@:
            xMin = yMin = 0

        # Calculate the transformed distance for each pairwise distance.  Keep
        # track of the minimum transformed distance, so that the corresponding
        # nodes can be joined.
        #
        # This is by far the most time-consuming portion of NJ.
        transMin = HUGE_VALF
        i = 0
        d = self.d
        rScaled = self.rScaled
        n = self.n
        for 0 <= x < n:
            rScaledX = rScaled[x]
            for x + 1 <= y < n:
                transCur = d[i] - (rScaledX + rScaled[y])
                i += 1

                # Since an arbitrary tie-breaking decision is being made anyway,
                # don't bother using distEq() here.
                if transCur < transMin:
                    xMin = x
                    yMin = y
                    transMin = transCur

        rXMin[0] = xMin
        rYMin[0] = yMin

    cdef Node _njNodesJoin(self, size_t xMin, size_t yMin,
      float *rDistX, float *rDistY):
        cdef Node node
        cdef Edge edgeX, edgeY
        cdef size_t iMin
        cdef float distX, distY

        node = Node(self.tree)

        edgeX = Edge(self.tree)
        edgeX.attach(node, self.nodes[xMin])
        iMin = nxy2i(self.n, xMin, yMin)
        distX = (self.d[iMin] + self.rScaled[xMin] - self.rScaled[yMin]) / 2
        edgeX.length = distX

        edgeY = Edge(self.tree)
        edgeY.attach(node, self.nodes[yMin])
        distY = self.d[iMin] - distX
        edgeY.length = distY

        rDistX[0] = distX
        rDistY[0] = distY
        return node

    cdef void _njRSubtract(self, size_t xMin, size_t yMin):
        cdef size_t n, x, iX, iY
        cdef float *d, *r
        cdef float dist

        # Subtract old distances from r.
        d = self.d
        r = self.r
        n = self.n
        iX = xMin - 1
        iY = yMin - 1
        for 0 <= x < xMin:
            dist = d[iX]
            iX += n - 2 - x
            r[x] -= dist

            dist = d[iY]
            iY += n - 2 - x
            r[x] -= dist

        # (x == xMin)
        iY += n - 2 - x
        x += 1

        for x <= x < yMin:
            iX += 1
            dist = d[iX]
            r[x] -= dist

            dist = d[iY]
            iY += n - 2 - x
            r[x] -= dist

        # (x == yMin)
        iX += 1
        dist = d[iX]
        r[x] -= dist
        x += 1

        for x <= x < n:
            iX += 1
            dist = d[iX]
            r[x] -= dist

            iY += 1
            dist = d[iY]
            r[x] -= dist

        # Rather than repeatedly subtracting distances from aR[aXMin] and
        # aR[aYMin] (and accumulating floating point error), simply clear these
        # two elements of r.
        r[xMin] = 0.0
        r[yMin] = 0.0

    cdef void _njCompact(self, size_t xMin, size_t yMin, Node node,
      float distX, float distY) except *:
        cdef size_t n, x, iX, iY
        cdef float *d, *r
        cdef float dist

        # Insert the new node, such that it overwrites one of its children.
        self.nodes[xMin] = node

        # Calculate distances to the new node, and add them to r.  This clobbers
        # old distances, just after the last time they are needed.
        d = self.d
        r = self.r
        n = self.n
        iX = xMin - 1
        iY = yMin - 1
        for 0 <= x < xMin:
            dist = ((d[iX] - distX) + (d[iY] - distY)) / 2
            d[iX] = dist
            iX += n - 2 - x
            iY += n - 2 - x
            r[x] += dist
            r[xMin] += dist

        # (x == xMin)
        iY += n - 2 - x
        x += 1

        for x <= x < yMin:
            iX += 1
            dist = ((d[iX] - distX) + (d[iY] - distY)) / 2
            d[iX] = dist
            iY += n - 2 - x
            r[x] += dist
            r[xMin] += dist

        # (x == yMin)
        iX += 1
        x += 1

        for x <= x < n:
            iX += 1
            iY += 1
            dist = ((d[iX] - distX) + (d[iY] - distY)) / 2
            d[iX] = dist
            r[x] += dist
            r[xMin] += dist

        # Fill in the remaining gap (yMin row/column), by moving the first row
        # into the gap.  The first row can be removed from the matrix in
        # constant time, whereas collapsing the gap directly would require a
        # series of memmove() calls, and leaving the gap would result in
        # increased cache misses.
        iX = 0
        iY = n + yMin - 3
        for 1 <= x < yMin:
            d[iY] = d[iX]
            iY += n - 2 - x
            iX += 1

        # (x == yMin)
        iX += 1
        x += 1

        for x <= x < n:
            iY += 1
            d[iY] = d[iX]
            iX += 1

        # Fill in the gap in r, and nodes.  rScaled is re-calculated from
        # scratch, so there is no need to touch it here.
        r[yMin] = r[0]
        self.nodes[yMin] = self.nodes[0]

    cdef void _njDiscard(self):
        # Move pointers forward, which removes the first row.
        self.d = &self.d[self.n - 1]
        self.r = &self.r[1]
        self.rScaled = &self.rScaled[1]
        self.nodes.pop(0)

    cdef void _njFinalJoin(self) except *:
        cdef Edge edge

        # Join the remaining two nodes.
        edge = Edge(self.tree)
        edge.attach(self.nodes[0], self.nodes[1])
        edge.length = self.d[0]

        self.tree.base = self.nodes[0]

    cdef void prepare(self, float *d, size_t n, Taxa.Map taxaMap) \
      except *:
        assert self.dBase == NULL # Was self.prepare() already called?

        self.dBase = d
        self.d = d
        self.nBase = n
        self.n = self.nBase

        self._rInit()
        self._rScaledInit()
        self._rScaledUpdate()
        self.tree = Tree(rooted=False)
        self._nodesInit(taxaMap)

    cdef Tree nj(self, bint random):
        cdef size_t xMin, yMin
        cdef float distX, distY
        cdef Node node

        assert self.tree is not None # Was self.prepare() called?

        # Iteratively join two nodes in the matrix, until only two are left.
        while self.n > 2:
            # Standard neighbor joining.
            self._rScaledUpdate()
            IF NjDebug:
                self._njDump()
            if random:
                self._njRandomMinFind(&xMin, &yMin)
            else:
                self._njDeterministicMinFind(&xMin, &yMin)
            node = self._njNodesJoin(xMin, yMin, &distX, &distY)
            self._njRSubtract(xMin, yMin)
            self._njCompact(xMin, yMin, node, distX, distY)
            self._njDiscard()
            self.n -= 1

        # Join last two nodes.
        IF NjDebug:
            self._njDump()
        self._njFinalJoin()

        return self.tree

cdef class Rnj(Nj):
    cdef size_t _rnjRowAllMinFind(self, size_t x, float *rDist):
        cdef size_t ret, n, y, i, nmins
        cdef float *d, *rScaled
        cdef float dist, minDist

        IF @enable_cc_silence@:
            ret = <size_t>-1
            nmins = 0

        d = self.d
        rScaled = self.rScaled
        n = self.n

        minDist = HUGE_VALF

        # Find the minimum distance from the node on row x to any other node
        # that comes before it in the matrix.
        if x != 0:
            i = nxy2i(n, 0, x)
            for 0 <= y < x:
                dist = d[i] - (rScaled[y] + rScaled[x])
                i += (n - 2 - y)

                if distEq(dist, minDist):
                    # Choose y such that all tied distances have an equal
                    # probability of being chosen.
                    nmins += 1
                    if gen_rand64_range(self.prng, nmins) == 0:
                        ret = y
                elif dist < minDist:
                    nmins = 1
                    minDist = dist
                    ret = y
            assert minDist != HUGE_VALF

        # Find the minimum distance from the node on row x to any other node
        # that comes after it in the matrix.
        if x < n - 1:
            i = nxy2i(n, x, x + 1)
            for x + 1 <= y < n:
                dist = d[i] - (rScaled[x] + rScaled[y])
                i += 1

                if distEq(dist, minDist):
                    # Choose y such that all tied distances have an equal
                    # probability of being chosen.
                    nmins += 1
                    if gen_rand64_range(self.prng, nmins) == 0:
                        ret = y
                elif dist < minDist:
                    nmins = 1
                    minDist = dist
                    ret = y
            assert minDist != HUGE_VALF

        rDist[0] = minDist
        return ret

    cdef bint _rnjRowAllMinOk(self, size_t x, float minDist):
        cdef size_t n, y, i
        cdef float *d, *rScaled
        cdef float dist

        d = self.d
        rScaled = self.rScaled
        n = self.n

        # Make sure that minDist is <= any transformed distance in the row
        # portion of row x.
        if x + 1 < n:
            i = nxy2i(n, x, x + 1)
            for x + 1 <= y < n:
                dist = d[i] - (rScaled[x] + rScaled[y])
                i += 1

                if dist < minDist and not distEq(dist, minDist):
                    return False

        # Make sure that minDist is <= any transformed distance in the column
        # portion of row x.
        if x != 0:
            i = nxy2i(n, 0, x)
            for 0 <= y < x:
                dist = d[i] - (rScaled[y] + rScaled[x])
                i += (n - 2 - y)

                if dist < minDist and not distEq(dist, minDist):
                    return False

        return True

    cdef size_t _rnjRowMinFind(self, size_t x):
        cdef size_t ret, n, y, i
        cdef float *d, *rScaled
        cdef float dist, minDist

        IF @enable_cc_silence@:
            ret = <size_t>-1

        d = self.d
        rScaled = self.rScaled
        n = self.n

        # Find the minimum distance from the node on row x to any other node
        # that comes after it in the matrix.
        i = nxy2i(n, x, x + 1)
        minDist = HUGE_VALF
        for x + 1 <= y < n:
            dist = d[i] - (rScaled[x] + rScaled[y])
            i += 1

            # Don't bother using distEq() here, since this method is only
            # used for deterministic RNJ.
            if dist < minDist:
                minDist = dist
                ret = y
        assert minDist != HUGE_VALF

        return ret

    # Make sure that clustering a and b would not change the distances
    # between nodes.  This must be done in order to make sure that we get the
    # true tree, in the case that the distance matrix corresponds to precisely
    # one tree (distances are additive).  If the distances are non-additive
    # though, there is no need to do this check.
    cdef bint _rnjPairClusterAdditive(self, size_t a, size_t b):
        cdef size_t n, iAB, iA, iB
        cdef float *d, *rScaled
        cdef float distA, distB, dist

        d = self.d
        rScaled = self.rScaled
        n = self.n

        # Calculate distances from {a,b} to the new node.
        iAB = nxy2i(n, a, b)
        distA = (d[iAB] + rScaled[a] - rScaled[b]) / 2
        distB = d[iAB] - distA

        # Calculate one distance to the new node, and make sure that it is
        # consistent with current distances.  It is not necessary to check all
        # new distances, since an additivity violation will cause all new
        # distances to be inconsistent (ignoring variation in accumulated
        # floating point error).

        if (b + 1 < n):
            iA = nxy2i(n, a, b+1)
            iB = nxy2i(n, b, b+1)
            dist = ((d[iA] - distA) + (d[iB] - distB)) / 2
            return distEq(dist + distA, d[nxy2i(n, a, b+1)])
        elif a > 0:
            dist = ((d[a-1] - distA) + (d[b-1] - distB)) / 2
            return distEq(dist + distA, d[nxy2i(n, 0, a)])
        else:
            assert b > 1
            iA = a
            iB = b + n - 3
            dist = ((d[iA] - distA) + (d[iB] - distB)) / 2
            return distEq(dist + distA, d[nxy2i(n, a, 1)])

    # Finish checking whether it is okay to cluster rows a and b;
    # _rnjRowMinFind() or _rnjRowAllMinFind() has already done some of the work
    # by the time this method is called.
    #
    # Two nodes, a and b, can be clustered if the transformed distance between
    # them is less than or equal to the transformed distances from a or b to
    # any other node.
    cdef bint _rnjPairClusterOk(self, size_t a, size_t b):
        cdef size_t n, x, iA, iB
        cdef float *d, *rScaled
        cdef float distAB, dist

        assert a < b

        d = self.d
        rScaled = self.rScaled
        n = self.n

        # Calculate the transformed distance between a and b.
        distAB = d[nxy2i(n, a, b)] - (rScaled[a] + rScaled[b])

        # Iterate over the row portion of distances for b.  Distances for a were
        # already checked before this method was called.
        if b < n - 1:
            iB = nxy2i(n, b, b + 1)
            for b + 1 <= x < n:
                dist = d[iB] - (rScaled[x] + rScaled[b])
                # Don't bother using distEq() here, since this method is only
                # used for deterministic RNJ.
                if dist < distAB:
                    return False
                iB += 1

        # Iterate over the first column portion of distances for a and b.
        iA = a - 1
        iB = b - 1
        for 0 <= x < a:
            dist = d[iA] - (rScaled[x] + rScaled[a])
            # Don't bother using distEq() here, since this method is only used
            # for deterministic RNJ.
            if dist < distAB:
                return False

            dist = d[iB] - (rScaled[x] + rScaled[b])
            # Don't bother using distEq() here, since this method is only used
            # for deterministic RNJ.
            if dist < distAB:
                return False

            iA += (n - 2 - x)
            iB += (n - 2 - x)

        # (x == a)
        iB += (n - 2 - x)
        x += 1

        # Iterate over the second column portion of distances for b.  Distances
        # for a were already checked before this method was called.
        for x <= x < b:
            dist = d[iB] - (rScaled[x] + rScaled[b])
            # Don't bother using distEq() here, since this method is only
            # used for deterministic RNJ.
            if dist < distAB:
                return False
            iB += (n - 2 - x)

        return True

    # Try all clusterings of two rows in the matrix, in a random order.  Do
    # this in as cache-friendly a manner as possible (keeping in mind that the
    # matrix is stored in row-major form).  This means:
    #
    # 1) For each randomly chosen row (x), find the row which is the closest
    #    (y), according to transformed distances, by calling
    #    _rnjRowAllMinFind().
    #
    # 2) If the additivity constraint is enabled, check whether clustering x
    #    and y would violate additivity, by calling _rnjPairClusterAdditive().
    #
    # 2) Check whether it is okay to cluster x and y, by calling _rnjAllMinOk().
    cdef bint _rnjRandomCluster(self, bint additive) except -1:
        cdef CxtRi ri
        cdef bint clustered, done
        cdef size_t randomRow, closestRow, x, y
        cdef float dist, distX, distY
        cdef Node node

        CxRiNew(&ri, self.prng)
        if CxRiInit(&ri, self.n):
            raise MemoryError("Error in CxRiInit(..., %d)" % self.n)

        try:
            clustered = True
            done = False
            while not done:
                if not clustered:
                    additive = False
                    if CxRiInit(&ri, self.n):
                        raise MemoryError("Error in CxRiInit(..., %d)" % self.n)
                clustered = False
                while CxRiIndGet(&ri) < CxRiNintsGet(&ri):
                    # Randomly sample a row, without replacement.
                    randomRow = CxRiRandomGet(&ri)

                    # Find a row that is closest to randomRow.
                    closestRow = self._rnjRowAllMinFind(randomRow, &dist)

                    if randomRow < closestRow:
                        x = randomRow
                        y = closestRow
                    else:
                        x = closestRow
                        y = randomRow

                    # Make sure that no row is closer to y than x is.

                    if ((not additive) or self._rnjPairClusterAdditive(x, y)) \
                      and self._rnjRowAllMinOk(closestRow, dist):
                        clustered = True
                        IF NjDebug:
                            self._njDump()
                        node = self._njNodesJoin(x, y, &distX, &distY)
                        self._njRSubtract(x, y)
                        self._njCompact(x, y, node, distX, distY)
                        self._njDiscard()
                        self.n -= 1
                        self._rScaledUpdate()

                        # Shrinking the matrix may have reduced it to the point
                        # that the enclosing loop will no longer function
                        # correctly.  Check this condition here, in order to
                        # reduce branch overhead for the case where no join is
                        # done.
                        if self.n == 2:
                            done = True
                            break

                        if CxRiInit(&ri, self.n):
                            raise MemoryError("Error in CxRiInit(..., %d)" % \
                              self.n)
        finally:
            CxRiDelete(&ri)
        IF NjDebug:
            self._njDump()
        return additive

    # Iteratively try all clusterings of two rows in the matrix.  Do this in a
    # cache-friendly manner (keeping in mind that the matrix is stored in
    # row-major form).  This means:
    #
    # 1) For each row (x), find the row after it which is the closest (y),
    #    according to transformed distances, by calling _rnjRowMinFind().  This
    #    operation scans the row portion of the distances for x, which is a
    #    fast operation.
    #
    # 2) If the additivity constraint is enabled, check whether clustering x
    #    and y would violate additivity, by calling _rnjPairClusterAdditive().
    #
    # 2) Check whether it is okay to cluster x and y, by calling
    #    _rnjPairClusterOk().
    #
    # 3) If x and y can be clustered, do so, then immediately try to cluster
    #    with x again (as long as collapsing the matrix didn't move row x).
    cdef bint _rnjDeterministicCluster(self, bint additive) except -1:
        cdef bint clustered, done
        cdef size_t x, y
        cdef float distX, distY
        cdef Node node

        clustered = True
        done = False
        while not done:
            if not clustered:
                additive = False
            clustered = False
            x = 0
            while x < self.n - 1:
                y = self._rnjRowMinFind(x)

                if ((not additive) or self._rnjPairClusterAdditive(x, y)) \
                  and self._rnjPairClusterOk(x, y):
                    clustered = True
                    IF NjDebug:
                        self._njDump()
                    node = self._njNodesJoin(x, y, &distX, &distY)
                    self._njRSubtract(x, y)
                    self._njCompact(x, y, node, distX, distY)
                    self._njDiscard()
                    self.n -= 1
                    self._rScaledUpdate()

                    # Shrinking the matrix may have reduced it to the point
                    # that the enclosing loop will no longer function
                    # correctly.  Check this condition here, in order to reduce
                    # branch overhead for the case where no join is done.
                    if self.n == 2:
                        done = True
                        break

                    # The indexing of the matrix is shifted as a result of
                    # having removed the first row.  Set x such that joining
                    # with this row is immediately tried again.
                    #
                    # Note that if x is 0, then the row is now at (y - 1); in
                    # that case, stay on row 0.
                    if x > 0:
                        x -= 1
                else:
                    x += 1
        IF NjDebug:
            self._njDump()
        return additive

    cdef rnj(self, bint random, bint additive):
        assert self.tree is not None # Was self.prepare() called?

        if self.n > 2:
            if random:
                additive = self._rnjRandomCluster(additive)
            else:
                additive = self._rnjDeterministicCluster(additive)

        # Join last two nodes.
        IF NjDebug:
            self._njDump()
        self._njFinalJoin()

        return (self.tree, additive)
