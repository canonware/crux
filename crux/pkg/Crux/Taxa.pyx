from Cx cimport CxCmp2Richcmp

import exceptions

class ValueError(exceptions.ValueError):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message

import weakref

cdef class Taxon:
    def __init__(self, str label):
        self.label = label

    def __cmp__(self, Taxon other):
        return cmp(self.label, other.label)

    def __hash__(self):
#        return self.label.__hash__()
        return <int>(id(self) & 0x7fffffff) # Fragile, but faster than above.

    def __repr__(self):
        return "Crux.Taxon.taxon('%s')" % self.label

cdef _taxa
_taxa = weakref.WeakValueDictionary()

cpdef Taxon get(str label):
    """
        Return a globally unique Taxon instance associated with a label string.
    """
    global _taxa
    cdef Taxon ret

    ret = _taxa.get(label)
    if ret is None:
        ret = Taxon(label)
        _taxa[label] = ret
    return ret

cdef class Map:
    """
        Map between Taxon instances and indices.  This class is useful when
        taxa ordering matters, as for CTMatrix.Alignment and Tree.Bipart.Vec.
    """
    def __init__(self, list taxa=None):
        cdef int i
        cdef Taxon taxon

        self._taxon2ind = {}
        self._ind2taxon = {}

        # Construct dictionaries for lookups in both directions.
        if taxa is not None:
            for 0 <= i < len(taxa):
                taxon = taxa[i]
                self._taxon2ind[taxon] = i
                self._ind2taxon[i] = taxon

        self.ntaxa = len(self._ind2taxon)

    def __richcmp__(Map self, Map other, int op):
        if other is self:
            return CxCmp2Richcmp(0, op)
        elif self.ntaxa != other.ntaxa:
            return CxCmp2Richcmp(1, op)
        else:
            try:
                for 0 <= i < self.ntaxa:
                    if self.taxonGet(i) is not self.taxonGet(i):
                        return CxCmp2Richcmp(1, op)
                return CxCmp2Richcmp(0, op)
            except:
                raise ValueError("Map indices must be contiguous")

    cpdef Taxon taxonGet(self, int ind):
        """
            Get taxon associated with index.
        """
        try:
            return <Taxon>self._ind2taxon[ind]
        except:
            return None

    cpdef int indGet(self, Taxon taxon):
        """
            Get index associated with taxon.
        """
        try:
            return self._taxon2ind[taxon]
        except:
            return -1

    cdef list taxaGet(self):
        cdef list ret
        cdef int i

        ret = []
        try:
            for 0 <= i < len(self._ind2taxon):
                ret.append(self._ind2taxon[i])
        except:
            raise ValueError("Map indices must be contiguous")

        return ret
    property taxa:
        """
            Ordered list of taxa in map.
        """
        def __get__(self):
            return self.taxaGet()

    cpdef map(self, Taxon taxon, int ind, bint replace=False):
        """
            Map a taxon to an index.  Typical usage is something like:

              m.map(taxon, m.ntaxa)
        """
        if replace == False:
            # Make sure that taxon hasn't already been mapped to an index.
            if self._taxon2ind.has_key(taxon):
                raise ValueError("Taxon %r already mapped" % taxon)

        self._taxon2ind[taxon] = ind
        self._ind2taxon[ind] = taxon
        self.ntaxa = len(self._ind2taxon)
