################################################################################
#
# <Copyright = jasone>
# <License>
#
################################################################################
#
# Version: Crux <Version = crux>
#
################################################################################

import crux

class Exception(crux.Exception):
    pass

class ValueError(Exception, ValueError):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message

class TaxonMap(object):
    def __init__(self, taxa=None):
        self._label2ind = {}
        self._ind2label = {}

        # Construct dictionaries for lookups in both directions.
        if taxa != None:
            for i in forints(len(taxa)):
                self._label2ind[taxa[i]] = i
                self._ind2label[i] = taxa[i]

    def ntaxaGet(self):
        return len(self._label2ind)

    def labelGet(self, ind):
        if self._ind2label.has_key(ind):
            rVal = self._ind2label[ind]
        else:
            rVal = None

        return rVal

    def indGet(self, label):
        if self._label2ind.has_key(label):
            rVal = self._label2ind[label]
        else:
            rVal = None

        return rVal

    def taxaGet(self):
        rVal = self._ind2label.keys()
        rVal.sort()
        i = 0
        for ind in rVal:
            # Make sure that taxon indices are contiguous.
            if ind != i:
                raise crux.TaxonMap.ValueError(
                    "Taxon indices must be contiguous")

            rVal[ind] = self._ind2label[ind]

            i += 1

        return rVal

    def equal(self, other):
        if other is self:
            rVal = True
        elif self.ntaxaGet() != other.ntaxaGet():
            rVal = False
        else:
            keys = self._ind2label.keys()
            keys.sort()
            for key in keys:
                if other.labelGet(key) == None:
                    rVal = False
                    break
            rVal = True

        return rVal

    # Map a label to an index.  Typical usage is something like:
    #
    #   m.map('Label', m.ntaxaGet())
    def map(self, label, ind, replace=False):
        if replace == False:
            # Make sure that label hasn't already been mapped to an index.
            if self._label2ind.has_key(label):
                raise crux.TaxonMap.ValueError(
                    "Label '%s' already mapped" % label)

        self._label2ind[label] = ind
        self._ind2label[ind] = label
#EOF
