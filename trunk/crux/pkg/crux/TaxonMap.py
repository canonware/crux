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
    def __init__(self, taxa=[]):
        self._label2ind = {}
        self._ind2label = {}

        # Construct dictionaries for lookups in both directions.
        i = 0
        for label in taxa:
            self._label2ind[label] = i
            self._ind2label[i] = label

            i += 1

    def ntaxaGet(self):
        return len(self._label2ind)

    def labelGet(self, ind):
        if self._ind2label.has_key(ind):
            retval = self._ind2label[ind]
        else:
            retval = None

        return retval

    def indGet(self, label):
        if self._label2ind.has_key(label):
            retval = self._label2ind[label]
        else:
            retval = None

        return retval

    def taxaGet(self):
        retval = self._ind2label.keys()
        retval.sort()
        i = 0
        for ind in retval:
            # Make sure that taxon indices are contiguous.
            if ind != i:
                raise crux.TaxonMap.ValueError(
                    "Taxon indices must be contiguous")

            retval[ind] = self._ind2label[ind]

            i += 1

        return retval

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
