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

class taxon_map(object):
    def __init__(self, taxa=[]):
        self._label2ind = {}
        self._ind2label = {}

        # Construct dictionaries for lookups in both directions.
        i = 0
        for label in taxa:
            self._label2ind[label] = i
            self._ind2label[i] = label

            i += 1

    def ntaxa_get(self):
        return len(self._label2ind)

    def label_get(self, ind):
        if self._ind2label.has_key(ind):
            retval = self._ind2label[ind]
        else:
            retval = None

        return retval

    def ind_get(self, label):
        if self._label2ind.has_key(label):
            retval = self._label2ind[label]
        else:
            retval = None

        return retval

    def taxa_get(self):
        retval = self._ind2label.keys()
        retval.sort()
        i = 0
        for ind in retval:
            # Make sure that taxon indices are contiguous.
            if ind != i:
                raise ValueError

            retval[ind] = self._ind2label[ind]

            i += 1

        return retval

    def map(self, label, ind):
        self._label2ind[label] = ind
        self._ind2label[ind] = label
#EOF
