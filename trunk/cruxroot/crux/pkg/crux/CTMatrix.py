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

import TaxonMap

import crux.Exception

class Exception(crux.Exception):
    pass

class ValueError(Exception, ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

class CTMatrix(object):
    def __init__(self, map=TaxonMap.TaxonMap()):
        # This is an array of CharacterType objects, and stores the character
        # types for the columns in _taxonData.
        self._chars = []

        # This is used to manage the rows of the data matrix (_taxonData).
        self._taxonMap = map

        # Row-major character data.  Each element in _taxonData is a string of
        # characters that belong to the corresponding taxon in _taxonMap.  The
        # keys are the integer indexes, as reported by self._taxonMap.indGet().
        self._taxonData = {}

    # Append CharacterType objects to _chars.
    def charsAppend(self, chars):
        self._chars += chars

    # Return the array of character types.
    def charsGet(self):
        return self._chars

    # Return the taxon map.
    def taxonMapGet(self):
        return self._taxonMap

    # Return the character data for a taxon.
    def dataGet(self, taxon):
        if not self._taxonData.has_key(self._taxonMap.indGet(taxon)):
            retval = None
        else:
            retval = self._taxonData[self._taxonMap.indGet(taxon)]

        return retval

    # Set the character data for a taxon.
    def dataSet(self, taxon, data):
        if len(data) != len(self._chars):
            raise crux.CTMatrix\
                  .ValueError("Wrong number of characters (%d, expected %d)"
                              % (len(data), len(self._chars)))
        if (self._taxonMap.indGet(taxon) == None):
            raise crux.CTMatrix\
                  .ValueError("Taxon %r not in taxon map" % taxon)

        self._taxonData[self._taxonMap.indGet(taxon)] = data
#EOF
