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

import crux.Exception

class Exception(crux.Exception):
    pass

class ValueError(Exception, ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

class CTMatrix(object):
    def __init__(self):
        pass
