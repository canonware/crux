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

import _Ring

class Exception(_Ring.Exception):
    pass

class ValueError(Exception, _Ring.ValueError):
    pass

class Ring(_Ring.Ring):
    def __init__(self, edge, end):
        pass
#EOF
