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

import _Edge

class Exception(_Edge.Exception):
    pass

class ValueError(Exception, _Edge.ValueError):
    pass

class Edge(_Edge.Edge):
    def __init__(self, tree):
        pass
#EOF
