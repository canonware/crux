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

from C_NewickParser import *

class NewickParser(C_NewickParser):
    def __init__(self):
        pass

    # Parse input and call the appropriate *Accept method for each token that is
    # accepted.
    def parse(self, input):
        C_NewickParser._parse(self, input)

    def openParenAccept(self):
        # Virtual method.
        pass

    def closeParenAccept(self):
        # Virtual method.
        pass

    def rootLabelAccept(self):
        # Virtual method.
        pass

    def internalLabelAccept(self):
        # Virtual method.
        pass

    def leafLabelAccept(self):
        # Virtual method.
        pass

    def colonAccept(self):
        # Virtual method.
        pass

    def lengthAccept(self):
        # Virtual method.
        pass

    def commaAccept(self):
        # Virtual method.
        pass

    def semicolonAccept(self):
        # Virtual method.
        pass

    def commentAccept(self):
        # Virtual method.
        pass

    def whitespaceAccept(self):
        # Virtual method.
        pass
#EOF
