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

import sys
import tempfile

class testclass(crux.FastaParser.FastaParser):
    def __init__(self):
        pass

    def vprint(self, str):
        print "--%s--: %r" % (str, self.token())
        sys.stdout.flush()

    def labelAccept(self):
        self.vprint("label")

    def commentAccept(self):
        self.vprint("comment")

    def charsAccept(self):
        self.vprint("chars")

dnaTestStrs = (
    """>Label
NXVHMDRWABSYCKGT-
""",

    """>Label
A C
GT""",

    """
>Taxon_A
ACGT-
>Taxon_B
ACG-T
""",

    """
>Taxon_A This is a comment.
ACGT-
>Taxon_B	This is another comment.
ACG-T
""",

    # Error cases.
    """
>
""",

    """
> Comment
""",

    """
>Label
""",

    """
>Label
ACGTZ
""",

    """
>Taxon_A
>Taxon_B
ACGT
""",

    """
>Taxon_A
ACGT>Taxon_B
ACGT
""",
    
    """
>Taxon_A
ACGT
>Taxon_B
""",
    
    )

proteinTestStrs = (
    """>Label
ABCDEFGHIKLMNPQRSTUVWXYZ-
""",

    # Error cases.
    """
>Label
ABCDEFGHIJKLMN
""",

    """
>Label
ABCDE789FGHIK
""",

    )

def testStrs(strs, charType):
    for str in strs:
        print "\n=== %r ===" % str
        try:
            test.parse(str, charType)
        except crux.Exception:
            error = sys.exc_info()
            print "Exception %s: %s" % (error[0], error[1])

    for str in strs:
        print "\n=== %r ===" % str
        try:
            f = tempfile.TemporaryFile()
            f.write(str)
            f.seek(0, 0)
            test.parse(f, charType)
        except crux.Exception:
            error = sys.exc_info()
            print "Exception %s: %s" % (error[0], error[1])

test = testclass()

print "Test begin"

testStrs(dnaTestStrs, "DNA")
testStrs(proteinTestStrs, "protein")

print "Test end"