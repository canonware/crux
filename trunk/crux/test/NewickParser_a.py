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

class testclass(crux.NewickParser.NewickParser):
    def __init__(self):
        pass

    def vprint(self, str):
        print "--%s--: %r" % (str, self.token())
        sys.stdout.flush()

    def openParenAccept(self):
        self.vprint("(")

    def closeParenAccept(self):
        self.vprint(")")

    def rootLabelAccept(self):
        self.vprint("rootLabel")

    def internalLabelAccept(self):
        self.vprint("internalLabel")

    def leafLabelAccept(self):
        self.vprint("leafLabel")

    def colonAccept(self):
        self.vprint(":")

    def lengthAccept(self):
        self.vprint("length")

    def commaAccept(self):
        self.vprint(",")

    def semicolonAccept(self):
        self.vprint(";")

    def commentAccept(self):
        self.vprint("[]")

    def whitespaceAccept(self):
        self.vprint(" ")

teststrs = (
    # Basic trees.
    "A;",
    "(A,B);",
    "(A,B,C);",
    "(A,(B,C));",
    "((A,B),(C,D));",

    # Branch lengths.
    "A:42;",
    "A:42.3;",
    "(A,B):42.3;",
    "(1,2,3:4.2):4;",

    # Comments.
    "[hi](A,B);",
    "([hi]A,B);",
    "(A,[hi]B);",
    "(A,B[hi]);",
    "(A,B)[hi];",

    "[[hi][bye]](A,B);",
    "[[hi] [bye]](A,B);",
    "[a[b[c[d[e[f]g]h]i]j]k](A,B);",

    # Whitespace.
    " (A,B);",
    "( A,B);",
    "(A ,B);",
    "(A, B);",
    "(A,B );",
    "(A,B :4.2);",
    "(A,B: 4.2);",
    "(A,B:4.2 );",
    "(A,B) ;",

    "(A,B [hi]);",
    "(A,B[hi] );",
    "(A,B [hi] );",
    "(A,B[hi] [bye]);",
    "(A,B[hi] [bye] );",
    "(A,B [hi] [bye] );",

    "  (A,B);",
    " \t (A,B);",
    "\t\r\n (A,B);",

    # Labels.
    ";",
    "A;",
    "(,);",
    "(,)label;",
    "(,);",
    "(,,);",
    "(,,,);",
    "(,(,));",
    "(,(,),);",
    "(,(,)label,);",
    "(A,B,);",
    "(A,B,,C,);",
    "(A,B,,C,,,);",

    # Branch lengths.
    "(A:0,B);",
    "(A:0123456789,B);",
    "(A:0.0,B);",
    "(B,A:0123456789.0123456789);",
    "(B,A:+1);",
    "(B,A:-1);",
    "(B,A:+3.42);",
    "(,A:-3.42);",
    "(,A:3e1);",
    "(,A:3.52e1);",
    "(,A:3.52e+1);",
    "(A:3.52e-1,);",
    "(A:3.52E1,);",
    "(A:3.52E+1,);",
    "(A:3.52E-1,);",

    # Unquoted labels.
    """(SomeCharacters!@#$%^*&.<>/?"\|-_=+`~{}whee,B);""",

    # Quoted labels.
    "('A quoted label.  Let''s embed a single quote.',B);",
    "(A,'',C);",
    
    # Test error conditions.
    "(A,B));",
    "((A,B);",
    "(A,B)",
    "(A,B:42.43.44);",
    "(A:42.);",

    "(A);",
    "(A,(B));",
    "(A,((B,C)));",
    "(A,((B,(C,D))));",

    # General tests.
    "( [hi]1_1[bye[hi again]bye for real]: 4.2,3,('4_4',5));",
    "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);",
    "(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;",
    "((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);",
    "(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460);",
    "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);",
    "(1,(((((((((((((((((((((2,(377,(378,(379,380)))),381),((353,(354,371)),(370,376))),(((345,(394,414)),(352,369)),351)),((372,409),(400,410))),(((350,361),((395,397),398)),(452,453))),(((((170,(346,347)),349),348),(359,360)),((((171,((386,387),((388,(389,390)),391))),(384,385)),(((182,((183,185),184)),((186,187),(188,189))),190)),((((((181,191),195),(192,(194,199))),356),(357,375)),(((196,198),383),197))))),((((((((((((3,442),441),444),443),445),438),((433,439),440)),((434,437),435)),(432,436)),358),(((((362,367),374),363),((364,366),365)),(((425,(426,(427,(428,429)))),(430,431)),(((454,464),((((((((456,466),457),459),461),462),(463,469)),(467,468)),458)),((455,465),460))))),(401,(((404,406),(446,450)),((447,(448,449)),451))))),((((325,(335,(336,337))),((((326,(330,333)),327),((328,(331,332)),329)),334)),(((338,339),340),(((341,343),342),344))),((417,420),419))),(((((((((((((((4,110),(6,109)),(111,112)),124),((145,146),147)),(154,(155,156))),(((131,132),137),((133,139),(134,(135,136))))),(((143,153),(149,152)),((144,151),(148,150)))),(((((((((5,(294,295)),(286,287)),((((240,((244,293),(254,(290,(311,312))))),((243,285),(289,(291,296)))),292),((313,314),316))),(((((((((233,(234,238)),((235,236),(237,241))),256),249),251),((((239,((245,250),260)),255),((246,((247,252),253)),248)),(281,282))),(272,(273,(275,276)))),((((257,(262,263)),(259,271)),(((((258,270),(268,269)),(266,267)),(261,264)),265)),(274,283))),((277,(278,279)),280))),((((242,319),320),318),(((304,(305,(308,315))),((306,307),309)),(310,317)))),(297,303)),288),((284,((298,(301,302)),300)),(((321,323),322),324))),299)),(((114,115),((128,129),130)),((((((((116,126),(118,127)),117),119),(120,121)),(123,142)),(122,((138,140),141))),125))),(113,(((((470,((((471,(472,473)),476),475),474)),478),((477,((479,480),((481,((483,484),(485,489))),482))),491)),((((486,492),(((488,494),(495,(496,497))),493)),487),490)),((498,500),499)))),(((158,((160,(161,((162,164),163))),168)),159),((165,167),166))),((157,392),(172,175))),((169,174),193)),((((((27,407),((68,(412,422)),(((69,(403,411)),(355,(((382,421),413),393))),402))),(((177,178),180),179)),173),176),70))),(232,368)),(((((48,61),62),405),66),((((60,65),63),67),((373,418),((408,423),416))))),((((((((200,(201,209)),(208,(229,231))),(203,219)),206),(202,(204,205))),(((((((210,((216,(220,222)),223)),((((211,212),226),214),(213,218))),217),215),225),224),(221,((227,228),424)))),230),207)),(((((((19,52),51),(40,(41,42))),(((49,(76,399)),((50,74),75)),72)),92),(((((((43,94),(((((95,97),(99,107)),(96,(98,104))),(((100,(102,105)),101),(103,106))),108)),73),45),44),(46,93)),396)),((64,71),90))),((59,415),91)),((((((((((((15,17),29),((((((77,88),82),((81,84),85)),(78,86)),79),80)),(24,26)),21),(((((((20,37),39),(((31,33),36),32)),(34,38)),30),(23,87)),83)),16),22),25),(18,35)),28),((47,53),((((54,55),57),56),58)))),(11,(12,(13,14)))),89),(9,10)),7),8));"
    
    )

test = testclass()

print "Test begin"

for str in teststrs:
    print "\n=== %r ===" % str
    try:
        test.parse(str)
    except:
        error = sys.exc_info()
        print "Exception %s: %s" % (error[0], error[1])

print "Test end"
