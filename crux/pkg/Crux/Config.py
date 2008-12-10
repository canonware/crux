import random
import time

import Crux

version = "@crux_version@"

prefix = "@PREFIX@"
bindir = "@BINDIR@"
datadir = "@DATADIR@"
mandir = "@MANDIR@"

debug = (True if @enable_debug@ else False)

# Defaults for command line settings.
batch = False
quiet = False
Crux.seed(int(time.time())) # Sets Config.seed.
infile = None
scriptargs = []
