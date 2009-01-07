import random
import time

import Crux

version = "@crux_version@"

prefix = "@PREFIX@"
bindir = "@BINDIR@"
datadir = "@DATADIR@"
mandir = "@MANDIR@"

debug = (True if @enable_debug@ else False)
mpi = (True if @enable_mpi@ else False)

# Defaults for when Crux is imported by an application other than the crux
# front end.
batch = True
verbose = True
Crux.seed(int(time.time() * 1000000) & 0xffffffff) # Sets Config.seed.
infile = None
scriptargs = []
