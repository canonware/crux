import random
import time

version = "@crux_version@"

prefix = "@PREFIX@"
bindir = "@BINDIR@"
datadir = "@DATADIR@"
mandir = "@MANDIR@"

debug = (True if @enable_debug@ else False)

# Defaults for command line settings.
batch = False
quiet = False
seed = int(time.time())
random.seed(seed)
infile = None
scriptargs = []
