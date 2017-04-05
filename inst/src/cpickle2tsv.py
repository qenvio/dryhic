#!/usr/bin/env python

import sys
import cPickle as pickle

args = sys.argv
infile = args[1]
outfile = args[2]

dat = pickle.load(open(infile, "r"))

with open(outfile, 'w') as f:
    for i in dat.keys():
        f.write(str(i) + "\t" + str(dat[i]) + "\n")

