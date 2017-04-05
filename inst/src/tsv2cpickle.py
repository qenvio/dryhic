#!/usr/bin/env python

import sys
import cPickle as pickle

args = sys.argv
infile = args[1]
outfile = args[2]

with open(infile, 'r') as f:
    read_data = f.read().splitlines()

dat = {}

for i in read_data:
    a = i.split("\t")
    dat[int(a[0])] = float(a[1])

pickle.dump(dat, open(outfile, "w"))
