#!/usr/bin/python
# retrieve zsim's cache stats
#

import h5py # presents HDF5 files as numpy arrays
import numpy as np

# Open stats file
f = h5py.File('zsim-ev.h5', 'r')

# Get the single dataset in the file
dset = f["stats"]["root"]

print "=== L2 Cache Stats ==="

# Each dataset is first indexed by record. A record is a snapshot of all the
# stats taken at a specific time.  All stats files have at least two records,
# at beginning (dest[0])and end of simulation (dset[-1]).  Inside each record,
# the format follows the structure of the simulated objects. A few examples:

# Phase count at end of simulation
endPhase = dset[-1]['phase']
print "Phase count: " + str(endPhase)

# Hits into all L2s, which actually is a single one
l2_hits = np.sum(dset[-1]['l2']['hGETS'] + dset[-1]['l2']['hGETX'])

# Misses into all L2s
l2_misses = np.sum(dset[-1]['l2']['mGETS'] + dset[-1]['l2']['mGETXIM'] + dset[-1]['l2']['mGETXSM'])

# Total number of instructions executed, counted by adding per-core counts
# (you could also look at procInstrs)
totalInstrs = np.sum(dset[-1]['simpleCore']['instrs'])

l2_accesses = l2_hits + l2_misses

print "Total instructions: " + str(totalInstrs)
print "Accesses: " + str(l2_accesses) + "; Hits: " + str(l2_hits) + "; Misses: " + str(l2_misses)
print "Miss rate: " + str(1.0 * l2_misses / l2_accesses)
