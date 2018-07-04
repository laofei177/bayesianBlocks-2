# bayesianBlocks
This is a naive translation into Matlab
from the Julia Implementation by Michael P.H. Stumpf and T. Chan
Based on the Python Code of Jake Vanderplas.

This version implements Bayesian blocks for histograms, where the
data are sorted, then treated as event data (see Scargle 2012).
Defaults as suggested in Scargle 2012 are used.

References:
Scargle 2012: http://adsabs.harvard.edu/abs/2012arXiv1207.5578S
Python implementation: https://github.com/astroML/astroML/blob/master/astroML/density_estimation/bayesian_blocks.py
Julia implementation: https://github.com/sisl/Discretizers.jl/blob/master/src/disc_bayesianblocks.jl