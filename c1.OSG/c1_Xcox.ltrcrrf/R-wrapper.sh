#!/bin/bash

# Uncompress the tarball
tar -xzf R-packages.tar.gz

# Set the library location
export R_LIBS="$PWD/R-packages"
# set TMPDIR variable
export TMPDIR=$_CONDOR_SCRATCH_DIR

# run the R program
Rscript main_boot.ltrcrrf.R
