#!/bin/bash

tar -xzf R*.tar.gz
tar -xzf fits*.tar.gz

mkdir -p data/
spec="$(ls * | grep "[0-9]\{4\}\..*")"

tar -C data/ --strip 1 -xzf "${spec}"

# export paths and run R
export PATH="$PWD/R/bin:$PATH"
export RHOME="$PWD/R"
export R_LIBS="$PWD/packages"
Rscript ./lyman_sub.R

# tar -czf "plots-${spec%%.*}.tgz" plots/
mv -t . plots/* 2>/dev/null
