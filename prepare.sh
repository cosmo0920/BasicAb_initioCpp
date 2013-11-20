#!/usr/bin/bash
git submodule init
cd include
ln -s ../clapack/INCLUDE/clapack.h  .
cd ..
