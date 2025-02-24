#!/bin/bash
set -eux

#export PREFIX="${PREFIX:-${CONDA_PREFIX}}"

export CMAKE_PREFIX_PATH=${PREFIX}

#Set library path
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="${LIBRARY_PATH:-}"

# Set compilers
export FC=$(which mpif90)
export CC=$(which mpicc)
export CXX=$(which mpicxx)

# Clone and build SciFortran
git clone https://github.com/SciFortran/SciFortran.git scifor
cd scifor
mkdir build
cd build
cmake .. -DPREFIX=${PREFIX}
make
make install
cd ../../

for d in ${PREFIX}/scifor/gnu/*/etc; do
    if [ -d "$d" ]; then
        export PKG_CONFIG_PATH=${d}:${PKG_CONFIG_PATH}
    fi
done


export GLOB_INC=$( pkg-config --cflags scifor )
export GLOB_LIB=$( pkg-config --libs   scifor  | sed  "s/;/ /g"  | sed 's/\\/  /g' )

# Build EDIpack
git clone https://github.com/edipack/edipack2.0.git edipack2
cd edipack2
mkdir build
cd build
cmake .. -DPREFIX=${PREFIX}
make
make install
cd ../


for d in ${PREFIX}/edipack2/gnu/*/etc; do
    if [ -d "$d" ]; then
        export PKG_CONFIG_PATH=${d}:${PKG_CONFIG_PATH}
    fi
done
export GLOB_INC=$( pkg-config --cflags scifor edipack2)
export GLOB_LIB=$( pkg-config --libs   scifor edipack2 | sed  "s/;/ /g"  | sed 's/\\/  /g' )

$PYTHON -m pip install . --prefix=${PREFIX} --no-deps --ignore-installed
