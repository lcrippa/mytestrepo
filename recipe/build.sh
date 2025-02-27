#!/bin/bash
set -eux

# Set compilers
export FC=$(which mpif90)
export CC=$(which mpicc)
export CXX=$(which mpicxx)

#Create pkg-config directory, if it doesn't exist
mkdir -p ${PREFIX}/lib/pkgconfig

# Ensure the unlink.d directory exists and copy post-uninstall script
mkdir -p ${PREFIX}/etc/conda/unlink.d
chmod +x ${RECIPE_DIR}/post-unlink.sh
cp ${RECIPE_DIR}/post-unlink.sh ${PREFIX}/etc/conda/unlink.d/

# Clone and build SciFortran
git clone https://github.com/SciFortran/SciFortran.git scifor
cd scifor
mkdir build
cd build
cmake .. -DPREFIX=${PREFIX}/opt
make -j
make install
cd ../../

for d_scifor in ${PREFIX}/opt/scifor/gnu/*/etc; do
    if [ -d "$d_scifor" ]; then
        export PKG_CONFIG_PATH=${d_scifor}:${PKG_CONFIG_PATH}
        cp ${d_scifor}/scifor.pc ${PREFIX}/lib/pkgconfig/
    fi
done

export GLOB_INC=$( pkg-config --cflags scifor )
export GLOB_LIB=$( pkg-config --libs   scifor  | sed  "s/;/ /g"  | sed 's/\\/  /g' )

# Clone and build EDIpack
git clone https://github.com/edipack/edipack2.0.git edipack2
cd edipack2
mkdir build
cd build
cmake .. -DPREFIX=${PREFIX}/opt
make -j
make install
cd ../


for d_edipack in ${PREFIX}/opt/edipack2/gnu/*/etc; do
    if [ -d "$d_edipack" ]; then
        export PKG_CONFIG_PATH=${d_edipack}:${PKG_CONFIG_PATH}
        cp ${d_edipack}/edipack2.pc ${PREFIX}/lib/pkgconfig/
    fi
done
export GLOB_INC=$( pkg-config --cflags scifor edipack2)
export GLOB_LIB=$( pkg-config --libs   scifor edipack2 | sed  "s/;/ /g"  | sed 's/\\/  /g' )


#Replace ambigouos variable name in the .pc files
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS (BSD sed) requires an explicit '' argument for in-place editing
    sed -i '' "1s|prefix|edipack_prefix|g" "${PREFIX}/lib/pkgconfig/edipack2.pc"
    sed -i '' "s|\${prefix}|\${edipack_prefix}|g" "${PREFIX}/lib/pkgconfig/edipack2.pc"
    sed -i '' "1s|prefix|scifor_prefix|g" "${PREFIX}/lib/pkgconfig/scifor.pc"
    sed -i '' "s|\${prefix}|\${scifor_prefix}|g" "${PREFIX}/lib/pkgconfig/scifor.pc"
else
    # Linux (GNU sed)
    sed -i "1s|prefix|edipack_prefix|g" "${PREFIX}/lib/pkgconfig/edipack2.pc"
    sed -i "s|\${prefix}|\${edipack_prefix}|g" "${PREFIX}/lib/pkgconfig/edipack2.pc"
    sed -i "1s|prefix|scifor_prefix|g" "${PREFIX}/lib/pkgconfig/scifor.pc"
    sed -i "s|\${prefix}|\${scifor_prefix}|g" "${PREFIX}/lib/pkgconfig/scifor.pc"
    # Fix for mkl linking
    sed -i 's|-L/usr/lib/x86_64-linux-gnu||g' "${PREFIX}/lib/pkgconfig/scifor.pc"
fi


#Install edipy
$PYTHON -m pip install . --prefix=${PREFIX} --no-deps --ignore-installed
