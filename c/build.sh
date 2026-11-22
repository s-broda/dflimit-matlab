#!/bin/sh
# builds the C engine and its tests (unix); CFLAGS as in .github/workflows/c.yml
set -e
cd "$(dirname "$0")"
CFLAGS="-std=c11 -O2 -Wall -Wno-unused-function -Wno-logical-op-parentheses ${EXTRA_CFLAGS:-}"
# (the classical engine, then the general engine, dfl_det_*)
for f in dfl_special dfl_mgf dfl_single dfl_moments dfl_double dfl_gilpelaez dfl_api dfl_threads tables \
         dfl_det_la dfl_det_design dfl_det_mgf dfl_det_path dfl_det_inner dfl_det_api dfl_det_paper2 dfl_det_route dfl_det_line dfl_det_tables; do cc $CFLAGS -c $f.c -o $f.o; done
cc $CFLAGS -O3 -c dfl_gen.c -o dfl_gen.o            # -O3: the l-rule of the noncentral kernel is written to be vectorized
cc $CFLAGS -O3 -c dfl_batch.c -o dfl_batch.o        # (and the sums over the nodes of the fixed layouts)
c++ -std=c++11 -O2 -c faddeeva/Faddeeva.cc -o Faddeeva.o
c++ -std=c++11 -O2 -c faddeeva/faddeeva_shim.cpp -o faddeeva_shim.o
OBJ="dfl_special.o dfl_mgf.o dfl_single.o dfl_batch.o dfl_gen.o dfl_moments.o dfl_double.o dfl_gilpelaez.o dfl_api.o dfl_threads.o tables.o Faddeeva.o faddeeva_shim.o"
OBJ="$OBJ dfl_det_la.o dfl_det_design.o dfl_det_mgf.o dfl_det_path.o dfl_det_inner.o dfl_det_api.o dfl_det_paper2.o dfl_det_route.o dfl_det_line.o dfl_det_tables.o"
for t in ${TESTS:-prim single batch gen fixed vmath mom dbl app bounds det}; do cc $CFLAGS -c dfl_test_$t.c -o dfl_test_$t.o; c++ -o dfl_test_$t dfl_test_$t.o $OBJ -lm -lpthread; done
