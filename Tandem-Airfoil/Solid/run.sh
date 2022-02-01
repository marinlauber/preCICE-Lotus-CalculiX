#!/bin/sh
set -e -u

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1

# making sure the solver is uptodate
# make -C $CCXHOME

# run the exe
echo running
echo

if [ $# -eq 0 ]; then
    ccx_preCICE -i calculix -bg -precice-participant Calculix &
else
    mv pressure.rout restart.rin &
    ccx_preCICE -i $1 -bg -precice-participant Calculix &
fi

wait
exit 0
