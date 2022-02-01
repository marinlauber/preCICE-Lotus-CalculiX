#
# set-up the FSI flag
export FSI=true
#
# make the exe
make -C $MGLHOME/src/geom libgeom.a
make -C $MGLHOME/src/oop  libfluid.a
make -f $MGLHOME/src/oop/Makefile lotus
echo made Lotus executable
echo

# run the exe
echo running
echo
if [ $1 -eq 0 ]; then
    if [ $# -eq 2 ]; then
        time ./lotus $3 &
    else
        time ./lotus &
    fi
else
    if [ $# -eq 2 ]; then
        time mpirun -n $1 ./lotus $2 &
    else
	    time mpirun -n $1 ./lotus &
    fi
fi

wait
exit 0