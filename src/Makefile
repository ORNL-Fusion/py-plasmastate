# FOR CORI
# module load python swig cray-netcdf 

SWIG = swig
SWIG_FLAGS = -python -c++
CXX = icc
CFLAGS = -g -c -fPIC
LIBBLAS = 

NTCCDIR = /global/common/software/atom/cori/share/ntcc_intel

NETCDFDIR = ${NETCDF_DIR}/lib
LDFLAGS = -shared -fPIC
LIBS = -L. \
    -L${NTCCDIR}/lib \
    -L${NETCDFDIR} \
	-lccps_lib -lplasma_state -lps_xplasma2 -lplasma_state_kernel -lxplasma2 \
    -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch -lfluxav -lr8bloat -lpspline -lezcdf \
    -llsode -llsode_linpack -lsmlib -lcomput -lportlib \
    -lnetcdff -lnetcdf \
    -L${PYTHON_DIR}/lib -lpython3.7m \
    ${LIBBLAS}

INCDIR = -I${NTCCDIR}/include -I/usr/include/python3.6m
INCDIR = -I${NTCCDIR}/include_ntcc -I${PYTHON_DIR}/include/python3.7m

all: _pyplasmastate.so

pyplasmastate_wrap.cxx: pyplasmastate.i 
	${SWIG} ${SWIG_FLAGS} $? 

pyplasmastate_wrap.o: pyplasmastate_wrap.cxx 
	${CXX} ${CFLAGS} ${INCDIR} $?

_pyplasmastate.so: pyplasmastate_wrap.o
	${CXX} ${LDFLAGS} $? ${LIBS} -o $@

clean:
	rm -rf *.o *.so __pycache__ pyplasmastate.py pyplasmastate_wrap.cxx
