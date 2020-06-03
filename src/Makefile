SWIG = swig
SWIG_FLAGS = -python -c++
CXX_GCC = g++
CXX = g++
CFLAGS = -g -c -fPIC 

ATOM = /project/projectdirs/atom/atom-install-edison
LOCAL = ${ATOM}/cesol

#NTCCDIR =  ${LOCAL}/local/ntcc_intel
#NTCCDIR =  /home/eiw/Projects/Docker/PS/install
#NTCCDIR : Parent directory for PPPL Plasma State and NTCC Libraries
NTCCDIR =  ../../install

NETCDFDIR = ${NETCDF_DIR}/lib
LDFLAGS = -shared -fPIC
LIBS = -L. \
    -L${NTCCDIR}/lib \
    -L${NETCDFDIR} \
	-lccps_lib -lplasma_state -lps_xplasma2 -lplasma_state_kernel -lxplasma2 \
    -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch -lfluxav -lr8bloat -lpspline -lezcdf \
    -llsode -llsode_linpack -lsmlib -lcomput -lportlib \
    -lnetcdff -lnetcdf \
    -lpython3.6m -lopenblas

INCDIR = -I${NTCCDIR}/include -I/usr/include/python3.6m

all: _pyplasmastate.so

pyplasmastate_wrap.cxx: pyplasmastate.i 
	${SWIG} ${SWIG_FLAGS} $? 

pyplasmastate_wrap.o: pyplasmastate_wrap.cxx 
	${CXX} ${CFLAGS} ${INCDIR} $?

_pyplasmastate.so: pyplasmastate_wrap.o
	${CXX} ${LDFLAGS} $? ${LIBS} -o $@

clean:
	rm -rf *.o *.so __pycache__ pyplasmastate.py pyplasmastate_wrap.cxx
