PROJ_ROOT=${WPDY2_ROOT}

# -- Directories --
SRC=${PROJ_ROOT}/src
BUILD=${PROJ_ROOT}/build
EXTERNAL=${PROJ_ROOT}/external
VPATH=${BUILD}:${SRC}

# -- common --
ifeq (${FC},gfortran)
	FF=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace -g -cpp -ffree-line-length-512 -fopenmp
	LDFLAGS=-llapack -lblas
endif
INCLUDE=-I${BUILD} -I${SRC}
MODS0=const sys math timer strutil

# -- command --
clean_all:
	rm -f ${BUILD}/*.o
	rm -f ${BUILD}/*.mod
	rm -f ${BUILD}/*.x

check_%: ${BUILD}/utest_%.x
	$<

# -- utility function
mod2obj = $(addprefix ${BUILD}/, $(addsuffix .o, $(1)))

# -- compile --
%.x:
	${FC} ${FF} $^ -o $@ -cpp ${LIBS} ${LDFLAGS} ${INCLUDE}

${BUILD}/%.o: ${SRC}/%.f90
	@if [ ! -d ${BUILD} ]; \
	   then mkdir -p ${BUILD}; \
	fi
	cd ${BUILD}; ${FC} ${FF} ${INCLUDE} -c $< -o $@

${BUILD}/fft4g.o:
	${FC} ${EXTERNAL}/fft/fft4g.f -c -o $@ 
${BUILD}/fftsg2d.o:
	${FC} ${EXTERNAL}/fft2d/fftsg2d.f -c -o $@

${BUILD}/con.o:
	cd ${BUILD}; ${FC} ${MAN4_ROOT}/src/con.f ${FF} ${INCLUDE} -c -o $@


# -- exe --
WPDYMODS=$(call mod2obj, ${MODS0} fft4g fftsg2d fft con wpdy wpdy_sop lambda)


