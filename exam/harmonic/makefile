include ${HOME}/src/wpdy2/local.mk
include ${WPDY2_ROOT}/conf.mk

TARGET=run.x

${TARGET}: ${WPDYMODS} ${TARGET:%.x=%.f90}
all: ${TARGET}
	./$< > run.log
	python plot.py
	imgcat rho.pdf

clean:
	rm -rf *.mod
	rm -rf _con
	rm -rf run.x
	rm -rf *.log
	rm -rf *.pdf
