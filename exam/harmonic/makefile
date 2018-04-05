include ../../local.mk
include ../../conf.mk

TARGET=run.x

${TARGET}: ${WPDYMODS} ${TARGET:%.x=%.f90}
all: ${TARGET}
	./$< > run.log

clean:
	rm -rf *.mod
	rm -rf _con
	rm -rf run.x
	rm -rf *.log

