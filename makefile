help:
	@echo "\"make lib\" to compile library."
	@echo "\"make test\" test."

ext/pfapack2/fortran/make.inc: make.inc
	cd ext/pfapack2/fortran && ln -s ../../../make.inc . && cd ../../..

ALL_INC = make.inc ext/pfapack2/fortran/make.inc

lib: $(ALL_INC)
	make -C ext/pfapack2/fortran
	make -C ltl2inv

clean: $(ALL_INC)
	make -C ltl2inv clean
	make -C test clean

tests: $(ALL_INC)
	make -C test
