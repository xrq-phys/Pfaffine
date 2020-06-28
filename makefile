help:
	@echo "\"make lib\" to compile library."
	@echo "\"make install\" to compile library and install to subdirectory \"dest\"."

lib: make.inc
	cd src; make -f makefile

install: make.inc
	cd src; make -f makefile install
	# TODO: Install rule for tests?

clean: make.inc
	cd src; make -f makefile clean
	cd test; make -f makefile clean

tests: make.inc
	cd test; make -f makefile
