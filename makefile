help:
	@echo "\"make lib\" to compile library."
	@echo "\"make install\" to compile library and install to subdirectory \"dest\"."

lib: make.inc
	cd src; make -f makefile

sandy: make.inc
	cd src; make -f makefile sandy

install: make.inc
	cd src; make -f makefile install

clean: make.inc
	cd src; make -f makefile clean

