help:
	@echo "\"make lib\" to compile library."
	@echo "\"make install\" to compile library and install to subdirectory \"dest\"."

lib: make.inc
	make -C src

install: make.inc
	make -C src install
	# TODO: Install rule for tests?

clean: make.inc
	make -C src clean
	make -C test clean

tests: make.inc
	make -C test
