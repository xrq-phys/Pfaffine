INC=-I $(BLIS_ROOT)/include -I/usr/local/include
LIB=-L $(BLIS_ROOT)/lib -lblis
CXX=g++ -std=c++11 $(INC) -g

simple_bli_ker_call.x: simple_bli_ker_call.cc
	$(CXX) $(LIB) simple_bli_ker_call.cc -o simple_bli_ker_call.x

small_dskr2k_call.x: small_dskr2k_call.cc ../src/skr2k.tcc
	$(CXX) $(LIB) small_dskr2k_call.cc -o small_dskr2k_call.x

skr2k.o: ../src/skr2k.tcc
	$(CXX) ../src/skr2k.tcc -o skr2k.o

exec: simple_bli_ker_call.x small_dskr2k_call.x
	LD_LIBRARY_PATH=$(BLIS_ROOT)/lib ./simple_bli_ker_call.x
	LD_LIBRARY_PATH=$(BLIS_ROOT)/lib ./small_dskr2k_call.x
