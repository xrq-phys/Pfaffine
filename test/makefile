INC=-I $(BLIS_ROOT)/include -I/usr/local/include
LIB=$(BLIS_ROOT)/lib/libblis.a
CXX=clang++ -std=c++11 $(INC)

simple_bli_ker_call.x: simple_bli_ker_call.cc
	$(CXX) $(LIB) simple_bli_ker_call.cc -o simple_bli_ker_call.x

small_dskr2k_call.x: small_dskr2k_call.cc ../src/skr2k.tcc
	$(CXX) $(LIB) small_dskr2k_call.cc -o small_dskr2k_call.x

skr2k.o: ../src/skr2k.tcc
	$(CXX) ../src/skr2k.tcc -o skr2k.o

exec: simple_bli_ker_call.x small_dskr2k_call.x
	./simple_bli_ker_call.x
	./small_dskr2k_call.x
