simple_bli_ker_call.x: simple_bli_ker_call.cc
	clang++ simple_bli_ker_call.cc -o simple_bli_ker_call.x -std=c++11 -I ~/opt/BLIS/include -I/usr/local/include ~/opt/BLIS/lib/libblis.a

exec: simple_bli_ker_call.x
	./simple_bli_ker_call.x
