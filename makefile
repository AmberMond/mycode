compile:
	g++ -O2 -std=c++17 src/main.cpp src/Mesh.cpp src/Solver.cpp -Iinclude -o 2d_ns

run:
	.\2d_ns.exe

clean: 
	rm -f 2d_ns *.exe *.o

