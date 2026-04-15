compile:
	g++ -O2 -std=c++17 main.cpp Mesh.cpp Solver.cpp -o 2d_ns

#      ./dg_quad_nodal [p] [Nx] [Ny] [t_end] [CFL] [benchmark] [flux] [limiter_policy]
#      ./dg_quad_nodal 2 40 40 0.05 0.05 densitywave rusanov new
#      ./dg_quad_nodal 2 40 40 0.05 0.05 densitywave hllc new
#      ./dg_quad_nodal 1 64 64 0.15 0.02 quadrant rusanov new
#      ./dg_quad_nodal 1 64 64 0.15 0.02 quadrant hllc new

run:
	.\2d_ns.exe

