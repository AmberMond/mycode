#include <iostream>
#include <fstream>
#include <string>
#include "Mesh.hpp"
#include "Solver.hpp"

// vtk output
void saveResults(const Mesh& mesh, const std::vector<State>& U, int step) {
    std::string filename = "sol_" + std::to_string(step) + ".vtk";
    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "# vtk DataFile Version 3.0\n2D Bump\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    
    // Points
    file << "POINTS " << mesh.nodes.size()-1 << " double\n";
    for(size_t i=1; i<mesh.nodes.size(); ++i) 
        file << mesh.nodes[i].x << " " << mesh.nodes[i].y << " 0.0\n";
    
    // Cells
    int listSize = 0;
    for(const auto& c : mesh.cells) listSize += (c.nodes.size() + 1);
    file << "CELLS " << mesh.cells.size() << " " << listSize << "\n";
    for(const auto& c : mesh.cells) {
        file << c.nodes.size();
        for(int n : c.nodes) file << " " << n-1;
        file << "\n";
    }
    
    file << "CELL_TYPES " << mesh.cells.size() << "\n";
    for(const auto& c : mesh.cells) {
        file << (c.nodes.size() == 3 ? "5\n" : "9\n");
    }

    // Data
    file << "CELL_DATA " << mesh.cells.size() << "\n";
    
    file << "SCALARS Mach double 1\nLOOKUP_TABLE default\n";
    for(const auto& u : U) {
        double p = getPressure(u);
        double c = getSoundSpeed(p, u.rho);
        double v = std::sqrt(u.rhou*u.rhou + u.rhov*u.rhov)/u.rho;
        file << v/c << "\n";
    }

    file << "SCALARS Pressure double 1\nLOOKUP_TABLE default\n";
    for(const auto& u : U) file << getPressure(u) << "\n";

    file << "SCALARS Density double 1\nLOOKUP_TABLE default\n";
    for(const auto& u : U) file << u.rho << "\n";
}
// for 2d bump
// int main() {
//     Mesh mesh;
//     // mesh
//     mesh.load("bump_1230_t8"); 

//     Solver solver(mesh);

//     // Ma = 0.5, AoA = 0
//     // P_static = 101325 Pa, T_static = 300 K

//     solver.setFlowConditions(0.5, 101325.0, 300.0);
    
//     solver.initialize();

//     // explicit time stepping
//     // L and c for calculate CFL number (explicit)
//     // dt ~ L / c ~ 0.01 / 340 ~ 3e-5
    
//     // double dt = 1e-5; 
//     // int maxSteps = 5000;

//     // std::cout << "Starting simulation..." << std::endl;
//     // for (int i = 0; i <= maxSteps; ++i) {
//     //     solver.timeStepExplicit(dt);
      
//     //     if (i % 100 == 0) {
//     //         std::cout << "Step " << i << std::endl;
//     //         saveResults(mesh, solver.getSolution(), i);
//     //     }
//     // }

//     std::cout << "JFNK" << std::endl;
//     solver.solveImplicit(1e-5, 2500);
//     saveResults(mesh, solver.getSolution(), 9999);

//     std::cout << " Grs " << std::endl;
//     solver.calculateEntropyError();

//     std::cout << " Drag " << std::endl;
//     solver.calculateDrag();

//     std::cout << "Done" << std::endl;



//     return 0;
//     }



// for wind tunnel

// int main() {
//     Mesh mesh;
//     // read mesh
//     mesh.load("wind_tunnel_4t"); 

//     Solver solver(mesh);

//     // 2. initialize
//     solver.initialize();

//     // 3. explicit time stepping
    
//     double dt = 5e-4;
//     int maxSteps = 20000; 

//     std::cout << "Starting Explicit Mach 3 Simulation..." << std::endl;
//     for (int i = 0; i <= maxSteps; ++i) {
//         solver.timeStepExplicit(dt);
      
        
//         if (i % 100 == 0) {
//             std::cout << "Step " << i << " | Physical Time: " << i * dt << std::endl;
//             saveResults(mesh, solver.getSolution(), i);
//         }
//     }

//     std::cout << "Done" << std::endl;

//     return 0;
// }


// for TGV

// int main() {
//     Mesh mesh;
//     mesh.load("tgv2x");

//     Solver solver(mesh);
//     solver.initializeTGV(); 

//     double dt_phy = 5e-4; 
//     int maxSteps = 300;

//     std::cout << "Starting Unsteady JFNK TGV Simulation..." << std::endl;
    
//     for (int i = 0; i < maxSteps; ++i) {
       
//         solver.solveUnsteady(dt_phy, 1);
//         std::cout << "Step: " << i + 1 << " Done." << std::endl;
//         saveResults(mesh, solver.getSolution(), i + 1);
//     }
    
//     std::cout << "Done" << std::endl;
//     return 0;
// }

// for NASA supersonic test case

// int main() {
//     Mesh mesh;
//     mesh.load("supersonic_wedge_flow");

//     Solver solver(mesh);
//     solver.setFlowConditions(2.5, 101325.0, 300.0); // Ma=2.5, P=101325 Pa, T=300 K
//     solver.initialize();

//     std::cout << "NASA supersonic wedge case with JFNK" << std::endl;

//     solver.solveImplicit(1e-2, 150); 
    
//     saveResults(mesh, solver.getSolution(), 99999);

//     std::cout << "Done" << std::endl;
//     return 0;
// }


// for NASA laminar flat plate


// int main() {
//     Mesh mesh;
//     // Load your flat plate unstructured mesh
//     mesh.load("new_laminar_flat_plate"); 

//     Solver solver(mesh);

//     solver.setNumericalMethod(FluxScheme::HLLC, SpatialScheme::SECOND_ORDER_CENTRAL);
//     // Convert Flow conditions to SI if needed: 
//     // Mach 0.1, P_static = 41368.5 Pa, T_static = 388.89 K
//     solver.setFlowConditions(0.1, 41368.5, 388.89);
    
//     // Add Viscosity initialization for flat plate (Re=200k, L=0.3048m)
//     solver.setKinematicViscosity(200000.0, 0.3048);
    
//     solver.initialize();

//     std::cout << "Starting simulation JFNK..." << std::endl;
//     // You might need a smaller CFL or local timestep since viscous cells near the wall are tiny
//     solver.solveImplicit(1e-4, 50000); 
//     double L = 0.3048;
//     solver.extractBoundaryLayerProfile(L, "bl_profile.csv");
//     solver.extractSkinFriction("cf_profile.csv");
//     saveResults(mesh, solver.getSolution(), 9999);
//     return 0;
// }


// SWBLI mach 5 shock wave:
// int main() {
//     Mesh mesh;
//     mesh.load("swbli_8x"); 

//     Solver solver(mesh);
//     solver.setFreestreamStatic(5.0, 4006.88, 68.33);
   
//     solver.setKinematicViscosity(1.88e7, 0.5);
//     solver.initialize();
    
//     std::cout<< "First part of JFNK" << std::endl;
//     solver.setNumericalMethod(FluxScheme::HLLC, SpatialScheme::FIRST_ORDER);
//     solver.solveImplicit(1e3, 100000);

//     std::cout << "Starting NASA Mach 5 SWBLI simulation..." << std::endl;
//     solver.setNumericalMethod(FluxScheme::HLLC, SpatialScheme::SECOND_ORDER_LIMITED);
//     solver.solveImplicit(1e-4, 200000);
//     solver.extractSkinFriction("swbli_cf_profile.csv");
//     saveResults(mesh, solver.getSolution(), 9999);
//     return 0;   
// }


// pm

int main(){
    Mesh mesh;
    mesh.load("pm_1x");

    Solver solver(mesh);
    solver.setNumericalMethod(FluxScheme::HLLC, SpatialScheme::SECOND_ORDER_LIMITED);
    solver.initialize();

    std::cout<< "Start implicit" << std::endl;

    solver.solveImplicit(1e-6, 500000);
    std::cout << "Done implicit" << std::endl;
    saveResults(mesh, solver.getSolution(), 9999);
    solver.extractWallPressure("pm_centerline.csv");
    std::cout<< "done" << std::endl;
    return 0;
}