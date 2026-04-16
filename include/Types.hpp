#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>
#include <cmath>

// geometry and mesh structures

// Node
struct Node {
    double x, y;
};

// Face
// connect relationship
struct Face {
    int node1, node2; // two nodes that form the face
    int owner;        // Owner Cell
    int neighbor;     // Neighbor Cell, -1 for physical boundary
    
    double nx, ny;    // unit normal vector (pointing owner -> neighbor)
    double area;      // face length 
    double x_mid, y_mid; // face center coordinates (for higher-order reconstruction)
    
    int bcTag;        // boundary tag (corresponding to Gmsh Physical ID)
                      // 0=internal face, 1=Inlet, 2=Outlet, 3=Wall
};

// Cell
struct Cell {
    int id;                 // Cell global index
    std::vector<int> nodes; // list of node indices that form the cell
    double volume;          // cell area
    double x_c, y_c;        // centroid coordinates
};


enum class FluxScheme{
    LLF,
    ROE,
    HLLC,
    AUSM
};

enum class SpatialScheme{
    FIRST_ORDER,
    SECOND_ORDER_CENTRAL,
    SECOND_ORDER_LIMITED
};

// 2. Physical State Structures


// U = [rho, rhou, rhov, rhoE]
struct State {
    double rho;   // density
    double rhou;  // x-momentum
    double rhov;  // y-momentum
    double rhoE;  // total energy

    // add variables for SST

    double rhok;
    double rhow;

   
    State(double r=0, double ru=0, double rv=0, double re=0, double rk=0, double rw=0) 
        : rho(r), rhou(ru), rhov(rv), rhoE(re), rhok(rk), rhow(rw) {}


    State operator+(const State& other) const {
        return {rho + other.rho, rhou + other.rhou, rhov + other.rhov, rhoE + other.rhoE, rhok+ other.rhok, rhow + other.rhow};
    }
    
    State operator-(const State& other) const {
        return {rho - other.rho, rhou - other.rhou, rhov - other.rhov, rhoE - other.rhoE, rhok - other.rhok, rhow - other.rhow};
    }

    State operator*(double scalar) const {
        return {rho * scalar, rhou * scalar, rhov * scalar, rhoE * scalar, rhok * scalar, rhow * scalar};
    }
};

#endif // TYPES_HPP