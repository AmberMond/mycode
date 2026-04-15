#include "Mesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

// Main loading function

void Mesh::load(const std::string& filename) {
    std::cout << "Loading mesh: " << filename << "..." << std::endl;
    
    // 1. read data
    readGmsh(filename);
    
    // 2. build connectivity
    buildConnectivity();
    
    // 3. compute geometric properties (area, normal)
    computeGeometry();
    
    std::cout << "Mesh initialization complete:" << std::endl;
    std::cout << "  - Nodes: " << nodes.size() - 1 << std::endl; // Subtract placeholder node 0
    std::cout << "  - Cells: " << cells.size() << std::endl;
    std::cout << "  - Faces: " << faces.size() << std::endl;
}

// 1. Read Gmsh file

void Mesh::readGmsh(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    int nNodes = 0, nElems = 0;

    // --- Read Nodes ---
    while (std::getline(file, line)) {
        if (line == "$Nodes") {
            file >> nNodes;
            nodes.resize(nNodes + 1); // Gmsh indices start at 1
            for (int i = 1; i <= nNodes; ++i) {
                int id; double z;
                file >> id >> nodes[i].x >> nodes[i].y >> z;
            }
            break;
        }
    }

    // --- Read Elements ---
    while (std::getline(file, line)) {
        if (line == "$Elements") {
            file >> nElems;
            for (int i = 0; i < nElems; ++i) {
                int id, type, nTags, tag, physTag = 0;
                file >> id >> type >> nTags;
                
                // Read Tags
                for (int j = 0; j < nTags; ++j) {
                    file >> tag;
                    if (j == 0) physTag = tag;
                }

                // Type 1: 2-node Line (boundary edge)
                if (type == 1) {
                    int n1, n2; file >> n1 >> n2;
                    // Store boundary tag for use in buildConnectivity
                    // Key is sorted node pair to ensure (1,2) and (2,1) are the same
                    boundaryEdges[{std::min(n1,n2), std::max(n1,n2)}] = physTag;
                }
                // Type 2: 3-node Triangle (internal cell)
                else if (type == 2) {
                    int n1, n2, n3; file >> n1 >> n2 >> n3;
                    Cell c;
                    c.id = cells.size(); // 0-based index
                    c.nodes = {n1, n2, n3};
                    cells.push_back(c);
                }
                // Type 3: 4-node Quad (if exists)
                else if (type == 3) {
                    int n1, n2, n3, n4; file >> n1 >> n2 >> n3 >> n4;
                    Cell c;
                    c.id = cells.size();
                    c.nodes = {n1, n2, n3, n4};
                    cells.push_back(c);
                }
                
                // After reading a line, consume the remaining characters
                std::getline(file, line); 
            }
            break;
        }
    }
}


// 2. Build Connectivity

void Mesh::buildConnectivity() {
    // Temporary map: Key = {min(n1,n2), max(n1,n2)}, Value = Face Index
    // Used to identify shared edges
    std::map<std::pair<int,int>, int> edgeToFaceIndex;

    for (int i = 0; i < cells.size(); ++i) {
        const auto& cell = cells[i];
        int n = cell.nodes.size(); 

        // Iterate over each edge of the cell
        for (int j = 0; j < n; ++j) {
            int n1 = cell.nodes[j];
            int n2 = cell.nodes[(j + 1) % n]; // Wrap around to the first node
            
            std::pair<int, int> edgeKey = {std::min(n1, n2), std::max(n1, n2)};

            if (edgeToFaceIndex.count(edgeKey)) {
                // Case A: Edge already exists in the map
                // This means it's the second time encountering this edge -> it's an internal face
                int faceIdx = edgeToFaceIndex[edgeKey];
                faces[faceIdx].neighbor = i; // Record right-side cell (Neighbor)
            } else {
                // Case B: Edge appears for the first time
                // Create a new face, current cell is the left-side cell (Owner)
                Face f;
                f.node1 = edgeKey.first;
                f.node2 = edgeKey.second;
                f.owner = i;
                f.neighbor = -1; // Temporarily mark as boundary, if not visited later, it's a true boundary
                f.bcTag = 0;     // Default to internal

                // Check if marked as physical boundary in readGmsh
                if (boundaryEdges.count(edgeKey)) {
                    f.bcTag = boundaryEdges[edgeKey];
                }

                faces.push_back(f);
                edgeToFaceIndex[edgeKey] = faces.size() - 1; // Record index
            }
        }
    }
}


// 3. Compute Geometry

void Mesh::computeGeometry() {
    // A. Compute cell geometry (area and centroid)
    for (auto& c : cells) {
        // Simple handling for triangles (if quadrilaterals exist, need to split or use shoelace formula)
        if (c.nodes.size() == 3) {
            double x1 = nodes[c.nodes[0]].x, y1 = nodes[c.nodes[0]].y;
            double x2 = nodes[c.nodes[1]].x, y2 = nodes[c.nodes[1]].y;
            double x3 = nodes[c.nodes[2]].x, y3 = nodes[c.nodes[2]].y;

            // Cross product's half = area
            c.volume = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
            c.x_c = (x1 + x2 + x3) / 3.0;
            c.y_c = (y1 + y2 + y3) / 3.0;
        }
        else if (c.nodes.size() == 4) {
             // Simple quadrilateral centroid approximation
             double sumX=0, sumY=0;
             for(int nid : c.nodes) { sumX += nodes[nid].x; sumY += nodes[nid].y; }
             c.x_c = sumX/4.0; c.y_c = sumY/4.0;
             // Shoelace formula for area
             double area = 0.0;
             for(int i=0; i<4; ++i) {
                 int n1 = c.nodes[i];
                 int n2 = c.nodes[(i+1)%4];
                 area += (nodes[n1].x * nodes[n2].y - nodes[n2].x * nodes[n1].y);
             }
             c.volume = 0.5 * std::abs(area);
        }
    }

    // B. Compute face geometry (normal vector)
    for (auto& f : faces) {
        double dx = nodes[f.node2].x - nodes[f.node1].x;
        double dy = nodes[f.node2].y - nodes[f.node1].y;
        
        f.area = std::sqrt(dx*dx + dy*dy);
        f.x_mid = 0.5 * (nodes[f.node1].x + nodes[f.node2].x);
        f.y_mid = 0.5 * (nodes[f.node1].y + nodes[f.node2].y);

        // Initial normal vector (rotate 90 degrees: -dy, dx)
        f.nx = dy / f.area;
        f.ny = -dx / f.area;

        // Critical check: Ensure normal vector points from Owner to Neighbor (or outside)
        // Vector V = Face_mid - Cell_Centroid_Owner
        double vec_x = f.x_mid - cells[f.owner].x_c;
        double vec_y = f.y_mid - cells[f.owner].y_c;

        // If dot product < 0, normal vector points back into Owner, needs flipping
        if (vec_x * f.nx + vec_y * f.ny < 0) {
            f.nx = -f.nx;
            f.ny = -f.ny;
        }
    }
}