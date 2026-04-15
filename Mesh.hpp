#ifndef MESH_HPP
#define MESH_HPP

#include <string>
#include <vector>
#include <map>
#include "Types.hpp" // 

class Mesh {
public:
    
    std::vector<Node> nodes; // node list
    std::vector<Cell> cells; // cell list
    std::vector<Face> faces; // face list (this is the object the solver loops over)

    // --- interface functions ---
    
    // constructor
    Mesh() = default;

    // main loading function: read file -> build topology -> compute geometry
    void load(const std::string& filename);

    // get mesh statistics
    int getNumCells() const { return cells.size(); }
    int getNumFaces() const { return faces.size(); }
    int getNumNodes() const { return nodes.size(); }

private:
   

    // 1. Read Gmsh .msh file 
    void readGmsh(const std::string& filename);

    // 2. Build face connectivity
    // Convert cell list to face list and find owner/neighbor relationships
    void buildConnectivity();

    // 3. Compute geometric properties (area, normal vector, centroid)
    void computeGeometry();

    // Temporary storage for boundary markers, used to tag Faces in buildConnectivity
    // Key: pair<min_node, max_node>, Value: Physical Tag
    std::map<std::pair<int,int>, int> boundaryEdges;
};

#endif // MESH_HPP