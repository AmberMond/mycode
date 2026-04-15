#ifndef SOLVER_HPP
#define SOLVER_HPP





#include "Mesh.hpp"
#include "Types.hpp"
#include "Flux.hpp"


struct CellGradient{
    double drhodx = 0.0, drhody = 0.0;
    double dudx = 0.0, dudy = 0.0;
    double dvdx = 0.0, dvdy = 0.0;
    double dpdx = 0.0, dpdy = 0.0;
    double dTdx = 0.0, dTdy = 0.0;

    // turbulence variable gradients

    double dkdx = 0.0, dkdy = 0.0;
    double dwdx = 0.0, dwdy = 0.0;
};

#include <vector>

class Solver {
public:
    FluxScheme iFluxScheme = FluxScheme::HLLC;
    SpatialScheme iSpatialScheme = SpatialScheme::SECOND_ORDER_LIMITED;

    void setNumericalMethod(FluxScheme flux, SpatialScheme spatial) {
        iFluxScheme = flux;
        iSpatialScheme = spatial;
    }
    double M_inf = 0.0;
    double physical_time = 0.0;
    double nu = 0.01;

    Solver(const Mesh& mesh);
    void setFreestreamStatic(double M_inf_in, double P_static, double T_static);
    void solveUnsteady(double dt_phy, int maxTimeSteps);
    void extractBoundaryLayerProfile(double x_target, const std::string& filename);
    void extractSkinFriction(const std::string& filename);
    // set flow conditions: Ma, P_ref, T_ref
    void setFlowConditions(double M_inf, double P_ref, double T_ref);

    // uniform flow
    void initialize();

    // residual
    void computeResiduals(const std::vector<State>& U_in, std::vector<State>& R_out, bool update_limiters = true);

    // explicit time stepping
    void timeStepExplicit(double dt);
    void extractWallPressure(const std::string& filename);
    const std::vector<State>& getSolution() const { return U_current; }

    // Implicit solver JFNK
    // Runs 'maxNewtonsteps' to solve steady state
    void solveImplicit(double tolerance = 1e-6, int maxNewtonstep = 50);

    void calculateEntropyError();
    void calculateDrag();
    void initializeTGV();

    void setKinematicViscosity(double Re, double L_ref);
    
private:
    const Mesh& mesh;
    std::vector<State> U_current; 
    std::vector<State> Residuals; 
    std::vector<double> dt_local; // local time step for each cell
    std::vector<State> U_old;
    std::vector<double> wallDistances;
    std::vector<double> limiters;
    
    std::vector<CellGradient> cellGradients;
    

    void computeLimiters(const std::vector<State>& U_in);
    void computeGradients(const std::vector<State>& U_in);
    void computeWallDistances();
    void computeLocalTimeStep(double CFL);
    double getTotalVolume() const;
    double norm_L2_VolumeWeighted(const std::vector<State>& R) const;
    double gamma = 1.4;
    double R_gas = 287.0;
    
    // BC
    double P_total_inlet; // inlet total pressure
    double T_total_inlet; // inlet total temprature
    double P_static_exit; // outlet static pressure

    // set boundary ghost state
    State getBoundaryGhostState(const State& U_inner, const Face& face);
    
    // like set others function in 2D explicit bump
    State primitiveToConservative(double rho, double u, double v, double p, double k = 0.0, double omega = 0.0);

    void computeDiagonalPreconditioner(std::vector<double>& D);
    // [New] Linear solver GMRES
    // Solves J * delta_U = -R
    // Returns the number of iterations used
    int solveGMRES(const std::vector<State>& U_base, 
                    const std::vector<State>& R_spatial_base,
                   const std::vector<State>& RHS, 
                   std::vector<State>& delta_U);

    // [New] JFNK Matrix-Vector Product
    // Computes J * v approx. (R(U + eps*v) - R(U)) / eps
    void applyJacobian(const std::vector<State>& U_base, 
                       const std::vector<State>& R_base, // pass R(U_base) to save cost
                       const std::vector<State>& v, 
                       std::vector<State>& Jv);

    // [New] Linear Algebra Helpers (Inline implementations or move to LinAlg.hpp)
    double norm_L2(const std::vector<State>& v) const;
    double dot_product(const std::vector<State>& a, const std::vector<State>& b) const;
};

#endif