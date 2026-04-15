#ifndef FLUX_HPP
#define FLUX_HPP

#include <cmath>
#include <algorithm>
#include "Types.hpp" //

// constant
const double GAMMA = 1.4;

// p
inline double getPressure(const State& U) {
    double V2 = (U.rhou * U.rhou + U.rhov * U.rhov) / (U.rho * U.rho);
    return (GAMMA - 1.0) * (U.rhoE - 0.5 * U.rho * V2);
}

// soundspeed
inline double getSoundSpeed(double p, double rho) {
    return std::sqrt(GAMMA * p / rho);
}

// normal vector: F_n(U)
inline State computeNormalFlux(const State& U, double nx, double ny) {
    double p = getPressure(U);
    double u = U.rhou / U.rho;
    double v = U.rhov / U.rho;
    double un = u * nx + v * ny; // normal velocity

    return State(
        U.rho * un,
        U.rhou * un + p * nx,
        U.rhov * un + p * ny,
        (U.rhoE + p) * un
    );
}
inline State computeFluxHLLC(const State& UL, const State& UR, double nx, double ny) {
    double pL = getPressure(UL), pR = getPressure(UR);
    double uL = UL.rhou / UL.rho, uR = UR.rhou / UR.rho;
    double vL = UL.rhov / UL.rho, vR = UR.rhov / UR.rho;
    double unL = uL * nx + vL * ny, unR = uR * nx + vR * ny; 
    double cL = std::sqrt(GAMMA * pL / UL.rho), cR = std::sqrt(GAMMA * pR / UR.rho);

    double hL = (UL.rhoE + pL) / UL.rho, hR = (UR.rhoE + pR) / UR.rho;

    // Roe average
    double sq_rhoL = std::sqrt(UL.rho), sq_rhoR = std::sqrt(UR.rho);
    double u_t = (sq_rhoL * uL + sq_rhoR * uR) / (sq_rhoL + sq_rhoR);
    double v_t = (sq_rhoL * vL + sq_rhoR * vR) / (sq_rhoL + sq_rhoR);
    double h_t = (sq_rhoL * hL + sq_rhoR * hR) / (sq_rhoL + sq_rhoR);
    double un_t = u_t * nx + v_t * ny;
    double c_t = std::sqrt((GAMMA - 1.0) * (h_t - 0.5 * (u_t * u_t + v_t * v_t)));

    double SL = std::min(unL - cL, un_t - c_t);
    double SR = std::max(unR + cR, un_t + c_t);
    double S_star = (pR - pL + UL.rho * unL * (SL - unL) - UR.rho * unR * (SR - unR)) / 
                    (UL.rho * (SL - unL) - UR.rho * (SR - unR));

    if (SL >= 0.0) return computeNormalFlux(UL, nx, ny);
    if (SR <= 0.0) return computeNormalFlux(UR, nx, ny);

    State F_star;
    if (S_star >= 0.0) { 
        double coef = UL.rho * (SL - unL) / (SL - S_star);
        State U_star(coef, coef * (uL + (S_star - unL) * nx), coef * (vL + (S_star - unL) * ny),
                     coef * (UL.rhoE / UL.rho + (S_star - unL) * (S_star + pL / (UL.rho * (SL - unL)))));
        F_star = computeNormalFlux(UL, nx, ny) + (U_star - UL) * SL;
    } else { 
        double coef = UR.rho * (SR - unR) / (SR - S_star);
        State U_star(coef, coef * (uR + (S_star - unR) * nx), coef * (vR + (S_star - unR) * ny),
                     coef * (UR.rhoE / UR.rho + (S_star - unR) * (S_star + pR / (UR.rho * (SR - unR)))));
        F_star = computeNormalFlux(UR, nx, ny) + (U_star - UR) * SR;
    }
    return F_star;
}


// LLF numerical flux function

// Inputs: left state UL, right state UR, face normal vector (nx, ny)
// Output: flux vector
inline State computeFluxLLF(const State& UL, const State& UR, double nx, double ny) {
    // 1. Compute left and right physical fluxes
    State FL = computeNormalFlux(UL, nx, ny);
    State FR = computeNormalFlux(UR, nx, ny);

    // 2. Compute wave speeds (eigenvalues)
    double pL = getPressure(UL);
    double uL = UL.rhou / UL.rho; 
    double vL = UL.rhov / UL.rho;
    double unL = uL * nx + vL * ny;
    double cL = getSoundSpeed(pL, UL.rho);

    double pR = getPressure(UR);
    double uR = UR.rhou / UR.rho; 
    double vR = UR.rhov / UR.rho;
    double unR = uR * nx + vR * ny;
    double cR = getSoundSpeed(pR, UR.rho);

    // alpha = max(|un| + c)
    double alphaL = std::abs(unL) + cL;
    double alphaR = std::abs(unR) + cR;
    double alpha = std::max(alphaL, alphaR);

    // 3. LLF 
    // F = 0.5*(FL + FR) - 0.5*alpha*(UR - UL)
    State Dissipation = (UR - UL) * (0.5 * alpha);
    State AverageFlux = (FL + FR) * 0.5;

    return AverageFlux - Dissipation;
}
// Roe numerical flux function (with Harten's Entropy Fix)
inline State computeFluxROE(const State& UL, const State& UR, double nx, double ny) {
    
    double pL = getPressure(UL);
    double uL = UL.rhou / UL.rho; 
    double vL = UL.rhov / UL.rho;
    double hL = (UL.rhoE + pL) / UL.rho; 

    double pR = getPressure(UR);
    double uR = UR.rhou / UR.rho; 
    double vR = UR.rhov / UR.rho;
    double hR = (UR.rhoE + pR) / UR.rho;

 
    double sqrt_rhoL = std::sqrt(UL.rho);
    double sqrt_rhoR = std::sqrt(UR.rho);
    double inv_sum = 1.0 / (sqrt_rhoL + sqrt_rhoR);

    double u_tilde = (sqrt_rhoL * uL + sqrt_rhoR * uR) * inv_sum;
    double v_tilde = (sqrt_rhoL * vL + sqrt_rhoR * vR) * inv_sum;
    double h_tilde = (sqrt_rhoL * hL + sqrt_rhoR * hR) * inv_sum;

    double V2_tilde = u_tilde * u_tilde + v_tilde * v_tilde;
    double c2_tilde = (GAMMA - 1.0) * (h_tilde - 0.5 * V2_tilde);
    double c_tilde = std::sqrt(std::max(c2_tilde, 1e-8));

  
    double un_tilde = u_tilde * nx + v_tilde * ny;
    double ut_tilde = -u_tilde * ny + v_tilde * nx;


    double lambda[4];
    lambda[0] = un_tilde - c_tilde;
    lambda[1] = un_tilde;
    lambda[2] = un_tilde;
    lambda[3] = un_tilde + c_tilde;


    double epsilon = 0.1 * c_tilde; 
    for (int i = 0; i < 4; ++i) {
        if (std::abs(lambda[i]) < epsilon) {
            lambda[i] = 0.5 * (lambda[i] * lambda[i] / epsilon + epsilon);
        }
        lambda[i] = std::abs(lambda[i]); 
    }


    double delta_rho = UR.rho - UL.rho;
    double delta_u   = uR - uL;
    double delta_v   = vR - vL;
    double delta_p   = pR - pL;

    double delta_un  = delta_u * nx + delta_v * ny;
    double delta_ut  = -delta_u * ny + delta_v * nx;


    double rho_tilde = sqrt_rhoL * sqrt_rhoR;
    double alpha[4];
    alpha[0] = 0.5 / c2_tilde * (delta_p - rho_tilde * c_tilde * delta_un);
    alpha[1] = delta_rho - delta_p / c2_tilde;
    alpha[2] = rho_tilde * delta_ut;
    alpha[3] = 0.5 / c2_tilde * (delta_p + rho_tilde * c_tilde * delta_un);


    // dW = sum( |lambda_k| * alpha_k * K_k )
    double dW1 = lambda[0] * alpha[0];
    double dW2 = lambda[1] * alpha[1];
    double dW3 = lambda[2] * alpha[2];
    double dW4 = lambda[3] * alpha[3];

    State Dissipation;

    Dissipation.rho = dW1 * 1.0 + dW2 * 1.0 + dW3 * 0.0 + dW4 * 1.0;
    

    Dissipation.rhou = dW1 * (u_tilde - c_tilde * nx) +
                       dW2 * u_tilde +
                       dW3 * (-ny) +
                       dW4 * (u_tilde + c_tilde * nx);


    Dissipation.rhov = dW1 * (v_tilde - c_tilde * ny) +
                       dW2 * v_tilde +
                       dW3 * nx +
                       dW4 * (v_tilde + c_tilde * ny);


    Dissipation.rhoE = dW1 * (h_tilde - un_tilde * c_tilde) +
                       dW2 * (0.5 * V2_tilde) +
                       dW3 * ut_tilde +
                       dW4 * (h_tilde + un_tilde * c_tilde);

  
    double sgn_un = (un_tilde > 0.0) ? 1.0 : -1.0;
    Dissipation.rhok = std::abs(un_tilde) * (UR.rhok / UR.rho - UL.rhok / UL.rho) * rho_tilde;
    Dissipation.rhow = std::abs(un_tilde) * (UR.rhow / UR.rho - UL.rhow / UL.rho) * rho_tilde;


    State FL = computeNormalFlux(UL, nx, ny);
    State FR = computeNormalFlux(UR, nx, ny);

    return (FL + FR) * 0.5 - (Dissipation * 0.5);
}

// Different numerical flux functions
inline State computeInviscidFlux(const State& UL, const State& UR, double nx, double ny, FluxScheme scheme){
    switch (scheme){
        case FluxScheme::LLF:
            return computeFluxLLF(UL, UR, nx, ny);
        case FluxScheme::ROE:
            return computeFluxROE(UL, UR, nx, ny);
        case FluxScheme::HLLC:
            return computeFluxHLLC(UL, UR, nx, ny);
        default:
            return computeFluxLLF(UL, UR, nx, ny);

    }
}

#endif // FLUX_HPP