#include "Solver.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>

Solver::Solver(const Mesh& m) : mesh(m) {
    U_current.resize(mesh.cells.size());
    Residuals.resize(mesh.cells.size());
}
//test
void Solver::setFlowConditions(double M_inf_in, double P_ref, double T_ref) {
    // calculate total pressure and temperature at inlet
    // T0 = T * (1 + (gamma-1)/2 * M^2)
    // P0 = P * (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1))
    
    M_inf = M_inf_in;
    double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
    double factor_P = std::pow(factor_T, gamma / (gamma - 1.0));

    // inlet 
    T_total_inlet = T_ref * factor_T;
    P_total_inlet = P_ref * factor_P;

    // outlet
    P_static_exit = P_ref;

    std::cout << "--- Flow Conditions Set ---" << std::endl;
    std::cout << "Mach: " << M_inf << std::endl;
    std::cout << "Inlet P0: " << P_total_inlet << ", T0: " << T_total_inlet << std::endl;
    std::cout << "Outlet P_static: " << P_static_exit << std::endl;
}

void Solver::setKinematicViscosity(double Re, double L_ref) {
    double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
    double T_static = T_total_inlet / factor_T;
    double c = std::sqrt(gamma * R_gas * T_static);
    double u_inf = M_inf * c;
    
    // nu = U * L / Re
    nu = (u_inf * L_ref) / Re;
    std::cout << "Kinematic Viscosity (nu) set to: " << nu << std::endl;
}

// for 2d bump
// void Solver::initialize() {
//     // inititalize : uniform flow
    
//     double p = P_static_exit;
//     // T = T0 / (1 + (g-1)/2 * M^2) 
//     double factor_T = std::pow(P_total_inlet / p, (gamma - 1.0) / gamma); // Isentropic relation
//     double T = T_total_inlet / factor_T;
    
//     double rho = p / (R_gas * T);
//     double c = std::sqrt(gamma * R_gas * T);
    
//     // u^2 = M^2 * c^2 
//     double V = std::sqrt(2.0 * (gamma * R_gas / (gamma - 1.0)) * (T_total_inlet - T));
    
    
//     State init = primitiveToConservative(rho, V, 0.0, p);
//     std::fill(U_current.begin(), U_current.end(), init);
    
//     std::cout << "Field initialized with M=" << V/c << ", P=" << p << std::endl;
// }


// for wind tunnel
// void Solver::initialize(){
//     double rho = 1.4;
//     double u = 3.0;
//     double v = 0.0;
//     double p = 1.0;

//     State init = primitiveToConservative(rho, u, v, p);
//     std::fill(U_current.begin(), U_current.end(), init);

//     std::cout << "Field initialized with Mach 3 forward step condition." << std::endl;
// }

void Solver::computeGradients(const std::vector<State>& U_in) {
    if (cellGradients.size() != mesh.cells.size()) {
        cellGradients.resize(mesh.cells.size());
    }
   
    std::fill(cellGradients.begin(), cellGradients.end(), CellGradient());

    for (const auto& f : mesh.faces) {
        State UL = U_in[f.owner];
        State UR = (f.neighbor >= 0) ? U_in[f.neighbor] : getBoundaryGhostState(UL, f);

        // left cell primitive variables
        double uL = UL.rhou / UL.rho;
        double vL = UL.rhov / UL.rho;
        double pL = getPressure(UL);
        double TL = pL / (UL.rho * R_gas);

        // right cell primitive variables
        double uR = UR.rhou / UR.rho;
        double vR = UR.rhov / UR.rho;
        double pR = getPressure(UR);
        double TR = pR / (UR.rho * R_gas);

        // face-centered values (simple average)
        double rho_f = 0.5 * (UL.rho + UR.rho);
        double p_f = 0.5 * (pL + pR);
        double u_f = 0.5 * (uL + uR);
        double v_f = 0.5 * (vL + vR);
        double T_f = 0.5 * (TL + TR);

        
        // owner
        cellGradients[f.owner].drhodx += rho_f * f.nx * f.area;
        cellGradients[f.owner].drhody += rho_f * f.ny * f.area;
        cellGradients[f.owner].dpdx += p_f * f.nx * f.area;
        cellGradients[f.owner].dpdy += p_f * f.ny * f.area;
        cellGradients[f.owner].dudx += u_f * f.nx * f.area;
        cellGradients[f.owner].dudy += u_f * f.ny * f.area;
        cellGradients[f.owner].dvdx += v_f * f.nx * f.area;
        cellGradients[f.owner].dvdy += v_f * f.ny * f.area;
        cellGradients[f.owner].dTdx += T_f * f.nx * f.area;
        cellGradients[f.owner].dTdy += T_f * f.ny * f.area;

        // neighbor
        if (f.neighbor >= 0) {
            cellGradients[f.neighbor].drhodx -= rho_f * f.nx * f.area;
            cellGradients[f.neighbor].drhody -= rho_f * f.ny * f.area;
            cellGradients[f.neighbor].dpdx -= p_f * f.nx * f.area;
            cellGradients[f.neighbor].dpdy -= p_f * f.ny * f.area;
            cellGradients[f.neighbor].dudx -= u_f * f.nx * f.area;
            cellGradients[f.neighbor].dudy -= u_f * f.ny * f.area;
            cellGradients[f.neighbor].dvdx -= v_f * f.nx * f.area;
            cellGradients[f.neighbor].dvdy -= v_f * f.ny * f.area;
            cellGradients[f.neighbor].dTdx -= T_f * f.nx * f.area;
            cellGradients[f.neighbor].dTdy -= T_f * f.ny * f.area;
        }
    }

    // divide by cell volume to get final gradients
    for (int i = 0; i < mesh.cells.size(); ++i) {
        double vol = mesh.cells[i].volume;
        cellGradients[i].drhodx /= vol;  cellGradients[i].drhody /= vol;
        cellGradients[i].dpdx /= vol;  cellGradients[i].dpdy /= vol;

        cellGradients[i].dudx /= vol;  cellGradients[i].dudy /= vol;
        cellGradients[i].dvdx /= vol;  cellGradients[i].dvdy /= vol;
        cellGradients[i].dTdx /= vol;  cellGradients[i].dTdy /= vol;
    }
}

// B-J limiter
void Solver::computeLimiters(const std::vector<State>& U_in) {
    int N = mesh.cells.size();
    if (limiters.size() != N) limiters.resize(N);
    
    // 初始化所有单元限值为 1.0 (无限制)
    std::fill(limiters.begin(), limiters.end(), 1.0);

    // 对于结构松散的网格，我们需要用到节点的邻接信息
    // 这里我们假设你网格中的面只提供了 f.owner 和 f.neighbor
    
    // 采用局部标量：例如使用压力 (Pressure) 作为限制的对象
    // 如果你愿意，也可以对 u, v, T 分别计算 phi 然后取最小值
    std::vector<double> P_cents(N);
    for (int i = 0; i < N; ++i) {
        P_cents[i] = getPressure(U_in[i]);
    }

    std::vector<double> P_max(N), P_min(N);
    for (int i = 0; i < N; ++i) {
        P_max[i] = P_cents[i];
        P_min[i] = P_cents[i];
    }
    
    // 1. 扫描所有面，找出每个单元的局部最大压力和最小压力
    for (const auto& f : mesh.faces) {
        int o = f.owner;
        int n = f.neighbor;
        
        double p_ghost = 0.0;
        if (n >= 0) {
            p_ghost = P_cents[n];
        } else {
            // 如果是边界，获取边界幽灵单元压力
            State U_ghost = getBoundaryGhostState(U_in[o], f);
            p_ghost = getPressure(U_ghost);
        }
        
        P_max[o] = std::max(P_max[o], p_ghost);
        P_min[o] = std::min(P_min[o], p_ghost);
        
        if (n >= 0) {
            P_max[n] = std::max(P_max[n], P_cents[o]);
            P_min[n] = std::min(P_min[n], P_cents[o]);
        }
    }
    
    // 2. 再次扫描所有面，计算每个面上的无限制重构压力，推导限制系数
    // 我们需要能够累积向各个单元施加最严格（最小）的限制系数
    const double eps = 1e-12; // 防止除以零的小量
    
    for (const auto& f : mesh.faces) {
        int o = f.owner;
        int n = f.neighbor;
        
        // 计算面向所有者单元中心的向量
        double dx_o = f.x_mid - mesh.cells[o].x_c;
        double dy_o = f.y_mid - mesh.cells[o].y_c;
        
        // 计算左侧无限制重构压力 P_L_f
        // 假设你要通过当前代码结构里的梯度来算。对于压力，你需要提取它的梯度。
        // 这里提供一个通过状态量重构出来的近似压力增量的方法：
        // deltaP \approx P_cent * (dT/dx * dx + ...) ... 这种依赖已有的 u,v,T 梯度
        // 或者比较稳妥的是在 computeGradients 里多存一个 dPdx，dPdy。
        
        // 简便法：利用我们已有的 T 梯度和已有的 u, v, rho 假设它反映了类似压力梯度形变
        CellGradient gradO = cellGradients[o];
        double u_o = U_in[o].rhou / U_in[o].rho;
        double v_o = U_in[o].rhov / U_in[o].rho;
        double T_o = P_cents[o] / (U_in[o].rho * R_gas);
        
        double T_o_f = T_o + 1.0 * (gradO.dTdx * dx_o + gradO.dTdy * dy_o); // 假设只限制温度主导的参数
        double P_o_f = U_in[o].rho * R_gas * T_o_f; 
        
        // 计算 alpha
        double delta_2 = P_o_f - P_cents[o];
        double alpha_o = 1.0;
        if (delta_2 > 0.0) {
            alpha_o = (P_max[o] - P_cents[o]) / (delta_2 + eps);
            alpha_o = std::min(1.0, alpha_o);
        } else if (delta_2 < 0.0) {
            alpha_o = (P_min[o] - P_cents[o]) / (delta_2 - eps);
            alpha_o = std::min(1.0, alpha_o);
        }
        
        limiters[o] = std::min(limiters[o], alpha_o);
        
        // --- 对于右侧单元(neighbor) ---
        if (n >= 0) {
            double dx_n = f.x_mid - mesh.cells[n].x_c;
            double dy_n = f.y_mid - mesh.cells[n].y_c;
            
            CellGradient gradN = cellGradients[n];
            double T_n = P_cents[n] / (U_in[n].rho * R_gas);
            double T_n_f = T_n + 1.0 * (gradN.dTdx * dx_n + gradN.dTdy * dy_n);
            double P_n_f = U_in[n].rho * R_gas * T_n_f;
            
            double delta_n = P_n_f - P_cents[n];
            double alpha_n = 1.0;
            if (delta_n > 0.0) {
                alpha_n = (P_max[n] - P_cents[n]) / (delta_n + eps);
                alpha_n = std::min(1.0, alpha_n);
            } else if (delta_n < 0.0) {
                alpha_n = (P_min[n] - P_cents[n]) / (delta_n - eps);
                alpha_n = std::min(1.0, alpha_n);
            }
            
            limiters[n] = std::min(limiters[n], alpha_n);
        }
    }
}






// Venkat Limiter
// void Solver::computeLimiters(const std::vector<State>& U_in) {
//     int N = mesh.cells.size();
//     if (limiters.size() != N) limiters.resize(N);
    
//     // 初始化所有单元限值为 1.0 (无限制)
//     std::fill(limiters.begin(), limiters.end(), 1.0);

//     std::vector<double> P_cents(N);
//     for (int i = 0; i < N; ++i) {
//         P_cents[i] = getPressure(U_in[i]);
//     }

//     std::vector<double> P_max(N), P_min(N);
//     for (int i = 0; i < N; ++i) {
//         P_max[i] = P_cents[i];
//         P_min[i] = P_cents[i];
//     }
    
//     // 1. 扫描所有面，找出每个单元的局部最大压力和最小压力
//     for (const auto& f : mesh.faces) {
//         int o = f.owner;
//         int n = f.neighbor;
        
//         double p_ghost = 0.0;
//         if (n >= 0) {
//             p_ghost = P_cents[n];
//         } else {
//             State U_ghost = getBoundaryGhostState(U_in[o], f);
//             p_ghost = getPressure(U_ghost);
//         }
        
//         P_max[o] = std::max(P_max[o], p_ghost);
//         P_min[o] = std::min(P_min[o], p_ghost);
        
//         if (n >= 0) {
//             P_max[n] = std::max(P_max[n], P_cents[o]);
//             P_min[n] = std::min(P_min[n], P_cents[o]);
//         }
//     }
    
//     // 2. 再次扫描所有面，利用 Venkatakrishnan 光滑限制器计算 alpha
//     // Venkatakrishnan 光滑系数 K，通常取 0.5 到 5.0 之间。
//     // K 越大，允许的光滑波动越大，收敛性越好，但激波捕捉可能会变宽一点。默认先用 5.0。
//     double K_venkat = 5.0; 
    
//     for (const auto& f : mesh.faces) {
//         int o = f.owner;
//         int n = f.neighbor;
        
//         // ------------- Owner 侧 -------------
//         double dx_o = f.x_mid - mesh.cells[o].x_c;
//         double dy_o = f.y_mid - mesh.cells[o].y_c;
        
//         CellGradient gradO = cellGradients[o];
//         double T_o = P_cents[o] / (U_in[o].rho * R_gas);
//         double T_o_f = T_o + 1.0 * (gradO.dTdx * dx_o + gradO.dTdy * dy_o);
//         double P_o_f = U_in[o].rho * R_gas * T_o_f; 
        
//         double delta_2_o = P_o_f - P_cents[o];
        
//         // 计算 Venkat 保护微小量 eps^2 = (K * dx)^3 
//         double h_o = std::sqrt(mesh.cells[o].volume);
//         double P_ref_o = P_cents[o];
//         double eps_venkat_o = K_venkat * K_venkat * h_o * h_o * h_o * P_ref_o * P_ref_o; // 也可以不乘 P_ref，直接用物理量的立方
//         // double eps_venkat_o = K_venkat * K_venkat * h_o * h_o * h_o;

//         double alpha_o = 1.0;
//         if (delta_2_o > 1e-12) {
//             double D_max = P_max[o] - P_cents[o];
//             double num = D_max * D_max + eps_venkat_o + 2.0 * delta_2_o * D_max;
//             double den = D_max * D_max + 2.0 * delta_2_o * delta_2_o + D_max * delta_2_o + eps_venkat_o;
//             alpha_o = num / den;
//         } else if (delta_2_o < -1e-12) {
//             double D_min = P_min[o] - P_cents[o];
//             double num = D_min * D_min + eps_venkat_o + 2.0 * delta_2_o * D_min;
//             double den = D_min * D_min + 2.0 * delta_2_o * delta_2_o + D_min * delta_2_o + eps_venkat_o;
//             alpha_o = num / den;
//         }
//         limiters[o] = std::min(limiters[o], alpha_o);
        
//         // ------------- Neighbor 侧 -------------
//         if (n >= 0) {
//             double dx_n = f.x_mid - mesh.cells[n].x_c;
//             double dy_n = f.y_mid - mesh.cells[n].y_c;
            
//             CellGradient gradN = cellGradients[n];
//             double T_n = P_cents[n] / (U_in[n].rho * R_gas);
//             double T_n_f = T_n + 1.0 * (gradN.dTdx * dx_n + gradN.dTdy * dy_n);
//             double P_n_f = U_in[n].rho * R_gas * T_n_f;
            
//             double delta_2_n = P_n_f - P_cents[n];
            
//             double h_n = std::sqrt(mesh.cells[n].volume);
//             double P_ref_n = P_cents[n];
//             double eps_venkat_n = K_venkat * K_venkat * h_n * h_n * h_n * P_ref_n * P_ref_n;

//             double alpha_n = 1.0;
//             if (delta_2_n > 1e-12) {
//                 double D_max = P_max[n] - P_cents[n];
//                 double num = D_max * D_max + eps_venkat_n + 2.0 * delta_2_n * D_max;
//                 double den = D_max * D_max + 2.0 * delta_2_n * delta_2_n + D_max * delta_2_n + eps_venkat_n;
//                 alpha_n = num / den;
//             } else if (delta_2_n < -1e-12) {
//                 double D_min = P_min[n] - P_cents[n];
//                 double num = D_min * D_min + eps_venkat_n + 2.0 * delta_2_n * D_min;
//                 double den = D_min * D_min + 2.0 * delta_2_n * delta_2_n + D_min * delta_2_n + eps_venkat_n;
//                 alpha_n = num / den;
//             }
//             limiters[n] = std::min(limiters[n], alpha_n);
//         }
//     }
// }

// void Solver::computeLimiters(const std::vector<State>& U_in) {
//     int N = mesh.cells.size();
//     if (limiters.size() != N) limiters.resize(N);
    
//     std::fill(limiters.begin(), limiters.end(), 1.0);

//     // Primitive variable arrays
//     std::vector<double> P_cents(N), U_cents(N), V_cents(N), Rho_cents(N);
//     for (int i = 0; i < N; ++i) {
//         P_cents[i] = getPressure(U_in[i]);
//         Rho_cents[i] = U_in[i].rho;
//         U_cents[i] = U_in[i].rhou / U_in[i].rho;
//         V_cents[i] = U_in[i].rhov / U_in[i].rho;
//     }

//     std::vector<double> P_max(N), P_min(N);
//     std::vector<double> U_max(N), U_min(N);
//     std::vector<double> V_max(N), V_min(N);
//     std::vector<double> Rho_max(N), Rho_min(N);

//     for (int i = 0; i < N; ++i) {
//         P_max[i] = P_cents[i]; P_min[i] = P_cents[i];
//         U_max[i] = U_cents[i]; U_min[i] = U_cents[i];
//         V_max[i] = V_cents[i]; V_min[i] = V_cents[i];
//         Rho_max[i] = Rho_cents[i]; Rho_min[i] = Rho_cents[i];
//     }
    
//     // 1. Gather local max/min bounds for all primitive variables
//     for (const auto& f : mesh.faces) {
//         int o = f.owner;
//         int n = f.neighbor;
        
//         double p_ghost = 0.0, rho_ghost = 0.0, u_ghost = 0.0, v_ghost = 0.0;
//         if (n >= 0) {
//             p_ghost = P_cents[n]; rho_ghost = Rho_cents[n];
//             u_ghost = U_cents[n]; v_ghost = V_cents[n];
//         } else {
//             State U_ghost = getBoundaryGhostState(U_in[o], f);
//             p_ghost = getPressure(U_ghost); rho_ghost = U_ghost.rho;
//             u_ghost = U_ghost.rhou / U_ghost.rho; v_ghost = U_ghost.rhov / U_ghost.rho;
//         }
        
//         P_max[o] = std::max(P_max[o], p_ghost); P_min[o] = std::min(P_min[o], p_ghost);
//         U_max[o] = std::max(U_max[o], u_ghost); U_min[o] = std::min(U_min[o], u_ghost);
//         V_max[o] = std::max(V_max[o], v_ghost); V_min[o] = std::min(V_min[o], v_ghost);
//         Rho_max[o] = std::max(Rho_max[o], rho_ghost); Rho_min[o] = std::min(Rho_min[o], rho_ghost);
        
//         if (n >= 0) {
//             P_max[n] = std::max(P_max[n], P_cents[o]); P_min[n] = std::min(P_min[n], P_cents[o]);
//             U_max[n] = std::max(U_max[n], U_cents[o]); U_min[n] = std::min(U_min[n], U_cents[o]);
//             V_max[n] = std::max(V_max[n], V_cents[o]); V_min[n] = std::min(V_min[n], V_cents[o]);
//             Rho_max[n] = std::max(Rho_max[n], Rho_cents[o]); Rho_min[n] = std::min(Rho_min[n], Rho_cents[o]);
//         }
//     }
    
//     // 2. Loop over faces and calculate alphas for each variable
//     double K_venkat = 5.0; 
    
//     // Helper lambda to calculate alpha
//     auto calc_alpha = [K_venkat](double v_c, double v_max, double v_min, double delta, double h, double v_ref) -> double {
//         if (std::abs(delta) < 1e-12) return 1.0;
//         double eps = K_venkat * K_venkat * h * h * h * v_ref * v_ref;
//         if (v_ref == 0.0) eps = K_venkat * K_venkat * h * h * h; // fallback

//         if (delta > 0) {
//             double D_max = v_max - v_c;
//             double num = D_max * D_max + eps + 2.0 * delta * D_max;
//             double den = D_max * D_max + 2.0 * delta * delta + D_max * delta + eps;
//             return std::max(0.0, num/den);
//         } else {
//             double D_min = v_min - v_c;
//             double num = D_min * D_min + eps + 2.0 * delta * D_min;
//             double den = D_min * D_min + 2.0 * delta * delta + D_min * delta + eps;
//             return std::max(0.0, num/den);
//         }
//     };

//     for (const auto& f : mesh.faces) {
//         int o = f.owner;
//         int n = f.neighbor;
        
//         // ------------- Owner side limits -------------
//         double dx_o = f.x_mid - mesh.cells[o].x_c;
//         double dy_o = f.y_mid - mesh.cells[o].y_c;
//         CellGradient gradO = cellGradients[o];
//         double h_o = std::sqrt(mesh.cells[o].volume);

//         double delta_p_o   = gradO.dpdx * dx_o + gradO.dpdy * dy_o;
//         double delta_rho_o = gradO.drhodx * dx_o + gradO.drhody * dy_o;
//         double delta_u_o   = gradO.dudx * dx_o + gradO.dudy * dy_o;
//         double delta_v_o   = gradO.dvdx * dx_o + gradO.dvdy * dy_o;

//         double a_p = calc_alpha(P_cents[o], P_max[o], P_min[o], delta_p_o, h_o, P_cents[o]);
//         double a_rho = calc_alpha(Rho_cents[o], Rho_max[o], Rho_min[o], delta_rho_o, h_o, Rho_cents[o]);
//         double a_u = calc_alpha(U_cents[o], U_max[o], U_min[o], delta_u_o, h_o, std::max(1.0, std::abs(U_cents[o]))); // u ref is 1 to avoid zeroing
//         double a_v = calc_alpha(V_cents[o], V_max[o], V_min[o], delta_v_o, h_o, std::max(1.0, std::abs(V_cents[o]))); // v ref is 1 to avoid zeroing

//         limiters[o] = std::min({limiters[o], a_p, a_rho, a_u, a_v});
        
//         // ------------- Neighbor side limits -------------
//         if (n >= 0) {
//             double dx_n = f.x_mid - mesh.cells[n].x_c;
//             double dy_n = f.y_mid - mesh.cells[n].y_c;
//             CellGradient gradN = cellGradients[n];
//             double h_n = std::sqrt(mesh.cells[n].volume);

//             double delta_p_n   = gradN.dpdx * dx_n + gradN.dpdy * dy_n;
//             double delta_rho_n = gradN.drhodx * dx_n + gradN.drhody * dy_n;
//             double delta_u_n   = gradN.dudx * dx_n + gradN.dudy * dy_n;
//             double delta_v_n   = gradN.dvdx * dx_n + gradN.dvdy * dy_n;

//             double a_p_n = calc_alpha(P_cents[n], P_max[n], P_min[n], delta_p_n, h_n, P_cents[n]);
//             double a_rho_n = calc_alpha(Rho_cents[n], Rho_max[n], Rho_min[n], delta_rho_n, h_n, Rho_cents[n]);
//             double a_u_n = calc_alpha(U_cents[n], U_max[n], U_min[n], delta_u_n, h_n, std::max(1.0, std::abs(U_cents[n])));
//             double a_v_n = calc_alpha(V_cents[n], V_max[n], V_min[n], delta_v_n, h_n, std::max(1.0, std::abs(V_cents[n])));

//             limiters[n] = std::min({limiters[n], a_p_n, a_rho_n, a_u_n, a_v_n});
//         }
//     }
// }


State Solver::primitiveToConservative(double rho, double u, double v, double p, double k, double omega) {
    double rhoE = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
    return State(rho, rho * u, rho * v, rhoE, rho * k, rho * omega);
}
// for Euler
// void Solver::computeResiduals(const std::vector<State>& U_in, std::vector<State>& R_out) {
//     std::fill(R_out.begin(), R_out.end(), State(0,0,0,0));

//     for (const auto& f : mesh.faces) {
//         State UL = U_in[f.owner];
//         State Flux;

//         if (f.neighbor >= 0) {
//             // internal face
//             State UR = U_in[f.neighbor];
//             Flux = computeFluxLLF(UL, UR, f.nx, f.ny);
//             double mag = f.area;
//             R_out[f.owner] = R_out[f.owner] + Flux * mag;
//             R_out[f.neighbor] = R_out[f.neighbor] - Flux * mag;
//         } else {
//             // boundary face
//             State U_ghost = getBoundaryGhostState(UL, f);
//             Flux = computeFluxLLF(UL, U_ghost, f.nx, f.ny);
//             R_out[f.owner] = R_out[f.owner] + Flux * f.area;
//         }
//     }
// }




void Solver::computeResiduals(const std::vector<State>& U_in, std::vector<State>& R_out, bool update_limiters) {
    std::fill(R_out.begin(), R_out.end(), State(0,0,0,0));

    // need gradients first for viscous flux
    computeGradients(U_in);
    if (iSpatialScheme == SpatialScheme::SECOND_ORDER_LIMITED) {
        if (update_limiters){
            computeLimiters(U_in);
        }
    }
    // some consts
    double Cp = gamma * R_gas / (gamma - 1.0);
    double Pr = 0.72; // Prandtl number

    for (const auto& f : mesh.faces) {
        State UL = U_in[f.owner];
        State UR = (f.neighbor >= 0) ? U_in[f.neighbor] : getBoundaryGhostState(UL, f);
        

        //  A. LLF for Euler part (Inviscid Flux)
        // when testing the flat plate, there is an error, so update
        // State FluxEuler = computeFluxLLF(UL, UR, f.nx, f.ny);

        State UL_recon;
        State UR_recon;

        switch (iSpatialScheme){
            case SpatialScheme::FIRST_ORDER:
                UL_recon = UL;
                UR_recon = UR;
                break;

            case SpatialScheme::SECOND_ORDER_CENTRAL:
            case SpatialScheme::SECOND_ORDER_LIMITED:{
                if (f.neighbor >= 0){
                    double dx_L = f.x_mid - mesh.cells[f.owner].x_c;
                    double dy_L = f.y_mid - mesh.cells[f.owner].y_c;
                    double dx_R = f.x_mid - mesh.cells[f.neighbor].x_c;
                    double dy_R = f.y_mid - mesh.cells[f.neighbor].y_c;

                    CellGradient gradL = cellGradients[f.owner];
                    CellGradient gradR = cellGradients[f.neighbor];

                    //double phi = (iSpatialScheme == SpatialScheme::SECOND_ORDER_LIMITED) ? 0.3 : 0.5;
                    double phi_L = (iSpatialScheme == SpatialScheme::SECOND_ORDER_LIMITED) ? limiters[f.owner] : 0.5;
                    double phi_R = (iSpatialScheme == SpatialScheme::SECOND_ORDER_LIMITED) ? limiters[f.neighbor] : 0.5;
                    double uL = UL.rhou / UL.rho;
                    double vL = UL.rhov / UL.rho;

                    // double TL = getPressure(UL) / (UL.rho * R_gas);

                    // double uL_f = uL + phi_L * (gradL.dudx * dx_L + gradL.dudy * dy_L);
                    // double vL_f = vL + phi_L * (gradL.dvdx * dx_L + gradL.dvdy * dy_L);
                    // double TL_f = TL + phi_L * (gradL.dTdx * dx_L + gradL.dTdy * dy_L);

                    // double pL_f = UL.rho * R_gas * TL_f;
                    // UL_recon = primitiveToConservative(UL.rho, uL_f, vL_f, pL_f);

                    // double uR = UR.rhou / UR.rho;
                    // double vR = UR.rhov / UR.rho;
                    // double TR = getPressure(UR) / (UR.rho * R_gas);

                    // double uR_f = uR + phi_R * (gradR.dudx * dx_R + gradR.dudy * dy_R);
                    // double vR_f = vR + phi_R * (gradR.dvdx * dx_R + gradR.dvdy * dy_R);
                    // double TR_f = TR + phi_R * (gradR.dTdx * dx_R + gradR.dTdy * dy_R);
                    // double pR_f = UR.rho * R_gas * TR_f;
                    // --- Left State Reconstruction ---
                    double pL = getPressure(UL);
                    
                    double rhoL_f = UL.rho + phi_L * (gradL.drhodx * dx_L + gradL.drhody * dy_L);
                    double pL_f   = pL     + phi_L * (gradL.dpdx   * dx_L + gradL.dpdy   * dy_L);
                    double uL_f   = uL     + phi_L * (gradL.dudx   * dx_L + gradL.dudy   * dy_L);
                    double vL_f   = vL     + phi_L * (gradL.dvdx   * dx_L + gradL.dvdy   * dy_L);
                    
                    // Safety clip: Prevent unphysical negative pressure/density interpolation
                    if (rhoL_f < 1e-4) rhoL_f = UL.rho; 
                    if (pL_f < 1e-4)   pL_f = pL;

                    UL_recon = primitiveToConservative(rhoL_f, uL_f, vL_f, pL_f);

                    // --- Right State Reconstruction ---
                    double uR = UR.rhou / UR.rho;
                    double vR = UR.rhov / UR.rho;
                    double pR = getPressure(UR);
                    
                    double rhoR_f = UR.rho + phi_R * (gradR.drhodx * dx_R + gradR.drhody * dy_R);
                    double pR_f   = pR     + phi_R * (gradR.dpdx   * dx_R + gradR.dpdy   * dy_R);
                    double uR_f   = uR     + phi_R * (gradR.dudx   * dx_R + gradR.dudy   * dy_R);
                    double vR_f   = vR     + phi_R * (gradR.dvdx   * dx_R + gradR.dvdy   * dy_R);

                    // Safety clip
                    if (rhoR_f < 1e-4) rhoR_f = UR.rho;
                    if (pR_f < 1e-4)   pR_f = pR;

                    UR_recon = primitiveToConservative(rhoR_f, uR_f, vR_f, pR_f);
                    // UR_recon = primitiveToConservative(UR.rho, uR_f, vR_f, pR_f);
                } else {
                    // 如果是物理边界边界，采用一阶（即直接等于 Ghost State 或 Cell Center）
                    UL_recon = UL;
                    UR_recon = UR;
                }
                break;
            }
        }
        // if (f.neighbor >=0){
        //     double dx_L = f.x_mid - mesh.cells[f.owner].x_c;
        //     double dy_L = f.y_mid - mesh.cells[f.owner].y_c;

        //     double dx_R = f.x_mid - mesh.cells[f.neighbor].x_c;
        //     double dy_R = f.y_mid - mesh.cells[f.neighbor].y_c;

        //     CellGradient gradL = cellGradients[f.owner];
        //     CellGradient gradR = cellGradients[f.neighbor];

        //     double phi = 0.5;

        //     double uL = UL.rhou / UL.rho;
        //     double vL = UL.rhov / UL.rho;
        //     double TL = getPressure(UL) / (UL.rho * R_gas);

        //     double uL_f = uL + phi * (gradL.dudx * dx_L + gradL.dudy * dy_L);
        //     double vL_f = vL + phi * (gradL.dvdx * dx_L + gradL.dvdy * dy_L);
        //     double TL_f = TL + phi * (gradL.dTdx * dx_L + gradL.dTdy * dy_L);

        //     double pL_f = UL.rho * R_gas * TL_f;
        //     UL_recon = primitiveToConservative(UL.rho, uL_f, vL_f, pL_f);

        //     double uR = UR.rhou / UR.rho;
        //     double vR = UR.rhov / UR.rho;
        //     double TR = getPressure(UR) / (UR.rho * R_gas);

        //     double uR_f = uR + phi * (gradR.dudx * dx_R + gradR.dudy * dy_R);
        //     double vR_f = vR + phi * (gradR.dvdx * dx_R + gradR.dvdy * dy_R);
        //     double TR_f = TR + phi * (gradR.dTdx * dx_R + gradR.dTdy * dy_R);
        //     double pR_f = UR.rho * R_gas * TR_f;
        //     UR_recon = primitiveToConservative(UR.rho, uR_f, vR_f, pR_f);
        // }
        State FluxEuler;
        if (f.neighbor < 0 && (f.bcTag == 3 || f.bcTag == 4)) {
            // EXTREMELY IMPORTANT FOR NAVIER-STOKES:
            // Do not use a Riemann solver at a solid wall boundary. It creates massive
            // artificial numerical viscosity that destroys the boundary layer.
            // The exact inviscid flux at a solid wall is purely the static pressure.
            double p_wall = getPressure(UL_recon);
            FluxEuler.rho = 0.0;
            FluxEuler.rhou = p_wall * f.nx;
            FluxEuler.rhov = p_wall * f.ny;
            FluxEuler.rhoE = 0.0;
        } else {
            // For interior cells or farfield boundaries, standard LLF is safe.
            FluxEuler = computeInviscidFlux(UL_recon, UR_recon, f.nx, f.ny, iFluxScheme);
        }
        // State FluxEuler = computeFluxLLF(UL_recon, UR_recon, f.nx, f.ny);

        //  B. Viscous Flux (2D Navier-Stokes)
       
        double rho_f = 0.5 * (UL.rho + UR.rho);
        double u_f   = 0.5 * (UL.rhou / UL.rho + UR.rhou / UR.rho);
        double v_f   = 0.5 * (UL.rhov / UL.rho + UR.rhov / UR.rho);

        // using cell-centered gradients to approximate face gradients (simple averaging)
        CellGradient gradL = cellGradients[f.owner];
        CellGradient gradR = (f.neighbor >= 0) ? cellGradients[f.neighbor] : gradL; 

        double dudx_f, dudy_f, dvdx_f, dvdy_f, dTdx_f, dTdy_f;
        if (f.neighbor < 0 && f.bcTag == 4) {
            // NO-SLIP WALL: Compute strict face normal gradients (velocity is 0 at wall).
            // Using a cell-average gradient drastically underestimates wall shear stress.
            double dx = f.x_mid - mesh.cells[f.owner].x_c;
            double dy = f.y_mid - mesh.cells[f.owner].y_c;
            double dist2 = dx * dx + dy * dy;

            double uL_c = UL.rhou / UL.rho;
            double vL_c = UL.rhov / UL.rho;
            
            // gradient = (0 - U_center) / distance 
            dudx_f = -uL_c * dx / dist2;
            dudy_f = -uL_c * dy / dist2;
            dvdx_f = -vL_c * dx / dist2;
            dvdy_f = -vL_c * dy / dist2;

            // Adiabatic wall keeps temperature gradient 0 (as defined by your ghost state)
            // dTdx_f = 0.0;
            // dTdy_f = 0.0;
            double T_wall = 300.0;
            double pL_c = getPressure(UL);
            double TL_c = pL_c / (UL.rho * R_gas);
            dTdx_f = 0.0;
            dTdy_f = 0.0;
            // dTdx_f = (T_wall - TL_c) * dx / dist2;
            // dTdy_f = (T_wall - TL_c) * dy / dist2;
        } else {
            // Simple averaging for interior faces
            dudx_f = 0.5 * (gradL.dudx + gradR.dudx);
            dudy_f = 0.5 * (gradL.dudy + gradR.dudy);
            dvdx_f = 0.5 * (gradL.dvdx + gradR.dvdx);
            dvdy_f = 0.5 * (gradL.dvdy + gradR.dvdy);
            dTdx_f = 0.5 * (gradL.dTdx + gradR.dTdx);
            dTdy_f = 0.5 * (gradL.dTdy + gradR.dTdy);

            double dx_LR = mesh.cells[f.neighbor].x_c - mesh.cells[f.owner].x_c;
            double dy_LR = mesh.cells[f.neighbor].y_c - mesh.cells[f.owner].y_c;
            double dist2 = dx_LR * dx_LR + dy_LR * dy_LR;


            if (dist2 > 1e-14){
                double du = (UR.rhou / UR.rho) - (UL.rhou / UL.rho);
                double dv = (UR.rhov / UR.rho) - (UL.rhov / UL.rho);

                double dudx_0 = dudx_f; double dudy_0 = dudy_f;
                double dvdx_0 = dvdx_f; double dvdy_0 = dvdy_f;

                dudx_f = dudx_0 - (dudx_0 * dx_LR + dudy_0 * dy_LR - du) * dx_LR / dist2;
                dudy_f = dudy_0 - (dudx_0 * dx_LR + dudy_0 * dy_LR - du) * dy_LR / dist2;
                dvdx_f = dvdx_0 - (dvdx_0 * dx_LR + dvdy_0 * dy_LR - dv) * dx_LR / dist2;
                dvdy_f = dvdy_0 - (dvdx_0 * dx_LR + dvdy_0 * dy_LR - dv) * dy_LR / dist2;
            }
        }

        // double dudx_f = 0.5 * (gradL.dudx + gradR.dudx);
        // double dudy_f = 0.5 * (gradL.dudy + gradR.dudy);
        // double dvdx_f = 0.5 * (gradL.dvdx + gradR.dvdx);
        // double dvdy_f = 0.5 * (gradL.dvdy + gradR.dvdy);
        // double dTdx_f = 0.5 * (gradL.dTdx + gradR.dTdx);
        // double dTdy_f = 0.5 * (gradL.dTdy + gradR.dTdy);

        // calculate viscosity and thermal conductivity
        // mu
        double factor_T_inf = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
        double T_static_inf = T_total_inlet / factor_T_inf;
        double rho_inf = P_static_exit / (R_gas * T_static_inf);

        double mu = rho_f * nu;  // dynamic viscosity 
        double k  = mu * Cp / Pr; // thermal conductivity 

        // divergence of V
        double div_V = dudx_f + dvdy_f;

        // Stress Tensor
        double tau_xx = mu * (2.0 * dudx_f - (2.0/3.0) * div_V);
        double tau_yy = mu * (2.0 * dvdy_f - (2.0/3.0) * div_V);
        double tau_xy = mu * (dudy_f + dvdx_f);

        // Heat Flux
        double q_x = -k * dTdx_f;
        double q_y = -k * dTdy_f;

        // f_visc * nx + g_visc * ny
        State FluxVisc;
        FluxVisc.rho = 0.0; // no mass diffusion
        FluxVisc.rhou = tau_xx * f.nx + tau_xy * f.ny;
        FluxVisc.rhov = tau_xy * f.nx + tau_yy * f.ny;

        double work_x = u_f * tau_xx + v_f * tau_xy - q_x;
        double work_y = u_f * tau_xy + v_f * tau_yy - q_y;
        FluxVisc.rhoE = work_x * f.nx + work_y * f.ny;

        //  C. Total Flux = F_Euler - F_Viscous  _______________ for viscous term
        State TotalFlux = FluxEuler - FluxVisc;


        // for inviscous term:

        // State TotalFlux = FluxEuler;

        // accumulate residuals 
        double mag = f.area;
        R_out[f.owner] = R_out[f.owner] + TotalFlux * mag;
        if (f.neighbor >= 0) {
            R_out[f.neighbor] = R_out[f.neighbor] - TotalFlux * mag;
        }
    }
}
// void Solver::initializeTGV(){
//     double V0 = 1.0;
//     double p0 = 100000.0;
//     double rho0 = 1.225;

//     for (int i = 0; i < mesh.cells.size(); ++i){
//         double x = mesh.cells[i].x_c;
//         double y = mesh.cells[i].y_c;

//         //theorytical TGV solution

//         double u = V0 * std::sin(x) * std::cos(y);
//         double v = -V0 * std::cos(x) * std::sin(y);

//         // p
//         double p = p0 + (rho0 * V0 * V0 / 4.0) * (std::cos(2.0 * x) + std::cos(2.0 * y));

//         U_current[i] = primitiveToConservative(rho0, u, v, p);
//     }
//     std::cout << "Field initialized with 2D TGV exact solution" << std::endl;
// }

// for 2d bump
// // set boundary ghost state
// State Solver::getBoundaryGhostState(const State& U_in, const Face& f) {
//     State U_ghost = U_in;

//     // --- BC Tag 1: Subsonic Inlet (Given P0, T0) ---
//     if (f.bcTag == 1) {
        
        
//         double p_in = getPressure(U_in);
        
//         // avoid NaN value
//         double p_use = std::min(p_in, P_total_inlet - 1e-4); 

//         //Isentropic relations:
//         // T = T0 * (P/P0)^((g-1)/g)
//         double T_ghost = T_total_inlet * std::pow(p_use / P_total_inlet, (gamma - 1.0) / gamma);
//         double rho_ghost = p_use / (R_gas * T_ghost);
        
//         // Velocity magnitude from Energy equation: V^2 = 2Cp(T0 - T)
//         double Cp = gamma * R_gas / (gamma - 1.0);
//         double V2 = 2.0 * Cp * (T_total_inlet - T_ghost);
//         double V = (V2 > 0) ? std::sqrt(V2) : 0.0;

       
//         U_ghost = primitiveToConservative(rho_ghost, V, 0.0, p_use);
//     }
//     // --- BC Tag 2: Outlet (Given Static Pressure) ---
//     else if (f.bcTag == 2) {
//         // subsonic outlet P_exit
        
//         double p_target = P_static_exit;
        
//         double rho = U_in.rho;
//         double u = U_in.rhou / U_in.rho;
//         double v = U_in.rhov / U_in.rho;
        
//         // set others
//         U_ghost = primitiveToConservative(rho, u, v, p_target);
//     }
//     // --- BC Tag 3: Slip Wall ---
//     else if (f.bcTag == 3) {
        
//         double u = U_in.rhou / U_in.rho;
//         double v = U_in.rhov / U_in.rho;
        
//         double vn = u * f.nx + v * f.ny;
        
//         // V_ghost
//         double u_g = u - 2.0 * vn * f.nx;
//         double v_g = v - 2.0 * vn * f.ny;
        
//         U_ghost.rho = U_in.rho;
//         U_ghost.rhou = U_in.rho * u_g;
//         U_ghost.rhov = U_in.rho * v_g;
//         U_ghost.rhoE = U_in.rhoE; 
//     }

//     return U_ghost;
// }


// // for wind tunnel
// State Solver::getBoundaryGhostState(const State& U_in, const Face& f) {
//     State U_ghost = U_in;

//     // --- BC Tag 1: Supersonic Inlet ---
//     if (f.bcTag == 1) {
//         double rho_ghost = 1.4;
//         double u_ghost   = 3.0;
//         double v_ghost   = 0.0;
//         double p_ghost   = 1.0;
//         U_ghost = primitiveToConservative(rho_ghost, u_ghost, v_ghost, p_ghost);
//     }
//     // --- BC Tag 2: Supersonic Outlet ---
//     else if (f.bcTag == 2) {
        
//         U_ghost = U_in; 
//     }
//     // --- BC Tag 3: Slip Wall ---
//     else if (f.bcTag == 3) {
//         double u = U_in.rhou / U_in.rho;
//         double v = U_in.rhov / U_in.rho;
        
//         double vn = u * f.nx + v * f.ny;
        
//         
//         double u_g = u - 2.0 * vn * f.nx;
//         double v_g = v - 2.0 * vn * f.ny;
        
//         U_ghost.rho = U_in.rho;
//         U_ghost.rhou = U_in.rho * u_g;
//         U_ghost.rhov = U_in.rho * v_g;
//         U_ghost.rhoE = U_in.rhoE; 
//     }

//     return U_ghost;
// }

// for NASA supersonic test case

// void Solver::initialize() {
//     // 从总温和总压反算静温和静压（或者直接使用设定的 exit 压力）
//     double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
//     double T_static = T_total_inlet / factor_T;
//     double p_static = P_static_exit; 
    
//     double rho = p_static / (R_gas * T_static);
//     double c = std::sqrt(gamma * R_gas * T_static);
//     double u = M_inf * c; // 来流仅有水平速度 (零攻角)
//     double v = 0.0;
    
//     State init = primitiveToConservative(rho, u, v, p_static);
//     std::fill(U_current.begin(), U_current.end(), init);
    
//     std::cout << "Field initialized with Ma = " << M_inf << ", P = " << p_static << std::endl;
// }

// corner case
// void Solver::initialize() {
//     // Freestream conditions for Prandtl-Meyer Expansion
//     this->M_inf = 2.5;
//     double p = 12.0;    // psia
//     double T = 550.0;   // R
//     // double M = 2.5;

//     // Use the solver's existing gamma and R_gas
//     double c = std::sqrt(gamma * R_gas * T); 
//     double u = this->M_inf * c;
//     double v = 0.0;
//     double rho = p / (R_gas * T);
    
//     State init = primitiveToConservative(rho, u, v, p);
//     std::fill(U_current.begin(), U_current.end(), init);
    
//     std::cout << "Field initialized: PM Expansion (Mach 2.5)" << std::endl;
// }
// 在 Solver::initialize() 的最开头强制统一英制常数
// void Solver::initialize() {
//     this->gamma = 1.4;
//     this->R_gas = 1716.48;  // <--- 修改为英制气体常数

//     this->M_inf = 2.5;
//     double p = 12.0 * 144.0; // <--- 注意：12 psia = 12 * 144 = 1728 lb/ft^2 (重要：压力必须使用 psf 来计算密度)
//     double T = 550.0;        // Rankine

//     double c = std::sqrt(gamma * R_gas * T); 
//     double u = this->M_inf * c; 
//     double v = 0.0;
//     double rho = p / (R_gas * T); // 如果 p 是 psia，这算出来的会错，必须换成 psf
    
//     State init = primitiveToConservative(rho, u, v, p);
//     std::fill(U_current.begin(), U_current.end(), init);
    
//     std::cout << "Field initialized: PM Expansion (Mach " << this->M_inf << ")" << std::endl;
// }
void Solver::initialize() {
    this->gamma = 1.4;
    this->R_gas = 1716.48;  

    this->M_inf = 2.5;
    
    // --- 重点：NASA 给的是总压和总温 ---
    double P0_psia = 12.0;    
    double T0_R = 550.0;   
    
    // 根据马赫数反算静压(psia)和静温(Rankine)
    double factor_T = 1.0 + 0.5 * (this->gamma - 1.0) * this->M_inf * this->M_inf;
    double p_psia = P0_psia / std::pow(factor_T, this->gamma / (this->gamma - 1.0));
    double T = T0_R / factor_T;

    // 转化为欧拉方程使用的基础单位 psf (1 psia = 144 psf)
    double p = p_psia * 144.0;
    
    double c = std::sqrt(this->gamma * this->R_gas * T); 
    double u = this->M_inf * c; 
    double v = 0.0;
    double rho = p / (this->R_gas * T); 
    
    State init = primitiveToConservative(rho, u, v, p);
    std::fill(U_current.begin(), U_current.end(), init);
    
    std::cout << "Field initialized: PM Expansion (P_static_psia = " << p_psia << ")" << std::endl;
}

// // For TGV
// State Solver::getBoundaryGhostState(const State& U_in, const Face& f){
//     double x = f.x_mid;
//     double y = f.y_mid;

//     double V0 = 1.0;
//     double p0 = 100000.0;
//     double rho0 = 1.225;

//     double F_decay = std::exp(-2.0 * nu * physical_time);
//     double P_decay = std::exp(-4.0 * nu * physical_time);

//     double u = V0 * std::sin(x) * std::cos(y) * F_decay;
//     double v = -V0 * std::cos(x) * std::sin(y) * F_decay;

//     double p = p0 + (rho0 * V0 * V0 / 4.0) * (std::cos(2.0 * x) + std::cos(2.0 * y)) * P_decay;

//     return primitiveToConservative(rho0, u, v, p);

// }

// for NASA supersonic test case

// State Solver::getBoundaryGhostState(const State& U_in, const Face& f) {
//     State U_ghost = U_in;

//     // --- BC Tag 1: Supersonic Inlet / Frozen Farfield ---
//     if (f.bcTag == 1) {
        
//         double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
//         double T_static = T_total_inlet / factor_T;
//         double rho_inf = P_static_exit / (R_gas * T_static);
//         double c_inf = std::sqrt(gamma * R_gas * T_static);
//         double u_inf = M_inf * c_inf;
        
//         U_ghost = primitiveToConservative(rho_inf, u_inf, 0.0, P_static_exit);
//     }
//     // --- BC Tag 2: Supersonic Extrapolating Outflow ---
//     else if (f.bcTag == 2) {
//         U_ghost = U_in; 
//     }
//     // --- BC Tag 3: Inviscid Slip Wall ---
//     else if (f.bcTag == 3) {
//         double u = U_in.rhou / U_in.rho;
//         double v = U_in.rhov / U_in.rho;
        
        
//         double vn = u * f.nx + v * f.ny;
        
        
//         double u_g = u - 2.0 * vn * f.nx;
//         double v_g = v - 2.0 * vn * f.ny;
        
//         U_ghost.rho = U_in.rho;
//         U_ghost.rhou = U_in.rho * u_g;
//         U_ghost.rhov = U_in.rho * v_g;
//         U_ghost.rhoE = U_in.rhoE; 
//     }
       

//     return U_ghost;
// }


// for laminar flat plate
// State Solver::getBoundaryGhostState(const State& U_in, const Face& f) {
//     State U_ghost = U_in;

//     // --- BC Tag 1 & 5: Subsonic Inlet & Farfield ---
//     if (f.bcTag == 1 || f.bcTag == 5) {
//         // Enforce incoming Mach, T, and P. 
//         double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
//         double T_static = T_total_inlet / factor_T;
//         double rho_inf = P_static_exit / (R_gas * T_static);
//         double c_inf = std::sqrt(gamma * R_gas * T_static);
//         double u_inf = M_inf * c_inf;
        
//         U_ghost = primitiveToConservative(rho_inf, u_inf, 0.0, P_static_exit);
//     }
//     // --- BC Tag 2: Subsonic Outlet ---
//     else if (f.bcTag == 2) {
//         // Enforce fixed pressure, extrapolate velocity and density
//         // double rho = U_in.rho;
//         // double u = U_in.rhou / U_in.rho;
//         // double v = U_in.rhov / U_in.rho;
//         // U_ghost = primitiveToConservative(rho, u, v, P_static_exit);


//         // for mach 5 case:
//         U_ghost = U_in;
//     }
//     // --- BC Tag 3: Slip Wall (Upstream of leading edge) ---
//     else if (f.bcTag == 3) {
//         double u = U_in.rhou / U_in.rho;
//         double v = U_in.rhov / U_in.rho;
//         double vn = u * f.nx + v * f.ny;
        
//         double u_g = u - 2.0 * vn * f.nx;
//         double v_g = v - 2.0 * vn * f.ny;
        
//         U_ghost.rho = U_in.rho;
//         U_ghost.rhou = U_in.rho * u_g;
//         U_ghost.rhov = U_in.rho * v_g;
//         U_ghost.rhoE = U_in.rhoE; 
//     }
//     // --- BC Tag 4: Flat Plate (No-Slip Wall) ---
//     else if (f.bcTag == 4) {
//         // Reverse velocity to enforce zero flow
//         // double u_g = - (U_in.rhou / U_in.rho);
//         // double v_g = - (U_in.rhov / U_in.rho);
        
//         // U_ghost.rho = U_in.rho;
//         // U_ghost.rhou = U_in.rho * u_g;
//         // U_ghost.rhov = U_in.rho * v_g;
//         // U_ghost.rhoE = U_in.rhoE; 

//         // for mach 5 case:
//         double u_g = - (U_in.rhou / U_in.rho);
//         double v_g = - (U_in.rhov / U_in.rho);

//         double ek_in = 0.5 * (U_in.rhou * U_in.rhou + U_in.rhov * U_in.rhov) / U_in.rho;
//         double p_in = (gamma - 1.0) * (U_in.rhoE - ek_in);
//         double T_in = p_in / (U_in.rho * R_gas);

//         // double T_wall = 300.0;
//         // double T_ghost = 2.0 * T_wall - T_in; // mirror temperature for adiabatic wall
//         // if (T_ghost <= 10.0) T_ghost = T_wall;
//         double T_ghost = T_in;
//         double rho_g = p_in / (R_gas * T_ghost);
//         double k_in = 0.0, omega_in = 0.0;

//         U_ghost = primitiveToConservative(rho_g, u_g, v_g, p_in, k_in, omega_in);
//     }

//     return U_ghost;
// }


//corner case

State Solver::getBoundaryGhostState(const State& U_in, const Face& f) {
    State U_ghost = U_in;

    // // --- BC Tag 1: Supersonic Inlet (Frozen Dirichlet) ---
    // if (f.bcTag == 1) {
    //     double p = 12.0;
    //     double T = 550.0;
    //     double M = 2.5;
        
    //     double c = std::sqrt(gamma * R_gas * T);
    //     double rho = p / (R_gas * T);
    //     double u = M * c;
    //     double v = 0.0;
        
    //     U_ghost = primitiveToConservative(rho, u, v, p);
    // }
    // --- BC Tag 1: Supersonic Inlet (Frozen Dirichlet) ---
    // if (f.bcTag == 1) {
    //     this->R_gas = 1716.48;     // 保持一致
    //     double p = 12.0 * 144.0;   // 1728.0 lb/ft^2
    //     double T = 550.0;
    //     double M = 2.5;
        
    //     double c = std::sqrt(gamma * R_gas * T);
    //     double rho = p / (R_gas * T);
    //     double u = M * c;
    //     double v = 0.0;
        
    //     U_ghost = primitiveToConservative(rho, u, v, p);
    // }
    // --- BC Tag 1: Supersonic Inlet (Frozen Dirichlet) ---
    if (f.bcTag == 1) {
        this->R_gas = 1716.48;     
        double P0_psia = 12.0;   
        double T0_R = 550.0;
        double M = 2.5;
        
        double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M * M;
        double p_psia = P0_psia / std::pow(factor_T, gamma / (gamma - 1.0));
        double T = T0_R / factor_T;
        
        double p = p_psia * 144.0;
        
        double c = std::sqrt(gamma * R_gas * T);
        double rho = p / (R_gas * T);
        double u = M * c;
        double v = 0.0;
        
        U_ghost = primitiveToConservative(rho, u, v, p);
    }
    // --- BC Tag 2: Supersonic Outlet (Extrapolation) ---
    else if (f.bcTag == 2) {
        // Zero-gradient for supersonic outflow: ghost state perfectly mirrors inside
        U_ghost = U_in; 
    }
    // --- BC Tag 3: Inviscid Wall (Slip Wall) ---
    else if (f.bcTag == 3) {
        double u = U_in.rhou / U_in.rho;
        double v = U_in.rhov / U_in.rho;
        
        // Wall normal velocity
        double vn = u * f.nx + v * f.ny;
        
        // Reflect velocity across the face normal
        double u_g = u - 2.0 * vn * f.nx;
        double v_g = v - 2.0 * vn * f.ny;
        
        U_ghost.rho = U_in.rho;
        U_ghost.rhou = U_in.rho * u_g;
        U_ghost.rhov = U_in.rho * v_g;
        U_ghost.rhoE = U_in.rhoE; 
    }

    return U_ghost;
}

void Solver::timeStepExplicit(double dt) {
    computeResiduals(U_current, Residuals);
    static int exp_step = 0;
    exp_step++;
    if (exp_step % 100 == 0){
        double currentNorm = norm_L2(Residuals);
        std::cout << "Explict Step" << exp_step << "| Residual Norm:" << currentNorm << std::endl;
    }
    for (int i = 0; i < mesh.cells.size(); ++i) {
        double vol = mesh.cells[i].volume;
        U_current[i] = U_current[i] - Residuals[i] * (dt / vol);
    


        if (U_current[i].rho < 1e-3){
            U_current[i].rho = 1e-3;
        }

        double kenetic_energy = 0.5 * (U_current[i].rhou * U_current[i].rhou + U_current[i].rhov * U_current[i].rhov) / U_current[i].rho;
        double internal_energy = U_current[i].rhoE - kenetic_energy;
        if (internal_energy < 1e-3){
            U_current[i].rhoE = kenetic_energy + 1e-3;
    }
}
}



// Linear Algebra


double Solver::norm_L2(const std::vector<State>& v) const {
    double sum = 0.0;
    for (const auto& s : v) {
        sum += s.rho*s.rho + s.rhou*s.rhou + s.rhov*s.rhov + s.rhoE*s.rhoE + s.rhok*s.rhok + s.rhow*s.rhow;
    }
    return std::sqrt(sum);
}

double Solver::dot_product(const std::vector<State>& a, const std::vector<State>& b) const {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i].rho*b[i].rho + a[i].rhou*b[i].rhou + a[i].rhov*b[i].rhov + a[i].rhoE*b[i].rhoE + a[i].rhok*b[i].rhok + a[i].rhow*b[i].rhow;
    }
    return sum;
}


// residual convergence

double Solver::getTotalVolume() const{
    double totalVol = 0.0;
    for (const auto& c : mesh.cells) {
        totalVol += c.volume;
    }
    return totalVol;
}


// L2 norm

double Solver::norm_L2_VolumeWeighted(const std::vector<State>& R) const {
    double sum_sq = 0.0;
    double total_vol = 0.0;

    for (size_t i = 0; i < R.size(); ++i){
        double vol = mesh.cells[i].volume;
        const auto& r = R[i];

        double rho_res = r.rho / vol;
        double u_res = r.rhou / vol;
        double v_res = r.rhov / vol;
        double E_res = r.rhoE / vol;
        double k_res = r.rhok / vol;
        double w_res = r.rhow / vol;

        double local_sq = rho_res*rho_res + u_res*u_res + v_res*v_res + E_res*E_res + k_res*k_res + w_res*w_res;

        sum_sq += local_sq * vol;
        total_vol += vol;
    }
    if (total_vol == 0) return 0.0;
    return std::sqrt(sum_sq / total_vol);

}

// drag effiecient

void Solver::calculateDrag(){
    double drag = 0.0;

    int wallTag = 3;

    for (const auto& f : mesh.faces){
        if (f.bcTag == wallTag){
            State U = U_current[f.owner];
            double p = getPressure(U);

            double Fx = -p * f.area * f.nx;
            
            drag += Fx;
           
        }
    }

    std::cout << " Drag :" << drag<< std::endl;
}

// ------------------------------------------------------------------
// JFNK Core: Jacobian-Vector Product
// ------------------------------------------------------------------

void Solver::applyJacobian(const std::vector<State>& U_base, 
                           const std::vector<State>& R_base,
                           const std::vector<State>& v, 
                           std::vector<State>& Jv) {
    // 1. Determine perturbation epsilon
    // Simple strategy: eps = 1e-6
    // Origin version: solid value: double epsilon = 1e-5; 
    
    // New version: automatically determine
    // scaled

    double unorm = norm_L2(U_base);
    double vnorm = norm_L2(v);

    double sqrt_machine_eps = 1e-7;

    double epsilon = sqrt_machine_eps * (1.0 + unorm) / (vnorm + 1e-30);
    // double epsilon = 1e-6;




    // Ideally, eps should be scaled: eps = b / ||v||
    // double v_norm = norm_L2(v); 
    // if(v_norm > 1e-12) epsilon = 1e-6 / v_norm;

    // 2. Perturb State: U_pert = U + eps * v
    int N = U_base.size();
    std::vector<State> U_pert(N);
    for(int i=0; i<N; ++i) {
        U_pert[i] = U_base[i] + v[i] * epsilon;
    }

    // 3. Compute Residual at Perturbed State: R(U + eps*v)
    std::vector<State> R_pert(N);
    computeResiduals(U_pert, R_pert, false);

    // 4. Finite Difference: Jv = (R_pert - R_base) / eps
    // Note: We are solving J * dU = -R. If Residual R is defined as RHS... 
    // Usually Implicit Eq: (U_new - U_old)/dt + R(U_new) = 0
    // Jacobian J = I/dt + dR/dU
    // Here we compute dR/dU * v first.
    
    // For pure steady state JFNK (infinite time step, dt -> inf):
    // J = dR/dU.
    // Jv approx [R(U+ev) - R(U)] / e
    
    for(int i=0; i<N; ++i) {
        State term1 = (R_pert[i] - R_base[i]) * (1.0 / epsilon);

        double time_coeff = mesh.cells[i].volume / (dt_local[i] + 1e-30);
        State term2 = v[i] * time_coeff;

        Jv[i] = term1 + term2;
    }
}

void Solver::computeLocalTimeStep(double CFL){
    if (dt_local.size() != mesh.cells.size()){
        dt_local.resize(mesh.cells.size());
    }
    std::vector<double> lambda_sum(mesh.cells.size(), 1e-12);
    for (const auto& f:mesh.faces){
        State UL = U_current[f.owner];
        double pL = getPressure(UL);
        double cL = std::sqrt(gamma * pL / UL.rho);
        double unL = (UL.rhou * f.nx + UL.rhov * f.ny) / UL.rho;

        lambda_sum[f.owner] += (std::abs(unL) + cL) * f.area;

        if (f.neighbor >=0 ){
            State UR = U_current[f.neighbor];
            double pR = getPressure(UR);
            double cR = std::sqrt(gamma * pR / UR.rho);
            double unR = (UR.rhou * (-f.nx) + UR.rhov * (-f.ny)) / UR.rho;

            lambda_sum[f.neighbor] += (std::abs(unR) + cR) * f.area;
        }
    }

    for (size_t i = 0; i < mesh.cells.size(); ++i){
        // State U = U_current[i];
        // double p = getPressure(U);
        // double c = std::sqrt(gamma * p / U.rho);
        // double vel = std::sqrt(U.rhou * U.rhou + U.rhov * U.rhov) / U.rho;

        // double h = std::sqrt(mesh.cells[i].volume);

        // dt_local[i] = CFL * h / (vel + c + 1e-12);
        dt_local[i] = CFL * mesh.cells[i].volume / lambda_sum[i];
    }
}

void Solver::computeDiagonalPreconditioner(std::vector<double>& D) {
    int N = mesh.cells.size();
    D.assign(N, 0.0);

    for (const auto& f : mesh.faces) {
        State UL = U_current[f.owner];
        State UR = (f.neighbor >= 0) ? U_current[f.neighbor] : getBoundaryGhostState(UL, f);
        
        double cL = std::sqrt(gamma * getPressure(UL) / UL.rho);
        double cR = std::sqrt(gamma * (getPressure(UR) > 0 ? getPressure(UR) : 1e-3) / (UR.rho > 0 ? UR.rho : 1e-3));
        double c_f = 0.5 * (cL + cR);

        double u_f = 0.5 * (UL.rhou/UL.rho + UR.rhou/UR.rho);
        double v_f = 0.5 * (UL.rhov/UL.rho + UR.rhov/UR.rho);
        double vn = std::abs(u_f * f.nx + v_f * f.ny);
        
        // 1. 无粘谱半径 (Inviscid Spectral Radius)
        double lambda_inv = vn + c_f;
        
        // 2. 粘性谱半径估算 (Viscous Spectral Radius)
        double dx = std::sqrt(mesh.cells[f.owner].volume); 
        double mu = 0.5 * (UL.rho + UR.rho) * nu;
        double lambda_visc = 2.0 * mu / (0.5 * (UL.rho + UR.rho) * dx + 1e-12);
        
        // LLF 和 粘性通量的对角线贡献
        double lambda_total = (lambda_inv + lambda_visc) * f.area;
        
        // 累加到主对角线
        D[f.owner] += lambda_total * 0.5;
        if (f.neighbor >= 0) {
            D[f.neighbor] += lambda_total * 0.5;
        }
    }

    // 添加源项对角线：伪时间步 (Vol / dt)
    for (int i = 0; i < N; ++i) {
        double vol = mesh.cells[i].volume;
        D[i] = (vol / (dt_local[i] + 1e-30)) + D[i];
        
        // 保底防止除零
        if (D[i] < 1e-12) D[i] = 1e-12;
    }
}

// ------------------------------------------------------------------
// GMRES Solver
// Solves A * x = b, where A is Jacobian matrix operator
// ------------------------------------------------------------------

int Solver::solveGMRES(const std::vector<State>& U_base,
                        const std::vector<State>& R_spatial_base, 
                       const std::vector<State>& RHS, // This is usually -Residual(U_base)
                       
                       std::vector<State>& x) { // x is delta_U (output)
    
    // int m = 15; // Krylov subspace size (Restart parameter)
    int m = 40;
    int max_iter = 2; // Number of restarts (keep small for Newton steps)
    int N = U_base.size();
    
    std::vector<double> diag_M;
    computeDiagonalPreconditioner(diag_M);


    //*****for TGV******** */
    // // Pre-calculate R(U_base) for JFNK usage (it usually equals -RHS)
    // std::vector<State> R_base(N);
    // // Since RHS = -R(U), then R(U) = -RHS
    // for(int i=0; i<N; ++i) R_base[i] = RHS[i] * (-1.0); 
    // ***** end *****



    // Vectors for Arnoldi process
    // V[i] is the i-th basis vector
    std::vector<std::vector<State>> V(m + 1, std::vector<State>(N)); 
    
    // Hessenberg matrix H (m+1 x m) stored in 1D or 2D. 
    // Small matrix, use simple 2D vector
    std::vector<std::vector<double>> H(m + 1, std::vector<double>(m, 0.0));

    // Rotations for QR decomposition
    std::vector<double> sn(m, 0.0);
    std::vector<double> cs(m, 0.0);
    
    // Algorithm start
    // Initial guess x = 0 (delta_U starts at 0)
    std::fill(x.begin(), x.end(), State(0,0,0,0));

    // Residual vector r = b - A*x0 -> since x0=0, r = b (= RHS)
    std::vector<State> r = RHS;
    double beta = norm_L2(r);
    
    if (beta < 1e-12) return 0; // Already converged

    V[0] = r;
    // Normalize V[0] = r / beta
    double invBeta = 1.0 / beta;
    for(int i=0; i<N; ++i) V[0][i] = V[0][i] * invBeta;

    // RHS of the small least squares problem g = [beta, 0, 0, ...]^T
    std::vector<double> g(m + 1, 0.0);
    g[0] = beta;

    // Arnoldi Iteration
    int k = 0;
    for (int j = 0; j < m; ++j) {
        k = j;
        
        // 1. Compute w = A * V[j] (Matrix-Vector Product)
        // Here A is the Jacobian J approx

        std::vector<State> z(N);
        for (int idx = 0; idx < N; ++idx){
            z[idx] = V[j][idx] * (1.0 / diag_M[idx]); // Preconditioning step: z = M^-1 * v
        }
        applyJacobian(U_base, R_spatial_base, z, V[j+1]); // Store result temporarily in V[j+1]

        // 2. Gram-Schmidt Orthogonalization
        for (int i = 0; i <= j; ++i) {
            H[i][j] = dot_product(V[j+1], V[i]); // h_ij = (Av_j, v_i)
            // w = w - h_ij * v_i
            for(int idx=0; idx<N; ++idx) {
                V[j+1][idx] = V[j+1][idx] - V[i][idx] * H[i][j];
            }
        }
        
        H[j+1][j] = norm_L2(V[j+1]); // h_{j+1, j} = ||w||

        // 3. Normalize next basis vector
        if (H[j+1][j] > 1e-12) { // breakdown check
            double invH = 1.0 / H[j+1][j];
            for(int idx=0; idx<N; ++idx) {
                V[j+1][idx] = V[j+1][idx] * invH;
            }
        }

        // 4. Apply Givens Rotation to preserve Hessenberg structure
        // Previous rotations
        for (int i = 0; i < j; ++i) {
            double temp = cs[i] * H[i][j] + sn[i] * H[i+1][j];
            H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
            H[i][j] = temp;
        }
        // New rotation
        double a = H[j][j];
        double b = H[j+1][j];
        // Calculate rotation parameters
        if (std::abs(b) < 1e-14) {
            cs[j] = 1.0; sn[j] = 0.0;
        } else {
            double t = std::sqrt(a*a + b*b);
            cs[j] = a / t;
            sn[j] = b / t;
        }
        // Apply new rotation
        H[j][j] = cs[j] * a + sn[j] * b;
        H[j+1][j] = 0.0;
        
        // Update RHS 'g'
        g[j+1] = -sn[j] * g[j];
        g[j]   = cs[j] * g[j];

        // Check convergence based on residual of small system (g[j+1])
        double GMRESTOL = 1e-2;
        if (std::abs(g[j+1]) < GMRESTOL * beta) break; 
    }

    // Solve upper triangular system Hy = g
    std::vector<double> y(k + 1);
    for (int i = k; i >= 0; --i) {
        y[i] = g[i];
        for (int l = i + 1; l <= k; ++l) {
            y[i] -= H[i][l] * y[l];
        }
        y[i] /= H[i][i];
    }

    // Form solution x = x0 + V_m * y
    // x += sum(y[i] * V[i])
    // for (int i = 0; i <= k; ++i) {
    //     for(int idx=0; idx<N; ++idx) {
    //         x[idx] = x[idx] + V[i][idx] * y[i];
    //     }
    // }

    // add pre-conditioner
    std::fill(x.begin(), x.end(), State(0,0,0,0)); // Start from zero
    for (int i = 0; i <= k; ++i) {
        for(int idx=0; idx<N; ++idx) {
            x[idx] = x[idx] + V[i][idx] * y[i] * (1.0 / diag_M[idx]); // Apply M^-1 to the solution
        }
    }
    
    return k + 1;
}

void Solver::calculateEntropyError() {
    double T_inf = T_total_inlet / (1.0 + 0.5*(gamma-1.0)*M_inf*M_inf);
    double P_inf = P_total_inlet / std::pow(1.0 + 0.5*(gamma-1.0)*M_inf*M_inf, gamma/(gamma-1.0));

    double rho_inf = P_inf / (R_gas * T_inf);
    double S_ref = P_inf / std::pow(rho_inf, gamma);
    double sum_error_sq_vol = 0.0;
    double total_vol = 0.0;

    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        State U = U_current[i];
        double p = getPressure(U);
        double rho = U.rho;
        if (p <= 0 || rho <= 0) continue;

        double S_local = p/std::pow(rho, gamma);
        double error = (S_local - S_ref) / S_ref;
        sum_error_sq_vol += (error * error) * mesh.cells[i].volume;
        total_vol += mesh.cells[i].volume;
    }

    double L2_Entropy_Error = std::sqrt(sum_error_sq_vol / total_vol);

    std::cout << " GRS " << std::endl;
    std::cout << " Cells: " << mesh.cells.size() << std::endl;
    std::cout << " L2 Entropy Error: " << L2_Entropy_Error << std::endl; 
}

// ------------------------------------------------------------------
// Newton Solver Interface
// ------------------------------------------------------------------

// original version
// void Solver::solveImplicit(double tolerance, int maxSteps) {
//     std::cout << "--- Starting JFNK Implicit Solver ---" << std::endl;
//     int N = mesh.cells.size();
//     double CFL = 50.0;
//     // Newton Loop
//     for (int step = 0; step < maxSteps; ++step) {

//         computeLocalTimeStep(CFL);
//         // 1. Calculate Residual R(U^k)
//         computeResiduals(U_current, Residuals);
        
//         double resVM = norm_L2_VolumeWeighted(Residuals);
//         // Check convergence (L2 norm of all density residuals)
//         double resNorm = norm_L2(Residuals);
//         double normalizedRes = resNorm / std::sqrt(static_cast<double>(mesh.cells.size()));
//         std::cout << "Newton Step " << step << " | Abs Residual Norm: " << resNorm << " VM ||R||:" << resVM
//                   << " | Normalized Residual Norm: " << normalizedRes << std::endl;
        
//         if (resVM < tolerance) {
//             std::cout << "Converged!" << std::endl;
//             break;
//         }

//         // 2. Solve Linear System J * dU = -R
//         // RHS = -R
//         std::vector<State> RHS(N);
//         for(int i=0; i<N; ++i) RHS[i] = Residuals[i] * (-1.0);
        
//         std::vector<State> delta_U(N); // Update vector
        
//         // Call GMRES
//         solveGMRES(U_current, Residuals, RHS, delta_U);

//         // 3. Update U^{k+1} = U^k + delta_U
//         // (Optional: perform Line Search here to ensure residual decreases)
//         // Simple full step:
//         double relaxation = 0.5; // Under-relaxation often helps stability initially
//         for(int i=0; i<N; ++i) {
//             U_current[i] = U_current[i] + delta_U[i] * relaxation;
//         }
//     }
// }

void Solver::solveImplicit(double tolerance, int maxSteps) {
    std::cout << "--- Starting JFNK Implicit Solver ---" << std::endl;
    int N = mesh.cells.size();
    
    // give the initial CFL with a reasonable number
    double CFL = 0.1; 
    std::ofstream resFile("residual_history.csv");
    resFile << "Step, Residual_VM\n";
    // Newton Loop
    for (int step = 0; step < maxSteps; ++step) {

        // calculate local time step based on current flow field and CFL condition
        // calculate local time step based on current flow field and CFL condition
        computeLocalTimeStep(CFL);
        
        bool should_update_limiter = true;
        
        // --- 重点：冻结限制器打破 Limiter Ringing ---
        // 在前几十步完成流场大结构建立后，果断停止更新 Limiter 向量
        // 允许牛顿本轮迭代精确下降到机器零
        if (step > 150 || CFL > 10.0) { 
            should_update_limiter = false;
        }

        // 1. Calculate Residual R(U^k)
        computeResiduals(U_current, Residuals, should_update_limiter);
        
        double resVM = norm_L2_VolumeWeighted(Residuals);
        double resNorm = norm_L2(Residuals);
        double normalizedRes = resNorm / std::sqrt(static_cast<double>(N));
        
        std::cout << "Newton Step " << step << " | Abs Residual Norm: " << resNorm 
                  << " | VM ||R||: " << resVM << " | CFL: " << CFL << std::endl;
        resFile << step << ", " << resVM << "\n";

        if (resVM < tolerance) {
            std::cout << "Converged!" << std::endl;
            break;
        }

        // 2. Solve Linear System J * dU = -R
        std::vector<State> RHS(N);
        for(int i = 0; i < N; ++i) {
            // 牛顿方程：(Vol/dt * I + dR/dU) * dU = -R
            // 但我们的 applyJacobian 中已经加了 Vol/dt 项，所以这里精确是 -R
            RHS[i] = Residuals[i] * (-1.0);
        }
        
        std::vector<State> delta_U(N); 
        
        // Call GMRES
        int gmres_iters = solveGMRES(U_current, Residuals, RHS, delta_U);
        std::cout << "  GMRES Iters: " << gmres_iters << std::endl;
        if (gmres_iters >= 29){
            CFL *= 0.5;
            if (CFL <0.05) CFL = 0.05;
        }
        // ---------------- 替换从第1682行开始的 Relaxation 逻辑 ----------------
        // 3. Update U^{k+1} = U^k + delta_U
        // 对于平滑的膨胀波，直接使用全步长牛顿更新，只有在残差极大时稍微截断
        double relaxation = 1.0; 
        if (resVM > 1e6) {
            relaxation = 0.5; // 只在启动初期缓冲一下
        } else if (resVM > 2e3 && resVM < 1e2){
            relaxation = 0.9;
        }
        
        for(int i = 0; i < N; ++i) {
            State U_new = U_current[i] + delta_U[i] * relaxation;
            
            // 强下溢保护
            double k_energy = 0.5 * (U_new.rhou * U_new.rhou + U_new.rhov * U_new.rhov) / (U_new.rho + 1e-12);
            double p_new = (gamma - 1.0) * (U_new.rhoE - k_energy);
            
            if (U_new.rho < 1e-5) U_new.rho = 1e-5;
            // if (U_new.rho < 1e-2 || p_new < 1e-2) {
            //     U_new = U_current[i] + delta_U[i] * 0.05; 
            //     if (U_new.rho < 1e-2) U_new.rho = 1e-2;
            //     double k_e2 = 0.5 * (U_new.rhou * U_new.rhou + U_new.rhov * U_new.rhov) / U_new.rho;
            //     if ((U_new.rhoE - k_e2) < 1e-2) U_new.rhoE = k_e2 + 1e-2;
            // }
            if (p_new < 1e-5){
                p_new = 1e-5;
                U_new.rhoE = p_new / (gamma - 1.0) + k_energy;
            }
            U_current[i] = U_new;
        }

        // ---------------- 替换从第1732行开始的 CFL Ramping 逻辑 ----------------
        if (step > 5) {
            // 平滑流场非常适合大步长激进增长
            if (resVM > 1e6) {
                CFL *= 1.2; 
            } else {
                CFL *= 1.5; // 残差一旦回落，指数级攀升 CFL 寻找稳态
            }
            if (gmres_iters >= 25){
                CFL *= 0.5;
            }
            
            // 取消之前硬编码的 50 的限制，允许拉高到 10000.0 以上
            double max_cfl = 10000; 
            
            if (CFL > max_cfl) {
                CFL = max_cfl;
            }
            if (CFL < 0.1) CFL = 0.1;
        }
        // 3. Update U^{k+1} = U^k + delta_U
        // 对于 2.5 马赫，激波刚形成时更新步很大，我们需要保守一点：
        // double relaxation = 1.0; 
        // if (resVM > 1e7){
        //     relaxation = 0.1;
        // } else if (resVM > 2e6){
        //     relaxation = 0.3;
        // } else if (resVM > 5e5) {
        //     relaxation = 0.5;
        // } else if (resVM > 1e5) {
        //     relaxation = 0.8;
        // } else {
        //     relaxation = 1.0; // FULL NEWTON UPDATE
        // }
        // // if (step > 100) relaxation = 1.0;
        // // if (step > 50) relaxation = 1.0; // 后期激波稳定后允许全步长更新

        // for(int i = 0; i < N; ++i) {
        //     State U_new = U_current[i] + delta_U[i] * relaxation;
            
        //     // --- 强下溢保护：防止密度和压力跨越零 ---
        //     double k_energy = 0.5 * (U_new.rhou * U_new.rhou + U_new.rhov * U_new.rhov) / (U_new.rho + 1e-12);
        //     double p_new = (gamma - 1.0) * (U_new.rhoE - k_energy);
            
        //     // 如果更新后导致物理量变得荒谬，我们就对增量进行截断 (或者放弃该网格的面更新)
        //     if (U_new.rho < 1e-2 || p_new < 1e-2) {
        //         // 非常极端的阻尼，几乎只让他随显式残差缓慢蠕动
        //         U_new = U_current[i] + delta_U[i] * 0.05; 
                
        //         // 再兜底一次
        //         if (U_new.rho < 1e-2) U_new.rho = 1e-2;
        //         double k_e2 = 0.5 * (U_new.rhou * U_new.rhou + U_new.rhov * U_new.rhov) / U_new.rho;
        //         if ((U_new.rhoE - k_e2) < 1e-2) U_new.rhoE = k_e2 + 1e-2;
        //     }
            
        //     U_current[i] = U_new;
        // }

        // // --- CFL Ramping ---
        // // if (step > 5){
        // //     if (resVM < 1e-3 && CFL < 10){
        // //         CFL *= 1.1;
        // //     } else if (CFL < 5.0){
        // //         CFL *= 1.05;
        // //         double max_cfl = (iSpatialScheme == SpatialScheme::FIRST_ORDER)? 5.0 : 0.8;
        // //         if (CFL > max_cfl){
        // //             CFL = max_cfl;
        // //         }
        // //         if (CFL > 1.0) CFL = 1.0;
        // //     }
        // // }
        // // --- CFL Ramping ---
        // // --- CFL Ramping ---
        // if (step > 5) {
        //     if (resVM > 1e7) {
        //         // EXTREMELY HIGH RESIDUAL: Drop CFL drastically to recover physics!
        //         CFL *= 0.5;
        //     } else if (resVM > 2e6) {
        //         // High but stable residual: hold steady 
        //         CFL *= 1.0;
        //     } else if (resVM > 5e5) {
        //         CFL *= 1.05;
        //     } else if (resVM > 1e5) {
        //         CFL *= 1.2;
        //     } else {
        //         CFL *= 2.0;    // Once transient initial shock passes, ramp wildly!
        //     }
            
        //     double max_cfl = (iSpatialScheme == SpatialScheme::FIRST_ORDER) ? 500.0 : 50.0;
        //     if (CFL > max_cfl) {
        //         CFL = max_cfl;
        //     }
        //     if (CFL < 0.05) CFL = 0.05; // Prevent CFL from collapsing completely
        // }
        // if (step > 500){
        //     if (gmres_iters >= 25 || resVM > 1e7){
        //         CFL *= 0.5;
        //     }else if (gmres_iters < 15 && resVM < 1e5){
        //         CFL *= 1.1;
        //     }

        //     double max_cfl = (iSpatialScheme == SpatialScheme::FIRST_ORDER)? 3.0 : 0.8;
        //     if (CFL > max_cfl) CFL = max_cfl;
        //     if (CFL < 0.01) CFL = 0.01;
        // }
        // if (step > 50){
        //     if (resVM > 2e6 || gmres_iters >= 29){
        //         CFL *= 0.5;
        //         relaxation = 0.3;
        //         if (CFL < 0.05) CFL = 0.05;
        //     } else if (gmres_iters < 20 && resVM < 5e5 ) {
        //         CFL *= 1.1;
        //         relaxation = 1.0;
        //     }
        //     double max_cfl = (iSpatialScheme == SpatialScheme::FIRST_ORDER)? 1000.0 : 100.0;
        //     if (CFL > max_cfl) CFL = max_cfl;
        //     if (CFL < 0.1) CFL = 0.1;
        // }
        // if (resVM < 2e-6){
        //     std::cout << "Judged as successful" << std::endl;
        //     break;but
        // }
    }
}

void Solver::solveUnsteady(double dt_phy, int maxTimeSteps) {
    std::cout << "--- Starting JFNK Unsteady Solver (Dual Time Stepping) ---" << std::endl;
    int N = mesh.cells.size();
    
    
    if (dt_local.size() != N) dt_local.resize(N);
    std::fill(dt_local.begin(), dt_local.end(), dt_phy); 
    
    U_old.resize(N);

    
    for (int t = 1; t <= maxTimeSteps; ++t) {
        U_old = U_current; 
        physical_time += dt_phy;
        int maxNewton = 15;
        double tolerance = 1e-6;
        
       
        for (int step = 0; step < maxNewton; ++step) {
            
            computeResiduals(U_current, Residuals); 
            
            std::vector<State> RHS(N);
            double norm_sq = 0.0;
            double total_vol = 0.0;
            
            for(int i = 0; i < N; ++i) {
                double vol = mesh.cells[i].volume;
                
                State dU_dt = (U_current[i] - U_old[i]) * (vol / dt_phy);
                State F_unsteady = dU_dt + Residuals[i];
                
                RHS[i] = F_unsteady * (-1.0); 
                
                double rho_res = F_unsteady.rho / vol;
                double u_res = F_unsteady.rhou / vol;
                double v_res = F_unsteady.rhov / vol;
                double E_res = F_unsteady.rhoE / vol;
                norm_sq += (rho_res*rho_res + u_res*u_res + v_res*v_res + E_res*E_res) * vol;
                total_vol += vol;
            }
            
            double resVM = std::sqrt(norm_sq / total_vol);
            
            if (resVM < tolerance && step > 0) {
                std::cout << "  [Time " << t * dt_phy << "] Converged in " << step << " Newton steps." << std::endl;
                break;
            }
            if (step == maxNewton - 1) {
                std::cout << "  [Time " << t * dt_phy << "] Reached max Newton steps: VM ||R|| = " << resVM << std::endl;
            }

            std::vector<State> delta_U(N);
            
            
            solveGMRES(U_current, Residuals, RHS, delta_U);

            for(int i = 0; i < N; ++i) {
                 
                U_current[i] = U_current[i] + delta_U[i] * 0.7;
            }
        }
    }
}

void Solver::extractBoundaryLayerProfile(double x_target, const std::string& filename) {
    std::ofstream outFile(filename);
    outFile << "y,u,v,eta\n";

    double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
    double T_static = T_total_inlet / factor_T;
    double c_inf = std::sqrt(gamma * R_gas * T_static);
    double u_inf = M_inf * c_inf;
    double rho_inf = P_static_exit / (R_gas * T_static);

    double x_tol = 0.005; 

    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        double x_c = mesh.cells[i].x_c;
        double y_c = mesh.cells[i].y_c;

        if (std::abs(x_c - x_target) < x_tol) {
            State U = U_current[i];
            double u = U.rhou / U.rho;
            double v = U.rhov / U.rho;

            
            // eta = y * sqrt(U_inf / (nu * x))
            double x_safe = std::max(x_target, 1e-6); 
            double eta = y_c * std::sqrt(u_inf / (nu * x_safe));

            
            outFile << y_c << "," << u / u_inf << "," << v / u_inf << "," << eta << "\n";
        }
    }
    std::cout << "Boundary layer profile at x = " << x_target << " exported to " << filename << std::endl;
}

void Solver::extractSkinFriction(const std::string& filename) {
    std::ofstream outFile(filename);
    outFile << "x,Cf,Cf_theory\n";

    // get freestream conditions for Cf calculation
    double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
    double T_static = T_total_inlet / factor_T;
    double c_inf = std::sqrt(gamma * R_gas * T_static);
    double u_inf = M_inf * c_inf;
    double rho_inf = P_static_exit / (R_gas * T_static);
    double q_inf = 0.5 * rho_inf * u_inf * u_inf; // inlet dynamic pressure
    double mu_inf = rho_inf * nu;                 // reference dynamic viscosity

    struct WallData {
        double x;
        double Cf;
        double Cf_theory;
    };
    std::vector<WallData> cfCurve;

    for (const auto& f : mesh.faces) {
        // bcTag == 4 
        if (f.neighbor < 0 && f.bcTag == 4) {
            State UL = U_current[f.owner];            
            
            // 1. 获取紧贴壁面那层网格中心的流速
            double u_c = UL.rhou / UL.rho;
            double v_c = UL.rhov / UL.rho;
            
            // 2. 计算流向与壁面平行方向（切向）的速度大小
            // 定义切向量 t = (-ny, nx)
            double tx = -f.ny;
            double ty =  f.nx;
            if (tx < 0) { // 强行让切向沿着 x 的正方向
                tx = -tx;
                ty = -ty;
            }
            
            // 这个速度带有正负号，正值表示纯流向，负值表示回流（分离）
            double u_tangent = u_c * tx + v_c * ty; 
            
            // 3. 精确计算从网格形心 (x_c, y_c) 到壁面 (面中心 x_mid, y_mid) 的纯法向距离 dn
            double dx = f.x_mid - mesh.cells[f.owner].x_c;
            double dy = f.y_mid - mesh.cells[f.owner].y_c;
            
            // 使用平面几何中的点到直线距离公式: dn = |dx * nx + dy * ny|
            double dn = std::abs(dx * f.nx + dy * f.ny);
            if (dn < 1e-12) dn = 1e-12; // 防止浮点数分母过小除零
            
            // 4. 使用标准的一阶差分计算壁面法向梯度: dudy = (u_tangent - u_wall) / dn
            // 因为不可滑移壁面 u_wall = 0
            double dudn = u_tangent / dn;
            
            // 5. 计算表面剪切应力和 Cf
            double mu = UL.rho * nu; 
            double tau_w = mu * dudn; // 这里如果出现回流，dudn 会是负的，刚好得到负的 Cf
            
            double Cf = tau_w / q_inf;

            // calculate Blasius Cf
            double x_local = f.x_mid;
            double Re_x = (u_inf * x_local) / nu;
            double Cf_theory = 0.0;
            if (Re_x > 0) {
                Cf_theory = 0.664 / std::sqrt(Re_x);
            }

            cfCurve.push_back({x_local, Cf, Cf_theory});
        }
    }

    std::sort(cfCurve.begin(), cfCurve.end(), [](const WallData& a, const WallData& b) {
        return a.x < b.x;
    });

    for (const auto& data : cfCurve) {
        outFile << data.x << "," << data.Cf << "," << data.Cf_theory << "\n";
    }

    std::cout << "Skin friction coefficient distribution exported to " << filename << std::endl;
}

void Solver::computeWallDistances() {
    int N = mesh.cells.size();
    wallDistances.assign(N, 1e10); // 初始化为非常大的数

    // 收集所有的壁面 (比如 Tag 4 是粘性壁面)
    std::vector<Face> viscousWalls;
    for (const auto& f : mesh.faces) {
        if (f.neighbor < 0 && f.bcTag == 4) { // 请根据湍流网格的设定调整Tag
            viscousWalls.push_back(f);
        }
    }

    if (viscousWalls.empty()) {
        std::cout << "Warning: No viscous walls found for distance calculation." << std::endl;
        return;
    }

    // 暴力法遍历计算每个 cell 到壁面的最短距离 
    // (对于 2D 问题且网格只有几万的话，O(N*N_walls) 初始化大概瞬间就跑完了)
    for (int i = 0; i < N; ++i) {
        double min_d2 = 1e20;
        double cx = mesh.cells[i].x_c;
        double cy = mesh.cells[i].y_c;
        
        for (const auto& f : viscousWalls) {
            // 这里用面中心作为近似，严格来说应该计算点到线段的垂距
            double dx = cx - f.x_mid;
            double dy = cy - f.y_mid;
            double d2 = dx*dx + dy*dy;
            if (d2 < min_d2) min_d2 = d2;
        }
        wallDistances[i] = std::sqrt(min_d2);
    }
    std::cout << "Wall distances computed for SST model." << std::endl;
}

void Solver::setFreestreamStatic(double M_inf_in, double P_static, double T_static) {
    M_inf = M_inf_in;
    double factor_T = 1.0 + 0.5 * (gamma - 1.0) * M_inf * M_inf;
    double factor_P = std::pow(factor_T, gamma / (gamma - 1.0));

    T_total_inlet = T_static * factor_T;
    P_total_inlet = P_static * factor_P;
    P_static_exit = P_static;
}

void Solver::extractWallPressure(const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open wall pressure file for writing." << std::endl;
        return;
    }

    // 表头：包含 x(inch), p(psi) 为了直接和图7对照
    outFile << "x_in, p_psi\n";

    // 用于收集和排序各个被采样点
    struct SamplePoint {
        double x;
        double p;
        bool operator<(const SamplePoint& other) const {
            return x < other.x;
        }
    };
    std::vector<SamplePoint> samples;

    for (const auto& f : mesh.faces) {
        if (f.bcTag == 3) {
            int c_idx = f.owner;
            double y_c = mesh.cells[c_idx].y_c;
            
            // 过滤出真正代表 "Expansion fan centerline / wall" 的底部边界单元
            // 考虑到上边界的 y_c 接近 1.0 (英尺)，底部接近 0.0 到 -0.32
            if (y_c < 0.2) {  
                double x_ft = f.x_mid; // 网格自带的通常是 ft
                double x_in = x_ft * 12.0; // 转换为 inch
                
                // 获取物面邻近网格的绝对压力 (psf)
                double p_psf = getPressure(U_current[c_idx]);
                
                // 将 psf (lb/ft^2) 转换回 psia (lb/in^2) 
                double p_psi = p_psf / 144.0;
                
                samples.push_back({x_in, p_psi});
            }
        }
    }

    // 按照 X 坐标从入口到出口进行排序
    std::sort(samples.begin(), samples.end());

    // 写入文件
    for (const auto& sp : samples) {
        outFile << sp.x << ", " << sp.p << "\n";
    }

    std::cout << "Centerline/Wall pressures extracted to: " << filename << std::endl;
}