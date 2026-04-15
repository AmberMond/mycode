#ifndef LINALG_HPP
#define LINALG_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include "Types.hpp" 

// 1. result = a + b
inline std::vector<State> vecAdd(const std::vector<State>& a, const std::vector<State>& b) {
    std::vector<State> res(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + b[i]; 
    }
    return res;
}

// 2. (AXPY): y = alpha * x + y
inline void axpy(double alpha, const std::vector<State>& x, std::vector<State>& y) {
    for (size_t i = 0; i < x.size(); ++i) {
        y[i] = y[i] + x[i] * alpha;
    }
}

// 3. x = alpha * x
inline void scale(double alpha, std::vector<State>& x) {
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = x[i] * alpha;
    }
}

// 4. (Dot Product): result = sum(a_i * b_i)

inline double dotProduct(const std::vector<State>& a, const std::vector<State>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        // State dot product: rho*rho + rhou*rhou + ...
        sum += a[i].rho * b[i].rho + 
               a[i].rhou * b[i].rhou + 
               a[i].rhov * b[i].rhov + 
               a[i].rhoE * b[i].rhoE;
    }
    return sum;
}

// 5. L2 norm: ||x||_2
inline double norm2(const std::vector<State>& x) {
    return std::sqrt(dotProduct(x, x));
}

#endif // LINALG_HPP