#include <iostream>
#include <cassert>
#include <cmath>
#include "matrix.h"
#include "RREF.h"

void testSymmetricDiagonalization() {
    // 2x2 Symmetric Matrix: A = {{2, 1}, {1, 2}}
    // Eigenvalues should be 3 and 1
    std::vector<std::vector<double>> data = {{2, 1}, {1, 2}};
    Matrix<double> A(data);
    
    assert(A.isSymmetric());
    assert(A.isDiagonalizable());
    
    auto res = A.diagonalize();
    
    // Check D
    // Values might be 3 and 1 or 1 and 3.
    double d1 = res.D.at(0,0);
    double d2 = res.D.at(1,1);
    assert((std::abs(d1 - 3.0) < 1e-7 && std::abs(d2 - 1.0) < 1e-7) ||
           (std::abs(d1 - 1.0) < 1e-7 && std::abs(d2 - 3.0) < 1e-7));
    
    // Check A = P * D * P^T
    Matrix<double> reconstructed = res.P * res.D * res.P.transpose();
    for(size_t i=0; i<2; i++) {
        for(size_t j=0; j<2; j++) {
            assert(std::abs(reconstructed.at(i, j) - A.at(i, j)) < 1e-7);
        }
    }
    std::cout << "Symmetric diagonalization test passed!" << std::endl;
}

void testNonSymmetricDiagonalization() {
    // 2x2 Non-symmetric but diagonalizable: A = {{1, 1}, {0, 2}}
    // Eigenvalues: 1, 2
    std::vector<std::vector<double>> data = {{1, 1}, {0, 2}};
    Matrix<double> A(data);
    
    assert(!A.isSymmetric());
    assert(A.isDiagonalizable());
    
    auto res = A.diagonalize();
    
    // Check A = P * D * P_inv
    Matrix<double> reconstructed = res.P * res.D * res.P.getInverseMatrix();
    for(size_t i=0; i<2; i++) {
        for(size_t j=0; j<2; j++) {
            assert(std::abs(reconstructed.at(i, j) - A.at(i, j)) < 1e-7);
        }
    }
    std::cout << "Non-symmetric diagonalization test passed!" << std::endl;
}

int main() {
    try {
        testSymmetricDiagonalization();
        testNonSymmetricDiagonalization();
    } catch(const std::exception& e) {
        std::cerr << "Test failed with error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
