#include <iostream>
#include "matrix.h"
#include "RREF.h"

int main() {
    // Identity matrix has repeated eigenvalue 1
    std::vector<std::vector<double>> data = {{1, 0}, {0, 1}};
    Matrix<double> A(data);
    
    auto eig = A.eigen();
    
    std::cout << "Eigenvalues: ";
    for(auto v : eig.eigenvalues) std::cout << v << " ";
    std::cout << std::endl;
    
    std::cout << "Number of eigenvectors: " << eig.eigenvectors.size() << std::endl;
    for(size_t i=0; i<eig.eigenvectors.size(); i++){
        std::cout << "v" << i << ": ";
        eig.eigenvectors[i].print();
    }
    
    return 0;
}
