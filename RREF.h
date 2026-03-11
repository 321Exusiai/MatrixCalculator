#pragma once

#include "matrix.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> 

template <typename T>
class RREF {
private:
    Matrix<T> mat;
    size_t rank = 0;
    std::vector<size_t> pivotCols;
    std::vector<size_t> pivotRows;
    bool isREF = false;
    bool isRREF = false;

public:
    RREF(const Matrix<T>& inputMat) : mat(inputMat), rank(0) {}

    void toREF(T eps = static_cast<T>(1e-9)) {
        size_t rows = mat.getRows();
        size_t cols = mat.getCols();
        size_t pivotRow = 0;
        rank = 0;
        pivotCols.clear();
        pivotRows.clear();

        for (size_t col = 0; col < cols && pivotRow < rows; col++) {
            // 鲁棒的主元选择：寻找当前列绝对值最大的行
            size_t max_index = pivotRow;
            T max_val = std::abs(mat.at(pivotRow, col));
            for (size_t row = pivotRow + 1; row < rows; row++) {
                T current_val = std::abs(mat.at(row, col));
                if (current_val > max_val) {
                    max_val = current_val;
                    max_index = row;
                }
            }

            if (max_val < eps) continue;

            if (max_index != pivotRow) {
                mat.exchangeRows(max_index, pivotRow);
            }

            rank++;
            pivotCols.push_back(col);
            pivotRows.push_back(pivotRow);

            for (size_t row = pivotRow + 1; row < rows; row++) {
                if (std::abs(mat.at(row, col)) < eps) {
                    mat.at(row, col) = 0;
                    continue;
                }
                T factor = -mat.at(row, col) / mat.at(pivotRow, col);
                mat.addScaledRow(row, pivotRow, factor);
                mat.at(row, col) = 0; // 强制置零避免数值误差
            }
            pivotRow++;
        }
        isREF = true;
    }

    void toRREF(T eps = static_cast<T>(1e-9)) {
        size_t rows = mat.getRows();
        size_t cols = mat.getCols();
        if (!isREF) toREF(eps);

        for (size_t i = 0; i < rank; i++) {
            size_t row = pivotRows[i];
            size_t col = pivotCols[i];
            T scaleFactor = static_cast<T>(1) / mat.at(row, col);
            mat.scaleRow(row, scaleFactor);
        }

        for (size_t i = rank; i > 0; i--) {
            size_t row = pivotRows[i - 1];
            size_t col = pivotCols[i - 1];
            for (size_t upperRow = row; upperRow > 0; upperRow--) {
                size_t actualUpperRow = upperRow - 1;
                if (std::abs(mat.at(actualUpperRow, col)) < eps) {
                    mat.at(actualUpperRow, col) = 0;
                    continue;
                }
                T factor = -mat.at(actualUpperRow, col);
                mat.addScaledRow(actualUpperRow, row, factor);
                mat.at(actualUpperRow, col) = 0;
            }
        }

        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                if (std::abs(mat.at(i, j)) < eps) mat.at(i, j) = 0;
            }
        }
        isRREF = true;
    }

    void displayResult() const {
        mat.display();
        std::cout << "Rank of the matrix is: " << rank << std::endl;
    }

    size_t getRank() const noexcept { return rank; }
    const Matrix<T>& getMatrix() const noexcept { return mat; }
    const std::vector<size_t>& getPivotCols() const noexcept { return pivotCols; }
    const std::vector<size_t>& getPivotRows() const noexcept { return pivotRows; }

    void setMatrix(const Matrix<T>& m) {
        mat = m;
        rank = 0;
        pivotCols.clear();
        pivotRows.clear();
        isREF = false;
        isRREF = false;
    }

    std::vector<Vector<T>> getKernel(T eps = static_cast<T>(1e-9)) {
        if(!isRREF) toRREF(eps);

        size_t n = mat.getCols();
        std::vector<Vector<T>> basis;
        std::vector<size_t> freeCols;
 
        for(size_t j = 0; j < n; j++){
            bool isPivot = false;
            for(size_t pc : pivotCols){
                if(pc == j){ isPivot = true; break; }
            }
            if(!isPivot) freeCols.push_back(j);
        }
        for (size_t freeCol : freeCols) {
            std::vector<T> v_vec(n, 0);
            v_vec[freeCol] = 1;
            for (size_t i = 0; i < rank; i++) {
                size_t pCol = pivotCols[i];
                size_t pRow = pivotRows[i];
                v_vec[pCol] = -mat.at(pRow, freeCol);
            }
            basis.push_back(Vector<T>(std::move(v_vec)));
        }
        return basis;
    }
};

// =========================================================
// 延迟实现的 Matrix<T> 成员函数
// ---------------------------------------------------------
// 架构原因: Matrix::rank() 和 Matrix::eigen() 需要完整的
// RREF<T> 定义，而 RREF<T> 又依赖完整的 Matrix<T> 定义，
// 形成循环依赖。解决方案:
//   1. matrix.h 中前置声明 RREF<T>，仅声明这两个方法
//   2. 本文件在 RREF<T> 定义之后提供实现
//   3. matrix.h 末尾 #include "RREF.h" 确保调用方自动获得实现
// =========================================================

template <typename T>
int Matrix<T>::rank() const {
    RREF<T> rrefSolver(*this);
    rrefSolver.toRREF();
    return static_cast<int>(rrefSolver.getRank());
}

template <typename T>
typename Matrix<T>::EigenDecomposition Matrix<T>::eigen(int max_iter) const {
    if (!isSquare()) throw std::logic_error("Eigen decomposition only for square matrices");
    
    Matrix<T> A_iter = *this;
    for(size_t k=0; k<static_cast<size_t>(max_iter); k++) {
        std::pair<Matrix<T>, Matrix<T>> qr = A_iter.qr_decomposition();
        Matrix<T> Q = qr.first;
        Matrix<T> R = qr.second;
        A_iter = R * Q;
    }

    EigenDecomposition result;
    for(size_t i=0; i<rows; i++) result.eigenvalues.push_back(A_iter.at(i, i));

    for (T lam : result.eigenvalues) {
        Matrix<T> LambdaI = Matrix<T>::identity(rows) * lam;
        Matrix<T> CharacteristicMat = (*this) - LambdaI;
        
        RREF<T> solver(CharacteristicMat);
        std::vector<Vector<T>> basis = solver.getKernel();

        if (!basis.empty()) {
            for(auto& v : basis) {
                try { result.eigenvectors.push_back(v.normalized()); } 
                catch (...) { result.eigenvectors.push_back(v); }
            }
        } else {
            result.eigenvectors.push_back(Vector<T>(rows, 0)); 
        }
    }
    return result;
}