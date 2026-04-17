#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <type_traits>
#include "vector.h"

// 前置声明 RREF 类，解决循环依赖
template <typename T> class RREF;

template <typename T>
class Matrix {
private:
    size_t rows;
    size_t cols;
    std::vector<std::vector<T>> data;

public:
    struct EigenDecomposition {
        std::vector<T> eigenvalues;
        std::vector<Vector<T>> eigenvectors;
    };

    struct DiagonalizationResult {
        Matrix<T> P;
        Matrix<T> D;
    };

    template <typename U>
    friend class RREF;

    // -------- Constructors --------
    Matrix(size_t r, size_t r3)
        : rows(r), cols(r3), data(r, std::vector<T>(r3, T())) {
        if (r == 0 || r3 == 0) {
            throw std::invalid_argument("Matrix dimensions must be positive");
        }
    }

    Matrix() : rows(0), cols(0) {}

    Matrix(const std::vector<std::vector<T>>& v) {
        if (v.empty() || v[0].empty())
            throw std::invalid_argument("Input vector must be non-empty");

        rows = static_cast<int>(v.size());
        cols = static_cast<int>(v[0].size());
        for (auto& row : v) {
            if (row.size() != cols)
                throw std::invalid_argument("All rows must have the same length");
        }
        data = v;
    }

    Matrix(std::vector<std::vector<T>>&& v) noexcept {
        if (v.empty() || v[0].empty()) {
            rows = 0; cols = 0; return;
        }
        rows = static_cast<int>(v.size());
        cols = static_cast<int>(v[0].size());
        data = std::move(v);
    }

    // 从 Vector 构造列矩阵 (n x 1)
    Matrix(const Vector<T>& v) : rows(v.size()), cols(1), data(v.size(), std::vector<T>(1)) {
        for (size_t i = 0; i < v.size(); ++i) {
            data[i][0] = v[i];
        }
    }

    // 移动语义
    Matrix(Matrix&& other) noexcept 
        : rows(other.rows), cols(other.cols), data(std::move(other.data)) {
        other.rows = 0; other.cols = 0;
    }

    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            data = std::move(other.data);
            rows = other.rows;
            cols = other.cols;
            other.rows = 0;
            other.cols = 0;
        }
        return *this;
    }

    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default;

    // -------- Static Constructors --------
    static Matrix identity(int n) {
        Matrix<T> mat(n, n);
        for (int i = 0; i < n; i++) {
            mat.at(i, i) = 1;
        }
        return mat;
    }

    static Matrix zero(int n) {
        Matrix<T> mat(n, n);
        return mat;
    }

    // -------- Basic Accessors --------
    size_t getRows() const noexcept { return rows; }
    size_t getCols() const noexcept { return cols; }

    T& at(size_t r, size_t c) {
        if (r >= rows || c >= cols)
            throw std::out_of_range("Matrix index out of bounds");
        return data[r][c];
    }

    const T& at(size_t r, size_t c) const {
        if (r >= rows || c >= cols)
            throw std::out_of_range("Matrix index out of bounds");
        return data[r][c];
    }

    // -------- Printing --------
    void display() const {
        std::cout << "\033[36m" << "Matrix (" << rows << "x" << cols << "):" << "\033[0m" << std::endl;
        // 顶部边框
        std::cout << "  \u250c";
        for (size_t j = 0; j < cols; j++) std::cout << "           ";
        std::cout << " \u2510" << std::endl;

        for (size_t i = 0; i < rows; i++) {
            std::cout << "  \u2502 ";
            for (size_t j = 0; j < cols; j++) {
                T val = data[i][j];
                if constexpr (std::is_floating_point_v<T>) {
                    if (std::abs(val) < 1e-9) val = static_cast<T>(0);
                }
                if (val == 0) std::cout << "\033[90m"; // 灰色表示0
                else std::cout << "\033[37m";
                std::cout << std::setw(10) << val << " ";
                std::cout << "\033[0m";
            }
            std::cout << " \u2502" << std::endl;
        }

        // 底部边框
        std::cout << "  \u2514";
        for (size_t j = 0; j < cols; j++) std::cout << "           ";
        std::cout << " \u2518" << std::endl;
    }

    // -------- Row Operations --------
    void exchangeRows(size_t r1, size_t r2) {
        if (r1 >= rows || r2 >= rows)
            throw std::out_of_range("Row index out of bounds");
        std::swap(data[r1], data[r2]);
    }

    void scaleRow(size_t r, T scalar) {
        if (r >= rows)
            throw std::out_of_range("Row index out of bounds");
        if (std::is_floating_point<T>::value) {
            if (std::abs(static_cast<double>(scalar)) < 1e-9)
                throw std::invalid_argument("Scaling factor too small");
        }
        for (size_t j = 0; j < cols; j++) {
            data[r][j] *= scalar;
        }
    }

    void addScaledRow(size_t targetRow, size_t sourceRow, T scalar) {
        if (targetRow >= rows || sourceRow >= rows)
            throw std::out_of_range("Row index out of bounds");

        if (std::is_floating_point<T>::value) {
            if (std::abs(static_cast<double>(scalar)) < 1e-9) return;
        }

        for (size_t j = 0; j < cols; j++) {
            data[targetRow][j] += data[sourceRow][j] * scalar;
        }
    }

    // -------- Matrix Operations --------
    Matrix<T> transpose() const {
        Matrix<T> result(cols, rows);
        for(size_t col = 0; col < cols; col++)
            for(size_t row = 0; row < rows; row++)
                result.at(col,row) = data[row][col];
        return result;
    }

    Matrix<T> operator+(const Matrix<T>& other) const {
        if(rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition");
        Matrix<T> result(rows, cols);
        for(size_t i = 0; i < rows; i++)
            for(size_t j = 0; j < cols; j++)
                result.at(i,j) = data[i][j] + other.data[i][j];
        return result;
    }

    Matrix<T> operator-(const Matrix<T>& other) const {
        if(rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for subtraction");
        Matrix<T> result(rows, cols);
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                result.at(i,j) = data[i][j] - other.data[i][j];
        return result;
    }

    Matrix<T> operator*(const Matrix<T>& other) const {
        if(cols != other.rows)
            throw std::invalid_argument("Matrix dimensions must match for multiplication");
        Matrix<T> result(rows, other.cols);
        for(size_t i = 0; i < rows; i++)
            for(size_t j = 0; j < other.cols; j++)
                for(size_t k = 0; k < cols; k++)
                    result.at(i,j) += data[i][k] * other.data[k][j];
        return result;
    }

    // 矩阵乘以向量 -> 向量
    Vector<T> operator*(const Vector<T>& v) const {
        if (cols != v.size())
            throw std::invalid_argument("Matrix columns must match vector size for multiplication");
        std::vector<T> result_vec(rows, T());
        for (size_t i = 0; i < rows; ++i) {
            for (size_t k = 0; k < cols; ++k) {
                result_vec[i] += data[i][k] * v[k];
            }
        }
        return Vector<T>(std::move(result_vec));
    }

    Matrix<T> operator*(T scalar) const {
        Matrix<T> result(rows, cols);
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                result.at(i,j) = data[i][j] * scalar;
        return result;
    }

    Matrix<T> operator/(T scalar) const {
        if(std::fabs(scalar) < 1e-9)
            throw std::invalid_argument("Scalar cannot be zero");
        Matrix<T> result(rows, cols);
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                result.at(i,j) = data[i][j] / scalar;
        return result;
    }

    Matrix<T> operator-() const {
        Matrix<T> result(rows, cols);
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                result.at(i,j) = -data[i][j];
        return result;
    }

    Matrix<T>& operator+=(const Matrix<T>& other) {
        if(rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition");
        for(size_t i = 0; i < rows; i++)
            for(size_t j = 0; j < cols; j++)
                data[i][j] += other.data[i][j];
        return *this;
    }

    Matrix<T>& operator-=(const Matrix<T>& other) {
        if(rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for subtraction");
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                data[i][j] -= other.data[i][j];
        return *this;
    }

    Matrix<T>& operator*=(const Matrix<T>& other) {
        if(cols != other.rows)
            throw std::invalid_argument("Matrix dimensions must match for multiplication");
        Matrix<T> result(rows, other.cols);
        for(size_t i = 0; i < rows; i++)
            for(size_t j = 0; j < other.cols; j++)
                for(size_t k = 0; k < cols; k++)
                    result.at(i,j) += data[i][k] * other.data[k][j];
        *this = result;
        return *this;
    }

    Matrix<T>& operator*=(T scalar) {
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                data[i][j] *= scalar;
        return *this;
    }

    Matrix<T>& operator/=(T scalar) {
        if(std::fabs(scalar) < 1e-9)
            throw std::invalid_argument("Scalar cannot be zero");
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                data[i][j] /= scalar;
        return *this;
    }

    // -------- Helpers --------
    static Matrix<T> rowSwap(int n, int i, int j) {
        Matrix<T> mat = identity(n);
        mat.exchangeRows(i, j);
        return mat;
    }

    static Matrix<T> rowScale(int n, int i, T c) {
        Matrix<T> mat = identity(n);
        mat.scaleRow(i, c);
        return mat;
    }

    static Matrix<T> rowadd(int n, int i, int j, T k) {
        Matrix<T> mat = identity(n);
        mat.addScaledRow(i, j, k);
        return mat;
    }

    Vector<T> getRow(size_t r) const {
        if (r >= rows) throw std::out_of_range("Row index out of bounds");
        return Vector<T>(data[r]);
    }

    Vector<T> getCol(size_t c) const {
        if (c >= cols) throw std::out_of_range("Col index out of bounds");
        std::vector<T> col(rows);
        for (size_t i = 0; i < rows; i++) col[i] = data[i][c];
        return Vector<T>(col);
    }

    Matrix<T> augment(const Matrix<T>& other) const {
        if (rows != other.rows) throw std::invalid_argument("Row count must match for augment");
        Matrix<T> result(rows, cols + other.cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) result.data[i][j] = data[i][j];
            for (size_t j = 0; j < other.cols; j++) result.data[i][cols + j] = other.data[i][j];
        }
        return result;
    }

    bool isSymmetric(T eps = static_cast<T>(1e-9)) const {
        if (rows != cols) return false;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = i + 1; j < cols; j++) {
                if (std::abs(data[i][j] - data[j][i]) > eps)
                    return false;
            }
        }
        return true;
    }

    bool isSkewSymmetric(T eps = static_cast<T>(1e-9)) const {
        if(rows != cols) return false;
        for(int i = 0; i < rows; i++) {
            for(int j = i + 1; j < cols; j++) {
                if(std::abs(data[i][i]) > eps) return false;
                if(std::abs(data[i][j] + data[j][i]) > eps) return false;
            }
        }
        return true;
    }

    void setToIdentity() {
        if (rows != cols) throw std::invalid_argument("Matrix must be square");
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < rows; j++)
                this->at(i, j) = (i == j) ? T(1) : T(0);
        }
    }

    Matrix<T> getInverseMatrix(T eps = static_cast<T>(1e-9)) const {
        if (this->rows != this->cols) throw std::invalid_argument("Matrix not square");
        T det = this->determinant(eps);
        if (std::abs(det) < eps) throw std::invalid_argument("Matrix is singular");
        int n = static_cast<int>(this->getCols());
        Matrix<T> augmentedMatrix = this->augment(Matrix<T>::identity(n));
        for (int i = 0; i < n; i++) {
            if (std::abs(augmentedMatrix.at(i, i)) < eps) throw std::invalid_argument("Numerical singularity");
            T scale = augmentedMatrix.at(i, i);
            for (int j = 0; j < 2 * n; j++) augmentedMatrix.at(i, j) /= scale;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    T factor = augmentedMatrix.at(j, i);
                    for (int k = 0; k < 2 * n; k++)
                        augmentedMatrix.at(j, k) -= factor * augmentedMatrix.at(i, k);
                }
            }
        }
        Matrix<T> inverseMatrix(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                inverseMatrix.at(i, j) = augmentedMatrix.at(i, j + n);
        return inverseMatrix;
    }

    bool isOrthogonal(T eps = static_cast<T>(1e-9)) const {
        if(this->getRows() != this->getCols()) throw std::invalid_argument("Must be square");
        Matrix<T> qt = this->transpose();
        Matrix<T> res = qt * (*this);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < rows; j++) {
                T val = res.at(i, j);
                if (i == j) { if (std::abs(val - 1) > eps) return false; }
                else { if (std::abs(val) > eps) return false; }
            }
        }
        return true;
    }

    // 延迟定义：实现位于 RREF.h (需要 RREF<T> 完整定义)
    // matrix.h 末尾已自动 include RREF.h，调用方无需额外引入
    int rank() const; 

    static Matrix<T> getRankNormalForm(size_t rows, size_t cols, size_t rank){
        Matrix<T> result(rows, cols);
        for(size_t i = 0; i < rank; i++) result.at(i, i)= T(1);
        return result;
    }

    Matrix<T> getEquivalenceNormalForm() const {
        size_t r = static_cast<size_t>(this->rank());
        return getRankNormalForm(rows, cols, r);
    } 

    static T dotProduct(const Matrix<T>& a, const Matrix<T>& b) {
        size_t aLen = (a.getRows() == 1) ? a.getCols() : a.getRows();
        size_t bLen = (b.getRows() == 1) ? b.getCols() : b.getRows();
        if (aLen != bLen) throw std::invalid_argument("Length mismatch");
        T sum = 0;
        for (size_t i = 0; i < aLen; i++) {
            T va = (a.getRows() == 1) ? a.at(0, i) : a.at(i, 0);
            T vb = (b.getRows() == 1) ? b.at(0, i) : b.at(i, 0);
            sum += va * vb;
        }
        return sum;
    }

    // 延迟定义：实现位于 RREF.h (需要 RREF<T> 完整定义)
    EigenDecomposition eigen(int max_iter = 1000) const;

    bool isSquare() const { return rows == cols; }

    bool isDiagonalizable() const;

    DiagonalizationResult diagonalize() const;

    T determinant(T eps = static_cast<T>(1e-9)) const {
        if (rows != cols) throw std::domain_error("Must be square");
        Matrix<T> temp(*this);
        T det = 1;
        int sign = 1;
        for (size_t i = 0; i < rows; i++) {
            size_t maxindex = i;
            for (size_t row = i + 1; row < rows; row++) {
                if (std::abs(temp.data[row][i]) > std::abs(temp.data[maxindex][i]))
                    maxindex = row;
            }
            if (std::abs(temp.data[maxindex][i]) < eps) return 0;
            if (maxindex != i) {
                temp.exchangeRows(maxindex, i);
                sign *= -1;
            }
            for (size_t row = i + 1; row < rows; row++) {
                if (std::abs(temp.data[row][i]) < eps) continue;
                T factor = -temp.data[row][i] / temp.data[i][i];
                temp.addScaledRow(row, i, factor);
            }
        }
        det = static_cast<T>(sign);
        for (size_t i = 0; i < rows; i++) det *= temp.data[i][i];
        return (std::abs(det) < eps) ? 0 : det;
    }

    Matrix<T> similarityTransform(const Matrix<T>& P) const {
        if (rows != cols) throw std::invalid_argument("Must be square");
        Matrix<T> Pinv = P.getInverseMatrix();
        return Pinv * (*this) * P;
    }

    bool isPossiblySimilarTo(const Matrix<T>& other) const {
        if (rows != cols || other.rows != other.cols) return false;
        if (rows != other.rows) return false;
        if (this->rank() != other.rank()) return false;
        if (this->determinant() != other.determinant()) return false;
        return true;
    }

    std::pair<Matrix<T>, Matrix<T>> qr_decomposition() const {
        if (rows != cols) throw std::invalid_argument("Must be square");
        int n = rows;
        std::vector<Vector<T>> a_cols;
        for(int j=0; j<n; j++) a_cols.push_back(this->getCol(j));
        std::vector<Vector<T>> q_cols;
        
        for (size_t i = 0; i < n; i++) {
            Vector<T> u = a_cols[i];
            for (size_t j = 0; j < i; j++) {
                T r_ji = q_cols[j].dot(a_cols[i]);
                u -= q_cols[j] * r_ji;
            }
            T r_ii = u.norm();
            if (std::abs(r_ii) < 1e-9) q_cols.push_back(u);
            else q_cols.push_back(u / r_ii);
        }

        Matrix<T> Q(n, n);
        for(size_t j=0; j<n; j++) 
            for(size_t i=0; i<n; i++) 
                Q.at(i, j) = q_cols[j][i];

        Matrix<T> R = Q.transpose() * (*this);
        for(size_t i=0; i<n; i++)
            for(size_t j=0; j<i; j++)
                R.at(i, j) = 0;
        return {Q, R};
    }

    // 矩阵 1-范数：列模和的最大值
    T norm1() const {
        if (cols == 0) return 0;
        T maxColSum = 0;
        for (size_t j = 0; j < cols; j++) {
            T colSum = 0;
            for (size_t i = 0; i < rows; i++) {
                colSum += std::abs(data[i][j]);
            }
            if (j == 0 || colSum > maxColSum) maxColSum = colSum;
        }
        return maxColSum;
    }

    // 矩阵 ∞-范数：行模和的最大值
    T normInf() const {
        if (rows == 0) return 0;
        T maxRowSum = 0;
        for (size_t i = 0; i < rows; i++) {
            T rowSum = 0;
            for (size_t j = 0; j < cols; j++) {
                rowSum += std::abs(data[i][j]);
            }
            if (i == 0 || rowSum > maxRowSum) maxRowSum = rowSum;
        }
        return maxRowSum;
    }

    // Frobenius 范数：所有元素平方和的平方根
    T normFrobenius() const {
        T sumSq = 0;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                sumSq += data[i][j] * data[i][j];
            }
        }
        return std::sqrt(sumSq);
    }
};

// =========================================================
// 自动引入 RREF.h，使 Matrix::rank() 和 Matrix::eigen() 的
// 延迟定义在所有 include matrix.h 的翻译单元中可用。
// 由 #pragma once 防止重复包含。
// =========================================================

