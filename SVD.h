#pragma once

#include "matrix.h"
#include "RREF.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>

/**
 * SVD 分解实验室 (Singular Value Decomposition)
 * 目标：将矩阵 A (m x n) 分解为 U * Sigma * V^T
 */
template <typename T>
class SVDSolver {
public:
    struct SVDResult {
        Matrix<T> U;
        Matrix<T> Sigma;
        Matrix<T> V; // 注意：此处 V 不是 V^T，展示时通常需要转置
    };

private:
    Matrix<T> A;
    size_t m, n;
    T eps;
    // 存储计算结果，由构造函数填充
    Matrix<T> U, Sigma, V;

public:
    SVDSolver(const Matrix<T>& input, T epsilon = static_cast<T>(1e-9)) 
        : A(input), m(input.getRows()), n(input.getCols()), eps(epsilon),
          U(input.getRows(), input.getRows()), 
          Sigma(input.getRows(), input.getCols()), 
          V(input.getCols(), input.getCols()) 
    {
        Matrix<T> At = A.transpose();
        
        if (m >= n) {
            // --- 情形 A: m >= n ---
            Matrix<T> AtA = At * A; 
            auto eig = AtA.eigen();

            // 1. 使用方法 A 对特征值/向量进行降序排列
            struct SortEntry { T val; Vector<T> vec; };
            std::vector<SortEntry> entries;
            for (size_t i = 0; i < eig.eigenvalues.size(); i++) {
                entries.push_back({eig.eigenvalues[i], eig.eigenvectors[i]});
            }
            std::sort(entries.begin(), entries.end(), [](const SortEntry& a, const SortEntry& b) {
                return a.val > b.val;
            });

            // 2. 构造 V 和 Sigma (对角线为 sqrt(lambda))
            std::vector<Vector<T>> sorted_v;
            Sigma = Matrix<T>(m, n); // 确保全零
            for (size_t i = 0; i < entries.size(); i++) {
                sorted_v.push_back(entries[i].vec);
                if (entries[i].val > eps) {
                    Sigma.at(i, i) = std::sqrt(std::max(T(0), entries[i].val));
                }
            }
            V = Matrix<T>::fromColumns(sorted_v);

            // 3. 计算 U 及其基补全
            std::vector<Vector<T>> u_cols;
            for (size_t i = 0; i < n; i++) {
                T sigma = Sigma.at(i, i);
                if (sigma > eps) {
                    // 公式：u_i = (1/sigma_i) * A * v_i
                    u_cols.push_back((A * V.getCol(i)) * (T(1) / sigma));
                }
            }
            
            // 补全基向量：添加标准基辅助 Gram-Schmidt 寻找左零空间
            for (size_t j = 0; j < m; j++) {
                Vector<T> ei(m, T(0));
                ei[j] = T(1);
                u_cols.push_back(ei);
            }
            auto full_u_basis = VectorSet<T>::gramSchmidt(u_cols, true);
            U = Matrix<T>::fromColumns(full_u_basis);
        }
        else {
            // --- 情形 B: m < n (留给你大显身手！) ---
            Matrix<T> AAt = A * At; 
            
            // 【提示】你的开发路线：
            // 1. 对 AAt 求特征分解 eig = AAt.eigen()
            // 2. 同样的排序逻辑 (方法 A)
            // 3. 构造 U (此时 U 由 AAt 的特征向量构成)
            // 4. 利用公式 v_i = (1/sigma_i) * At * u_i 计算 V 的前 m 列
            // 5. 使用 Gram-Schmidt 补齐 V 的剩余 n-m 列
        }
    }

    /**
     * compute() 现在仅需返回构造中预计算好的结果
     */
    SVDResult compute() {
        return {U, Sigma, V};
    }

    /**
     * 赛博风格的结果展示
     */
    void display(const SVDResult& res) {
        std::cout << NEON_CYAN << "\n--- [ SVD 分解实验室系统报告 ] ---" << RESET << std::endl;
        
        std::cout << YELLOW << "\n[ 左奇异向量矩阵 U (" << m << "x" << m << ") ]" << RESET << std::endl;
        res.U.display();

        std::cout << YELLOW << "\n[ 奇异值矩阵 Sigma (" << m << "x" << n << ") ]" << RESET << std::endl;
        res.Sigma.display();

        std::cout << YELLOW << "\n[ 右奇异向量矩阵 V (" << n << "x" << n << ") ]" << RESET << std::endl;
        res.V.display();

        std::cout << GREEN << BOLD << "\n验证: A \u2248 U * Sigma * V^T" << RESET << std::endl;
    }
};
