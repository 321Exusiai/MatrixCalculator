#pragma once
#include "RREF.h"  // 已传递引入 matrix.h 和 vector.h
#include<iostream>
#include<vector>
#include<cmath>
#include<stdexcept>

template <typename T>
class SolvingEquation
{
    public:
        enum class SolutionType {
        NoSolution,
        UniqueSolution,
        InfiniteSolutions
    };
    private:
        Matrix<T> augmented;
        Matrix<T> rrefMatrix;
        RREF<T> rrefSolver;
        SolutionType type;
        
        Vector<T> particular;//特解
        std::vector<Vector<T>> nullspace;//齐次解空间的基
    public:
        SolvingEquation(const Matrix<T>& A, const Matrix<T>& b)
    : rrefSolver(A.augment(b)), augmented(A.augment(b))
    {
        if (b.getRows() != A.getRows() || b.getCols() != 1)
            throw std::invalid_argument("Invalid dimensions");

        rrefSolver.toRREF();
        rrefMatrix = rrefSolver.getMatrix();
        type = solve();  // 立即求解
    }



        SolutionType solve()
        {
            size_t n = augmented.getCols() - 1;
            const auto& pivotCols = rrefSolver.getPivotCols();
            // 无解判断: 增广列（最后一列）是否有主元
            for(size_t col : pivotCols)
            {
                if(col == n)
                {
                    type = SolutionType::NoSolution;
                    return type;
                }
            }
            if(pivotCols.size() == n)
                type = SolutionType::UniqueSolution;
            else 
                type = SolutionType::InfiniteSolutions;
            return type;
        }

        void computeSolution(T eps = static_cast<T>(1e-12)) 
        {
            if (type == SolutionType::NoSolution)
                solve();
            
            if (type == SolutionType::NoSolution)
                throw std::logic_error("无解：系统不一致。");

            size_t rows = augmented.getRows();
            size_t cols = augmented.getCols() - 1; // 矩阵 A 的列数
            
            // 提取矩阵 A 和 向量 b
            Matrix<T> A(rows, cols);
            Vector<T> b(rows);
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    A.at(i, j) = augmented.at(i, j);
                }
                b[i] = augmented.at(i, cols);
            }

            // --- 高效引擎尝试 (针对方阵且唯一解情况) ---
            if (rows == cols && type == SolutionType::UniqueSolution) {
                try {
                    // 1. 尝试 Cholesky (仅对正定矩阵)
                    if (A.isSymmetric(eps)) {
                        try {
                            auto decomp = A.cholesky(eps);
                            // 解 LL^T x = b
                            Vector<T> y = Matrix<T>::forwardSubstitution(decomp.L, b);
                            particular = Matrix<T>::backwardSubstitution(decomp.L.transpose(), y);
                            std::cout << "[引擎消息]: 检测到正定性，已启用 Cholesky 分解加速求解。" << std::endl;
                            return; 
                        } catch (...) { /* 非正定，降级到 LU */ }
                    }

                    // 2. 尝试 LU 分解 (PA = LU)
                    auto decomp = A.lu(eps);
                    // 处理置换: b -> Pb
                    Vector<T> Pb(rows);
                    for (size_t i = 0; i < rows; ++i) Pb[i] = b[decomp.P[i]];
                    
                    Vector<T> y = Matrix<T>::forwardSubstitution(decomp.L, Pb);
                    particular = Matrix<T>::backwardSubstitution(decomp.U, y);
                    std::cout << "[引擎消息]: 已启用 LU 分解 (Partial Pivoting) 加速求解。" << std::endl;
                    return;
                } catch (...) {
                    // 如果 LU 也失败（数值极不稳定），则回退到 RREF 保证结果
                    std::cout << "[系统警告]: LU 分析遭遇困难，正在回退到 RREF 兜底流程..." << std::endl;
                }
            }

            // --- 兜底流程: 基于 RREF 提取解 (适用于奇异矩阵、非方阵、无穷多解) ---
            const auto& pivotCols = rrefSolver.getPivotCols();
            if (type == SolutionType::UniqueSolution) {
                std::vector<T> part_vec(cols);
                for (size_t i = 0; i < cols; i++)
                    part_vec[i] = rrefMatrix.at(i, cols);
                particular = Vector<T>(std::move(part_vec));
            } else {
                std::vector<size_t> freeCols;
                for (size_t j = 0; j < cols; j++) {
                    bool isPivot = false;
                    for(size_t pc : pivotCols) if(pc == j) { isPivot = true; break; }
                    if (!isPivot) freeCols.push_back(j);
                }
                
                std::vector<T> part_vec(cols, 0);
                for (size_t i = 0; i < pivotCols.size(); i++) {
                    size_t col = pivotCols[i];
                    part_vec[col] = rrefMatrix.at(i, cols);
                }
                particular = Vector<T>(std::move(part_vec));

                nullspace.clear();
                for (auto freeCol : freeCols) {
                    std::vector<T> v_vec(cols, 0);
                    v_vec[freeCol] = 1;
                    for (size_t i = 0; i < pivotCols.size(); i++) {
                        size_t pcol = pivotCols[i];
                        v_vec[pcol] = -rrefMatrix.at(i, freeCol);
                    }
                    nullspace.push_back(Vector<T>(std::move(v_vec)));
                }
            }
        }

        void printSolution() const {
            size_t n = particular.size();
            if (type == SolutionType::NoSolution) {
                std::cout << "The system has NO solution" << std::endl;
                return;
            }
            if (type == SolutionType::UniqueSolution) 
            {
                std::cout << "Unique solution:" << std::endl;
                std::cout << "x = ( ";
                for (size_t i = 0; i < n; i++) {
                    std::cout << particular[i];
                    if (i != n - 1) std::cout << ", ";
                }
                std::cout << " )^T" << std::endl;
                return;
            }
            std::cout << "Infinite solutions:" << std::endl;
            //无穷解
            // 特解
            std::cout << "x = ( ";
            for (size_t i = 0; i < n; i++) {
                std::cout << particular[i];
                if (i != n - 1) std::cout << ", ";
            }
            std::cout << " )^T" << std::endl;

            // 齐次解空间
            for (size_t i = 0; i < nullspace.size(); i++) {
                std::cout << "  + t" << i + 1 << " * ( ";
                for (size_t j = 0; j < n; j++) {
                    std::cout << nullspace[i][j];
                    if (j != n - 1) std::cout << ", ";
                }
                std::cout << " )^T" << std::endl;
            }
        }
};