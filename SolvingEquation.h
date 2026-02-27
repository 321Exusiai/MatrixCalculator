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
            int n = augmented.getCols() - 1;
            const auto& pivotCols = rrefSolver.getPivotCols();
            // 无解判断: 增广列（最后一列）是否有主元
            for(int col : pivotCols)
            {
                if(col == n)
                {
                    type = SolutionType::NoSolution;
                    return type;
                }
            }
            if(static_cast<int>(pivotCols.size()) == n)
                type = SolutionType::UniqueSolution;
            else 
                type = SolutionType::InfiniteSolutions;
            return type;
        }

        void computeSolution(T eps = static_cast<T>(1e-9)) 
        {
            if (type == SolutionType::NoSolution)
                solve();
            
            if (type == SolutionType::NoSolution)
                throw std::logic_error("No solution exists");

            int n = augmented.getCols() - 1;
            const auto& pivotCols = rrefSolver.getPivotCols();
            
            if (type == SolutionType::UniqueSolution) {
                std::vector<T> part_vec(n);
                for (int i = 0; i < n; i++)
                    part_vec[i] = rrefMatrix.at(i, n);
                particular = Vector<T>(std::move(part_vec));
            } else {
                std::vector<int> freeCols;
                for (int j = 0; j < n; j++) {
                    bool isPivot = false;
                    for(int pc : pivotCols) if(pc == j) { isPivot = true; break; }
                    if (!isPivot) freeCols.push_back(j);
                }
                
                std::vector<T> part_vec(n, 0);
                for (int i = 0; i < static_cast<int>(pivotCols.size()); i++) {
                    int col = pivotCols[i];
                    part_vec[col] = rrefMatrix.at(i, n);
                }
                particular = Vector<T>(std::move(part_vec));

                nullspace.clear();
                for (auto freeCol : freeCols) {
                    std::vector<T> v_vec(n, 0);
                    v_vec[freeCol] = 1;
                    for (int i = 0; i < static_cast<int>(pivotCols.size()); i++) {
                        int pcol = pivotCols[i];
                        v_vec[pcol] = -rrefMatrix.at(i, freeCol);
                    }
                    nullspace.push_back(Vector<T>(std::move(v_vec)));
                }
            }
        }

        void printSolution() const {
            int n = particular.size();
            if (type == SolutionType::NoSolution) {
                std::cout << "The system has NO solution" << std::endl;
                return;
            }
            if (type == SolutionType::UniqueSolution) 
            {
                std::cout << "Unique solution:" << std::endl;
                std::cout << "x = ( ";
                for (int i = 0; i < n; i++) {
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
            for (int i = 0; i < n; i++) {
                std::cout << particular[i];
                if (i != n - 1) std::cout << ", ";
            }
            std::cout << " )^T" << std::endl;

            // 齐次解空间
            for (int i = 0; i < nullspace.size(); i++) {
                std::cout << "  + t" << i + 1 << " * ( ";
                for (int j = 0; j < n; j++) {
                    std::cout << nullspace[i][j];
                    if (j != n - 1) std::cout << ", ";
                }
                std::cout << " )^T" << std::endl;
            }
        }
};