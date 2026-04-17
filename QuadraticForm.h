#pragma once

#include "matrix.h"
#include "RREF.h"
#include "SolvingEquation.h"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

template <typename T>
class QuadraticForm {
private:
    size_t n;         // 未知数的个数
    Matrix<T> mat;    // 二次型对应的实对称矩阵

public:
    // 构造函数：接受维度 n 和长度为 n(n+1)/2 的系数一维数组
    // 数组顺序建议为：a11, a12, ..., a1n, a22, a23, ..., ann
    QuadraticForm(size_t dim, const std::vector<T>& coeffs) : n(dim), mat(dim, dim) {
        size_t expectedSize = n * (n + 1) / 2;
        if (coeffs.size() != expectedSize) {
            throw std::invalid_argument("二次型的系数个数必须恰好为 n(n+1)/2");
        }

        size_t index = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                if (i == j) {
                    mat.at(i, j) = coeffs[index];       // 对角元保留
                } else {
                    T halfVal = coeffs[index] / static_cast<T>(2);
                    mat.at(i, j) = halfVal;             // 非对角元除以 2 并对称分配
                    mat.at(j, i) = halfVal;
                }
                index++;
            }
        }
    }

    // 获取未知数个数
    size_t getN() const {
        return n;
    }

    // 获取二次型矩阵 (实对称矩阵)
    const Matrix<T>& getMatrix() const {
        return mat;
    }

    // 运用实对称矩阵对角化，化为标准型与规范型
    void orthogonalStandardize() const {
        std::cout << "\n--- [ 1. 二次型对应的实对称矩阵 A ] ---" << std::endl;
        mat.display();

        std::cout << "\n--- [ 2. 进行正交对角化 ] ---" << std::endl;
        auto res = mat.diagonalize();
        
        std::cout << "正交变换矩阵 P (特征向量列矩阵):" << std::endl;
        res.P.display();

        std::cout << "对角矩阵 D (特征值矩阵):" << std::endl;
        res.D.display();

        std::cout << "\n--- [ 3. 结果分析 ] ---" << std::endl;
        std::cout << "正交坐标变换为: X = PY" << std::endl;
        
        std::cout << "主轴标准型为: f = ";
        size_t p = 0; // 正惯性指数
        size_t q = 0; // 负惯性指数
        bool first = true;
        for (size_t i = 0; i < n; ++i) {
            T val = res.D.at(i, i);
            if (std::abs(val) < 1e-9) continue;

            if (val > 1e-9) p++;
            else if (val < -1e-9) q++;

            if (!first && val > 0) std::cout << " + ";
            if (val < 0) std::cout << " - ";
            
            std::cout << std::abs(val) << "y" << (i + 1) << "^2";
            first = false;
        }
        if (first) std::cout << "0";
        std::cout << std::endl;

        std::cout << "正惯性指数 p = " << p << ", 负惯性指数 q = " << q << std::endl;
        
        std::cout << "复规范型 (Canonical Form) 为: f = ";
        bool first_c = true;
        // 先打印 1 的项
        for (size_t i = 0; i < p; ++i) {
            if (!first_c) std::cout << " + ";
            std::cout << "z" << (i + 1) << "^2";
            first_c = false;
        }
        // 再打印 -1 的项
        for (size_t i = 0; i < q; ++i) {
            std::cout << " - z" << (p + i + 1) << "^2";
            first_c = false;
        }
        if (first_c) std::cout << "0";
        std::cout << std::endl;

        std::cout << "\n--- [ 4. 有定性判定 ] ---" << std::endl;
        if (p == n) {
            std::cout << "该二次型/矩阵是: " << GREEN << BOLD << "正定 (Positive Definite)" << RESET << std::endl;
        } else if (q == n) {
            std::cout << "该二次型/矩阵是: " << RED << BOLD << "负定 (Negative Definite)" << RESET << std::endl;
        } else if (p > 0 && q > 0) {
            std::cout << "该二次型/矩阵是: " << MAGENTA << BOLD << "不定 (Indefinite)" << RESET << std::endl;
        } else if (p > 0 && q == 0) {
            std::cout << "该二次型/矩阵是: " << YELLOW << BOLD << "半正定 (Positive Semi-definite)" << RESET << std::endl;
        } else if (q > 0 && p == 0) {
            std::cout << "该二次型/矩阵是: " << YELLOW << BOLD << "半负定 (Negative Semi-definite)" << RESET << std::endl;
        } else {
            std::cout << "该二次型/矩阵是: " << WHITE << BOLD << "零矩阵" << RESET << std::endl;
        }
    }

    // 寻找驻点（极值点）：解决 2Ax + b = 0
    void solveExtremumDirect(const Vector<T>& b_vec) const {
        std::cout << CYAN << BOLD << "\n--- [ 1. 直接解析法求解驻点 ] ---" << RESET << std::endl;
        std::cout << "目标函数: f(x) = x^T A x + b^T x" << std::endl;
        std::cout << "平稳条件: grad f = 2Ax + b = 0  =>  2Ax = -b" << std::endl;
        
        Matrix<T> twoA = mat * static_cast<T>(2);
        Vector<T> negB = b_vec * static_cast<T>(-1);
        
        try {
            // 将 Vector negB 转化为 n x 1 的 Matrix 传入方程组求解器
            SolvingEquation<T> solver(twoA, Matrix<T>(negB));
            solver.computeSolution();
            
            std::cout << GREEN << "求解结果:" << RESET << std::endl;
            solver.printSolution();
        } catch (const std::exception& e) {
            std::cerr << RED << "解析求解失败: " << e.what() << RESET << std::endl;
        }
    }

    // 牛顿法演练 (验证一步收敛)
    void solveExtremumNewton(const Vector<T>& b_vec, const Vector<T>& x0) const {
        std::cout << CYAN << BOLD << "\n--- [ 2. 牛顿法 (Newton's Method) 仿真 ] ---" << RESET << std::endl;
        std::cout << YELLOW << "初始猜测点 x0: " << RESET; x0.print();
        
        try {
            // 一阶导 (梯度): g = 2Ax + b
            // 二阶导 (Hessian): H = 2A
            // Newton step: x = x0 - H^-1 * g
            Matrix<T> twoA = mat * static_cast<T>(2);
            Vector<T> grad = (twoA * x0) + b_vec;
            
            std::cout << "在 x0 处的梯度 grad = "; grad.print();
            
            std::cout << "计算 Hessian (2A) 的逆矩阵并步进..." << std::endl;
            Matrix<T> invH = twoA.getInverseMatrix();
            Vector<T> step = invH * grad;
            Vector<T> x1 = x0 - step;
            
            std::cout << GREEN << BOLD << "牛顿迭代结果 (第 1 步): " << RESET;
            x1.print();
            std::cout << "(注：由于二次型是二阶函数，牛顿法理论上 1 步必达精确极值点)" << std::endl;
        } catch (const std::exception& e) {
            std::cerr << RED << "牛顿法执行失败: " << e.what() << RESET << std::endl;
        }
    }

    // 约束最值分析：在单位球面 ||x||=1 上的取值范围
    void reportRangeOnUnitSphere() const {
        std::cout << CYAN << BOLD << "\n--- [ 3. 约束最值分析 (||x||=1) ] ---" << RESET << std::endl;
        std::cout << "根据瑞利商 (Rayleigh Quotient) 定理，在单位球面上：" << std::endl;
        
        auto res = mat.diagonalize();
        T max_lambda = res.D.at(0, 0);
        T min_lambda = res.D.at(0, 0);
        
        for (size_t i = 1; i < n; ++i) {
            T val = res.D.at(i, i);
            if (val > max_lambda) max_lambda = val;
            if (val < min_lambda) min_lambda = val;
        }
        
        std::cout << ">>> " << YELLOW << "最大值 (Max): " << RESET << BOLD << max_lambda << RESET << std::endl;
        std::cout << ">>> " << YELLOW << "最小值 (Min): " << RESET << BOLD << min_lambda << RESET << std::endl;
        std::cout << "取值范围为: [" << min_lambda << ", " << max_lambda << "]" << std::endl;
    }
};
