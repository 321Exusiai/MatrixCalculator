#include <iostream>
#include <vector>
#include <limits>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

#include "vector.h"       // Layer 0: 核心向量
#include "matrix.h"       // Layer 1: 核心矩阵 (依赖 vector.h)
#include "RREF.h"         // Layer 2: 矩阵变换 (依赖 matrix.h)
#include "VectorSet.h"    // Layer 3: 向量组 (依赖 vector.h/matrix.h/RREF.h)
#include "SolvingEquation.h" // Layer 3: 方程组求解 (依赖 RREF.h)
#include "BlockMatrix.h"   // Layer 3: 分块矩阵 (依赖 matrix.h)

// --- ANSI 颜色代码 ---
#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define WHITE   "\033[37m"

void enableANSI() {
#ifdef _WIN32
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if (hOut == INVALID_HANDLE_VALUE) return;
    DWORD dwMode = 0;
    if (!GetConsoleMode(hOut, &dwMode)) return;
    dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    SetConsoleMode(hOut, dwMode);
#endif
}

// 辅助函数：输入矩阵
Matrix<double> inputMatrix(const std::string& name) {
    size_t r = 0, c = 0;
    while (true) {
        std::cout << YELLOW << "请输入矩阵 " << BOLD << name << RESET << YELLOW << " 的行数和列数: " << RESET;
        if (std::cin >> r >> c && r > 0 && c > 0) break;
        std::cout << RED << "输入无效，请输入两个正整数！" << RESET << std::endl;
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    Matrix<double> mat(r, c);
    std::cout << YELLOW << "请输入 " << BOLD << r * c << RESET << YELLOW << " 个元素 (按行输入):" << RESET << std::endl;
    for (size_t i = 0; i < r; i++) {
        for (size_t j = 0; j < c; j++) {
            while (!(std::cin >> mat.at(i, j))) {
                std::cout << RED << "输入无效，请重新输入第 (" << i+1 << "," << j+1 << ") 个数字: " << RESET;
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
    }
    return mat;
}

void demoLinearSystem() {
    std::cout << CYAN << BOLD << "\n--- [ 解线性方程组 Ax = b ] ---" << RESET << std::endl;
    Matrix<double> A = inputMatrix("A");
    Matrix<double> b = inputMatrix("b");

    try {
        SolvingEquation<double> solver(A, b);
        solver.computeSolution();
        std::cout << GREEN << BOLD << "\n[ 求解结果 ]" << RESET << std::endl;
        solver.printSolution();
    } catch (const std::exception& e) {
        std::cout << RED << "错误: " << e.what() << RESET << std::endl;
    }
}

void demoMatrixCalc() {
    std::cout << CYAN << BOLD << "\n--- [ 矩阵计算器 ] ---" << RESET << std::endl;
    Matrix<double> A = inputMatrix("A");

    std::cout << YELLOW << "A 的秩 (Rank): " << BOLD << A.rank() << RESET << std::endl;
    
    if (A.getRows() == A.getCols()) {
        std::cout << YELLOW << "A 的行列式 (Det): " << BOLD << A.determinant() << RESET << std::endl;
        try {
            Matrix<double> Inv = A.getInverseMatrix();
            std::cout << GREEN << BOLD << "A 的逆矩阵:" << RESET << std::endl;
            Inv.display();
        } catch (const std::exception& e) {
            std::cout << MAGENTA << "A 不可逆: " << e.what() << RESET << std::endl;
        }

        std::cout << GREEN << BOLD << "A 的特征值和特征向量:" << RESET << std::endl;
        try {
            auto eig = A.eigen();
            std::cout << YELLOW << "特征值: " << RESET;
            for (auto v : eig.eigenvalues) std::cout << CYAN << v << " " << RESET;
            std::cout << GREEN << "\n特征向量:" << RESET << std::endl;
            for (const auto& v : eig.eigenvectors) v.print();
        } catch (const std::exception& e) {
            std::cout << RED << "计算失败: " << e.what() << RESET << std::endl;
        }
    } else {
        std::cout << MAGENTA << "非方阵无法计算行列式/逆/特征值。" << RESET << std::endl;
    }
}

void demoVectorSet() {
    std::cout << CYAN << BOLD << "\n--- [ 向量组正交化 (Gram-Schmidt) ] ---" << RESET << std::endl;
    size_t n, dim;
    std::cout << YELLOW << "请输入向量个数: " << RESET;
    std::cin >> n;
    std::cout << YELLOW << "请输入向量维度: " << RESET;
    std::cin >> dim;

    std::vector<Vector<double>> vecs;
    for (size_t i = 0; i < n; i++) {
        std::cout << YELLOW << "请输入向量 " << BOLD << "v" << i + 1 << RESET << YELLOW << " (" << dim << "维): " << RESET;
        std::vector<double> raw(dim);
        for (size_t j = 0; j < dim; j++) {
            while (!(std::cin >> raw[j])) {
                std::cout << RED << "错误：请输入数字！请重新输入 v" << i + 1 << " 的第 " << j + 1 << " 个分量: " << RESET;
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
        vecs.emplace_back(std::move(raw));
    }

    try {
        VectorSet<double> vSet(vecs);
        std::cout << YELLOW << "向量组秩: " << BOLD << vSet.dimension() << RESET << std::endl;
        if (vSet.isLinearIndependent()) std::cout << GREEN << "向量组线性无关" << RESET << std::endl;
        else std::cout << MAGENTA << "向量组线性相关" << RESET << std::endl;

        std::cout << GREEN << BOLD << "\n正交化结果:" << RESET << std::endl;
        auto ortho = VectorSet<double>::gramSchmidt(vecs, true); // true for normalize
        for (const auto& v : ortho) v.print();

    } catch (const std::exception& e) {
        std::cout << RED << "错误: " << e.what() << RESET << std::endl;
    }
}

void demoBlockMatrix() {
    std::cout << CYAN << BOLD << "\n--- [ 分块矩阵演示 ] ---" << RESET << std::endl;
    std::cout << YELLOW << "正在构建一个 4x4 矩阵，由 4 个 2x2 的块组成..." << RESET << std::endl;
    
    BlockMatrix<double> BM(2, 2, 2);

    Matrix<double> I = Matrix<double>::identity(2);
    Matrix<double> Z = Matrix<double>::zero(2);
    Matrix<double> Ones(2, 2); 
    Ones.at(0,0)=1; Ones.at(0,1)=1; Ones.at(1,0)=1; Ones.at(1,1)=1;

    BM.getBlock(0, 0) = I * 2.0;    // 左上: 2I
    BM.getBlock(0, 1) = Ones;       // 右上: 全1
    BM.getBlock(1, 0) = Z;          // 左下: 0
    BM.getBlock(1, 1) = I;          // 右下: I

    std::cout << GREEN << "构建的分块矩阵 M:" << RESET << std::endl;
    BM.display();

    std::cout << GREEN << "分块矩阵转置 M^T:" << RESET << std::endl;
    BM.transposeBlockMatrix().display();

    std::cout << GREEN << "分块矩阵乘法 M * M:" << RESET << std::endl;
    (BM * BM).display();
}

void printBanner() {
    std::cout << CYAN << BOLD;
    std::cout << "==========================================================" << std::endl;
    std::cout << "       __  ___      __         _      __" << std::endl;
    std::cout << "      /  |/  /___ _/ /______  (_)_ __/ /____  ____" << std::endl;
    std::cout << "     / /|_/ / __ `/ __/ ___/ / / / / / / __ `/ ___/" << std::endl;
    std::cout << "    / /  / / /_/ / /_/ /    / / /_/ / / /_/ / /" << std::endl;
    std::cout << "   /_/  /_/\\__,_/\\__/_/    /_/\\__,_/_/\\__,_/_/" << std::endl;
    std::cout << "                                              " << std::endl;
    std::cout << "           线性代数工具库 (Matrix Toolkit)" << std::endl;
    std::cout << "==========================================================" << RESET << std::endl;
}

int main() {
    enableANSI();
    while (true) {
        printBanner();
        std::cout << WHITE << "  1. 解线性方程组 (Ax = b)" << RESET << std::endl;
        std::cout << WHITE << "  2. 矩阵运算 (秩, 行列式, 逆, 特征值)" << RESET << std::endl;
        std::cout << WHITE << "  3. 向量组操作 (秩, Schmidt正交化)" << RESET << std::endl;
        std::cout << WHITE << "  4. 分块矩阵演示 (Block Matrix)" << RESET << std::endl;
        std::cout << RED   << "  0. 退出" << RESET << std::endl;
        std::cout << CYAN << "请选择操作 [0-4]: " << RESET;

        int choice;
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "\n\n";
            continue;
        }

        switch (choice) {
            case 1: demoLinearSystem(); break;
            case 2: demoMatrixCalc(); break;
            case 3: demoVectorSet(); break;
            case 4: demoBlockMatrix(); break;
            case 0: 
                std::cout << GREEN << BOLD << "感谢使用，再见！" << RESET << std::endl; 
                return 0;
            default: 
                std::cout << RED << "无效选择，请重试。" << RESET << std::endl;
        }
        std::cout << YELLOW << "\n按回车键继续..." << RESET;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cin.get();
        std::cout << "\n\n";
    }
}