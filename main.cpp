#include <iostream>
#include <vector>
#include <limits>
#include "vector.h"       // Layer 0: 核心向量
#include "matrix.h"       // Layer 1: 核心矩阵 (依赖 vector.h)
#include "RREF.h"         // Layer 2: 矩阵变换 (依赖 matrix.h)
#include "VectorSet.h"    // Layer 3: 向量组 (依赖 vector.h/matrix.h/RREF.h)
#include "SolvingEquation.h" // Layer 3: 方程组求解 (依赖 RREF.h)
#include "BlockMatrix.h"   // Layer 3: 分块矩阵 (依赖 matrix.h)

// 辅助函数：输入矩阵
Matrix<double> inputMatrix(const std::string& name) {
    int r = 0, c = 0;
    while (true) {
        std::cout << "请输入矩阵 " << name << " 的行数和列数: ";
        if (std::cin >> r >> c && r > 0 && c > 0) break;
        std::cout << "输入无效，请输入两个正整数！" << std::endl;
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    Matrix<double> mat(r, c);
    std::cout << "请输入 " << r * c << " 个元素 (按行输入):" << std::endl;
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            while (!(std::cin >> mat.at(i, j))) {
                std::cout << "输入无效，请重新输入第 (" << i+1 << "," << j+1 << ") 个数字: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
    }
    return mat;
}

void demoLinearSystem() {
    std::cout << "\n--- 解线性方程组 Ax = b ---" << std::endl;
    Matrix<double> A = inputMatrix("A");
    Matrix<double> b = inputMatrix("b");

    try {
        SolvingEquation<double> solver(A, b);
        solver.computeSolution();
        solver.printSolution();
    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
}

void demoMatrixCalc() {
    std::cout << "\n--- 矩阵计算器 ---" << std::endl;
    Matrix<double> A = inputMatrix("A");

    std::cout << "A 的秩 (Rank): " << A.rank() << std::endl;
    
    if (A.getRows() == A.getCols()) {
        std::cout << "A 的行列式 (Det): " << A.determinant() << std::endl;
        try {
            Matrix<double> Inv = A.getInverseMatrix();
            std::cout << "A 的逆矩阵:" << std::endl;
            Inv.display();
        } catch (const std::exception& e) {
            std::cout << "A 不可逆: " << e.what() << std::endl;
        }

        std::cout << "A 的特征值和特征向量:" << std::endl;
        try {
            auto eig = A.eigen();
            std::cout << "特征值: ";
            for (auto v : eig.eigenvalues) std::cout << v << " ";
            std::cout << "\n特征向量:" << std::endl;
            for (const auto& v : eig.eigenvectors) v.print();
        } catch (const std::exception& e) {
            std::cout << "计算失败: " << e.what() << std::endl;
        }
    } else {
        std::cout << "非方阵无法计算行列式/逆/特征值。" << std::endl;
    }
}

void demoVectorSet() {
    std::cout << "\n--- 向量组正交化 (Gram-Schmidt) ---" << std::endl;
    int n, dim;
    std::cout << "请输入向量个数: ";
    std::cin >> n;
    std::cout << "请输入向量维度: ";
    std::cin >> dim;

    std::vector<Vector<double>> vecs;
    for (int i = 0; i < n; i++) {
        std::cout << "请输入向量 v" << i + 1 << " (" << dim << "维): ";
        std::vector<double> raw(dim);
        for (int j = 0; j < dim; j++) {
            while (!(std::cin >> raw[j])) {
                std::cout << "错误：请输入数字！请重新输入 v" << i + 1 << " 的第 " << j + 1 << " 个分量: ";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
        vecs.emplace_back(std::move(raw));
    }

    try {
        VectorSet<double> vSet(vecs);
        std::cout << "向量组秩: " << vSet.dimension() << std::endl;
        if (vSet.isLinearIndependent()) std::cout << "向量组线性无关" << std::endl;
        else std::cout << "向量组线性相关" << std::endl;

        std::cout << "\n正交化结果:" << std::endl;
        auto ortho = VectorSet<double>::gramSchmidt(vecs, true); // true for normalize
        for (const auto& v : ortho) v.print();

    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
}

void demoBlockMatrix() {
    std::cout << "\n--- 分块矩阵演示 ---" << std::endl;
    std::cout << "正在构建一个 4x4 矩阵，由 4 个 2x2 的块组成..." << std::endl;
    
    // 构造 2x2 块矩阵，每个块大小为 2
    BlockMatrix<double> BM(2, 2, 2);

    // 设置块
    Matrix<double> I = Matrix<double>::identity(2);
    Matrix<double> Z = Matrix<double>::zero(2);
    Matrix<double> Ones(2, 2); 
    Ones.at(0,0)=1; Ones.at(0,1)=1; Ones.at(1,0)=1; Ones.at(1,1)=1;

    BM.getBlock(0, 0) = I * 2.0;    // 左上: 2I
    BM.getBlock(0, 1) = Ones;       // 右上: 全1
    BM.getBlock(1, 0) = Z;          // 左下: 0
    BM.getBlock(1, 1) = I;          // 右下: I

    std::cout << "构建的分块矩阵 M:" << std::endl;
    BM.display();

    std::cout << "分块矩阵转置 M^T:" << std::endl;
    BM.transposeBlockMatrix().display();

    std::cout << "分块矩阵乘法 M * M:" << std::endl;
    (BM * BM).display();
}

int main() {
    while (true) {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "   线性代数工具库 (Interactive Console)" << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "1. 解线性方程组 (Ax = b)" << std::endl;
        std::cout << "2. 矩阵运算 (秩, 行列式, 逆, 特征值)" << std::endl;
        std::cout << "3. 向量组操作 (秩, Schmidt正交化)" << std::endl;
        std::cout << "4. 分块矩阵演示 (Block Matrix)" << std::endl;
        std::cout << "0. 退出" << std::endl;
        std::cout << "请选择操作 [0-4]: ";

        int choice;
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        switch (choice) {
            case 1: demoLinearSystem(); break;
            case 2: demoMatrixCalc(); break;
            case 3: demoVectorSet(); break;
            case 4: demoBlockMatrix(); break;
            case 0: std::cout << "再见！" << std::endl; return 0;
            default: std::cout << "无效选择，请重试。" << std::endl;
        }
    }
}