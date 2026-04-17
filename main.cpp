// --- ANSI 256色 & 霓虹色彩代码 ---
#define RESET       "\033[0m"
#define BOLD        "\033[1m"
#define RED         "\033[31m"
#define GREEN       "\033[32m"
#define YELLOW      "\033[33m"
#define BLUE        "\033[34m"
#define MAGENTA     "\033[35m"
#define CYAN        "\033[36m"
#define WHITE       "\033[37m"

// 赛博朋克扩展色 (256色)
#define NIXIE_ORANGE "\033[38;5;208m"
#define NEON_PINK    "\033[38;5;198m"
#define NEON_PURPLE  "\033[38;5;93m"
#define NEON_CYAN    "\033[38;5;51m"
#define GRID_GRAY    "\033[38;5;240m"

#include <iostream>
#include <vector>
#include <limits>
#include <string>
#include <thread>
#include <chrono>
#include <iomanip>

#ifdef _WIN32
#include <windows.h>
#endif

#include "vector.h"       // Layer 0: 核心向量
#include "matrix.h"       // Layer 1: 核心矩阵 (依赖 vector.h)
#include "RREF.h"         // Layer 2: 矩阵变换 (依赖 matrix.h)
#include "VectorSet.h"    // Layer 3: 向量组 (依赖 vector.h/matrix.h/RREF.h)
#include "SolvingEquation.h" // Layer 3: 方程组求解 (依赖 RREF.h)
#include "BlockMatrix.h"   // Layer 3: 分块矩阵 (依赖 matrix.h)
#include "QuadraticForm.h" // Layer 3: 二次型 (依赖 matrix.h/RREF.h)

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
    std::cout << CYAN << "A 的 1-范数: " << BOLD << A.norm1() << RESET << std::endl;
    std::cout << CYAN << "A 的 ∞-范数: " << BOLD << A.normInf() << RESET << std::endl;
    std::cout << CYAN << "A 的 F-范数: " << BOLD << A.normFrobenius() << RESET << std::endl;
    
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

void demoDiagonalization() {
    std::cout << CYAN << BOLD << "\n--- [ 矩阵对角化 ] ---" << RESET << std::endl;
    Matrix<double> A = inputMatrix("A");

    try {
        if (!A.isSquare()) {
            std::cout << RED << "错误：只有方阵可以对角化。" << RESET << std::endl;
            return;
        }

        if (A.isSymmetric()) {
            std::cout << GREEN << "检测到实对称矩阵，将进行正交对角化。" << RESET << std::endl;
        }

        if (A.isDiagonalizable()) {
            auto res = A.diagonalize();
            std::cout << GREEN << BOLD << "\n[ 对角化成功 ]" << RESET << std::endl;
            std::cout << YELLOW << "变换矩阵 P (特征向量构成):" << RESET << std::endl;
            res.P.display();
            std::cout << YELLOW << "对角矩阵 D (特征值):" << RESET << std::endl;
            res.D.display();

            std::cout << CYAN << "\n验证 P * D * P_inv:" << RESET << std::endl;
            Matrix<double> reconstructed;
            if (A.isSymmetric()) {
                reconstructed = res.P * res.D * res.P.transpose();
            } else {
                try {
                    reconstructed = res.P * res.D * res.P.getInverseMatrix();
                } catch (...) {
                    std::cout << RED << "无法计算 P 的逆矩阵，验证跳过。" << RESET << std::endl;
                    return;
                }
            }
            reconstructed.display();
        } else {
            std::cout << MAGENTA << "矩阵不可对角化（几何重数不足）。" << RESET << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << RED << "错误: " << e.what() << RESET << std::endl;
    }
}

void demoQuadraticForm() {
    std::cout << CYAN << BOLD << "\n--- [ 二次型化标准型 ] ---" << RESET << std::endl;
    
    int n;
    std::cout << YELLOW << "请输入未知数个数 n: " << RESET;
    while (!(std::cin >> n) || n <= 0) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << RED << "无效输入，请输入正整数: " << RESET;
    }

    size_t count = static_cast<size_t>(n) * (n + 1) / 2;
    std::vector<double> coeffs;
    std::cout << YELLOW << "请输入 " << count << " 个系数 (依次为 ";
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j <= n; ++j) {
            std::cout << "a" << i << j << " ";
        }
    }
    std::cout << "): " << RESET << std::endl;

    for (size_t k = 0; k < count; ++k) {
        double val;
        while (!(std::cin >> val)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << RED << "无效输入，请输入数字: " << RESET;
        }
        coeffs.push_back(val);
    }

    try {
        QuadraticForm<double> qf(static_cast<size_t>(n), coeffs);
        qf.orthogonalStandardize();

        // 交互：是否进入最优化实验室
        std::cout << NEON_CYAN << "\n>> 是否进入最优化实验室 (Optimization Lab) 进行极值分析? [y/n]: " << RESET;
        char optChoice;
        std::cin >> optChoice;
        if (optChoice == 'y' || optChoice == 'Y') {
            std::cout << YELLOW << "请输入线性项系数向量 b (长度为 " << n << "): " << RESET << std::endl;
            std::vector<double> b_data(n);
            for (int i = 0; i < n; ++i) {
                std::cout << "b" << (i + 1) << ": ";
                while (!(std::cin >> b_data[i])) {
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << RED << "无效输入，请重新输入数字: " << RESET;
                }
            }
            Vector<double> b_vec(std::move(b_data));

            // 1. 直接解析法
            qf.solveExtremumDirect(b_vec);

            // 2. 牛顿法演示（需要用户提供一个起点）
            std::cout << NEON_CYAN << "\n>> 是否进行牛顿法 (Newton's Method) 仿真演示? [y/n]: " << RESET;
            std::cin >> optChoice;
            if (optChoice == 'y' || optChoice == 'Y') {
                std::cout << YELLOW << "请输入模拟起点 x0 (长度为 " << n << "): " << RESET << std::endl;
                std::vector<double> x0_data(n);
                for (int i = 0; i < n; ++i) {
                    std::cout << "x0_" << (i + 1) << ": ";
                    std::cin >> x0_data[i];
                }
                qf.solveExtremumNewton(b_vec, Vector<double>(std::move(x0_data)));
            }

            // 3. 约束最值
            qf.reportRangeOnUnitSphere();
        }
    } catch (const std::exception& e) {
        std::cout << RED << "执行失败: " << e.what() << RESET << std::endl;
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

void printSystemBoot() {
    auto sleep_ms = [](int ms) {
        std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    };

    std::cout << NEON_CYAN << "[ 系统自检中... ]" << RESET << std::endl;
    const char* tasks[] = {
        "正在初始化核心矩阵引擎...",
        "正在建立 QR 迭代特征值求解器连接...",
        "正在加载二次型惯性指数分析模块...",
        "正在同步标准正交基计算单元...",
        "系统状态：[ 就绪 ]"
    };

    for (auto task : tasks) {
        std::cout << GRID_GRAY << "> " << task << RESET ;
        std::cout.flush();
        sleep_ms(200);
        std::cout << GREEN << " [OK]" << RESET << std::endl;
    }
    sleep_ms(300);
}

void printCyberBanner() {
    std::cout << NEON_PURPLE << "╔══════════════════════════════════════════════════════════╗" << RESET << std::endl;
    std::cout << NEON_PURPLE << "║" << NIXIE_ORANGE << "  __  __       _        _              ____ ___  ____  _____ " << NEON_PURPLE << "  ║" << RESET << std::endl;
    std::cout << NEON_PURPLE << "║" << NIXIE_ORANGE << " |  \\/  | __ _| |_ _ __(_)_  __  ____ / ___/ _ \\|  _ \\| ____|" << NEON_PURPLE << " ║" << RESET << std::endl;
    std::cout << NEON_PURPLE << "║" << NIXIE_ORANGE << " | |\\/| |/ _` | __| '__| \\ \\/ / / ___| |  | | | | |_) |  _|  " << NEON_PURPLE << " ║" << RESET << std::endl;
    std::cout << NEON_PURPLE << "║" << NIXIE_ORANGE << " | |  | | (_| | |_| |  | |>  <  |___ | |__| |_| |  _ <| |___ " << NEON_PURPLE << " ║" << RESET << std::endl;
    std::cout << NEON_PURPLE << "║" << NIXIE_ORANGE << " |_|  |_|\\__,_|\\__|_|  |_/_/\\_\\      \\____\\___/|_| \\_\\_____|" << NEON_PURPLE << " ║" << RESET << std::endl;
    std::cout << NEON_PURPLE << "║" << NEON_PINK << "                                            v1.2-CYBER" << NEON_PURPLE << "   ║" << RESET << std::endl;
    std::cout << NEON_PURPLE << "╚══════════════════════════════════════════════════════════╝" << RESET << std::endl;
    std::cout << NEON_CYAN << "          >> 线性代数工具库 | 霓虹终端演示版 <<" << RESET << std::endl;
}

void printDecoratedMenu() {
    std::cout << GRID_GRAY << "┌──────────────────────────────────────────────────────────┐" << RESET << std::endl;
    std::cout << GRID_GRAY << "│" << NEON_PINK << "  请选择操作指令:" << std::setw(42) << " " << GRID_GRAY << "│" << RESET << std::endl;
    std::cout << GRID_GRAY << "├──────────────────────────────────────────────────────────┤" << RESET << std::endl;
    
    auto printItem = [](const char* num, const char* text) {
        std::cout << GRID_GRAY << "│" << NIXIE_ORANGE << "  [" << num << "]" << RESET << " " << WHITE << std::left << std::setw(48) << text << GRID_GRAY << "│" << RESET << std::endl;
    };

    printItem("1", "解线性方程组 (Ax = b)");
    printItem("2", "矩阵运算 (秩, 行列式, 逆, 特征值)");
    printItem("3", "向量组操作 (秩, Schmidt 正交化)");
    printItem("4", "分块矩阵演示 (Block Matrix)");
    printItem("5", "矩阵对角化 (Diagonalization)");
    printItem("6", "二次型化标准型 (Quadratic Form)");
    
    std::cout << GRID_GRAY << "├──────────────────────────────────────────────────────────┤" << RESET << std::endl;
    std::cout << GRID_GRAY << "│" << RED << "  [0] 退出系统" << std::setw(45) << " " << GRID_GRAY << "│" << RESET << std::endl;
    std::cout << GRID_GRAY << "└──────────────────────────────────────────────────────────┘" << RESET << std::endl;
    std::cout << NEON_CYAN << "  >> " << RESET;
}

int main() {
    enableANSI();
    printSystemBoot();
    
    while (true) {
        printCyberBanner();
        printDecoratedMenu();

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
            case 5: demoDiagonalization(); break;
            case 6: demoQuadraticForm(); break;
            case 0: 
                std::cout << NEON_PINK << BOLD << "感谢使用，系统正在下线..." << RESET << std::endl; 
                return 0;
            default: 
                std::cout << RED << "无效指令，请重试。" << RESET << std::endl;
        }
        std::cout << YELLOW << "\n按回车键继续..." << RESET;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cin.get();
        std::cout << "\n\n";
    }
}