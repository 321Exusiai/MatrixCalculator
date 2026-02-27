直接复制以下全部内容到你的 README.md 中：Markdown# MatrixCalculator: A Layered Linear Algebra Toolkit

![C++](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)
![Build](https://img.shields.io/badge/Build-MinGW--w64-orange.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

本项目是一个基于 **C++17** 开发的高性能线性代数工具库。通过严谨的四层架构设计，实现了从底层向量运算到高层方程组求解、特征值计算及分块矩阵运算的完整功能。

---

## 🏗 项目架构 (System Architecture)

代码采用了分层设计（Layered Design），确保了极高的模块化程度和逻辑一致性：

* **Layer 0: `vector.h`** - 原子向量空间 $V^n$ 的基本操作。
* **Layer 1: `matrix.h`** - 基础矩阵层。支持内存管理、QR 分解及基础初等变换。
* **Layer 2: `RREF.h`** - 核心算法层。实现带主元选择的 Gauss-Jordan 消元逻辑。
* **Layer 3: 综合应用层**
    * `SolvingEquation.h`: 线性方程组全自动化求解器。
    * `VectorSet.h`: 向量组线性相关性分析与 Gram-Schmidt 正交化。
    * `BlockMatrix.h`: 高阶分块矩阵运算逻辑。

---

## 🧪 数学原理与代码实现

### 1. 向量归一化 (Normalization)
在 `vector.h` 中，利用 $L^2$ 范数实现向量标准化，为后续正交化打下基础：
$$
\|\mathbf{v}\| = \sqrt{\sum_{i=1}^n v_i^2} \implies \hat{\mathbf{v}} = \frac{\mathbf{v}}{\|\mathbf{v}\|}
$$

### 2. 矩阵消元与秩 (Rank)
`RREF.h` 实现了行最简形变换。通过寻找列主元（Partial Pivoting）降低数值误差：
$$
\text{rank}(A) = \#\{\text{pivot elements in RREF}(A)\}
$$

### 3. 方程组解的判定
在 `SolvingEquation.h` 中，程序通过对系数矩阵 $A$ 与增广矩阵 $(A|\mathbf{b})$ 的秩进行比较，自动识别解的状态：
- **唯一解**: $\text{rank}(A) = \text{rank}(A|\mathbf{b}) = n$
- **无穷解**: $\text{rank}(A) = \text{rank}(A|\mathbf{b}) < n$
- **无解**: $\text{rank}(A) < \text{rank}(A|\mathbf{b})$

---

## 📸 功能演示 (Demo)

### ① 线性方程组自动化求解
程序不仅能给出唯一解，还能处理包含自由变量的齐次与非齐次方程组。
**输入：**
$$
\begin{cases} x_1 + 2x_2 = 3 \\ 2x_1 + 4x_2 = 6 \end{cases}
$$
**程序输出演示：**
```text
[Status] Infinite Solutions Detected.
[Particular Solution] x_p = (3, 0)^T
[Basis of Null Space] η1 = (-2, 1)^T
② 矩阵特征值提取 (QR Algorithm)通过连续相似变换 $A_{k+1} = R_k Q_k$，矩阵将逐步收敛至上三角阵（Schur Form）：$$A_0 \xrightarrow{QR} A_1 \xrightarrow{QR} \dots \xrightarrow{QR} \begin{pmatrix} \lambda_1 & \dots & * \\ 0 & \lambda_2 & \dots \\ 0 & 0 & \lambda_n \end{pmatrix}$$🚀 核心优化细节移动语义 (Move Semantics): 在 matrix.h 构造函数与运算符重载中大量使用 std::move，将大型矩阵操作的复杂度从 $O(n^2)$ 的数据拷贝降至 $O(1)$ 的指针转移。数值稳定性: 消元算法在每次寻找主元时，会搜索当前列绝对值最大的元素：C++int max_row = find_max_pivot(current_col);
swap_rows(current_row, max_row);
鲁棒异常处理: 针对不可逆矩阵求逆、维度不匹配等错误，系统会抛出 std::invalid_argument 并提供详细的错误上下文。🛠 快速上手1. 编译 (Requires C++17)Bashg++ -std=c++17 main.cpp -o matrix_calc
2. 运行交互控制台Bash./matrix_calc
Developer: 321ExusiaiCourse: 线性代数
