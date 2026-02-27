# MatrixCalculator: A Layered Linear Algebra Toolkit

![C++](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)
![Build](https://img.shields.io/badge/Build-MinGW--w64-orange.svg)

本项目是一个基于 **C++17** 开发的高性能线性代数工具库。通过严谨的四层架构设计，实现了从底层向量运算到高层方程组求解、特征值计算及分块矩阵运算的完整功能。

---

## 🏗 项目架构 (System Architecture)

代码采用了分层设计（Layered Design），确保了极高的模块化程度和可维护性：

* **Layer 0: `vector.h`** - 原子向量操作。实现向量空间 $V^n$ 的基本定义。
* **Layer 1: `matrix.h`** - 基础矩阵层。支持矩阵存储、QR 分解及基础行列变换。
* **Layer 2: `RREF.h`** - 矩阵变换算法。核心实现高斯-约当消元逻辑。
* **Layer 3: 应用层** * `SolvingEquation.h`: 线性方程组全自动化求解。
    * `VectorSet.h`: 向量组线性相关性分析及正交化。
    * `BlockMatrix.h`: 分块矩阵的高阶运算逻辑。

---

## 🧪 数学原理与代码实现

### 1. 向量归一化 (Normalization)
在 `vector.h` 中，我们实现了基于 $L^2$ 范数的归一化，通过点积 (Dot Product) 计算向量长度：
$$
\|\mathbf{v}\| = \sqrt{\sum_{i=1}^n v_i^2} \implies \hat{\mathbf{v}} = \frac{\mathbf{v}}{\|\mathbf{v}\|}
$$

### 2. 矩阵消元与 RREF 变换
`RREF.h` 实现了带主元选择的消元算法。通过初等行变换将矩阵 $A$ 转化为行最简形，并自动识别矩阵的秩 (Rank)：
$$
\text{rank}(A) = \text{number of non-zero rows in RREF}(A)
$$

### 3. 线性方程组判定逻辑
在 `SolvingEquation.h` 中，程序通过比较系数矩阵 $A$ 和增广矩阵 $(A|\mathbf{b})$ 的秩来判定解的状态：
* **唯一解**: $\text{rank}(A) = \text{rank}(A|\mathbf{b}) = n$
* **无穷多解**: $\text{rank}(A) = \text{rank}(A|\mathbf{b}) < n$
* **无解**: $\text{rank}(A) < \text{rank}(A|\mathbf{b})$

### 4. 特征值迭代 (QR Algorithm)
`matrix.h` 中实现了 QR 分解，结合 `RREF.h` 中的迭代逻辑，通过相似变换寻找特征值：
$$
A_k = Q_k R_k \xrightarrow{\text{Iterate}} A_{k+1} = R_k Q_k
$$

---

## 🚀 性能优化特性

* **移动语义 (Move Semantics)**：在 `BlockMatrix` 与 `Matrix` 类中广泛使用 `std::move`，避免了大型矩阵拷贝带来的性能损耗。
* **数值稳定性**：在消元逻辑中加入了 `max_index` 寻找最大主元的步骤，有效抑制了浮点数运算中的舍入误差。
* **鲁棒性 (Robustness)**：利用 C++ 异常处理机制（`std::invalid_argument`），对维度不匹配、除零等非法数学操作进行严格校验。

---

## 🛠 快速上手

### 编译
本项目无需外部依赖，仅需支持 C++17 的编译器即可：
```bash
g++ -std=c++17 main.cpp -o matrix_calc****
