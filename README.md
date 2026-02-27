# MatrixCalculator: A Layered Linear Algebra Toolkit

![C++](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)
![Build](https://img.shields.io/badge/Build-MinGW--w64-orange.svg)

本项目是一个基于 **C++17** 开发的高性能线性代数工具库。通过严谨的四层架构设计，实现了从底层向量运算到高层方程组求解、特征值计算及分块矩阵运算的完整功能。

---

## 🏗 项目架构 (System Architecture)

代码采用了分层设计（Layered Design），确保了极高的模块化程度和可维护性：

* **Layer 0: `vector.h`** - 原子向量操作。实现向量空间 $V^n$ 的基本定义。
* **Layer 1: `matrix.h`** - 基础矩阵层。支持内存管理、QR 分解及基础行列变换。
* **Layer 2: `RREF.h`** - 矩阵变换算法。核心实现带主元选择的高斯-约当消元逻辑。
* **Layer 3: 综合应用层**
    * `SolvingEquation.h`: 线性方程组全自动化求解。
    * `VectorSet.h`: 向量组线性相关性分析及正交化。
    * `BlockMatrix.h`: 分块矩阵的高阶运算逻辑。

---

## 🧪 数学原理与代码实现

### 1. 向量归一化 (Normalization)

在 `vector.h` 中，我们实现了基于 $L_2$ 范数的归一化：

$$\|\mathbf{v}\| = \sqrt{\sum_{i=1}^n v_i^2} \implies \hat{\mathbf{v}} = \frac{\mathbf{v}}{\|\mathbf{v}\|}$$

### 2. 矩阵消元与 RREF 变换

`RREF.h` 实现了带主元选择的消元算法。通过初等行变换将矩阵 $A$ 转化为行最简形：

$$\text{rank}(A) = \#\{\text{pivot elements in RREF}(A)\}$$


### 3. 线性方程组判定逻辑

在 `SolvingEquation.h` 中，程序通过比较系数矩阵 $A$ 和增广矩阵 $(A|\mathbf{b})$ 的秩来判定解的状态：

* **唯一解**: $\text{rank}(A) = \text{rank}(A|\mathbf{b}) = n$
* **无穷多解**: $\text{rank}(A) = \text{rank}(A|\mathbf{b}) < n$
* **无解**: $\text{rank}(A) < \text{rank}(A|\mathbf{b})$

### 4. 特征值迭代 (QR Algorithm)

`matrix.h` 中通过连续相似变换寻找特征值，矩阵将逐步收敛至上三角阵（Schur Form）：

$$A_k = Q_k R_k \implies A_{k+1} = R_k Q_k$$


---

## 📸 功能演示 (Demo)

### 线性方程组全自动化求解

程序能自动处理包含自由变量的非齐次方程组。

**输入示例：**

$$\begin{cases} x_1 + 2x_2 = 3 \\ 2x_1 + 4x_2 = 6 \end{cases}$$

**程序输出：**

```text
[Status] Infinite Solutions Detected.
[Particular Solution] x_p = (3, 0)^T
[Basis of Null Space] η1 = (-2, 1)^T
```
---
## 🚀 性能优化特性移动语义 (Move Semantics)：在 Matrix 类中广泛使用 std::move，将复杂度从 $O(n^2)$ 的拷贝降至 $O(1)$ 的指针转移。
数值稳定性：在消元时搜索当前列绝对值最大的元素作为主元，抑制浮点误差：C++int max_row = find_max_pivot(current_col);
swap_rows(current_row, max_row);
鲁棒性 (Robustness)：利用 std::invalid_argument 对维度不匹配或不可逆矩阵进行严格校验。
## 🛠 快速上手 
1. 编译 (Requires C++17)Bashg++ -std=c++17 main.cpp -o matrix_calc
2. 运行Bash./matrix_calc
# Developer: 321Exusiai 
# Course: 高等代数 (Freshman Year WHU 2026 寒假)
