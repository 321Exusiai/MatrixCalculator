#pragma once

#include "matrix.h"
#include <vector>
#include <stdexcept>

template <typename T>
class BlockMatrix
{
private:
    std::vector<std::vector<Matrix<T> > > blocks;
    int numRows;
    int numCols;
    int blockSize;
public:
    BlockMatrix(int numRows,int numCols,int blockSize):numRows(numRows), numCols(numCols), blockSize(blockSize) 
    {
        blocks.resize(numRows);
        for(int i = 0; i < numRows; i++)
        {
            blocks[i].resize(numCols);
            for(int j = 0; j < numCols; j++)
            {
                blocks[i][j] = Matrix<T>(blockSize, blockSize);
            }
        }
    }

    // 移动构造
    BlockMatrix(BlockMatrix&& other) noexcept 
        : blocks(std::move(other.blocks)), numRows(other.numRows), numCols(other.numCols), blockSize(other.blockSize) {
        other.numRows = 0; other.numCols = 0; other.blockSize = 0;
    }

    BlockMatrix& operator=(BlockMatrix&& other) noexcept {
        if (this != &other) {
            blocks = std::move(other.blocks);
            numRows = other.numRows;
            numCols = other.numCols;
            blockSize = other.blockSize;
            other.numRows = 0; other.numCols = 0; other.blockSize = 0;
        }
        return *this;
    }

    BlockMatrix(const BlockMatrix&) = default;
    BlockMatrix& operator=(const BlockMatrix&) = default;

    // 获取某个块
    Matrix<T>& getBlock(int row, int col)
    {
        if (row < 0 || row >= numRows || col < 0 || col >= numCols)
            throw std::out_of_range("Block index out of bounds");
        return blocks[row][col];
    }

    const Matrix<T>& getBlock(int row, int col) const
    {
        if (row < 0 || row >= numRows || col < 0 || col >= numCols)
            throw std::out_of_range("Block index out of bounds");
        return blocks[row][col];
    }

    // 获取分块矩阵的行数和列数
    int getTotalRows() const noexcept { return numRows * blockSize; }
    int getTotalCols() const noexcept { return numCols * blockSize; }
    int getBlockRows() const noexcept { return numRows; }
    int getBlockCols() const noexcept { return numCols; }
    
    static BlockMatrix<T> identity(int numRows, int numCols, int blockSize)
    {
        BlockMatrix<T> res(numRows, numCols, blockSize);
        if(numRows != numCols)
            throw std::invalid_argument("Identity block matrix must be square matrix");
        for(int i = 0; i < numRows; i++)
        {
            res.getBlock(i, i).setToIdentity();
        }
        return res;
    }

    // 矩阵操作
    BlockMatrix<T> transposeBlockMatrix() const {
        BlockMatrix<T> res(numCols, numRows, blockSize);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                res.getBlock(j, i) = this->blocks[i][j].transpose();
            }
        }
        return res;
    }

    BlockMatrix<T> operator+(const BlockMatrix<T>& other) const {
        if (numRows != other.numRows || numCols != other.numCols || blockSize != other.blockSize) {
            throw std::invalid_argument("Dimensions must match.");
        }
        BlockMatrix<T> res(numRows, numCols, blockSize);
        for(int i = 0; i < numRows; i++) {
            for(int j = 0; j < numCols; j++) { // 修正为 numCols
                // 使用 . 访问引用，并返回计算好的 res
                res.getBlock(i, j) = this->blocks[i][j] + other.blocks[i][j];
            }
        }
    return res;
    }

    BlockMatrix<T> operator-(const BlockMatrix<T>& other) const {
        if (numRows != other.numRows || numCols != other.numCols || blockSize != other.blockSize) {
            throw std::invalid_argument("BlockMatrix dimensions and block sizes must match for subtraction.");
        }
        BlockMatrix<T> res(numRows, numCols, blockSize);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                res.getBlock(i, j) = this->blocks[i][j] - other.blocks[i][j];
            }
        }
        return res;
    }

    // 标量乘法
    BlockMatrix<T> operator*(T scalar) const {
        BlockMatrix<T> res(numRows, numCols, blockSize);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                res.blocks[i][j] = this->blocks[i][j] * scalar;
            }
        }
        return res;
    }

    // 标量除法
    BlockMatrix<T> operator/(T scalar) const {
        BlockMatrix<T> res(numRows, numCols, blockSize);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                res.blocks[i][j] = this->blocks[i][j] / scalar;
            }
        }
        return res;
    }

    // 负号运算符
    BlockMatrix<T> operator-() const {
        BlockMatrix<T> res(numRows, numCols, blockSize);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                res.blocks[i][j] = -this->blocks[i][j];
            }
        }
        return res;
    }

    // -------- 复合赋值运算符 (Compound Assignment Operators) --------

    BlockMatrix<T>& operator+=(const BlockMatrix<T>& other) {
        if (numRows != other.numRows || numCols != other.numCols || blockSize != other.blockSize) 
            throw std::invalid_argument("Dimensions mismatch");
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                this->blocks[i][j] += other.blocks[i][j];
            }
        }
        return *this;
    }

    BlockMatrix<T>& operator-=(const BlockMatrix<T>& other) {
        if (numRows != other.numRows || numCols != other.numCols || blockSize != other.blockSize) 
            throw std::invalid_argument("Dimensions mismatch");
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                this->blocks[i][j] -= other.blocks[i][j];
            }
        }
        return *this;
    }

    BlockMatrix<T>& operator*=(T scalar) {
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                this->blocks[i][j] *= scalar;
            }
        }
        return *this;
    }

    BlockMatrix<T>& operator/=(T scalar) {
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                this->blocks[i][j] /= scalar;
            }
        }
        return *this;
    }

   BlockMatrix<T>& operator*=(const BlockMatrix<T>& other) {
    if (numCols != other.numRows || blockSize != other.blockSize)
        throw std::invalid_argument("Dimensions mismatch for multiplication.");

    // 1. 预分配空间（仅此一个临时对象）
    BlockMatrix<T> res(numRows, other.numCols, blockSize);

    // 2. 执行矩阵乘法逻辑
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < other.numCols; j++) {
            // 这里可以进一步优化：先清除 res 的当前块（如果是新构造的则已经是零矩阵）
            for (int k = 0; k < numCols; k++) {
                res.blocks[i][j] += this->blocks[i][k] * other.blocks[k][j];
            }
        }
    }

    // 3. 核心优化：利用 std::swap 交换内部数据
    // 这样不需要进行数据的逐个拷贝，只是交换了指针（vector 的内部指针）
    std::swap(this->blocks, res.blocks);
    this->numRows = res.numRows;
    this->numCols = res.numCols;
    // res 在函数结束处析构，由于它现在持有的是旧的（被替换掉的）数据，析构它是正常的
    
    return *this;
    }

    // 相应的，为了避免重复代码，让 operator* 调用 operator*=
BlockMatrix<T> operator*(const BlockMatrix<T>& other) const {
    BlockMatrix<T> res = *this; // 拷贝一份
    res *= other;               // 调用优化后的 *=
    return res;
}

    // 显示
    void display(T eps = static_cast<T>(1e-9)) const {
        std::cout << "BlockMatrix (" << getTotalRows() << "x" << getTotalCols() << "):" << std::endl;
        for (int i = 0; i < numRows; i++) {
            for (int rowInBlock = 0; rowInBlock < blockSize; rowInBlock++) {
                for (int j = 0; j < numCols; j++) {
                    for (int colInBlock = 0; colInBlock < blockSize; colInBlock++) {
                        T val = blocks[i][j].at(rowInBlock, colInBlock); 
                        if (std::abs(val) < eps) val = 0;
                        std::cout << std::setw(10) << val << " ";
                    }
                    std::cout << " | "; 
                }
                std::cout << std::endl;
            }
            for(int k=0; k < numCols * blockSize * 11; k++) std::cout << "-";
            std::cout << std::endl;
        }
    }

    // 初等行变换
    void exchangeBlockRows(int i, int j)
    {
        if(i >= numRows || j >= numRows)
            throw std::out_of_range("Block row index out of range");
        std::swap(blocks[i],blocks[j]);
    }

    void scaleBlockRow(int i, const Matrix<T>& multiplier) {
        if (i < 0 || i >= numRows)
            throw std::out_of_range("Block row index out of range");
        if (multiplier.getCols() != blockSize)
            throw std::invalid_argument("Scaling factor dimensions don't match");
        if constexpr (std::is_floating_point_v<T>) {
            if (std::abs(multiplier.determinant()) < static_cast<T>(1e-9))
                throw std::invalid_argument("Scaling factor too small");
        }
        for (int j = 0; j < numCols; j++) {
            blocks[i][j] = multiplier * blocks[i][j];
        }
    }

    void addScaledBlockRow(int targetRow, int sourceRow, const Matrix<T>& multiplier){
        if (targetRow < 0 || targetRow >= numRows || sourceRow < 0 || sourceRow >= numRows || multiplier.getCols() != blockSize) {
            throw std::out_of_range("Block row index out of range");
        }
        for(int j = 0; j < numCols; j++){
            blocks[targetRow][j] += multiplier * blocks[sourceRow][j]; 
        }
    }

    // 分块初等矩阵
    static BlockMatrix<T> blockSwapMatrix(int totalBlockRows, int blockSize, int i, int j) {
        // 1. 先生成一个分块单位阵
        BlockMatrix<T> E = BlockMatrix<T>::identity(totalBlockRows, totalBlockRows, blockSize);
        // 2. 交换第 i 行和第 j 行的块（此时 E 的 i,i 和 j,j 变成了 0，i,j 和 j,i 变成了 I）
        std::swap(E.getBlock(i, i), E.getBlock(j, i)); // 注意这里逻辑，本质是把单位阵的行换了
        std::swap(E.getBlock(i, j), E.getBlock(j, j));
        return E;
    }

    static BlockMatrix<T> blockScalingMatrix(int totalBlockRows, int blockSize, int i, const Matrix<T>& M) {
        BlockMatrix<T> E = BlockMatrix<T>::identity(totalBlockRows, totalBlockRows, blockSize);
        // 将第 i 个对角块替换为 M
        E.getBlock(i, i) = M;
        return E;
    }

    static BlockMatrix<T> blockAdditionMatrix(int totalBlockRows, int blockSize, int i, int j, const Matrix<T>& M) {
        // E = I
        BlockMatrix<T> E = BlockMatrix<T>::identity(totalBlockRows, totalBlockRows, blockSize);
        // 在 (i, j) 位置放入 M
        // 意味着当 A 左乘 E 时：Row_i = Row_i + M * Row_j
        E.getBlock(i, j) = M;
        return E;
    }

    // 把分块矩阵转换为普通矩阵
    // 在BlockMatrix类中添加：
    Matrix<T> toMatrix() const {
        int totalRows = getTotalRows();
        int totalCols = getTotalCols();
        Matrix<T> result(totalRows, totalCols);
    
        for(int blockRow = 0; blockRow < numRows; blockRow++) {
            for(int blockCol = 0; blockCol < numCols; blockCol++) {
                const Matrix<T>& block = blocks[blockRow][blockCol];
                for(int i = 0; i < blockSize; i++) {
                    for(int j = 0; j < blockSize; j++) {
                        int globalRow = blockRow * blockSize + i;
                        int globalCol = blockCol * blockSize + j;
                        result.at(globalRow, globalCol) = block.at(i, j);
                    }
                }
            }
        }
    
        return result;
    }
};
