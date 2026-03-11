// =========================================================
// VectorSet.h — 向量组级操作 (Layer 3, 应用层)
// ---------------------------------------------------------
// 职责: 线性无关判定、基提取、秩、Gram-Schmidt 正交化
// 单个向量的原子操作 (加减乘除、点积、范数) 见 vector.h
// =========================================================
#pragma once

#include "vector.h"
#include "matrix.h"  // 已自动引入 RREF.h (见 matrix.h 末尾)
#include<iostream>
#include<vector>
#include<cmath>
#include<stdexcept>
#include<algorithm>

template <typename T>
class VectorSet {
public:
    enum class VectorOrientation {
        Column,  // 默认：教材风格
        Row
    };

private:
    std::vector<Vector<T>> vecs;
    size_t dim;
    size_t rank;
    std::vector<size_t> pivotIndices;
    VectorOrientation orientation;

public:
    VectorSet(const std::vector<Vector<T>>& inputVecs,
              VectorOrientation orient = VectorOrientation::Column)
        : vecs(inputVecs), orientation(orient)
    {
        if (vecs.empty())
            throw std::invalid_argument("Vector set cannot be empty");

        dim = vecs[0].size();
        for (const auto& v : vecs) {
            if (v.size() != dim)
                throw std::invalid_argument("Vectors have inconsistent dimensions");
        }

        Matrix<T> A(
            orientation == VectorOrientation::Column ? dim : vecs.size(),
            orientation == VectorOrientation::Column ? vecs.size() : dim
        );

        if (orientation == VectorOrientation::Column) {
            for (size_t c = 0; c < vecs.size(); c++)
                for (size_t r = 0; r < dim; r++)
                    A.at(r, c) = vecs[c][r];
        } else {
            for (size_t r = 0; r < vecs.size(); r++)
                for (size_t c = 0; c < dim; c++)
                    A.at(r, c) = vecs[r][c];
        }

        RREF<T> rref(std::move(A));
        rref.toRREF();

        rank = rref.getRank();
        pivotIndices = (orient == VectorOrientation::Column)
            ? rref.getPivotCols()
            : rref.getPivotRows();
    }

    // 辅助构造函数：支持原始 std::vector<std::vector<T>>
    VectorSet(const std::vector<std::vector<T>>& inputRaw,
              VectorOrientation orient = VectorOrientation::Column)
        : orientation(orient) {
        for (const auto& row : inputRaw) {
            vecs.emplace_back(row);
        }
        // 不能在构造函数体内赋值委派，直接进行初始化逻辑
        if (vecs.empty())
            throw std::invalid_argument("Vector set cannot be empty");

        dim = vecs[0].size();
        for (const auto& v : vecs) {
            if (v.size() != dim)
                throw std::invalid_argument("Vectors have inconsistent dimensions");
        }

        Matrix<T> A(
            orientation == VectorOrientation::Column ? dim : vecs.size(),
            orientation == VectorOrientation::Column ? vecs.size() : dim
        );

        if (orientation == VectorOrientation::Column) {
            for (size_t c = 0; c < vecs.size(); c++)
                for (size_t r = 0; r < dim; r++)
                    A.at(r, c) = vecs[c][r];
        } else {
            for (size_t r = 0; r < vecs.size(); r++)
                for (size_t c = 0; c < dim; c++)
                    A.at(r, c) = vecs[r][c];
        }

        RREF<T> rref(std::move(A));
        rref.toRREF();

        rank = rref.getRank();
        pivotIndices = (orient == VectorOrientation::Column)
            ? rref.getPivotCols()
            : rref.getPivotRows();
    }

    bool isLinearIndependent() const noexcept {
        return rank == vecs.size();
    }

    std::vector<Vector<T>> basis() const {
        std::vector<Vector<T>> b;
        for (size_t idx : pivotIndices)
            b.push_back(vecs[idx]);
        return b;
    }

    size_t dimension() const noexcept { return rank; }

    static std::vector<Vector<T>>
    gramSchmidt(const std::vector<Vector<T>>& vectors, bool normalize = false) {
        if (vectors.empty())
            throw std::invalid_argument("Gram-Schmidt: empty vector set");

        int dim = vectors[0].size();
        for (const auto& v : vectors)
            if (v.size() != dim)
                throw std::invalid_argument("Gram-Schmidt: dimension mismatch");

        std::vector<Vector<T>> orth;

        for (const auto& v : vectors) {
            Vector<T> u = v;

            for (const auto& uj : orth) {
                T ip_uj = uj.dot(uj);
                if (std::abs(ip_uj) > 1e-9) {
                    T coeff = v.dot(uj) / ip_uj;
                    u -= (uj * coeff);
                }
            }

            if (u.norm() < 1e-9)
                continue;

            if (normalize) {
                u = u.normalized();
            }

            orth.push_back(u);
        }

        return orth;
    }

};
