// =========================================================
// vector.h — 单个向量的原子操作 (Layer 0, 无项目内依赖)
// ---------------------------------------------------------
// 职责: 四则运算、点积、范数、归一化、正交判定
// 向量组级操作 (线性无关、基、Gram-Schmidt) 见 VectorSet.h
// =========================================================
#pragma once

#include<vector>
#include<cmath>
#include<stdexcept>

template<typename T>
class Vector{
    private:
        std::vector<T> data;
    public:
        Vector() = default;
        explicit Vector(int n, T val = T()) : data(n, val) {}
        explicit Vector(const std::vector<T>& v) : data(v) {}
        explicit Vector(std::vector<T>&& v) noexcept : data(std::move(v)) {}

        // 移动语义
        Vector(Vector&& other) noexcept : data(std::move(other.data)) {}
        Vector& operator=(Vector&& other) noexcept {
            if (this != &other) {
                data = std::move(other.data);
            }
            return *this;
        }
        Vector(const Vector&) = default;
        Vector& operator=(const Vector&) = default;

        int size() const noexcept { return static_cast<int>(data.size()); }
        
        T& operator[](int i) { return data.at(i); }
        const T& operator[](int i) const { return data.at(i); }

        const std::vector<T>& raw() const noexcept { return data; }

        // 四则运算
        Vector<T> operator+(const Vector<T>& other) const {
            if (size() != other.size())
                throw std::invalid_argument("Vector size mismatch");
            std::vector<T> res_vec(size());
            for (int i = 0; i < size(); i++)
                res_vec[i] = data[i] + other[i];
            return Vector<T>(res_vec);
        }

        Vector<T> operator-(const Vector<T>& other) const {
            if (size() != other.size())
                throw std::invalid_argument("Vector size mismatch");
            std::vector<T> res_vec(size());
            for (int i = 0; i < size(); i++)
                res_vec[i] = data[i] - other[i];
            return Vector<T>(res_vec);
        }

        Vector<T> operator*(T scalar) const {
            std::vector<T> res_vec(size());
            for (int i = 0; i < size(); i++)
                res_vec[i] = data[i] * scalar;
            return Vector<T>(res_vec);
        }

        Vector<T> operator/(T scalar) const {
            if (std::abs(scalar) < 1e-9)
                throw std::invalid_argument("Division by zero");
            return (*this) * (1.0 / scalar);
        }

        Vector<T>& operator+=(const Vector<T>& other) {
            if (size() != other.size())
                throw std::invalid_argument("Vector size mismatch");
            for (int i = 0; i < size(); i++)
                data[i] += other.data[i];
            return *this;
        }

        Vector<T>& operator-=(const Vector<T>& other) {
            if (size() != other.size())
                throw std::invalid_argument("Vector size mismatch");
            for (int i = 0; i < size(); i++)
                data[i] -= other.data[i];
            return *this;
        }

        Vector<T>& operator*=(T scalar) {
            for (auto& el : data) el *= scalar;
            return *this;
        }

        Vector<T>& operator/=(T scalar) {
            if (std::abs(scalar) < 1e-9)
                throw std::invalid_argument("Division by zero");
            for (auto& el : data) el /= scalar;
            return *this;
        }

        friend Vector<T> operator*(T scalar, const Vector<T>& v) {
            return v * scalar;
        }

        T dot(const Vector<T>& other) const {
            if (size() != other.size())
                throw std::invalid_argument("Dot product size mismatch");

            T sum = 0;
            for (int i = 0; i < size(); i++)
                sum += data[i] * other[i];
            return sum;
        }

        T norm() const {
            return std::sqrt(this->dot(*this));
        }

        Vector<T> normalized(T eps = 1e-9) const {
            T n = norm();
            if (n < eps)
                throw std::invalid_argument("Cannot normalize zero vector");
            return (*this) * (1 / n);
        }

        bool isOrthogonalTo(const Vector<T>& other, T eps = 1e-9) const {
            return std::fabs(dot(other)) < eps;
        }

        void print() const {
            for (const auto& el : data) {
                std::cout << el << " ";
            }
            std::cout << std::endl;
        }
};