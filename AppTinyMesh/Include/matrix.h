#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <mathematics.h>

class Matrix{

    public:
        Matrix();
        Matrix(int r, int c);
        Matrix(Vector, Vector, Vector);
        Matrix(const Matrix &m1);
        ~Matrix();

        void InsertValue(double value, int r, int c);
        double GetValue(int r, int c) const;
        int getRow();
        int getCol();

        Matrix identity();
        Matrix transpose();
        Matrix inverse();

        static const Matrix rotateX(const double&);
        static const Matrix rotateY(const double&);
        static const Matrix rotateZ(const double&);
        Matrix homothetie(const double&, const double&, const double&);

        Vector multiplyVector(const Vector& v) const;

        Matrix& operator=(const Matrix& b);
        Matrix& operator+=(const Matrix& b);
        Matrix& operator-=(const Matrix& b);
        Matrix& operator*=(const Matrix& b);
        Matrix& operator*=(const double& x);
        Vector operator*=(const Vector& v);
        Matrix& operator/=(const double& x);

    private:
        int col, row;
        double** m;
};

Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator-(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, const double& b);
Vector operator*(const Matrix& a, const Vector& v);
Matrix operator/(const Matrix& a, const double& b);


#endif // MATRIX_H
