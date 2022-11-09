#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

class Matrix{

    public:
        Matrix();
        Matrix(int r, int c);
        Matrix(const Matrix &m1);
        ~Matrix();

        void InsertValue(double value, int r, int c);
        double GetValue(int r, int c) const;
        int getRow();
        int getCol();

        Matrix& operator=(const Matrix& b);
        Matrix& operator+=(const Matrix& b);
        Matrix& operator-=(const Matrix& b);
        Matrix& operator*=(const Matrix& b);
        Matrix& operator*=(const double& x);

    private:
        int col, row;
        double** m;
};

Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator-(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, const double& b);


#endif // MATRIX_H
