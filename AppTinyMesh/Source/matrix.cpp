#include "matrix.h"


Matrix::Matrix()
{
    row = 0;
    col = 0;
    m = new double*[row];
}

Matrix::Matrix(int r, int c)
{
    row = r;
    col = c;
    m = new double*[row];
    for (int i = 0; i<row; ++i){
        m[i] = new double[col];
    }
    for(int i=0; i<row; ++i){
        for(int j=0; j<col; j++){
            m[i][j]=0;
        }
    }
}

Matrix::Matrix(const Matrix& m1) : row(m1.row), col(m1.col)
{
    m = new double*[row];
    for (int i = 0; i<row; ++i){
        m[i] = new double[col];
    }
    for(int i=0; i<row; ++i){
        for(int j=0; j<col; j++){
            m[i][j]=m1.m[i][j];
        }
    }
}

Matrix::~Matrix()
{
    for (int i = 0; i<row; ++i) {
        delete[] m[i];
    }
    delete[] m;
}


void Matrix::InsertValue(double value, int r, int c)
{
    m[r][c] = value;
}

double Matrix::GetValue(int r, int c) const
{
    return m[r][c];
}

int Matrix::getRow()
{
    return row;
}

int Matrix::getCol()
{
    return col;
}

//Member Operators
Matrix& Matrix::operator=(const Matrix& b)
{
    if (this == &b) {
        return *this;
    }

    if (row != b.row || col != b.col) {
        for (int i = 0; i < row; ++i) {
            delete[] m[i];
        }
        delete[] m;

        row = b.row;
        col = b.col;
        m = new double*[row];
        for (int i = 0; i<row; ++i){
            m[i] = new double[col];
        }
        for(int i=0; i<row; ++i){
            for(int j=0; j<col; j++){
                m[i][j]=0;
            }
        }
    }

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            m[i][j] = b.m[i][j];
        }
    }
    return *this;
}


Matrix& Matrix::operator+=(const Matrix& b)
{
    for (int i = 0; i<row; ++i) {
        for (int j = 0; j<col; ++j) {
            m[i][j] += b.m[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& b)
{
    for (int i = 0; i<row; ++i) {
        for (int j = 0; j<col; ++j) {
            m[i][j] -= b.m[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& b)
{
    Matrix t(row, b.col);
        for (int i = 0; i<t.row; ++i) {
            for (int j = 0; j<t.col; ++j) {
                for (int k = 0; k<col; ++k) {
                    t.m[i][j] += (m[i][k] * b.m[k][j]);
                }
            }
        }
        return (*this = t);
}

Matrix& Matrix::operator*=(const double& x)
{
    for (int i = 0; i<row; ++i) {
        for (int j = 0; j<col; ++j) {
            m[i][j] *= x;
        }
    }
    return *this;
}

//Operators


Matrix operator+(const Matrix& a, const Matrix& b)
{
    Matrix t(a);
    return (t+=b);
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
    Matrix t(a);
    return (t-=b);
}

Matrix operator*(const Matrix& a, const Matrix& b)
{
    Matrix t(a);
    return (t*=b);
}

Matrix operator*(const Matrix& a, const double& b)
{
    Matrix t(a);
    return (t*=b);
}
