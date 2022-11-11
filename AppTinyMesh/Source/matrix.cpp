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

Matrix::Matrix(Vector row1, Vector row2, Vector row3)
{
    row = 3;
    col = 3;
    m = new double*[3];
    for (int i = 0; i<row; ++i){
        m[i] = new double[col];
    }
    for (int i = 0; i<3; ++i){
        m[0][i] = row1[i];
        m[1][i] = row2[i];
        m[2][i] = row3[i];
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

Vector Matrix::operator*=(const Vector& v)
{
    return Vector(v[0]*m[0][0]+v[1]*m[0][1]+v[2]*m[0][2],
            v[0]*m[1][0]+v[1]*m[1][1]+v[2]*m[1][2],
            v[0]*m[2][0]+v[1]*m[2][1]+v[2]*m[2][2]);
}

Vector Matrix::multiplyVector(const Vector& v) const
{
    return Vector(v[0]*m[0][0]+v[1]*m[0][1]+v[2]*m[0][2],
            v[0]*m[1][0]+v[1]*m[1][1]+v[2]*m[1][2],
            v[0]*m[2][0]+v[1]*m[2][1]+v[2]*m[2][2]);
}


Matrix& Matrix::operator/=(const double& x)
{
    for (int i = 0; i<row; ++i) {
        for (int j = 0; j<col; ++j) {
            m[i][j] = m[i][j]/x;
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


Vector operator*=(const Matrix& a, const Vector& v)
{
    Matrix t(a);
    return (t*=v);
}

Matrix operator/(const Matrix& a, const double& b)
{
    Matrix t(a);
    return (t/=b);
}


Matrix Matrix::identity()
{
    return Matrix(Vector(1.0,0.0,0.0),Vector(0.0,1.0,0.0),Vector(0.0,0.0,1.0));
}

Matrix Matrix::transpose()
{
    return Matrix(Vector(m[0][0],m[1][0],m[2][0]),Vector(m[0][1],m[1][1],m[2][1]),Vector(m[0][2],m[1][2],m[2][2]));
}

Matrix Matrix::inverse()
{
    Vector row1, row2, row3;
    double det;

    row1 = Vector(m[1][1]*m[2][2]-m[1][2]*m[2][1],
            m[1][1]*m[2][2]-m[1][2]*m[2][0],
            m[1][1]*m[2][1]-m[1][1]*m[2][0]);
    row2 = Vector(m[0][1]*m[2][2]-m[0][2]*m[2][1],
            m[0][0]*m[2][2]-m[0][2]*m[2][0],
            m[0][0]*m[2][1]-m[0][1]*m[2][0]);
    row3 = Vector(m[0][1]*m[1][2]-m[0][2]*m[1][1],
            m[0][0]*m[1][2]-m[0][2]*m[1][0],
            m[0][0]*m[1][1]-m[0][1]*m[1][0]);

    Matrix comp(row1, row2, row3);

    comp = comp.transpose();

    det = m[0][0]*row1[0] - m[0][1] * row1[1] + m[0][2] * row1[2];

    if(det == 0)
    {
        std::cout<<"Matrix can't be inverted.";
        return *this;
    }
    else {
        return comp/det;
    }
}

const Matrix Matrix::rotateX(const double & a)
{
    return Matrix(Vector(1.0,0.0,0.0), Vector(0.0, cos(a), -sin(a)), Vector(0.0,sin(a),cos(a)));
}

const Matrix Matrix::rotateY(const double & a)
{
    return Matrix(Vector(cos(a),0.0,-sin(a)), Vector(0.0, 1.0,0.0), Vector(sin(a),0.0,cos(a)));
}

const Matrix Matrix::rotateZ(const double & a)
{
    return Matrix(Vector(cos(a),-sin(a),0.0), Vector(sin(a), cos(a),0.0), Vector(0.0,0.0,1.0));
}

Matrix Matrix::homothetie(const double & x, const double& y, const double& z)
{
    return Matrix(Vector(x,0.0,0.0), Vector(0.0, y,0.0), Vector(0.0,0.0,z));
}

