#ifndef CYLINDER_H
#define CYLINDER_H

#pragma once

#include <vector>
#include <iostream>
#include "mathematics.h"
#include "matrix.h"

class Cylinder
{
    public:
        Cylinder() {}
        explicit Cylinder(double,double);



    public:
        static const int nbr_edges;
        static const int* edge;
        static const int nbr_normals;
        Vector normal[20];
        Vector triangles[20][3];
        Matrix Normals[20];
        Matrix TabTriangles[20][3];
        double coordonnee(int, int);
};

#endif // CYLINDER_H
