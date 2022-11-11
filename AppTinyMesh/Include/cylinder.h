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
        explicit Cylinder(double,double,const Vector&);

        ~Cylinder() {};

        double Radius() const;
        double Height() const;
        Vector Center() const;
        Vector Vertex(int) const;

    private:
        Vector lower, upper;
        Vector center;
        double radius, height;
};

#endif // CYLINDER_H
