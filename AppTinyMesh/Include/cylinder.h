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
        Vector Vertex(int) const;

    private:
        Vector lower, upper;
        double radius;
};

#endif // CYLINDER_H
