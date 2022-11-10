#ifndef SPHERE_H
#define SPHERE_H
#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>

#include "mathematics.h"

class Sphere{

    public:
        Sphere() {}
        explicit Sphere(double, const Vector &);

        ~Sphere() {};

        double Radius() const;
        Vector Center() const;

    private:
        double radius;
        Vector center;
};

#endif // SPHERE_H
