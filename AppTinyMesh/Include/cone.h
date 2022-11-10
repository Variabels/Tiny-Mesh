#ifndef CONE_H
#define CONE_H
#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>

#include "mathematics.h"

class Cone{

    public:
        Cone() {}
        explicit Cone(double,const Vector&,const Vector&);

        ~Cone() {};

        Vector Vertex(int) const;
        double Radius() const;

    private:
        Vector lower, upper;
        double radius;
};

#endif // CONE_H
