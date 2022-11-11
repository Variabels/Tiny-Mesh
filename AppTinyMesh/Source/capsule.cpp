#include "capsule.h"

Capsule::Capsule(double r, double h, const Vector & c)
{
    center = c;
    height = h;
    radius = r;
}

double Capsule::Radius() const {
    return radius;
}

double Capsule::Height() const {
    return height;
}

Vector Capsule::Center() const {
    return center;
}
