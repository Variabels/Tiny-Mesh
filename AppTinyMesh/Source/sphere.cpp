#include "sphere.h"

Sphere::Sphere(double r, const Vector & c){
    center = c;
    radius = r;
}


double Sphere::Radius() const{
    return radius;
}

Vector Sphere::Center() const{
    return center;
}
