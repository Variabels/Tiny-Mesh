#ifndef CAPSULE_H
#define CAPSULE_H

#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"

class Capsule
{
    public:
        //! Empty.
        Capsule() {}
        explicit Capsule(const double, const double, const Vector&);

        ~Capsule() {};

        Vector Vertex(int) const;
        double Radius() const;

    private:
        double radius;
        Vector upper, lower;
   };

#endif // CAPSULE_H
