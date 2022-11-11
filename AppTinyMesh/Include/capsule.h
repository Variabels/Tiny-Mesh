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

        Vector Center() const;
        double Radius() const;
        double Height() const;

    private:
        double radius, height;
        Vector center;
   };

#endif // CAPSULE_H
