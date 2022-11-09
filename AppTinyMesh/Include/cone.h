#ifndef CONE_H
#define CONE_H

#include <vector>
#include <iostream>

#include "mathematics.h"

class Cone{

public:
    Cone();
    Vector Vertex(int) const;
    static const int edge[12];
    static const Vector normal[4];
};

#endif // CONE_H
