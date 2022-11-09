#include "cylinder.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

Cylinder::Cylinder(double r, double h){
    int div = 20;
    double step = 2.0 * M_PI / (div);
    double alpha;

    for(int i=0; i<div; i++)
    {
        alpha = i*step;
        triangles[i][0] = Vector(cos(alpha),0.0,sin(alpha));
        triangles[i][1] = Vector(cos((i+1)*step),0.0,sin((i+1)*step));
        triangles[i][2] = Vector(0.0,h,0.0);
        Matrix m = Matrix(3, 1);
        m.InsertValue(cos(alpha), 0, 0);
        m.InsertValue(0.0, 1, 0);
        m.InsertValue(sin(alpha), 2, 0);
        TabTriangles[i][0] = m;
        Normals[i] = m;
        m.InsertValue(cos(i+1)*step, 0, 0);
        m.InsertValue(0.0, 1, 0);
        m.InsertValue(sin(i+1)*step, 2, 0);
        TabTriangles[i][1] = m;
        m.InsertValue(0.0, 0, 0);
        m.InsertValue(h, 1, 0);
        m.InsertValue(0.0, 2, 0);
        TabTriangles[i][2] = m;


        normal[i] = Vector(cos(alpha), 0.0, sin(alpha));

    }
}

