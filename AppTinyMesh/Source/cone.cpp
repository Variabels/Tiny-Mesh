#include "cone.h"

Cone::Cone()
{

}

const int Cone::edge[] =
{
    0, 1, 2,
    2, 1, 3,
    3, 1, 4,
    4, 1, 0
};

const Vector Cone::normal[] =
{
  Vector(cos(0)/sqrt(2.f),1.f/sqrt(2.f),sin(0)/sqrt(2.f)),
  Vector(cos(M_PI/2)/sqrt(2.f),1.f/sqrt(2.f),sin(M_PI/2)/sqrt(2.f)),
  Vector(cos(M_PI)/sqrt(2.f),1.f/sqrt(2.f),sin(M_PI)/sqrt(2.f)),
  Vector(cos(3*M_PI/2)/sqrt(2.f),1.f/sqrt(2.f),sin(3*M_PI/2)/sqrt(2.f))
};
