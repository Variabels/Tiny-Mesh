#ifndef HEIGHTFIELD_H
#define HEIGHTFIELD_H
#pragma once
#include "mesh.h"
#include <QImage>
#include <QString>

class HeightField : public Mesh
{

    public:

        HeightField() {};
        HeightField(const QString&, double);

        ~HeightField(){};

};

#endif // HEIGHTFIELD_H
