#include "heightField.h"
#include <iostream>

HeightField::HeightField(const QString& name, double scale)
{
    QImage img(name);

    int width = img.width();
    int height = img.height();

    // Reservation espace
    varray.reserve(2*(width-1)*(height-1)*3);
    narray.reserve(2*(width-1)*(height-1)*3);
    vertices.reserve(width*height);
    normals.reserve(width*height);

    // Vertexes
    for (int i=0; i<height; i++)
    {
        for (int j=0; j<width; j++)
        {
          QRgb value_rgb(img.pixel(j, i));
          double value = qGray(value_rgb);
          vertices.push_back(Vector(10.0*double(j)/double(width-1)-5.0, 10*double(i)/double(height-1)-5.0, value/255.0*scale));
        }
    }

    //Normales
    Vector n1, n2, n3, n4, n;
    for (int i=0; i<height; i++)
    {
        for (int j=0; j<width; j++)
        {
          if (i == 0)
            n1 = vertices[i*width+j];
          else
            n1 = vertices[(i-1)*width+j];

          if (i == height-1)
            n4 = vertices[i*width+j];
          else
            n4 = vertices[(i+1)*width+j];

          if (j == 0)
            n2 = vertices[i*width+j];
          else
            n2 = vertices[i*width+j-1];

          if (j == width-1)
            n3 = vertices[i*width+j];
          else
            n3 = vertices[i*width+j+1];

          n = (n3 - n2) / (n4 - n1);
          normals.push_back(n/Norm(n));
        }
    }

    //Triangles
    for (int i=0; i<height-1; i++)
    {
        for (int j=0; j<width-1; j++)
        {
          AddSmoothTriangle(i*width+j, i*width+j, i*width+j+1, i*width+j+1, (i+1)*width+j+1, (i+1)*width+j+1);
          AddSmoothTriangle(i*width+j, i*width+j, (i+1)*width+j+1, (i+1)*width+j+1, (i+1)*width+j, (i+1)*width+j);
        }
    }
}
