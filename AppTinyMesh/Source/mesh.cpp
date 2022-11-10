#include "mesh.h"

/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
*/



/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Initialize the mesh from a list of vertices and a list of triangles.

Indices must have a size multiple of three (three for triangle vertices and three for triangle normals).

\param vertices List of geometry vertices.
\param indices List of indices wich represent the geometry triangles.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), varray(indices)
{
  normals.resize(vertices.size(), Vector::Z);
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& va, const std::vector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
  // Initialize 
  normals.resize(vertices.size(), Vector::Null);

  narray = varray;

  // Accumulate normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize 
  for (int i = 0; i < normals.size(); i++)
  {
    Normalize(normals[i]);
  }
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  varray.push_back(a);
  narray.push_back(na);
  varray.push_back(b);
  narray.push_back(nb);
  varray.push_back(c);
  narray.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.push_back(a);
  narray.push_back(n);
  varray.push_back(b);
  narray.push_back(n);
  varray.push_back(c);
  narray.push_back(n);
}

/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Index of the vertices.
\param na, nb, nc, nd Index of the normal for all vertices.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
  // First triangle
  AddSmoothTriangle(a, na, b, nb, c, nc);

  // Second triangle
  AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Index of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
  AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
  if (vertices.size() == 0)
  {
    return Box::Null;
  }
  return Box(vertices);
}

/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box)
{
  // Vertices
  vertices.resize(8);

  for (int i = 0; i < 8; i++)
  {
    vertices[i] = box.Vertex(i);
  }

  // Normals
  normals.push_back(Vector(-1, 0, 0));
  normals.push_back(Vector(1, 0, 0));
  normals.push_back(Vector(0, -1, 0));
  normals.push_back(Vector(0, 1, 0));
  normals.push_back(Vector(0, 0, -1));
  normals.push_back(Vector(0, 0, 1));

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 4);
  AddTriangle(1, 2, 3, 4);

  AddTriangle(4, 5, 6, 5);
  AddTriangle(5, 7, 6, 5);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 0);

  AddTriangle(1, 3, 5, 1);
  AddTriangle(3, 7, 5, 1);

  AddTriangle(0, 1, 5, 2);
  AddTriangle(0, 5, 4, 2);

  AddTriangle(3, 2, 7, 3);
  AddTriangle(6, 7, 2, 3);

}


Mesh::Mesh(const Cone& c, const int div)
{
    const Vector lower = c.Vertex(0);
    const Vector upper = c.Vertex(1);
    const double radius = c.Radius();

    const Vector n = Normalized(upper-lower);
    Vector x, y;
    n.Orthonormal(x, y);

    vertices.reserve(div+2);

    // Circle slice size
    const double theta = (2*M_PI)/div;

    //Base
    for (int i=0; i<div; i++)
    {
        Vector tmp(x*cos(theta*i)+y*sin(theta*i)+lower);
        tmp *= radius;
        vertices.push_back(tmp);
    }

    vertices.push_back(lower);
    normals.push_back(-n);
    for (int i=0; i<div; i++)
        AddTriangle(vertices.size()-1, i, (i+1)%div, normals.size()-1);

    //Face
    vertices.push_back(upper);
    for (int i=0; i<div; i++)
    {
        Vector normal = Normalized(vertices[i]-n);
        normals.push_back(normal);
        AddTriangle(vertices.size()-1, i, ((i+1)%div), normals.size()-1);
    }

}

Mesh::Mesh(const Sphere& s, int div)
{
    double r = s.Radius();
    Vector c = s.Center();

    int hStep = div;
    int vStep = div;
    double x, y, z;

    vertices.emplace_back(Vector(c[0], c[1], c[2]+r));
    normals.push_back(Normalized(vertices.back()));
    vertices.emplace_back(Vector(c[0], c[1], c[2]-r));
    normals.push_back(Normalized(vertices.back()));

    for(int i=1; i<hStep; i++){
        for(int j=0; j<vStep; j++){

            x = sin(M_PI*(double)i/(double)hStep)*cos(2*M_PI*(double)j/(double)vStep)*r;
            y = sin(M_PI*(double)i/(double)hStep)*sin(2*M_PI*(double)j/(double)vStep)*r;
            z = cos(M_PI*(double)i/(double)hStep)*r;

            vertices.emplace_back(Vector(x, y, z));
            normals.push_back(Normalized(vertices.back()));
        }
    }

    for(int i = 0; i<vStep-1; i++){
        AddSmoothTriangle(0, 0, i+2, i+2, i+3, i+3);
        AddSmoothTriangle(1, 1, Vertexes()-i-1, Vertexes()-i-1, Vertexes()-i-2, Vertexes()-i-2);
    }

    AddSmoothTriangle(2, 2, 0, 0, 2+vStep-1, 2+vStep-1);
    AddSmoothTriangle(1, 1, Vertexes()-vStep, Vertexes()-vStep, Vertexes()-1, Vertexes()-1);


    for(int i = 0; i<hStep-2; i++){
        for(int j = 0; j<vStep; j++){
            int v1 = i*vStep+j+ 2;
            int v4 = i*vStep+(j+1)%vStep+2;
            int v2 = (i+1)*vStep+j+2;
            int v3 = (i+1)*vStep+(j+1)%vStep+2;

            AddSmoothTriangle(v1, v1, v2, v2, v3, v3);
            AddSmoothTriangle(v4, v4, v1, v1, v3, v3);
        }
    }
}

Mesh::Mesh(const Capsule& c, const int div){

    double r = c.Radius();
    Vector lower = c.Vertex(0);
    Vector upper = c.Vertex(1);

    int hStep = div/2;
    int vStep = div;
    double x, y, z;


    //Sphere du bas
    vertices.emplace_back(Vector(lower[0], lower[1], lower[2]-r));
    normals.push_back(Normalized(vertices.back()));

    for(int i=1; i<hStep; i++){
        for(int j=0; j<vStep; j++){

            x = lower[0]+sin(M_PI/2.*(double)i/(double)hStep)*cos(2*M_PI*(double)j/(double)vStep)*r;
            y = lower[1]+sin(M_PI/2.*(double)i/(double)hStep)*sin(2*M_PI*(double)j/(double)vStep)*r;
            z = lower[2]-cos(M_PI/2.*(double)i/(double)hStep)*r;


            vertices.emplace_back(Vector(x, y, z));
            normals.push_back(Normalized(vertices.back()));
        }
    }

    for(int i=0; i<vStep-1; i++){
        AddSmoothTriangle( i+1, i+1, 0, 0, i+2, i+2);
    }

    AddSmoothTriangle( 0, 0, 1, 1, vStep, vStep);

    for(int i=0; i<hStep-2; i++){
        for(int j=0; j<vStep; j++){
            int v1 = i*vStep+j+1;
            int v4 = i*vStep+(j+1)%vStep+1;
            int v2 = (i+1)*vStep+j+1;
            int v3 = (i+1)*vStep+(j+1)%vStep+1;

            AddSmoothTriangle(v2,v2, v1, v1,v3, v3);
            AddSmoothTriangle(v1, v1, v4, v4, v3, v3);
        }
    }

    //Sphere du haut
    int nbVertex = Vertexes();

    vertices.emplace_back(Vector(upper[0], upper[1], upper[2]+r));
    normals.push_back(Normalized(vertices.back()));

    for(int i=1; i<hStep; i++){
        for(int j=0; j<vStep; j++){

            x = upper[0]+sin(M_PI/2.*(double)i/(double)hStep)*cos(2*M_PI*(double)j/(double)vStep)*r;
            y = upper[1]+sin(M_PI/2.*(double)i/(double)hStep)*sin(2*M_PI*(double)j/(double)vStep)*r;
            z = upper[2]+cos(M_PI/2.*(double)i/(double)hStep)*r;

            vertices.emplace_back(Vector(x, y, z));
            normals.push_back(Normalized(vertices.back()));
        }
    }


    for(int i=0; i<vStep-1; i++){
        AddSmoothTriangle(nbVertex, nbVertex, nbVertex+i+1, nbVertex+i+1, nbVertex+i+2, nbVertex+i+2);
    }

    AddSmoothTriangle(nbVertex+1, nbVertex+1, nbVertex, nbVertex, nbVertex+vStep, nbVertex+vStep);

    for(int i=0; i<hStep-2; i++){
        for(int j=0; j<vStep; j++){
            int v1 = i*vStep+j+1+nbVertex;
            int v4 = i*vStep+(j+1)%vStep+1+nbVertex;
            int v2 = (i+1)*vStep+j+1+nbVertex;
            int v3 = (i+1)*vStep+(j+1)%vStep+1+nbVertex;

            AddSmoothTriangle(v1, v1, v2, v2, v3, v3);
            AddSmoothTriangle(v4, v4, v1, v1, v3, v3);
        }
    }


    //todo: fix cylinder

    int v1 = Vertexes()/2-vStep;
    int v2 = Vertexes()-vStep;

    for(int i=0; i<vStep-1; i++){
        AddSmoothTriangle(v1+i, v1+i, v1+i+1, v1+i+1, v2+i, v2+i);
        AddSmoothTriangle(v2+i, v2+i, v1+i+1, v1+i+1, v2+i+1, v2+i+1);
    }

    AddSmoothTriangle(v1+vStep-1, v1+vStep-1, v1, v1, v2, v2);
    AddSmoothTriangle(v2, v2, v2+vStep-1, v2+vStep-1, v1+vStep-1, v1+vStep-1);

}

Mesh::Mesh(const Cylinder& c, int div)
{
    const Vector lower = c.Vertex(0);
    const Vector upper = c.Vertex(1);
    const double radius = c.Radius();

    const Vector n = Normalized(upper-lower);
    Vector x, y;
    n.Orthonormal(x, y);

    vertices.reserve((div*2)+2);

    const double theta = (2*M_PI)/div;

    //Cercle du haut
    for (int i=0; i<div; i++)
    {
        Vector tmp(x*cos(theta*i)+y*sin(theta*i)+lower);
        tmp *= radius;
        vertices.push_back(tmp);
    }

    vertices.push_back(lower);
    normals.push_back(-n);
    for (int i=0; i<div; i++)
        AddTriangle(vertices.size()-1, i, (i+1)%div, normals.size()-1);

    //Cercle du bas
    int offset = vertices.size();
    for (int i=0; i<div; i++)
    {
        Vector tmp(x*cos(theta*i)+y*sin(theta*i)+upper);
        tmp *= radius;
        vertices.push_back(tmp);
    }

    vertices.push_back(upper);
    normals.push_back(n);
    for (int i=0; i<div; i++)
        AddTriangle(vertices.size()-1, i+offset, ((i+1)%div)+offset, normals.size()-1);

    //Face
    for (int i=0; i<div; i++)
    {
        Vector normal = Normalized(vertices[i]-n);
        normals.push_back(normal);
        AddTriangle(i, (i+1)%div, offset+i, normals.size()-1);
        AddTriangle(offset+i, offset+((1+i)%div), (i+1)%div, normals.size()-1);
    }
}

Mesh::Mesh(const Torus& t, const int res)
{
    double r = t.Radius();
    double thickness = t.Thickness();
    double x, y, z;
    int hStep = res;
    int vStep = res;

    for (int i = 0; i<hStep; i++) {
        for (int j = 0; j<vStep; j++) {

            x = (r + thickness*cos(2*M_PI*(double)i/(double)hStep))*cos(2*M_PI*(double)j/(double)vStep);
            y = (r + thickness*cos(2*M_PI*(double)i/(double)hStep))*sin(2*M_PI*(double)j/(double)vStep);
            z = thickness*sin(2*M_PI*(double)i/(double)hStep);

            vertices.emplace_back(Vector(x,y,z));
            normals.push_back(Normalized(vertices.back()));
        }
    }

    for (int i = 0; i<hStep-1; i++) {
        for (int j = 0; j < vStep; j++) {
            int v1 = i*vStep+j;
            int v4 = i*vStep+(j+1)%vStep;
            int v2 = (i+1)*vStep+j;
            int v3 = (i+1)*vStep+(j+1)%vStep;

            AddSmoothQuadrangle(v1,v1,v2,v2,v3,v3,v4,v4);
        }
    }

    for (int i = 0; i<vStep-1; i++) {
        int v1 = (hStep-1)*vStep+i;
        int v2 = (hStep-1)*vStep+i+1;
        int v3 = i;
        int v4 = i+1;

        AddSmoothQuadrangle(v4,v4,v2,v2,v1,v1,v3,v3);
    }

    AddSmoothQuadrangle(0, 0,(hStep-1)*vStep,(hStep-1)*vStep,(hStep-1)*vStep+vStep-1,(hStep-1)*vStep+vStep-1,vStep-1,vStep-1);
}


/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
    // Vertexes
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] *= s;
    }

    if (s < 0.0)
    {
        // Normals
        for (int i = 0; i < normals.size(); i++)
        {
            normals[i] = -normals[i];
        }
    }
}



#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QRegularExpression>
#include <QtCore/qstring.h>

/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
  vertices.clear();
  normals.clear();
  varray.clear();
  narray.clear();

  QFile data(filename);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rexn("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rext("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    QRegularExpressionMatch matchN = rexn.match(line);
    QRegularExpressionMatch matchT = rext.match(line);
    if (match.hasMatch())//rexv.indexIn(line, 0) > -1)
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble()); vertices.push_back(q);
    }
    else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
    {
      Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());  normals.push_back(q);
    }
    else if (matchT.hasMatch())//rext.indexIn(line, 0) > -1)
    {
      varray.push_back(matchT.captured(1).toInt() - 1);
      varray.push_back(matchT.captured(3).toInt() - 1);
      varray.push_back(matchT.captured(5).toInt() - 1);
      narray.push_back(matchT.captured(2).toInt() - 1);
      narray.push_back(matchT.captured(4).toInt() - 1);
      narray.push_back(matchT.captured(6).toInt() - 1);
    }
  }
  data.close();
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
    return;
  QTextStream out(&data);
  out << "g " << meshName << Qt::endl;
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}

