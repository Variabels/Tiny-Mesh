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

/**
 * Creation d'un mesh pour un cone
 * c cone
 * div nombre de divisions
*/
Mesh::Mesh(const Cone& c, const int div)
{
    const Vector lower = c.Vertex(0);
    const Vector upper = c.Vertex(1);
    const double radius = c.Radius();
    const double theta = (2*M_PI)/div;

    const Vector n = Normalized(upper-lower);
    Vector x, y;
    n.Orthonormal(x, y);

    vertices.reserve(div+2);


    //Vertexes
    for (int i=0; i<div; i++)
    {
        Vector tmp(x*cos(theta*i)+y*sin(theta*i)+lower);
        tmp *= radius;
        vertices.push_back(tmp);
    }

    vertices.push_back(lower);
    normals.push_back(-n);

    //Ajout des triangles de la base
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

/**
 * Creation d'un mesh de sphere
 * s sphere
 * div nombre de divisions
*/
Mesh::Mesh(const Sphere& s, int div)
{
    double r = s.Radius();
    Vector c = s.Center();

    int hStep = div;
    int vStep = div;
    double x, y, z;

    //Vertexes et normales

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

    //Triangles de la sphere
    for(int i = 0; i<vStep-1; i++){
        AddSmoothTriangle(0, 0, i+2, i+2, i+3, i+3);
        AddSmoothTriangle(1, 1, Vertexes()-i-1, Vertexes()-i-1, Vertexes()-i-2, Vertexes()-i-2);
    }

    AddSmoothTriangle(2, 2, 0, 0, 2+vStep-1, 2+vStep-1);
    AddSmoothTriangle(1, 1, Vertexes()-vStep, Vertexes()-vStep, Vertexes()-1, Vertexes()-1);

    //Sommet
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

/**
 * Creation de la demi sphere
 * s sphere
 * div nombre de divisions
 * nb hauteur de la demi-sphere
*/
Mesh::Mesh(const Sphere& s, const int div, const int nb)
{
  double alpha,beta;
  int divBeta = div;

  vertices.resize(nb*divBeta+2);
  normals.resize(nb*divBeta+2);

  // Reservation espace
  varray.reserve((nb*divBeta)*6);
  narray.reserve((nb*divBeta)*6);

  //Vertexes
  for (int i=0; i<nb; i++)
  {

      alpha = -0.5*M_PI+ double(i+1)*M_PI/nb/2;

      for (int j=0; j<divBeta; j++)
      {
          beta = (double(j)*2.0*M_PI)/(double)divBeta;
          vertices[j+i*divBeta] = Vector(cos(alpha)*cos(beta), cos(alpha)*sin(beta), -sin(alpha))*s.Radius()+s.Center();
          normals[j+i*divBeta] = vertices[j+i*divBeta]-s.Center();
      }
  }

  //Vertexes sommet
  vertices[nb*divBeta] = Vector(0.0,0.0,1.0) * s.Radius()+s.Center();
  vertices[nb*divBeta+1] = s.Center();


  //normales
  normals[nb*divBeta] = Vector(0.0,0.0,1.0);
  normals[nb*divBeta+1] = Vector(0.0,0.0,-1.0);

  //Triangles de la demi-sphere
  for (int i=0; i<nb-1; i++)
  {
      for (int j=0; j<divBeta; j++)
      {
          AddSmoothTriangle(j+i*divBeta,j+i*divBeta, (j+1)%divBeta+i*divBeta,
                            (j+1)%divBeta+i*divBeta, j+divBeta+i*divBeta, j+divBeta + i*divBeta);
          AddSmoothTriangle((j+1)%divBeta+i*divBeta, (j+1)%divBeta+i*divBeta, (j+1)%divBeta+divBeta+i*divBeta,
                            (j+1)%divBeta+divBeta+i*divBeta, j+divBeta+i*divBeta, j+divBeta+i*divBeta);
      }
  }

  //Sommet de la demi-sphere
  for (int j=0; j<divBeta; j++)
  {
      AddSmoothTriangle(nb*divBeta,nb*divBeta, (j+1)%divBeta, (j+1)%divBeta, j, j);
      AddSmoothTriangle(nb*divBeta+1, nb*divBeta+1, (j+1)%divBeta+(nb-1)*divBeta,
                        (j+1)%divBeta+(nb-1)*divBeta, j+(nb-1)*divBeta, j+(nb-1)*divBeta);
  }
}


/**
 * Capsule faite grace a des merges d'un cylindre et deux demi-spheress
 * c une capsule
 * div nombre de divisions
*/
Mesh::Mesh(const Capsule& c, const int div){

    Sphere s(c.Radius(),c.Center()+Vector(0.0,0.0,c.Height()/2-c.Radius()));
    this->Merge(Mesh(s,div,10));
    Cylinder cyl(c.Radius(),c.Height()-c.Radius()*2,c.Center()+Vector(0.0,0.0,c.Height()/2+1));
    this->Merge(Mesh(cyl,div,10));
    Sphere cap(c.Radius(), c.Center()+Vector(0.0,0.0,c.Height()/2-c.Radius()-2.5));
    Mesh caps = Mesh(cap, div, 10);
    caps.Transform(Matrix::rotateX(M_PI));
    this->Merge(caps);

}

/**
 * Creation d'un cylindre a une face
 * cyl un cylindre
 * div le nombre de divisions du cercle
 * nb hauteur du cylindre
*/
Mesh::Mesh(const Cylinder& cyl, const int div, const int nb)
{
    double alpha;
    double step = 2.0*M_PI/(double)div;
    double step2 = cyl.Height()/(double)nb-1;
    Vector c = cyl.Center()-Vector(0.0,0.0,(cyl.Height()/2));


  //Reservation d'espace
  varray.reserve(div*nb*6);
  narray.reserve(div*nb*6);

  //Bouts/Disques
  for (int i=0; i<nb; ++i) {
      for (int j=0; j<div; j++)
      {
          alpha = (double)j*step;
          vertices.push_back(Vector(cos(alpha),sin(alpha),0.0)*cyl.Radius()+c+Vector(0.0,0.0,step2*i));
          normals.push_back(vertices[(j+div*i)]-(c+Vector(0.0,0.0,step2*i)));
      }
  }

  //La face
  for (int i=0; i<nb-1; ++i) {
      for(int j=0; j<div; j++){
          AddSmoothTriangle(j+div*i,j+div*i, j+div*(i+1),j+div*(i+1), (j+1)%div+div*(i+1),(j+1)%div+div*(i+1));
          AddSmoothTriangle(j+div*i,j+div*i, (j+1)%div+div*(i+1),(j+1)%div+div*(i+1), (j+1)%div+div*i,(j+1)%div+div*i);
      }
  }
}


/**
 * Creation d'un cylindre
 * c un cylindre
 * div le nombre de divisions
*/
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


/**
 * Creation d'un mesh pour un torus
 * t un torus
 * div un entier pour savoir la resolution
*/
Mesh::Mesh(const Torus& t, const int div)
{
    double r = t.Radius();
    double thickness = t.Thickness();
    double x, y, z;
    int hStep = div;
    int vStep = div;

    //Calcul des normales

    for (int i = 0; i<hStep; i++) {
        for (int j = 0; j<vStep; j++) {

            x = (r + thickness*cos(2*M_PI*(double)i/(double)hStep))*cos(2*M_PI*(double)j/(double)vStep);
            y = (r + thickness*cos(2*M_PI*(double)i/(double)hStep))*sin(2*M_PI*(double)j/(double)vStep);
            z = thickness*sin(2*M_PI*(double)i/(double)hStep);

            vertices.emplace_back(Vector(x,y,z));
            normals.push_back(Normalized(vertices.back()));
        }
    }

    //Ajout des triangles

    for (int i = 0; i<hStep-1; i++) {
        for (int j = 0; j < vStep; j++) {
            int v1 = i*vStep+j;
            int v4 = i*vStep+(j+1)%vStep;
            int v2 = (i+1)*vStep+j;
            int v3 = (i+1)*vStep+(j+1)%vStep;

            AddSmoothQuadrangle(v1,v1,v2,v2,v3,v3,v4,v4);
        }
    }

    //Derniere ligne de triangles

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


//Fusion de deux mesh
void Mesh::Merge(const Mesh& m){

    int v = vertices.size();
    int n = normals.size();

    for(int i=0; i<m.narray.size(); i++)
    {
        narray.push_back(m.narray[i]+n);
    }
    for(int i=0; i<m.varray.size(); i++)
    {
        varray.push_back(m.varray[i]+v);
    }
    for(int i=0; i<m.normals.size(); i++)
    {
        normals.push_back(m.normals[i]);
    }
    for(int i=0; i<m.vertices.size(); i++)
    {
        vertices.push_back(m.vertices[i]);
    }

}

//Permet d'appliquer des rotations grace a la matrice m
void Mesh::Transform(const Matrix& m){

    for(int i=0; i<vertices.size(); i++)
    {
        vertices[i] = m.multiplyVector(vertices[i]);
    }
    for(int i=0; i<normals.size(); i++)
    {
        normals[i]= m.multiplyVector(normals[i]);
    }
}


//Fonction de deformation
void Mesh::SphereWarp(const Sphere& s,const Vector& d){
    double n;
    Vector tmp;
    for(int i=0; i<vertices.size(); i++){
        n = (s.Center()[0]-vertices[i][0])*(s.Center()[0]-vertices[i][0])+
                (s.Center()[1]-vertices[i][1])*(s.Center()[1]-vertices[i][1])+
                (s.Center()[2]-vertices[i][2])*(s.Center()[2]-vertices[i][2]);
        if(n<(s.Radius()*s.Radius())){
            tmp=d*(1-(n/(s.Radius()*s.Radius())));
            vertices[i]=vertices[i]+tmp;
        }
    }

}
