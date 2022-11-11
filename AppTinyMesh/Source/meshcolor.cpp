#include "meshcolor.h"

/*!
\brief Create an empty mesh.
*/
MeshColor::MeshColor()
{
}

/*!
\brief Constructor from a Mesh with color array and indices.
\param m Base mesh.
\param cols Color array.
\param carr Color indexes, should be the same size as Mesh::varray and Mesh::narray.
*/
MeshColor::MeshColor(const Mesh& m, const std::vector<Color>& cols, const std::vector<int>& carr) : Mesh(m), colors(cols), carray(carr)
{
}

/*!
\brief Constructor from a Mesh.
\param m the base mesh
*/
MeshColor::MeshColor(const Mesh& m) : Mesh(m)
{
	colors.resize(vertices.size(), Color(1.0, 1.0, 1.0));
	carray = varray;
}

/*!
\brief Empty.
*/
MeshColor::~MeshColor()
{
}

/**
 * Couleurs pour les heightFields
 * m mesh
 * scale l'echelle
 * flatSea si vrai, applatie la mer
*/
MeshColor::MeshColor(const Mesh& m, double scale, bool flatSea) : Mesh(m)
{
    //Les niveaux
    double ocean = 0.0;
    double shore = 0.2;
    double plain = 0.3;
    double elevation = 0.4;
    double hill = 0.5;
    double mountain = 0.6;
    double peak = 0.9;

    //Couleurs
    Color oceanC = Color(82, 176, 234);
    Color shoreC = Color(255, 247, 223);
    Color plainC = Color(132, 180, 90);
    Color elevationC = Color(10, 141, 0);
    Color hillC = Color(0, 87,0);
    Color mountainC = Color(87, 99, 93);
    Color peakC = Color(255, 244, 255);

    //Reservation
    colors.reserve(vertices.size());
    carray = varray;

    //Distribution des couleurs selon le niveau du vertexe
    Color vColor;
    for (int i=0; i<vertices.size(); i++)
    {
        double altitude = vertices[i]*Vector::Z/scale; // Between 0 and 1
        if (altitude < shore)
        {
            vColor = Color::Lerp((altitude - ocean) / (shore - ocean), oceanC, shoreC);
            if(flatSea)
                vertices[i] = Vector(vertices[i]*Vector::X, vertices[i]*Vector::Y, shore*scale);
        }
        else if(altitude < plain)
        {
            vColor = Color::Lerp((altitude - shore) / (plain - shore), shoreC, plainC);
            if (flatSea)
            {
                vertices[i] = Vector(vertices[i]*Vector::X, vertices[i]*Vector::Y, plain * scale);
            }
        }
        else if(altitude < elevation)
        {
            vColor = Color::Lerp((altitude - plain) / (elevation - plain), plainC, elevationC);
        }
        else if(altitude < hill)
        {
            vColor = Color::Lerp((altitude - elevation) / (hill - elevation), elevationC, hillC);
        }
        else if(altitude < mountain)
        {
            vColor = Color::Lerp((altitude - hill) / (mountain - hill), hillC, mountainC);
        }
        else
            vColor = Color::Lerp((altitude - mountain) / (peak - mountain), mountainC, peakC);

        colors.push_back(vColor);
    }
}
