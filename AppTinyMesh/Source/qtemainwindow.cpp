#include "qte.h"
#include "implicits.h"
#include <chrono>

MainWindow::MainWindow()
{
	// Chargement de l'interface
	uiw.setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
	uiw.widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));

}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
    connect(uiw.boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw.coneMesh, SIGNAL(clicked()), this, SLOT(ConeMeshExample()));
    connect(uiw.sphereMesh, SIGNAL(clicked()), this, SLOT(SphereMeshExample()));
    connect(uiw.capsuleMesh, SIGNAL(clicked()), this, SLOT(CapsuleMeshExample()));
    connect(uiw.sphereWarpMesh, SIGNAL(clicked()), this, SLOT(SphereWarpExample()));
    connect(uiw.torusMesh, SIGNAL(clicked()), this, SLOT(TorusExample()));
    connect(uiw.cylinderMesh, SIGNAL(clicked()), this, SLOT(CylinderMeshExample()));
    connect(uiw.mushroomMesh, SIGNAL(clicked()), this, SLOT(MushroomExample()));
    connect(uiw.heightFieldMesh, SIGNAL(clicked()), this, SLOT(HeightFieldExample()));
	connect(uiw.sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
	connect(uiw.resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
	connect(uiw.wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
	connect(uiw.radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
	connect(uiw.radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BoxMeshExample()
{
    auto start = std::chrono::high_resolution_clock::now();

	Mesh boxMesh = Mesh(Box(1.0));

	std::vector<Color> cols;
	cols.resize(boxMesh.Vertexes());
	for (int i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
	UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for box: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::CylinderMeshExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    Mesh cylinderMesh = Mesh(Cylinder(1.0, 1.0, Vector(.0,0.0,0.0)), 64);

    std::vector<Color> cols;
    cols.resize(cylinderMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(cylinderMesh, cols, cylinderMesh.VertexIndexes());
    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for Cylinder: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::ConeMeshExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    Mesh coneMesh = Mesh(Cone(1.0,Vector(1.0, 0.0, 0.0), Vector(-1.0, 1.0, 0.0)), 64);

    std::vector<Color> cols;
    cols.resize(coneMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(coneMesh, cols, coneMesh.VertexIndexes());
    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for cone: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::SphereWarpExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    Mesh sphereMesh = Mesh(Sphere(5.0, Vector(0.0, 0.0, 0.0)), 64);
    sphereMesh.SphereWarp(Sphere(10.0, Vector(0.0,0.0,-7.0)), Vector(0.0, 0.0, 10.5));

    std::vector<Color> cols;
    cols.resize(sphereMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(sphereMesh, cols, sphereMesh.VertexIndexes());
    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for SphereWarp: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::MushroomExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    Mesh shroomMesh = Mesh(Sphere(5.0, Vector(0.0, 0.0, 0.0)), 64);
    shroomMesh.SphereWarp(Sphere(10.0, Vector(0.0,0.0,-7.0)), Vector(0.0, 0.0, 10.5));
    Mesh trunk = Mesh(Cylinder(1.0,10, Vector(-5.0,0.0,0.0)),64);
    trunk.Transform(Matrix::rotateY(M_PI/2));
    shroomMesh.Merge(trunk);
    Mesh sphere = Mesh(Sphere(5.0, Vector(0.0, 0.0, -10.0)), 64, 10);
    shroomMesh.Merge(sphere);
    Mesh sphere2 = Mesh(Sphere(10.0, Vector(0.0, 0.0, -17.0)), 64, 10);
    shroomMesh.Merge(sphere2);
    shroomMesh.SphereWarp(Sphere(2.0, Vector(5.0, 0.0, 0.0)), Vector(1.0,1.0,-1.0));
    Mesh sphere3 = Mesh(Sphere(2.0, Vector(0.0, 0.0, 0.0)), 64);
    sphere3.SphereWarp(Sphere(5.0, Vector(0.0,0.0, 0.0)), Vector(0.0, 0.0, 5.0));
    sphere3.Transform(Matrix::rotateX(M_PI));
    shroomMesh.Merge(sphere3);
    shroomMesh.Transform(Matrix::rotateY(M_PI/4));


    std::vector<Color> cols;
    cols.resize(shroomMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(shroomMesh, cols, shroomMesh.VertexIndexes());
    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for mushroom: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::CapsuleMeshExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    Mesh capsuleMesh = Mesh(Capsule(2.0, 10.0,Vector(0.0, 0.0, 0.0)), 64);

    std::vector<Color> cols;
    cols.resize(capsuleMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(capsuleMesh, cols, capsuleMesh.VertexIndexes());
    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for Capsule: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::SphereMeshExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    Mesh sphereMesh = Mesh(Sphere(1.0, Vector(0.0, 0.0, 0.0)), 64);

    std::vector<Color> cols;
    cols.resize(sphereMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(sphereMesh, cols, sphereMesh.VertexIndexes());
    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for Sphere: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::SphereImplicitExample()
{
    auto start = std::chrono::high_resolution_clock::now();

      AnalyticScalarField implicit;

      Mesh implicitMesh;
      implicit.Polygonize(31, implicitMesh, Box(2.0));

      std::vector<Color> cols;
      cols.resize(implicitMesh.Vertexes());
      for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(0.8, 0.8, 0.8);

      meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
      UpdateGeometry();

      auto end =std::chrono::high_resolution_clock::now();
      auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
      std::cout<<"Generation time for implicit Sphere: "<< duree.count() <<"ms"<<std::endl;

}

void MainWindow::TorusExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    Mesh torusMesh = Mesh(Torus(Vector(0.0, 0.0, 0.0), 1.0, 0.5), 32);

    std::vector<Color> cols;
    cols.resize(torusMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(torusMesh, cols, torusMesh.VertexIndexes());
    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for torus: "<< duree.count() <<"ms"<<std::endl;
}

void MainWindow::HeightFieldExample()
{
    auto start = std::chrono::high_resolution_clock::now();

    //Remplacer avec le path correspondant
    QString image = QString::fromStdString("D:/Master/Image/Tiny-Mesh/AppTinyMesh/Data/terrain2.png");

    HeightField field;
    field = HeightField(image, 10.0);
    field.Transform(Matrix::rotateZ(M_PI/2));
    meshColor = MeshColor(field, 10.0, true);


    UpdateGeometry();

    auto end =std::chrono::high_resolution_clock::now();
    auto duree = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<<"Generation time for HeightField: "<< duree.count() <<"ms"<<std::endl;

}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

	uiw.lineEdit->setText(QString::number(meshColor.Vertexes()));
	uiw.lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();

}

void MainWindow::UpdateMaterial()
{
	meshWidget->UseWireframeGlobal(uiw.wireframe->isChecked());

	if (uiw.radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
