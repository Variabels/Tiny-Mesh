#ifndef __Qte__
#define __Qte__

#include "ui-interface.h"
#include "realtime.h"
#include "meshcolor.h"
#include "heightField.h"

class MainWindow : public QMainWindow
{
  Q_OBJECT
private:
  Ui::Assets uiw;				//!< Interface

  MeshWidget* meshWidget;		//!< Viewer
  MeshColor meshColor;		//!< Mesh.

public:
  MainWindow();
  ~MainWindow();
  void CreateActions();
  void UpdateGeometry();

public slots:
  void editingSceneLeft(const Ray&);
  void editingSceneRight(const Ray&);
  void BoxMeshExample();
  void CylinderMeshExample();
  void SphereMeshExample();
  void CapsuleMeshExample();
  void ConeMeshExample();
  void MushroomExample();
  void SphereImplicitExample();
  void SphereWarpExample();
  void HeightFieldExample();
  void TorusExample();
  void ResetCamera();
  void UpdateMaterial();
};

#endif
