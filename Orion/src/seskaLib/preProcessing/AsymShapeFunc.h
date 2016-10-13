// Calculate for a given point a shapefunction set where the influence zones
// are different in negative and positive coordinate direction.

#ifndef AsymShapeFunc_h_
#define AsymShapeFunc_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "BasisPolynom.h"
#include "BasisPolyRegular.h"
#include "BasisPolyAsym.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "PETScFunctions.h"
#include "petscksp.h"
#include "WindowFuncAsym.h"



namespace AsymShapeFunc {

  // Calculate at a gauss point for all its supporting particles their 
  // shape functions.
  void calcShapes(InputFileData* InputData,
		  int& supportSize,
		  intVector& sPtcls,
		  std::vector<Particle>& ptcls,
		  double& x,double& y,double& z,
		  dbVector& shapeFuncs,
		  int& basisTermNum,
		  std::map<std::string,double>& modelData,
		  std::ofstream& logFile,
		  PetscViewer& viewerSEQ);

  // Calculate at a gauss point for all its supporting particles their 
  // shape functions and their first order derivations.
  void calcShapes(InputFileData* InputData,
		  int& supportSize,
		  intVector& sPtcls,
		  std::vector<Particle>& ptcls,
		  double& x,double& y,double& z,
		  dbVector& shapeFuncs,
		  dbVector& xDerivShapes,
		  dbVector& yDerivShapes,
		  dbVector& zDerivShapes,
		  std::map<std::string,double>& modelData,
		  std::ofstream& logFile,
		  PetscViewer& viewerSEQ);

  // Calculate at a point for all its supporting particles their 
  // shape functions and their first and second order derivations.
  void calcShapes(InputFileData* InputData,
		  int& supportSize,
		  intVector& sPtcls,
		  std::vector<Particle>& ptcls,
		  double& x,double& y,double& z,
		  dbVector& shapeFuncs,
		  dbVector& xDerivShapes,
		  dbVector& yDerivShapes,
		  dbVector& zDerivShapes,
		  dbVector& xxDerivShapes,
		  dbVector& yyDerivShapes,
		  dbVector& zzDerivShapes,
		  dbVector& xyDerivShapes,
		  dbVector& yzDerivShapes,
		  dbVector& zxDerivShapes,
		  std::map<std::string,double>& modelData,
		  std::ofstream& logFile,
		  PetscViewer& viewerSEQ);

}
#endif
