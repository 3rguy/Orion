// Calculate for a given point a shapefunction set with a regular
// basis (EFG).

#ifndef EFGShapeFunc_h_
#define EFGShapeFunc_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "BasisPolyRegular.h"
#include "BasisPolynom.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "PETScFunctions.h"
#include "petscksp.h"
#include "WindowFunctionSet.h"


namespace EFGShapeFunc {

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

  	// Calculate at a node for all its supporting particles their
    // shape functions (used in POD calculation)
    void calcShapes(InputFileData* InputData,
  		  int& supportSize,
  		  intVector& sPtcls,
  		  std::vector<Particle>& ptcls,
  		  double& x,double& y,double& z,
  		  dbVector& shapeFuncs,
  		  int& basisTermNum,
  		  std::ofstream& logFile);

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
