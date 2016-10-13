#ifndef MLSShapeFuncSet_h_
#define MLSShapeFuncSet_h_

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "AsymShapeFunc.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "EFGShapeFunc.h"
#include "GaussPoint.h"
#include "InputFileData.h"
#include "OrthoShapeFunc.h"
#include "Particle.h"
#include "petscsys.h"
#include "RKPMShapeFunc.h"
#include "ShepardShapeFunc.h"


class MLSShapeFuncSet : public ShepardShapeFunc {

  public:

    MLSShapeFuncSet(InputFileData* InputData,std::ofstream& logFile) {};

    ~MLSShapeFuncSet() {};

    // Calculate at a certain point for all its supporting particles their 
    // shape functions.
    void calcShapeFuncs(InputFileData* InputData,int& supportSize,
			intVector& sPtcls,std::vector<Particle>& ptcls,
			double& x,double& y,double& z,
			dbVector& shapeFuncs,int& basisTermNum,
			std::map<std::string,double>& modelData,
			std::ofstream& logFile,PetscViewer& viewerSEQ);

    // Calculate at a certain point for all its supporting particles their 
    // shape functions and their first order derivations.
    void calcShapeFuncs(InputFileData* InputData,int& supportSize,
			intVector& sPtcls,std::vector<Particle>& ptcls,
			double& x,double& y,double& z,
			dbVector& shapeFuncs,dbMatrix& firstDerivShapes,
			std::map<std::string,double>& modelData,
			std::ofstream& logFile,PetscViewer& viewerSEQ);

    // Calculate at a certain point for all its supporting particles their 
    // shape functions and their first order derivations.
    void calcShapeFuncs(InputFileData* InputData,
			int& supportSize,
			intVector& sPtcls,
			std::vector<Particle>& particles,
			double& x,double& y,double& z,
			dbMatrix& radiusDerivs,
			dbVector& shapeFuncs,
			dbMatrix& firstDerivShapes,
			std::map<std::string,double>& modelData,
			std::ofstream& logFile,
			PetscViewer& viewerSEQ);

    // Calculate at a certain point for all its supporting particles their 
    // shape functions and their first and second order derivations.
    void calcShapeFuncs(InputFileData* InputData,int& supportSize,
			intVector& sPtcls,std::vector<Particle>& ptcls,
			double& x,double& y,double& z,
			dbVector& shapeFuncs,dbMatrix& firstDerivShapes,
			dbMatrix& secondDerivShapes,
			std::map<std::string,double>& modelData,
			std::ofstream& logFile,PetscViewer& viewerSEQ);


};

#endif
