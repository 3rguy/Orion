#ifndef ModifiedMLSGauss_h_
#define ModifiedMLSGauss_h_

#include <float.h>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "BackgroundMesh.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "GaussPoint.h"
#include "InputFileData.h"
#include "MLSGaussIntegral.h"
#include "Particle.h"
#include "petscsys.h"
#include "petscvec.h"

class ModifiedMLSGauss : public virtual MLSGaussIntegral {

  public:

    ModifiedMLSGauss(InputFileData* InputData,std::ofstream& logFile);
    ~ModifiedMLSGauss() {};

    void clearClass();

    // Calculate for all Gauss points and all particles a vector containing 
    // the calculated shape functions and their derivatives for all 
    // supporting particles.
    void setAllShapeFuncs(InputFileData* InputData,
			  std::map<std::string,double>& calcData,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile,
			  PetscViewer& viewerMPI,
			  PetscViewer& viewerSEQ);

    // Create ghost particle on the deformation boundary (u=0/r=0 only!)
    void createBoundGhostPtcls(InputFileData* InputData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerMPI,
			       PetscViewer& viewerSEQ);

    // Determine for all Gauss points and particles their supporting 
    // boundary ghost particles.
    void setGhostPtcleConn(InputFileData* InputData,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile,
			   PetscViewer& viewerMPI,
			   PetscViewer& viewerSEQ);

    // Calculate for all gauss points a vector containing the calculated
    // shape functions and their derivations for all supporting particles.
    void setShapeFuncsOnGauss(InputFileData* InputData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,
			      PetscViewer& viewerMPI,
			      PetscViewer& viewerSEQ);

    // Calculate for all particles a vector containing the calculated
    // shape functions and their derivations for all supported particles.
    void setShapeFuncsOnPtcls(InputFileData* InputData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,
			      PetscViewer& viewerMPI,
			      PetscViewer& viewerSEQ);

    // Calculate for all boundary gauss points a vector containing the 
    // calculated shape functions for all supported particles.
    void setShapeFuncsOnBGauss(InputFileData* InputData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerMPI,
			       PetscViewer& viewerSEQ);
    

};

#endif
