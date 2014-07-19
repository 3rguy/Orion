#ifndef MLSGaussIntegral_h_
#define MLSGaussIntegral_h_

#include <float.h>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "BackgroundMesh.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "FEMGeometry.h"
#include "GaussPoint.h"
#include "GaussPointSets.h"
#include "InputFileData.h"
#include "MLSShapeFuncSet.h"
#include "Particle.h"
#include "ParticleDistribution.h"
#include "petscsys.h"
#include "petscvec.h"
#include "WindowFunctionSets.h"

class MLSGaussIntegral : public virtual MLSShapeFuncSet,
			 public virtual ParticleDistribution,
			 public virtual BackgroundMesh {

  public:

    MLSGaussIntegral(InputFileData* InputData,std::ofstream& logFile);
    ~MLSGaussIntegral() {};


    /*******************************************************************/
    // Calculating shape functions.

    // Calculate for all Gauss points and all particles a vector containing 
    // the calculated shape functions and their derivatives for all 
    // supporting particles.
    void setAllShapeFuncs(InputFileData* InputData,
			  std::map<std::string,double>& calcData,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile,
			  PetscViewer& viewerMPI,
			  PetscViewer& viewerSEQ);

//    // Compute derivatives of the influence radius field.
//    void setInflRadiusDerivsOnGauss(InputFileData* InputData,
//				    std::map<std::string,double>& calcData,
//				    std::map<std::string,double>& modelData,
//				    std::ofstream& logFile,
//				    PetscViewer& viewerMPI,
//				    PetscViewer& viewerSEQ);

    // Calculate for all gauss points  a vector containing the calculated
    // shape functions and their derivations for all supporting particles.
    void setShapeFuncsOnGauss(InputFileData* InputData,
			      std::map<std::string,double>& calcData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,PetscViewer& viewerMPI,
			      PetscViewer& viewerSEQ);

    // Calculate for all particles a vector containing the calculated
    // shape functions of its supporting neighbour particles.
    void setShapeFuncsOnPtcls(InputFileData* InputData,
			      std::map<std::string,double>& calcData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,PetscViewer& viewerMPI,
			      PetscViewer& viewerSEQ);

    // Calculate for all boundary gauss points to which a surface load is 
    // applied a vector containing the calculated shape functions for all 
    // supporting particles.
    void setShapeFuncsOnBGauss(InputFileData* InputData,
			       std::map<std::string,double>& calcData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerMPI,
			       PetscViewer& viewerSEQ);

    // Check the shape function and their derivatives on all Gauss points.
    void checkGaussShapes(InputFileData* InputData,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile);

   
  private:

    //intVector elementRootList;
};

#endif
