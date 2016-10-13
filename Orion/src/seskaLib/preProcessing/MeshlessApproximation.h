#ifndef MeshlessApproximation_h_
#define MeshlessApproximation_h_

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "FEMDiscretising.h"
#include "GaussPoint.h"
#include "InputFileData.h"
#include "MLSDiscretising.h"
#include "MaxEntDiscretising.h"
#include "Particle.h"
#include "petscsys.h"

class MeshlessApproximation : public virtual FEMDiscretising,
                              public virtual MLSDiscretising,
                              public virtual MaxEntDiscretising {

  public:

    MeshlessApproximation(FEMGeometry* FEData,
			  InputFileData* InputData,
			  std::map<std::string,double>& calcData,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile,
			  PetscViewer& viewerMPI,
			  PetscViewer& viewerSEQ);

    ~MeshlessApproximation() {};

    /******************************************************************/
    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    void copyFEData(FEMGeometry* FEData,InputFileData* InputData,
		    std::map<std::string,double>& modelData,
		    std::ofstream& logFile);

    // Determine the necessary influence radii for all particles.
    void setInfluenceRadii(InputFileData* InputData,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile,
			   PetscViewer& viewerMPI,
			   PetscViewer& viewerSEQ);

    // Calculation of the interpolanten.
    void calcInterpolants(InputFileData* InputData,
			  std::map<std::string,double>& calcData,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile,
			  PetscViewer& viewerMPI,
			  PetscViewer& viewerSEQ);


    /******************************************************************/
    // Get particles DOF

    // Assemble all deformation degrees of freedom. 
    void getAllPtcleDefDOF(InputFileData* InputData,
			   dbVector& allDOF,
			   std::map<std::string,double>& calcData,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile);


};

#endif
