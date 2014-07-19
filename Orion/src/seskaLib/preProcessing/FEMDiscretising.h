#ifndef FEMDiscretising_h_
#define FEMDiscretising_h_
#include <iostream>
#include <vector>
#include "math.h"

#include "BackgroundMesh.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "FEMElement.h"
#include "FEMGeometry.h"
#include "InputFileData.h"
#include "Particle.h"
#include "ParticleDistribution.h"
#include "mpi.h"


class FEMDiscretising : public virtual BackgroundMesh,
                        public virtual ParticleDistribution {

  public: 

    FEMDiscretising(InputFileData* InputData,
		    std::ofstream& logFile);
    ~FEMDiscretising() {};
    
    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    void copyFEData(FEMGeometry* FEData,
		    InputFileData* InputData,
		    std::map<std::string,double>& modelData,
		    std::ofstream& logFile);

    // Calculation of the shape functions and their derivatives with
    // respect to global coordinates.
    void calcInterpolants(InputFileData* InputData,
			  std::map<std::string,double>& calcData,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile,
                          PetscViewer& viewerMPI,
                          PetscViewer& viewerSEQ);

    // Determine for all nodes their connected neighbour nodes.
    void setPtclePtclsConn(InputFileData* InputData,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile);
    
    // Determine for each Gauss point its supporting nodes.
    void setGaussPtcleConn(InputFileData* InputData,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile);
    
    // Determine for each boundary gauss point its supporting nodes.
    void setBGaussPtcleConn(InputFileData* InputData,
			    std::map<std::string,double>& modelData,
			    std::ofstream& logFile);
    
    // Determine the Gauss point distribution among the processors.
    void setGaussDistribution(InputFileData* InputData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile);
    
    // Calculate for all Gauss points and all nodes a vector containing 
    // the calculated shape functions and their derivatives for all 
    // supporting nodes.
    void setAllShapeFuncs(InputFileData* InputData,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile,
			  PetscViewer& viewerMPI,
			  PetscViewer& viewerSEQ);
    
    // Calculate for all volume Gauss points a vector containing the 
    // shape functions and their derivations for their respective supporting
    // nodes.
    void setShapeFuncsOnGauss(InputFileData* InputData,
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
    
    // Calculate for all nodes a vector containing the calculated
    // shape functions and their derivations for all supported nodes.
    void setShapeFuncsOnPtcls(InputFileData* InputData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,
			      PetscViewer& viewerMPI,
			      PetscViewer& viewerSEQ);
    
  protected:

    std::vector<FEMElement> surfaceNodesElems;
    std::vector<FEMElement> lineNodesElems;

    
};
#endif
