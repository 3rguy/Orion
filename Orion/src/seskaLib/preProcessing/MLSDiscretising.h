#ifndef MLSDiscretising_h_
#define MLSDiscretising_h_

#include <float.h>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "FEMGeometry.h"
#include "InputFileData.h"
#include "MLSGaussIntegral.h"
#include "MLSPtcleIntegral.h"
#include "Particle.h"
#include "PETScFunctions.h"
#include "petscsys.h"
#include "petscvec.h"
#include "Point.h"


class MLSDiscretising : public virtual MLSGaussIntegral,
    public virtual MLSPtcleIntegral {

  public:

    MLSDiscretising(InputFileData* InputData,
		     std::ofstream& logFile);
    ~MLSDiscretising() {};

    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    void copyFEData(FEMGeometry* FEData,InputFileData* InputData,
                    std::map<std::string,double>& modelData,
                    std::ofstream& logFile);

    // Determine the necessary influence radii for all particles.
    void setInfluenceRadii(InputFileData* InputData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile,PetscViewer& viewerMPI,
                           PetscViewer& viewerSEQ);

    // Calculation of the RKPM interpolanten 
    void calcInterpolants(InputFileData* InputData,
                          std::map<std::string,double>& calcData,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile,PetscViewer& viewerMPI,
                          PetscViewer& viewerSEQ);

    // Calculation of the RKPM interpolanten with Gaussian quadrature.
    void calcGaussIntShapes(InputFileData* InputData,
                            std::map<std::string,double>& calcData,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile,PetscViewer& viewerMPI,
                            PetscViewer& viewerSEQ);

    // Calculation of the RKPM interpolanten with particle integration.
    void calcPtcleIntShapes(InputFileData* InputData,
                            std::map<std::string,double>& calcData,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile,PetscViewer& viewerMPI,
                            PetscViewer& viewerSEQ);

    // Rearrange Particle and Gauss point support lists according to the 
    // new particles vector ordering
    void rearrangeSupportLists(InputFileData* InputData,
                               intVector& newGlobalPtcls,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);

    // Compute spline coefficients and norming factor for all particles
    void setAllCustomPtcleSplines(InputFileData* InputData,
                                  std::map<std::string,double>& modelData,
                                  std::ofstream& logFile,PetscViewer& viewerMPI,
                                  PetscViewer& viewerSEQ);

    /*******************************************************************/
    // Incorporating of essential and natural boundary conditions.
    // Modify the disrete equation set's line and column identifiers being 
    // conform the equation set entries including essential and natural 
    // boundary conditions.
    void modifyEQIDs(InputFileData* InputData,std::ofstream& logFile,
                     PetscViewer& viewerMPI);

    // Set a special matrix containing the shape function values (for r=0)
    // of all particles essential boundary conditions applied to transform 
    // the discrete equation system to satisfy the essential boundary 
    // conditions.
    void setTransformingShapes(InputFileData* InputData,
                               intMatrix& boundaryPtcls,
                               intMatrix& localBoundPtcls,intVector& newDOFID,
                               std::map<std::string,double>& calcData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile,PetscViewer& viewerMPI,
                               PetscViewer& viewerSEQ);

    /******************************************************************/
    // Gather the degrees of freedom
    // Assemble all degrees of freedom. 
    void getAllPtcleDOF(InputFileData* InputData,dbVector& allDOF,
                        std::map<std::string,double>& calcData,
                        std::map<std::string,double>& modelData,
                        std::ofstream& logFile);

    // Print all degrees of freedom. 
    void printAllPtcleDOF(InputFileData* InputData,
                          std::map<std::string,double>& calcData,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile);

    // Assemble all step degrees of freedom. 
    void getAllPtcleStepDOF(InputFileData* InputData,dbVector& allStepDOF,
                            std::map<std::string,double>& calcData,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile);

    //    // Assemble all step degrees of freedom. 
    //    void getAllPtcleStepDOF(InputFileData* InputData,
    //			    dbVector& allStepDOF,
    //			    std::map<std::string,double>& calcData,
    //			    std::map<std::string,double>& modelData,
    //			    std::ofstream& logFile);

    // Assemble all deformation degrees of freedom. 
    void getAllPtcleDefDOF(InputFileData* InputData,dbVector& allDOF,
                           std::map<std::string,double>& calcData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Assemble all displacement degrees of freedom. 
    void getAllPtcleDispDOF(InputFileData* InputData,dbVector& allDOF,
                            std::map<std::string,double>& calcData,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile);

    /*******************************************************************/
    // Return the indices of all particles any loading is applied
    intMatrix& getPointForceBoundPtcleIdx(
        InputFileData* InputData,std::map<std::string,double>& modelData,
        std::ofstream& logFile);

    /*******************************************************************/
    // Return greatest number of supporting particles of a single
    // boundary particle (additional some of its boundary neighours).
    const int& getMaxBoundPtcleSupport() {
      return maxBoundPtcleSupport;
    };

    // at a specified distribution of points compute the MLS approximation 
    // of a given field of matrix values known at a distribution of particles
    void matrixFieldMLSApproximation(std::vector<Particle>& ptcls,
                                     std::vector<Point>& points,
                                     std::map<std::string,double>& data,
                                     std::ofstream& logFile,
                                     PetscViewer& viewerSEQ);
    
  private:

    // general particle's stuff
    int maxBoundPtcleSupport;


};

#endif
