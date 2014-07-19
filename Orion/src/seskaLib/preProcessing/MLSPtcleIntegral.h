#ifndef MLSPtcleIntegral_h_
#define MLSPtcleIntegral_h_

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "FEMGeometry.h"
#include "InputFileData.h"
#include "MLSShapeFuncSet.h"
#include "Particle.h"
#include "ParticleDistribution.h"
#include "petscsys.h"
#include "petscvec.h"


class MLSPtcleIntegral : public virtual MLSShapeFuncSet,
			 public virtual ParticleDistribution {

  public:

    MLSPtcleIntegral(InputFileData* InputData,std::ofstream& logFile);
    ~MLSPtcleIntegral() {};

    /********************************************************************/
    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    void copyFEData(FEMGeometry* FEData,InputFileData* InputData,
		    std::map<std::string,double>& modelData,
		    std::ofstream& logFile);

    /*******************************************************************/
    // Partitioning stuff.

    // Determine for each particle its supporting particles and set up
    // a particle connectivity list used for partitioning.
    //    void setPartitionGraph(InputFileData* InputData,std::ofstream& logFile);

    // Partition the particles to the processors.
    //    void particlePartitioning(InputFileData* InputData,
    //			      int* connList,int* connListStartPos,
    //			      int* globalGPtsRanges,
    //			      int* localGaussPtList,
    //			      std::ofstream& logFile);

    // Exchange the particle data between the processor that each of them 
    // the correct local data.
    //    void exchangeParticleData(InputFileData* InputData,int* gaussProcList,
    //			      std::ofstream& logFile);

    /*******************************************************************/

    // Rearrange Gauss point support lists according to the new particles
    // vector ordering
    void rearrangeSupportLists(InputFileData* InputData,
			       intVector& newGlobalPtcls,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile);

    /*******************************************************************/
    // Caluculating shape functions.

    // Calculate for all particles a vector containing the calculated
    // shape functions of its supporting neighbour particles.
    void setShapeFuncsOnPtcls(InputFileData* InputData,
			      std::map<std::string,double>& calcData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,PetscViewer& viewerMPI,
			      PetscViewer& viewerSEQ);

    /*******************************************************************/
      // Return the indices of all particles any body loading is applied.
    intMatrix& getBodyForcePtcleIdx() { return bodyForcePtcleIdx; };

    // Return the indices of all particles any boundary loading is applied
    intMatrix& getTractionBoundPtcleIdx() { 
      return tractionBoundPtcleIdx; };

    intMatrix& getLineForceBoundPtcleIdx() { 
      return lineForceBoundPtcleIdx; };

    intMatrix& getSurfacePressureBoundPtcleIdx() { 
      return surfacePressureBoundPtcleIdx; };

    intMatrix& getPointForceBoundPtcleIdx() { 
      return pointForceBoundPtcleIdx; };

    // Indices of particles deformation boundary conditions are 
    // applied.
    intMatrix& getPointDispBoundPtcleIdx() { 
      return pointDispBoundPtcleIdx; };

    intMatrix& getSurfaceDispBoundPtcleIdx() { 
      return surfaceDispBoundPtcleIdx; };

    // actually indices of all particles with Dirichlet conditions
   // (point,line,surface)
    intMatrix& getAllDisplacementBoundPtcleIdx(InputFileData* InputData,
					       std::map<std::string,double>& modelData,
					       std::ofstream& logFile);

    // actually indices of all particles with v. Neumann conditions 
    // (point,line,surface)
    intMatrix& getAllForceBoundPtcleIdx(InputFileData* InputData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile);

    /*******************************************************************/
    // index[0] = ptcleID, index[1] = weightID, index[2] = conditionID,
    // index[3] = surfaceID

    // Indices of particles body forces are applied.
    intMatrix bodyForcePtcleIdx;

    // Indices of all particles any boundary loading is applied
    intMatrix tractionBoundPtcleIdx;
    intMatrix lineForceBoundPtcleIdx;
    intMatrix surfacePressureBoundPtcleIdx;
    intMatrix pointForceBoundPtcleIdx;

    // Indices of particles deformation boundary conditions are 
    // applied.
    intMatrix pointDispBoundPtcleIdx;
    intMatrix lineDispBoundPtcleIdx;
    intMatrix surfaceDispBoundPtcleIdx;

    // indices of all particles with Dirichlet conditions 
    // (point,line,surface)
    intMatrix allDisplacementBoundPtcleIdx;

    // indices of all particles with v. Neumann conditions 
    // (point,line,surface)
    intMatrix allForceBoundPtcleIdx;
    
};

#endif
