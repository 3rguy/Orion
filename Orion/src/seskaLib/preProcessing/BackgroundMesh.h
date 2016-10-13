#ifndef BackgroundMesh_h_
#define BackgroundMesh_h_

#include <float.h>
#include <fstream>
#include <iostream>
#include <limits>
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
#include "ResultantReactionCondition.h"

class BackgroundMesh : public virtual ParticleDistribution {

public:

    BackgroundMesh(InputFileData* InputData,std::ofstream& logFile) :
      ParticleDistribution(InputData,logFile),
      globalGaussPtsNum(0),globalBGaussPtsNum(0),
      localGaussPtsNum(0),localBGaussPtsNum(0),
      maxGaussPtsPerVolumeElem(0),
      maxIndirectPtcleSupport(0),
      maxBGaussSupport(0),minBGaussSupport(0),
      maxGaussSupport(0),minGaussSupport(0)  {};

    ~BackgroundMesh() {};

    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    void copyFEData(FEMGeometry* FEData,InputFileData* InputData,
                    std::map<std::string,double>& modelData,
                    std::ofstream& logFile);
    
    void setPtcleRadsElemDepend(InputFileData* InputData,
                                dbVector& globalInfluenceRadii,
                                std::map<std::string,double>& modelData,
                                std::ofstream& logFile);
    
    void setPtcleRadsElemDepend2(InputFileData* InputData,
                                 std::map<std::string,double>& modelData,
                                 std::ofstream& logFile);

    // Return the maximal number of Gauss points per volume element.
    int getMaxGaussPtsPerVolumeElem();
    

    void clearClass();

    /******************************************************************/
    // Determine the necessary influence radii for all particles dependent
    // on the FEM background mesh.
    void setInfluenceRadii(InputFileData* InputData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile,
                           PetscViewer& viewerMPI,
                           PetscViewer& viewerSEQ);
    
    /*******************************************************************/
    // Connectivity lists between particles and gauss points, and
    // influencing spheres and influencing spheres.

    // Determine for all influence spheres their neighbour influence
    // spheres.
    void setInflSpheresConn(InputFileData* InputData,
                            std::ofstream& logFile);

    /*******************************************************************/
    // Connectivity lists between particles and particles, and point and
    // particles.

    // Determine for all particles their supporting neighbour particles.
    void setPtclePtclsConn(InputFileData* InputData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Determine for a particle its neighbour particles making use of the
    // finite element mesh.
    void setElemPtcleSupport(InputFileData* InputData,
                             int& betaPtcle,
                             blVector& includedPtcls,
                             blVector& includedElems,
                             int& supportSize,
                             intVector& suppPtcls,
                             dbVector& coords,
                             intVector& neededElemCounter,
                             intMatrix& neededElemIdx,
                             std::map<std::string,double>& modelData,
                             std::ofstream& logFile);

    // Determine for a point a portion of nonlocal supporting particles.
    bool getNonLocalSupport(InputFileData* InputData,
                            intVector& neededElemCounter,
                            intMatrix& neededElemIdx,
                            intVector& nonLocalNodes,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile);

    // Determine for a point a portion of nonlocal supporting particles.
    bool getNonLocalSupport(InputFileData* InputData,
			    int& alphaPtcle,
			    blVector& includedPtcls,
			    intVector& neededElemCounter,
			    intMatrix& neededElemIdx,
			    intVector& nonLocalNodes,
			    intVector& nonLocalBetaNodes,
			    std::map<std::string,double>& modelData,
			    std::ofstream& logFile);
    

    // Determine for each particle its included gauss points.
    void setPtcleGaussConn(InputFileData* InputData,
                           dbVector& allgCoords,
                           intVector& inclPtsStartPos,
                           intVector& inclPtsCounts,
                           intVector& inclPts,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Determine for each Gauss point its supporting particles.
    void setGaussPtcleConn(InputFileData* InputData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Determine for each boundary Gauss point its supporting particles.
    void setBGaussPtcleConn(InputFileData* InputData,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile);

    // Determine for a Gauss point its neighbour particles making use of the
    // finite element mesh.
    void setElemPtcleSupport(InputFileData* InputData,
                             int& betaPtcle,
                             blVector& includedPtcls,
                             int& supportSize,
                             intVector& suppPtcls,
                             dbVector& gaussCoords,
                             std::map<std::string,double>& modelData,
                             std::ofstream& logFile);

    // Rearrange Gauss point support lists according to the new particles
    // vector ordering
    void rearrangeSupportLists(InputFileData* InputData,
                               intVector& newGlobalPtcls,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);

    // Distribute the degrees of freedoms to the local gauss points.
    void assignGaussPtsDOF(InputFileData* InputData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);


    /*******************************************************************/
    // Partitioning stuff.

    // Determine the Gauss point distribution among the processors.
    void setGaussDistribution(InputFileData* InputData,
                              std::map<std::string,double>& modelData,
                              std::ofstream& logFile);

    // Determine for each gauss point its supporting particles and set up
    // a gauss point connectivity list used for partitioning.
    void setPartitionGraph(InputFileData* InputData,
                           intVector& gaussProcList,
                           dbVector& allGCoords,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Partition the gauss points to the processors.
    void partitionGaussPoints(InputFileData* InputData,
                              intVector& connList,intVector& weights,
                              intVector& connListStartPos,
                              intVector& globalGPtsRanges,
                              intVector& gaussProcList,
                              std::ofstream& logFile);

    // Exchange the volume Gauss data between the processor that each of
    // them possesses the correct local data.
    void exchangeGaussData(InputFileData* InputData,
                           intVector& gaussProcList,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Get the minimum distance between particles
    double getMinPtcleDistance(){ return minPtcleDistance;};

    // Determine the initial boundary Gauss point distribution.
    void getInitBGaussSplitting(InputFileData* InputData,
                                int& defaultPortion,
                                int& startIdx,int& endIdx,
                                std::map<std::string,double>& modelData,
                                std::ofstream& logFile);

    // Split the boundary Gauss points between the processor that each of
    // them possesses the correct local data.
    void splitBoundGaussPts(InputFileData* InputData,
                            intVector& gaussProcList,
                            intVector& globalLocalIdx,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile);


    // Determine the locally needed particles.
    void setLocalPtcls(InputFileData* InputData,
                       std::map<std::string,double>& modelData,
                       std::ofstream& logFile,
                       PetscViewer& viewerMPI);

    /******************************************************************/
    // FEM stuff

    intVector& getElementRootList() { return elementRootList; };
    intVector& getNewLocalElemIdx() { return newLocalElemIdx; };
    intMatrix& getElemGaussIdx() { return elemGaussIdx; };

    // Merge the local portions of element-Gauss point entries
    void mergeLocalElemGaussVectors(InputFileData* InputData,
                                    dbVector& localVec,
                                    int gaussEntries,
                                    dbMatrix3& globalMatrix,
                                    std::map<std::string,double>& calcData,
                                    std::map<std::string,double>& modelData,
                                    std::ofstream& logFile);


    /******************************************************************/
    // Gauss points' stuff.

    // Return number of local gauss points.
    int getLocalGaussNum() { return localGaussPtsNum; };

    // Return number of global gauss points.
    const int& getGlobalGaussNum() { return globalGaussPtsNum; };

    // Return the number of local boundary gauss points.
    const int& getLocalBGaussNum() { return localBGaussPtsNum; };

    // Return the number of global boundary gauss points.
    const int& getGlobalBGaussNum() { return globalBGaussPtsNum; };

    // Return greatest number of supporting particles of a single
    // gauss point.
    const int& getMaxGaussSupport() { return maxGaussSupport; };

    // Return greatest number of supporting particles of a single
    // boundary gauss point.
    const int& getMaxBoundGaussSupport() { return maxBGaussSupport; };



    // Create a vector containing all global gauss coordinates.
    dbVector getAllGaussCoords();

    // Create a matrix containing the penalty parameters of all boundary
    // Gauss points - local and otherwise.
    void getAllBoundPenaltyParameters(InputFileData* InputData,
                                      dbVector& allPenaltyParams,
                                      std::map<std::string,double>& modelData,
                                      std::ofstream& logFile);

    // Create a vector containing the history variables of all Gauss points.
    void getAllGPointHistoryVariables(InputFileData* InputData,
                                      dbVector& allHistoryVariables,
				      int& historySize,
                                      std::map<std::string,double>& modelData,
                                      std::ofstream& logFile);

    // Return vector containing all local gauss point data.
    std::vector<GaussPoint>& getGaussPointsVec() { return gaussPoints; };

    // Return vector containing all boundary gauss point data.
    std::vector<GaussPoint>& getBoundGaussPtsVec()
    { return boundGaussPoints; };
    
    // Return all boundary Gauss indices where a displacement boundary
    // condition is applied (point,line,surface).
    intMatrix& getAllDispBoundGaussPtsIdx(InputFileData* InputData,
                                          std::map<std::string,double>& modelData,
                                          std::ofstream& logFile);

    // Return all boundary Gauss indices where a rotation boundary
    // condition is applied (point,line,surface).
    intMatrix& getAllRotBoundGaussPtsIdx(InputFileData* InputData,
                                         std::map<std::string,double>& modelData,
                                         std::ofstream& logFile);

    // Return all boundary Gauss indices where an electric Dirichlet
    // boundary condition is applied (point,line,surface).
    intMatrix& getAllElectricBoundGaussPtsIdx(InputFileData* InputData,
                                              std::map<std::string,double>& modelData,
                                              std::ofstream& logFile);

    // Return all boundary Gauss indices where a TPM
    // boundary condition is applied (point,line,surface).
    intMatrix& getAllTPMBoundGaussPtsIdx(InputFileData* InputData,
                                              std::map<std::string,double>& modelData,
                                              std::ofstream& logFile);

    // Return all boundary Gauss indices where an depolarisation Dirichlet
    // boundary condition is applied (point,line,surface).
    intMatrix& getAllDepolarisationBoundGaussPtsIdx(InputFileData* InputData,
                                                    std::map<std::string,double>& modelData,
                                                    std::ofstream& logFile);

    // Return all boundary Gauss indices where a micro boundary
    // condition is applied (point,line,surface).
    intMatrix& getAllMicroBoundGaussPtsIdx(InputFileData* InputData,
                                           std::map<std::string,double>& modelData,
                                           std::ofstream& logFile);

    // Return all boundary Gauss indices where a deformation boundary
    // condition is applied (point,line,surface).
    intMatrix& getAllDefBoundGaussPtsIdx(InputFileData* InputData,
                                         std::map<std::string,double>& modelData,
                                         std::ofstream& logFile);

    // Return all boundary Gauss indices where a force boundary
    // condition is applied.
    intMatrix& getAllForceBoundGaussPtsIdx(InputFileData* InputData,
                                           std::map<std::string,double>& modelData,
                                           std::ofstream& logFile);

    // Return all surface integration Gauss point indices.
    intMatrix& getAllSurfaceGaussPtsIdx(std::ofstream& logFile);

    // Return all line integration Gauss point indices.
    intMatrix& getAllLineGaussPtsIdx(std::ofstream& logFile);

      
    // Return indices of boundary Gauss points with displacement boundary
    // conditions applied (point,line,surface).
    intMatrix& getSurfaceDispBoundGaussPtsIdx() {
        return surfaceDispBoundGaussPtsIdx; };

    intMatrix& getLineDispBoundGaussPtsIdx() {
        return lineDispBoundGaussPtsIdx; };

    // Return indices of boundary Gauss points with rotation boundary
    // conditions applied
    intMatrix& getSurfaceRotBoundGaussPtsIdx() {
        return surfaceRotBoundGaussPtsIdx; };

    intMatrix& getLineRotBoundGaussPtsIdx() {
        return lineRotBoundGaussPtsIdx; };

    // Return indices of boundary Gauss points with electric boundary
    // conditions applied
    intMatrix& getSurfaceElectricBoundGaussPtsIdx() {
        return surfaceElectricBoundGaussPtsIdx; };

    intMatrix& getLineElectricBoundGaussPtsIdx() {
        return lineElectricBoundGaussPtsIdx; };

    // TEMP: necessary?
    // Return indices of boundary Gauss points with depolarisation boundary
    // conditions applied
    intMatrix& getSurfaceDepolarisationBoundGaussPtsIdx() {
        return surfaceDepolarisationBoundGaussPtsIdx; };

    intMatrix& getLineDepolarisationBoundGaussPtsIdx() {
        return lineDepolarisationBoundGaussPtsIdx; };


    // Return indices of boundary Gauss points with a force loading
    // applied on
    intMatrix& getTractionBoundGaussPtsIdx() {
        return tractionBoundGaussPtsIdx; };

    intMatrix& getSurfacePressureBoundGaussPtsIdx() {
        return surfacePressureBoundGaussPtsIdx; };

    intMatrix& getLineForceBoundGaussPtsIdx() {
        return lineForceBoundGaussPtsIdx; };

    intMatrix& getPointForceBoundPtcleIdx() {
        return pointForceBoundPtcleIdx; };

    intMatrix& getElasticSurfaceForceBoundGaussPtsIdx() {
        return elasticSurfaceForceBoundGaussPtsIdx; };

    intMatrix& getElasticLineForceBoundGaussPtsIdx() {
        return elasticLineForceBoundGaussPtsIdx; };

    // Return indices of boundary Gauss points with a moment loading
    // applied on
    intMatrix& getSurfaceMomentBoundGaussPtsIdx() {
        return surfaceMomentBoundGaussPtsIdx; };

    intMatrix& getLineMomentBoundGaussPtsIdx() {
        return lineMomentBoundGaussPtsIdx; };

    intMatrix& getPointMomentBoundPtcleIdx() {
        return pointMomentBoundPtcleIdx; };

    // Return gauss points indices to that a volume load is applied.
    intMatrix& getBodyForceGaussPtsIdx()
    { return bodyForceGaussPtsIdx; };

    intMatrix& getBodyMomentGaussPtsIdx()
    { return bodyMomentGaussPtsIdx; };

    intMatrix& getBodyElectricChargeGaussPtsIdx()
    { return bodyElectricChargeGaussPtsIdx; };

    // Return gauss points indices which have electric charge loads
    // are applied.
    intMatrix& getSurfaceElectricChargeBoundGaussPtsIdx()
    { return surfaceElectricChargeBoundGaussPtsIdx; };

    // Return gauss points indices which have fluid volume flux
    // are applied.
    intMatrix& getFluidVolumeFluxBoundGaussPtsIdx()
    { return fluidVolumeFluxBoundGaussPtsIdx; };

    // Return the gauss points associated with cavity-volume-control
    // conditions
    intMatrix& getCavityVolumeControlBoundGaussPtsIdx()
    { return cavityVolumeControlBoundGaussPtsIdx; };

    // Return the gauss points associated with resultant reaction surfaces
    intMatrix& getResultantReactionBoundGaussPtsIdx()
    { return resultantReactionBoundGaussPtsIdx; }

    intVector& getGaussRootList()  { return gaussRootList; };
    intVector& getBoundGaussRootList() { return bGaussRootList; };

    /******************************************************************/
    // Get deformations on Gauss points

    // Calculate displacements and and rotations on all gauss points
    // from all processors.
    void getAllGaussDeforms(InputFileData* InputData,
                            dbVector& deformations,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile);

    // Calculate rotations of one integration point.
    void getIntPointRotation(InputFileData* InputData,
                             GaussPoint& gPoint,
                             dbVector& allDOF,
                             dbVector& rotation,
                             std::map<std::string,double>& modelData,
                             std::ofstream& logFile);

    /******************************************************************/
    // In order to build in the essential boundary conditions not only
    // the direct supporting particles have to take into account, but
    // also the indirect -> maxSphereSupport.

    // Return greatest number of direct and indirect supporting particles
    // of a single particle.
    const int& getMaxIndirectPtcleSupport()
    { return maxIndirectPtcleSupport; };

    // Set the greatest number of direct and indirect supporting particles
    // of a single particle.
    void setMaxIndirectPtcleSupport(int& size)
    { maxIndirectPtcleSupport = size; };

    /*******************************************************************/
    // Get all particle ID which need the rotation tensor as history
    // variable calculated.
    intVector& getRotTensPtcls() { return rotationTensorParticles; };

    /*******************************************************************/

    double minPtcleDistance;

    intVector elementRootList;
    intVector globalLocalElemIdx;
    std::vector<FEMElement> nodesElements;

    std::vector<GaussPoint> gaussPoints;
    std::vector<GaussPoint> boundGaussPoints;

    int globalBGaussPtsNum,localBGaussPtsNum;
    int globalGaussPtsNum,localGaussPtsNum;

    // background mesh data are used to calculuate influence radii yet.
    int maxGaussPtsPerVolumeElem;
    intVector newLocalElemIdx; // old global idx -> new local idx
    intMatrix elemGaussIdx; // new local elem idx -> all its local GPts

    // Indices of particles deformation boundary conditions are
    // applied.
    intMatrix pointDispBoundPtcleIdx;
    intMatrix pointRotBoundPtcleIdx;

    // Indices of local gauss points any kind of electric charge loads
    // are applied.
    intMatrix surfaceElectricChargeBoundGaussPtsIdx;
    intMatrix bodyElectricChargeGaussPtsIdx;

    // Indices of local gauss points which have a TPM
    // constraint applied
    intMatrix fluidVolumeFluxBoundGaussPtsIdx;
    intMatrix porePressureBoundGaussPtsIdx;

    // Indices of all boundary Gauss points any deformation boundary
    // conditions are applies
    intMatrix surfaceDispBoundGaussPtsIdx;
    intMatrix lineDispBoundGaussPtsIdx;

    intMatrix surfaceRotBoundGaussPtsIdx;
    intMatrix lineRotBoundGaussPtsIdx;

    intMatrix surfaceElectricBoundGaussPtsIdx;
    intMatrix lineElectricBoundGaussPtsIdx;

    intMatrix surfaceDepolarisationBoundGaussPtsIdx;
    intMatrix lineDepolarisationBoundGaussPtsIdx;

    intMatrix surfaceMicroBoundGaussPtsIdx;
    intMatrix lineMicroBoundGaussPtsIdx;

    // Indices of particles deformation boundary conditions are
    // applied.
    intMatrix pointElectricBoundPtcleIdx;
    intMatrix pointDepolarisationBoundPtcleIdx;

    // Indices of all boundary Gauss points any boundary loading is
    // applied
    intMatrix pointForceBoundPtcleIdx;
    intMatrix lineForceBoundGaussPtsIdx;
    intMatrix elasticLineForceBoundGaussPtsIdx;
    intMatrix tractionBoundGaussPtsIdx;
    intMatrix elasticSurfaceForceBoundGaussPtsIdx;
    intMatrix surfacePressureBoundGaussPtsIdx;

    // Indices of local gauss points any body loads are applied.
    intMatrix bodyForceGaussPtsIdx;
    intMatrix bodyMomentGaussPtsIdx;

    // Indices of local gauss points any body loads are applied.
    intMatrix surfaceMomentBoundGaussPtsIdx;

    int maxBGaussSupport,minBGaussSupport;
    int maxGaussSupport,minGaussSupport;

    /****************************/
    // index[0] = ptcleID, index[1] = weightID, index[2] = conditionID,
    // index[3] = surfaceID

    // Indices of local gauss points any body loads are applied.
    intMatrix pointMomentBoundPtcleIdx;
    intMatrix lineMomentBoundGaussPtsIdx;

    // indices of all gauss points associated with cavity-volume-control
    // conditions
    intMatrix cavityVolumeControlBoundGaussPtsIdx;

    // indices of all gauss points associated with resultant reaction surfaces
    intMatrix resultantReactionBoundGaussPtsIdx;


    // all Dirichlet boundary Gauss points indices (point,line,surface)
    intMatrix allDisplacementBoundGaussPtsIdx;
    intMatrix allRotationBoundGaussPtsIdx;
    intMatrix allElectricBoundGaussPtsIdx;
    intMatrix allDepolarisationBoundGaussPtsIdx;
    intMatrix allMicroBoundGaussPtsIdx;
    intMatrix allDeformationBoundGaussPtsIdx;

    // all loading boundary Gauss points indices (point,line,surface)
    intMatrix allForceBoundGaussPtsIdx;

    intMatrix allLineBoundGaussPtsIdx;
    intMatrix allSurfaceBoundGaussPtsIdx;

    intMatrix allTPMBoundGaussPtsIdx;


    /****************************/
    // root list for volume and boundary Gauss points.
    intVector gaussRootList;
    intVector bGaussRootList;

    /*******************************************************************/
    // background mesh data are used to calculuate influence radii yet.

    intMatrix globalUnsortedElemIdx;

    /****************************/
    // general particle's stuff
    int maxIndirectPtcleSupport;

    /****************************/
    intVector rotationTensorParticles;

};
#endif
