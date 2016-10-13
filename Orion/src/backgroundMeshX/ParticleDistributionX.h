/*
 * ParticleDistributionX.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef PARTICLEDISTRIBUTIONX_H_
#define PARTICLEDISTRIBUTIONX_H_


#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "FEMElementX.h"
#include "GaussPointX.h"
#include "InputFileData.h"
#include "ParticleX.h"
#include "petscsys.h"
//#include "WindowFuncAsym.h"

class ParticleDistributionX {

public:

  ParticleDistributionX(InputFileData* InputData,std::ofstream& logFile) :
    particlesNum(0),maxPtcleSupport(0) {};

    ~ParticleDistributionX() {};

    // Determine for each process a portion of particles and set a root list,
    // which process each particle belongs to.
    void setPtcleProcDistrib(InputFileData* InputData,std::ofstream& logFile);

    /*******************************************************************/
    // Particles' stuff.

    // Determine the necessary influence radii for all particles.
    void setInfluenceRadii(InputFileData* InputData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile,
                           PetscViewer& viewerMPI,
                           PetscViewer& viewerSEQ);

    // determine particle influence radii correspondingly to the distances
    // to their neighbouring particles
    void setPtcleRadsDistDepend(InputFileData* InputData,
                                std::map<std::string,double>& modelData,
                                std::ofstream& logFile);

    // determine particle influence radii correspondingly to the distances
    // to their neighbouring particles
    void setPtcleRadsDistDepend2(InputFileData* InputData,
                                 std::map<std::string,double>& modelData,
                                 std::ofstream& logFile);

    // Determine the influence radius of all particles support dependent.
    void setPtcleRadsSupportDepend(InputFileData* InputData,
                                   std::vector<GaussPointX>& gPoints,
                                   std::vector<GaussPointX>& bGPoints,
                                   std::map<std::string,double>& modelData,
                                   std::ofstream& logFile);

    // Post-process all particles' influence radii.
    void postProcInfRads(InputFileData* InputData,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile);

    // Return number of particles.
    int getParticlesNum() { return particlesNum; };

    // Return greatest number of supporting particles of a single
    // particle.
    const int& getMaxPtcleSupport() { return maxPtcleSupport; };

    // Return the vector containing all particle data.
    std::vector<ParticleX>& getParticlesVec() { return particles; };

    // Return the particle root list
    intVector& getPtcleRootList() { return ptcleRootList; };

    // Return vector containing connectivity sorted particles
    // identifiers to the original read ones.
    intVector& getOldIdx() { return oldIdx; };

    // Return vector containing connectivity original read particle
    // identifiers to the new sorted ones.
    intVector& getNewIdx() { return newIdx; };

    // Return vector containing those particle which are exclusively
    // local
    intVector& getExclusiveLocalPtcls() { return exclusiveLocalPtcls; };
    intVector& getExtendedLocalGlobalPtcls()
      { return extendedLocalGlobalPtcls; };

    // Return start position of processor's exclusively local particle
    // portion within allExclLocalPtcls
    intVector& getAllExclLocalPtclsStartIdx()
      { return allExclLocalPtclsStartIdx; };

    // Return consecutive ordering of all exclusively local particles
    intVector& getAllExclLocalPtcls() { return allExclLocalPtcls; };

    // Return mapping of global particle IDs to the all processor's
    // consecutive exclusively local particle ordering
    intVector& getAllGlobalExclLocalPtcleIdx()
      { return allGlobalExclLocalPtcleIdx; };

    // Return vector containing all locally needed particles
    intVector& getLocalGlobalPtcls() { return localGlobalPtcls; };

    // Return vector containing connectivity of global and locally needed
    // particle identifiers
    intVector& getGlobalLocalPtcleIdx() { return globalLocalPtcleIdx; };

    // Update the particle coordinates.
    void updateParticleCoords(InputFileData* InputData,
                              dbVector& allDOF,
                              std::map<std::string,double>& modelData,
                              std::ofstream& logFile);

    void clearAllPtcleDOF(InputFileData* InputData,
                          std::map<std::string,double>& calcData,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile);

    /*******************************************************************/
    // Connectivity lists between particles and particles, and point and
    // particles.

    // Determine for all particles their supporting neighbour particles.
    void setPtclePtclsConn(InputFileData* InputData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Determine for given point their neighbour particles.
    void setPointPtclsConn(InputFileData* InputData,
                           dbVector& coords,
                           int& supportSize,
                           intVector& supportingPtcls,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Determine for given point their neighbour particles and the minimum
    // distance with the supporting particle.
    void setPointPtclsConn(InputFileData* InputData,
                           dbVector& coords,
                           int& supportSize,
                           intVector& supportingPtcls,
                           double& localMinPtcleDist,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    // Determine for a particle its neighbour particles making use of the
    // finite element mesh.
    void setElemPtcleSuppList(InputFileData* InputData,
			      int& alphaPtcle,
			      int& betaPtcle,
			      int& supportSize,
			      intVector& suppPtcls,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile);

    /********************************************************************/
    // Create a list of the exclusively local particles from all processors
    // in consecutive ordering and establish the mapping of global ptcle ID to
    // to this sequential ordering.
    void setAllExclLocalPtcleOrder(InputFileData* InputData,
				   std::ofstream& logFile);

    // Get the exclusively local particles from all processors in consecutive
    // ordering.
    void getAllLocalPtcleOrdering(InputFileData* InputData,
				  intVector& localStartIdx,
				  intVector& allLocalPtcls,
				  intVector& allGlobalLocalIdx,
				  std::ofstream& logFile);

    // Merge the local portions of entries stored at particles
    void mergeLocalPtcleVectors(InputFileData* InputData,
				dbVector& localVec,int ptcleEntries,
				dbVector& globalVec,
				std::map<std::string,double>& calcData,
				std::map<std::string,double>& modelData,
				std::ofstream& logFile);

    // Deallocate all class arrays of particles which are not locally needed.
    void deleteNonlocalPtcls(InputFileData* InputData,
			     std::ofstream& logFile);

    // Create a vector containing all particle coordinates.
    void getAllPtcleCoords(InputFileData* InputData,
			   dbVector& allCoords,
			   std::map<std::string,double>& calcData,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile);


 protected:

    /********************************************************************/
    // general particle's stuff
    int particlesNum;
    int maxPtcleSupport;
    std::vector<ParticleX> particles;

    /****************************/
    // parallel particle distribution

    intVector ptcleRootList;
    intVector globalLocalPtcleIdx;

    // exclusively local particles (no overlap)
    intVector exclusiveLocalPtcls;

    // localGlobalPtcls: connectivity of local and global particles
    //                  (including particles in domain overlapping zones!)
    // extendedLocalGlobalPtcls: supporting particles of localGlobalPtcls
    //                           not included in localGlobalPtcls
    intVector localGlobalPtcls;
    intVector extendedLocalGlobalPtcls;

    // ordering and connectivity of exclusively local particles
    //
    // allExclLocalPtcls: consecutive ordering of all exclusively local
    //                    particles
    // allGlobalExclLocalPtcleIdx: mapping of global particle IDs to the all
    //                             processor's consecutive exclusively local
    //                             particle ordering
    // localPtclsStartIdx: start position of processor's exclusively local
    //                     particle portion within allExclLocalPtcls

    intVector allExclLocalPtclsStartIdx;
    intVector allExclLocalPtcls;
    intVector allGlobalExclLocalPtcleIdx;

    /****************************/
    // Particle indexing.

    // New particle number i has original read number oldIdx[i].
    intVector oldIdx;

    // Original read particle number i has new number newIdx[i].
    intVector newIdx;

};


#endif /* PARTICLEDISTRIBUTIONX_H_ */
