#ifndef MaxEntDiscretising_h_
#define MaxEntDiscretising_h_
#include <iostream>
#include <vector>
#include "math.h"

#include "BackgroundMesh.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "InputFileData.h"
#include "MLSDiscretising.h"
#include "MaxEntShapeFunc.h"
#include "Particle.h"
#include "ParticleDistribution.h"
#include "mpi.h"


class MaxEntDiscretising : public virtual BackgroundMesh,
                           public virtual MaxEntShapeFunc,
                           public virtual MLSDiscretising {

public: 

    MaxEntDiscretising(InputFileData* InputData,
                       std::ofstream& logFile);
    ~MaxEntDiscretising() {};

    // Calculation of the RKPM interpolanten
    void calcInterpolants(InputFileData* InputData,
			  std::map<std::string,double>& calcData,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile,
                          PetscViewer& viewerMPI,
                          PetscViewer& viewerSEQ);

    // Calculate the beta values for each particle
    void setBetaOnPtcls(InputFileData* InputData,
                        std::map<std::string,double>& modelData,
                        std::ofstream& logFile);

    // Calculate for all Gauss points and all particles a vector containing
    // the calculated shape functions and their derivatives for all
    // supporting particles.
    void setAllShapeFuncs(InputFileData* InputData,
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

    // Calculate at a certain point for all its supporting particles their
    // shape functions.
    void calcShapeFuncs(InputFileData* InputData,
                        intVector& sPtcls,
                        std::vector<Particle>& particles,
                        double& x,double& y,double& z,
                        dbVector& shapeFuncs,
                        double& MinPtcleDist,
                        std::map<std::string,double>& modelData,
                        std::ofstream& logFile,
                        PetscViewer& viewerSEQ);


    // Calculate at a certain point for all its supporting particles their
    // shape functions and their first order derivations.
    void calcShapeFuncs(InputFileData* InputData,
                        intVector& sPtcls,
                        std::vector<Particle>& particles,
                        double& x,double& y,double& z,
                        dbVector& shapeFuncs,
                        dbMatrix& firstDerivShapes,
                        double& MinPtcleDist,
                        std::map<std::string,double>& modelData,
                        std::ofstream& logFile,
                        PetscViewer& viewerSEQ);

    // Calculate the shape functions and derivatives of the shape functions
    // using MLS
    void calcMlsInterpolants(InputFileData* InputData,
			     std::map<std::string,double>& calcData,
                             std::map<std::string,double>& modelData,
                             std::ofstream& logFile,
                             PetscViewer& viewerMPI,
                             PetscViewer& viewerSEQ);

    // Approximate the value of Beta using MLS shape functions
    void approxBetaValues(InputFileData* InputData,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile);

    // Approximate the derivative values of Beta using MLS shape functions
    void approxBetaDerivValues(InputFileData* InputData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);

    // Clear shape Functions and derivations of shape functions of
    // each particle
    void clearShapeFuncs(InputFileData* InputData,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile);


    // Debugging Purposes
    void debugLogFile(InputFileData* InputData,
                       std::map<std::string,double>& modelData,
                       std::ofstream& logFile);

    // Shape Function Test
    void shapeFuncTest(InputFileData* InputData,
                       std::map<std::string,double>& modelData,
                       std::ofstream& logFile,
                       PetscViewer& viewerSEQ);

};
#endif
