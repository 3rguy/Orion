#ifndef MaxEntShapeFunc_h_
#define MaxEntShapeFunc_h_
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "math.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "InputFileData.h"
#include "Particle.h"
#include "ParticleDistribution.h"
#include "mpi.h"


class MaxEntShapeFunc {

public:

    MaxEntShapeFunc(InputFileData* InputData,std::ofstream& logFile){};
    ~MaxEntShapeFunc() {};

    // Calculate at a point for all its supporting particles their
    // shape functions.
    void calcShapes(InputFileData* InputData,
                    intVector& sPtcls,
                    std::vector<Particle>& Ptcls,
                    double& x, double& y,double& z,
                    dbVector& shapeFuncs,
                    std::map<std::string,double>& modelData,
                    std::ofstream& logFile,
                    PetscViewer& viewerSEQ);
    
    // For POD Calculation
    void calcShapes(InputFileData* InputData,
                        intVector& sPtcls,
                        std::vector<Particle>& Ptcls,
                        double& x, double& y,double& z,
                        dbVector& shapeFuncs,
                        std::ofstream& logFile);

    // Calculate at a point for all its supporting particles their
    // shape functions.
    void calcShapes(InputFileData* InputData,
                    intVector& sPtcls,
                    std::vector<Particle>& Ptcls,
                    double& x, double& y,double& z,
                    dbVector& shapeFuncs,
                    dbMatrix& firstDerivShapes,
                    std::map<std::string,double>& modelData,
                    std::ofstream& logFile,
                    PetscViewer& viewerSEQ);
    
    // Compute the Shape Function values of one element for
    // one Newton Iteration
    dbVector gamma_(InputFileData* InputData,
                    dbVector& mainPtcleCoord,
                    std::vector<Particle>&  Ptcls,
                    intVector& sPtcls,
                    int& sPtclsNum, int& usedDims,
                    dbMatrix& ngtv_inv_J,
                    std::ofstream& logFile);
    
    // Calculate the Partition function Z
    double CalPartitionFunc (InputFileData* InputData,
                             dbVector& mainPtcleCoord,
                             std::vector<Particle>& Ptcls,
                             intVector& sPtcls,
                             int& sPtclsNum, dbVector& lam,
                             dbVector& temp,
                             int& usedDims,
                             std::ofstream& logFile);
    
    // Calculate the First Derivative of Function
    // for Newton Method
    dbVector CalNewtonFirstDeriv(dbVector& mainPtcleCoord,
                                 std::vector<Particle>& Ptcls,
                                 intVector& sPtcls, int& sPtclsNum,
                                 dbVector& ShapeFuncs,
                                 int& usedDims,
                                 std::ofstream& logFile);
    
    // Calculate the Second Derivative of Function
    // for Newton Method
    dbMatrix CalNewtonSecDeriv(dbVector& mainPtcleCoord,
                               std::vector<Particle>& Ptcls,
                               intVector& sPtcls,
                               int& sPtclsNum,
                               dbVector& ShapeFuncs,
                               dbVector dgam, int& usedDims,
                               std::ofstream& logFile);

    // Calculate the First Derivative of Shape Function
    // Values
    dbMatrix CalFirstDerivFunc(dbVector& mainPtcleCoord,
                               std::vector<Particle>& Ptcls,
                               intVector& sPtcls,
                               int& sPtclsNum,
                               dbVector& shapeFuncs,
                               dbMatrix& ngtv_inv_J,
                               int& usedDims,
                               std::ofstream& logFile);

    // Calculate the First Derivatives of the shape functions
    // (including delta Beta)
    dbMatrix CalFirstDerivFunc_Beta(dbVector& mainPtcleCoord,
                                    std::vector<Particle>& Ptcls,
                                    intVector& sPtcls,int& sPtclsNum,
                                    dbVector& shapeFuncs,dbMatrix& ngtv_inv_J,
                                    int& usedDims,std::ofstream& logFile);

    // Calculate the First Derivatives of the shape functions using 2nd Method
    dbMatrix CalFirstDerivFunc_A(dbVector& mainPtcleCoord,
                                 std::vector<Particle>& Ptcls,
                                 intVector& sPtcls,int& sPtclsNum,
                                 dbVector& shapeFuncs,dbMatrix& ngtv_inv_J,
                                 int& usedDims,std::ofstream& logFile);

};
#endif
