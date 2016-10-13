// Calculate the window function ordinate of each particle supporting at a
// given point. The window functions have a prismatic shape.

#ifndef PrismaticWindowFunctionSet_h_
#define PrismaticWindowFunctionSet_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"
#include "PrismaticWindowFunction.h"
#include "WindowFunctionSet.h"

class PrismaticWindowFunctionSet : public virtual WindowFunctionSet {

  public:

    PrismaticWindowFunctionSet() {
    }
    ;
    PrismaticWindowFunctionSet(InputFileData* InputData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile) {
    }
    ;

    PrismaticWindowFunctionSet(InputFileData* InputData,
                               std::vector<Particle>& ptcls,intVector& sPtcls,
                               double& x,double& y,double& z,int& supportSize,
                               unsigned int derivationOrder,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);

    PrismaticWindowFunctionSet(InputFileData* InputData,
                               std::vector<Particle>& ptcls,intVector& sPtcls,
                               double& x,double& y,double& z,int& supportSize,
                               unsigned int derivationOrder,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile,PetscViewer& viewerSEQ) {
    }
    ;

    virtual
    ~PrismaticWindowFunctionSet() {
    }
    ;

    // Calculate window functions.
    void
    calcWinFunctions(InputFileData* InputData,std::vector<Particle>& ptcls,
                     intVector& sPtcls,double& x,double& y,double& z,
                     int& supportSize,std::ofstream& logFile);

    // Calculation of first order derivation of the window functions 
    void
    calcWinFunction1stDerivs(InputFileData* InputData,
                             std::vector<Particle>& ptcls,intVector& sPtcls,
                             double& x,double& y,double& z,int& supportSize,
                             std::ofstream& logFile);

    // Calculation of second order derivation of the window functions 
    void
    calcWinFunction2ndDerivs(InputFileData* InputData,
                             std::vector<Particle>& ptcls,intVector& sPtcls,
                             double& x,double& y,double& z,int& supportSize,
                             std::ofstream& logFile);

    dbVector&
    getWindowFuncs() {
      return windowFuncs;
    }
    ;

    dbVector&
    getXDerivWinFuncs() {
      return xDerivWinFuncs;
    }
    ;
    dbVector&
    getYDerivWinFuncs() {
      return yDerivWinFuncs;
    }
    ;
    dbVector&
    getZDerivWinFuncs() {
      return zDerivWinFuncs;
    }
    ;

    dbVector&
    getXXDerivWinFuncs() {
      return xxDerivWinFuncs;
    }
    ;
    dbVector&
    getYYDerivWinFuncs() {
      return yyDerivWinFuncs;
    }
    ;
    dbVector&
    getZZDerivWinFuncs() {
      return zzDerivWinFuncs;
    }
    ;

    dbVector&
    getXYDerivWinFuncs() {
      return xyDerivWinFuncs;
    }
    ;
    dbVector&
    getYZDerivWinFuncs() {
      return yzDerivWinFuncs;
    }
    ;
    dbVector&
    getZXDerivWinFuncs() {
      return zxDerivWinFuncs;
    }
    ;

  protected:

    dbVector windowFuncs;

    dbVector xDerivWinFuncs;
    dbVector yDerivWinFuncs;
    dbVector zDerivWinFuncs;

    dbVector xxDerivWinFuncs;
    dbVector yyDerivWinFuncs;
    dbVector zzDerivWinFuncs;

    dbVector xyDerivWinFuncs;
    dbVector yzDerivWinFuncs;
    dbVector zxDerivWinFuncs;

};

#endif
