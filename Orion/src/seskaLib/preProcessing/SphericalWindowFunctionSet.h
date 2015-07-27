// Calculate the window function ordinate of each particle supporting at a
// given point. The window functions have a spherical shape.

#ifndef SphericalWindowFunctionSet_h_
#define SphericalWindowFunctionSet_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"
#include "WindowFunctionSet.h"

class SphericalWindowFunctionSet : public virtual WindowFunctionSet {

  public:

    SphericalWindowFunctionSet() {
    }
    ;
    SphericalWindowFunctionSet(InputFileData* InputData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile) {
    }
    ;

    SphericalWindowFunctionSet(InputFileData* InputData,
                               std::vector<Particle>& ptcls,intVector& sPtcls,
                               dbVector& coords,int& supportSize,
                               unsigned int derivationOrder,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);

    SphericalWindowFunctionSet(InputFileData* InputData,
                               std::vector<Particle>& ptcls,intVector& sPtcls,
                               dbVector& coords,int& supportSize,
                               unsigned int derivationOrder,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile,PetscViewer& viewerSEQ) {
    }
    ;

    virtual
    ~SphericalWindowFunctionSet() {
    }
    ;

    // Calculate window functions.
    void
    calcWinFunctionSet(InputFileData* InputData,std::vector<Particle>& ptcls,
                       intVector& sPtcls,dbVector& coords,int& supportSize,
                       unsigned int derivationOrder,
                       std::map<std::string,double>& modelData,
                       std::ofstream& logFile);

    dbVector&
    getWindowFuncs() {
      return windowFuncs;
    }
    ;

    dbVector&
    getXDerivWinFuncs() {
      return firstDerivWinFuncs[0];
    }
    ;
    dbVector&
    getYDerivWinFuncs() {
      return firstDerivWinFuncs[1];
    }
    ;
    dbVector&
    getZDerivWinFuncs() {
      return firstDerivWinFuncs[2];
    }
    ;

    dbVector&
    getXXDerivWinFuncs() {
      return secondDerivWinFuncs[0];
    }
    ;
    dbVector&
    getYYDerivWinFuncs() {
      return secondDerivWinFuncs[1];
    }
    ;
    dbVector&
    getZZDerivWinFuncs() {
      return secondDerivWinFuncs[2];
    }
    ;
    dbVector&
    getXYDerivWinFuncs() {
      return secondDerivWinFuncs[3];
    }
    ;
    dbVector&
    getYZDerivWinFuncs() {
      return secondDerivWinFuncs[4];
    }
    ;
    dbVector&
    getZXDerivWinFuncs() {
      return secondDerivWinFuncs[5];
    }
    ;

  protected:

    dbVector windowFuncs;

    dbMatrix firstDerivWinFuncs;
    dbMatrix secondDerivWinFuncs;

  private:

    // Calculate the 1D cubic spline and its derivatives.
    double
    cubicSplineW1(double& rad);
    double
    cubicSplineW2(double& rad);
    double
    cubicSplineDW1(double& rad);
    double
    cubicSplineDW2(double& rad);
    double
    cubicSplineD2W1(double& rad);
    double
    cubicSplineD2W2(double& rad);

    // Calculate the 3D cubic window function and its derivatives.

    void
    calcCubicWinFunc(Particle& ptcle,dbVector& xnorm,double& wFunc,
                     std::map<std::string,double>& modelData,
                     std::ofstream& logFile);
    void
    calcCubicWinFunc1stDerivs(Particle& ptcle,dbVector& xnorm,dbVector& dWFuncs,
                              std::map<std::string,double>& modelData,
                              std::ofstream& logFile);
    void
    calcCubicWinFunc2ndDerivs(Particle& ptcle,dbVector& xnorm,
                              dbVector& d2WFuncs,
                              std::map<std::string,double>& modelData,
                              std::ofstream& logFile);

};

#endif
