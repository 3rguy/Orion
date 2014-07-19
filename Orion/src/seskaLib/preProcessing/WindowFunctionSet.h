// Calculate the window function of each particle supporting the 
// requested point.

#ifndef WindowFunctionSet_h_
#define WindowFunctionSet_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"
#include "WindowFunctionSets.h"

class WindowFunctionSet {

  public:

    WindowFunctionSet() {};
    WindowFunctionSet(InputFileData* InputData,
		      std::map<std::string,double>& modelData,
		      std::ofstream& logFile) {};

    WindowFunctionSet(InputFileData* InputData,
		      std::vector<Particle>& ptcls,
		      intVector& sPtcls,double& x,double& y,
		      double& z,int& supportSize,
		      unsigned int derivationOrder,
		      std::map<std::string,double>& modelData, 
		      std::ofstream& logFile);

    // For POD Calculation
    WindowFunctionSet(InputFileData* InputData,
    		      std::vector<Particle>& ptcls,
    		      intVector& sPtcls,double& x,double& y,
    		      double& z,int& supportSize,
    		      unsigned int derivationOrder,
    		      std::ofstream& logFile);

    WindowFunctionSet(InputFileData* InputData,
		      std::vector<Particle>& ptcls,
		      intVector& sPtcls,double& x,double& y,
		      double& z,int& supportSize,
		      unsigned int derivationOrder,
		      std::map<std::string,double>& modelData, 
		      std::ofstream& logFile,
		      PetscViewer& viewerSEQ) {};

    virtual ~WindowFunctionSet() {};

    // Calculate window functions.
    void calcWinFunctions(InputFileData* InputData,
			 std::vector<Particle>& ptcls,
			 intVector& sPtcls,double& x,double& y,
			 double& z,int& supportSize,
			 std::ofstream& logFile);

    // Calculation of first order derivation of the window functions 
    void calcWinFunction1stDerivs(InputFileData* InputData,
				  std::vector<Particle>& ptcls,
				  intVector& sPtcls,
				  double& x,double& y,
				  double& z,
				  int& supportSize,
				  std::ofstream& logFile);

    // Calculation of second order derivation of the window functions 
    void calcWinFunction2ndDerivs(InputFileData* InputData,
				  std::vector<Particle>& ptcls,
				  intVector& sPtcls,
				  double& x,double& y,
				  double& z,
				  int& supportSize,
				  std::ofstream& logFile);
    
    dbVector& getWindowFuncs() { return windowFuncs; };

    dbVector& getXDerivWinFuncs() { return xDerivWinFuncs; };
    dbVector& getYDerivWinFuncs() { return yDerivWinFuncs; };
    dbVector& getZDerivWinFuncs() { return zDerivWinFuncs; };

    dbVector& getXXDerivWinFuncs() { return xxDerivWinFuncs; };
    dbVector& getYYDerivWinFuncs() { return yyDerivWinFuncs; };
    dbVector& getZZDerivWinFuncs() { return zzDerivWinFuncs; };

    dbVector& getXYDerivWinFuncs() { return xyDerivWinFuncs; };
    dbVector& getYZDerivWinFuncs() { return yzDerivWinFuncs; };
    dbVector& getZXDerivWinFuncs() { return zxDerivWinFuncs; };

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
