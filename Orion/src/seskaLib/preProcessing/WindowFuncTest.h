// Calculate the window function of each particle supporting the 
// requested point.

#ifndef WindowFuncTest_h_
#define WindowFuncTest_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "defs.h"
#include "ElementTemplate.h"
#include "GaussPoint.h"
#include "GaussPointSets.h"
#include "InputFileData.h"
#include "LineElementTemplates.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"
#include "SurfaceElementTemplates.h"
#include "VolumeElementTemplates.h"
#include "WindowFunctionSets.h"



class WindowFuncTest {

  public:

    WindowFuncTest(InputFileData* InputData,
		   std::vector<Particle>& ptcls,
		   intVector& sPtcls,double& x,double& y,
		   double& z,int& supportSize,
		   unsigned int derivationOrder,
		   std::map<std::string,double>& modelData, 
		   std::ofstream& logFile);


    ~WindowFuncTest() {};


    // set-up brick-shaped Gauss quadrature cell
    double calcWeightFunctionIntegral(InputFileData* InputData,
				      dbVector& radii,
				      int& splineOrder,
				      std::map<std::string,double>& modelData, 
				      std::ofstream& logFile);

    // set-up brick-shaped Gauss quadrature cell
    void setBrickIntegrationCell(InputFileData* InputData,
				 dbVector& radii,
				 int& order,
				 std::vector<GaussPoint>& integrationPoints,
				 std::map<std::string,double>& modelData,
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

  private:
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

    // Calculate window functions.
    void calcWinFunctions(InputFileData* InputData,
			  std::vector<Particle>& ptcls,
			  intVector& sPtcls,double& x,double& y,
			  double& z,int& supportSize,
			  std::map<std::string,double>& modelData, 
			  std::ofstream& logFile);

    // Calculation of first order derivation of the window functions 
    void calcWinFunction1stDerivs(InputFileData* InputData,
				  std::vector<Particle>& ptcls,
				  intVector& sPtcls,
				  double& x,double& y,double& z,
				  int& supportSize,
				  std::map<std::string,double>& modelData, 
				  std::ofstream& logFile);

    // Calculation of second order derivation of the window functions 
    void calcWinFunction2ndDerivs(InputFileData* InputData,
				  std::vector<Particle>& ptcls,
				  intVector& sPtcls,
				  double& x,double& y,double& z,
				  int& supportSize,
				  std::map<std::string,double>& modelData, 
				  std::ofstream& logFile);

};

#endif
