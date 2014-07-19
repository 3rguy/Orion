// Calculate the window function of each particle supporting the 
// requested point considering the influence zones to be
// different in negative and positive coordinate direction.

#ifndef WindowFuncAsym_h_
#define WindowFuncAsym_h_

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
#include "PETScFunctions.h"
#include "petscsys.h"
#include "SurfaceElementTemplates.h"
#include "VolumeElementTemplates.h"
#include "WindowFunctionSet.h"
#include "WindowFunctionSets.h"


class WindowFuncAsym : public virtual WindowFunctionSet {

  public:

    WindowFuncAsym(): numOfSegments(0),integrationOrder(0) {};
    WindowFuncAsym(InputFileData* InputData,
		   std::vector<Particle>& ptcls,
		   intVector& sPtcls,double& x,double& y,
		   double& z,int& supportSize,
		   unsigned int derivationOrder,
		   std::map<std::string,double>& modelData, 
		   std::ofstream& logFile,
		   PetscViewer& viewerSEQ);


    ~WindowFuncAsym() {};


    // Calculate window functions.
    void calcWinFunctions(InputFileData* InputData,
			  std::vector<Particle>& ptcls,
			  intVector& sPtcls,double& x,double& y,
			  double& z,int& supportSize,
			  std::map<std::string,double>& modelData, 
			  std::ofstream& logFile,
			  PetscViewer& viewerSEQ);

    // Calculation of first order derivation of the window functions 
    void calcWinFunction1stDerivs(InputFileData* InputData,
				  std::vector<Particle>& ptcls,
				  intVector& sPtcls,
				  double& x,double& y,double& z,
				  int& supportSize,
				  std::map<std::string,double>& modelData, 
				  std::ofstream& logFile,
				  PetscViewer& viewerSEQ);

    // Calculation of second order derivation of the window functions 
    void calcWinFunction2ndDerivs(InputFileData* InputData,
				  std::vector<Particle>& ptcls,
				  intVector& sPtcls,
				  double& x,double& y,double& z,
				  int& supportSize,
				  std::map<std::string,double>& modelData, 
				  std::ofstream& logFile,
				  PetscViewer& viewerSEQ);

    // compute spline coefficients (currently only cubic)
    void setCustomPtcleSpline(InputFileData* InputData,
			      Particle& ptcle,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,
			      PetscViewer& viewerSEQ);

    // compute the window function integral
    void setCustomPtcleSplineIntegral(InputFileData* InputData,
				      Particle& ptcle,
				      std::map<std::string,double>& modelData,
				      std::ofstream& logFile,
				      PetscViewer& viewerSEQ);


    //protected:

    // compute the window function ordinate for local point x
    void calcSplineValue(InputFileData* InputData,
			 Particle& ptcle,dbVector& x,double& w,
			 std::map<std::string,double>& modelData, 
			 std::ofstream& logFile,
			 PetscViewer& viewerSEQ);

    // compute the window function ordinate and its for first order 
    // derivatives for local point x
    void calcSplineValue(InputFileData* InputData,
			 Particle& ptcle,dbVector& x,double& w,
			 double& dxW,double& dyW,double& dzW,
			 std::map<std::string,double>& modelData, 
			 std::ofstream& logFile,
			 PetscViewer& viewerSEQ);

    // compute the window function ordinate and its for first-and second-order 
    // derivatives for local point x
    void calcSplineValue(InputFileData* InputData,
			 Particle& ptcle,dbVector& x,double& w,
			 double& dxW,double& dyW,double& dzW,
			 double& dxxW,double& dyyW,double& dzzW,
			 double& dxyW,double& dyzW,double& dzxW,
			 std::map<std::string,double>& modelData, 
			 std::ofstream& logFile,
			 PetscViewer& viewerSEQ);

  private:

    int numOfSegments;
    int integrationOrder;

    double cubicSpline(double s);
    double dcubicSpline(double s);
    double d2cubicSpline(double s);
    double quarticSpline(double s);

    dbVector P(double x);
    dbVector dP(double x);
    dbVector d2P(double x);
    dbVector d3P(double x);
    dbVector d4P(double x);
    dbVector d5P(double x);

    dbVector expP(double x);
    dbVector dexpP(double x);
    dbVector d2expP(double x);


};


#endif
