// Calculate a regular basis polynom set for the requested point
// - one for each particle considering window-function with influence zones
// being different in negative and positive coordinate direction.

#ifndef BasisPolyAsym_h_
#define BasisPolyAsym_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "BasisPolynom.h"
#include "BasisPolynomSets.h"
#include "commonTypedefs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"

class BasisPolyAsym : public BasisPolynom {

  public:

    BasisPolyAsym(InputFileData* InputData,std::vector<Particle>& ptcls,
		  intVector& sPtcls,double& x,double& y,double& z,
		  int& supportSize,int& linEQSize,
		  unsigned int derivationOrder,
		  std::map<std::string,double>& modelData, 
		  std::ofstream& logFile);

    ~BasisPolyAsym() {};

  private:

    // Calculate basis polynom set.
    void calcPolynom(InputFileData* InputData,
		     std::vector<Particle>& ptcls,intVector& sPtcls,
		     double& x,double& y,double& z,int& supportSize,
		     int& linEQSize,std::ofstream& logFile);

    // Calculation of first order derivations of the basis polynom 
    void calcPoly1stDerivs(InputFileData* InputData,
			   std::vector<Particle>& ptcls,
			   intVector& sPtcls,double& x,double& y,double& z,
			   int& supportSize,int& linEQSize,
			   std::ofstream& logFile);

    // Calculation of second order derivations of the basis polynom 
    void calcPoly2ndDerivs(InputFileData* InputData,
			   std::vector<Particle>& ptcls,
			   intVector& sPtcls,double& x,double& y,double& z,
			   int& supportSize,int& linEQSize,
			   std::ofstream& logFile);

    // obtain the normalized coordinates
    void getNormedCoords(double& xnorm,double& ynorm,
			 double& znorm,dbVector& radii);

    // multiply 1st order polynomial derivative with dxnorm/dx, dynorm/dy
    // and dznorm/dz  
    void multiplyDXNormDerivative(double& xnorm,double& ynorm,
				  double& znorm,dbVector& radii,
				  dbVector& dPx,dbVector& dPy,
				  dbVector& dPz);

    // multiply 2nd order polynomial derivative with dxnorm^2/dx^2,
    // dynorm^2/dy^2,dznorm^2/dz^2,dxnorm/dx*dynorm/dy, dynorm/dy*dznorm/dz
    // and dznorm/dz*dxnorm/dx  
    void multiplyD2XNormDerivative(double& xnorm,double& ynorm,
				   double& znorm,dbVector& radii,
				   dbVector& dPxx,dbVector& dPyy,
				   dbVector& dPzz,dbVector& dPxy,
				   dbVector& dPyz,dbVector& dPzx);

};

#endif
