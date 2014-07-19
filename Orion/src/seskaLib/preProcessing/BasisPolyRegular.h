// Calculate a regular basis polynom set for the requested point
// - one for each particle.

#ifndef BasisPolyRegular_h_
#define BasisPolyRegular_h_

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

class BasisPolyRegular : public BasisPolynom {

  public:

    BasisPolyRegular(InputFileData* InputData,std::vector<Particle>& ptcls,
		     intVector& sPtcls,double& x,double& y,double& z,
		     int& supportSize,int& linEQSize,
		     unsigned int derivationOrder,
		     std::map<std::string,double>& modelData, 
		     std::ofstream& logFile);

    //For POD calculation
    BasisPolyRegular(InputFileData* InputData,std::vector<Particle>& ptcls,
    		     intVector& sPtcls,double& x,double& y,double& z,
    		     int& supportSize,int& linEQSize,
    		     unsigned int derivationOrder,std::ofstream& logFile);


    ~BasisPolyRegular() {};

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

};

#endif
