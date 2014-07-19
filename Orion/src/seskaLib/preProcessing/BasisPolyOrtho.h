// Calculate a orthogonal basis polynom set for the requested point
// and one for suporting each particle.

#ifndef BasisPolyOrtho_h_
#define BasisPolyOrtho_h_

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
#include "ShapefunctionSet.h"

class BasisPolyOrtho : public BasisPolynom {

  public:

    BasisPolyOrtho(InputFileData* InputData,
		   BasisPolynom* PolynomSet,
		   ShapefunctionSet* ShepardShapeSet,
		   std::vector<Particle>& ptcls,
		   intVector& sPtcls,double& x,double& y,double& z,
		   int& supportSize,unsigned int derivationOrder,
		   std::map<std::string,double>& modelData, 
		   std::ofstream& logFile,
		   PetscViewer& viewerSEQ);

    ~BasisPolyOrtho() {};

  private:

    // Calculate basis polynom set.
    void calcPolynom(InputFileData* InputData,BasisPolynom* PolynomSet,
		     ShapefunctionSet* ShepardShapeSet,
		     std::vector<Particle>& ptcls,intVector& sPtcls,
		     double& x,double& y,double& z,int& supportSize,
		     std::map<std::string,double>& modelData, 
		     std::ofstream& logFile,PetscViewer& viewerSEQ);

    // Calculation of first order derivations of the basis polynom 
    void calcPoly1stDerivs(InputFileData* InputData,BasisPolynom* PolynomSet,
			   ShapefunctionSet* ShepardShapeSet,
			   std::vector<Particle>& ptcls,intVector& sPtcls,
			   double& x,double& y,double& z,int& supportSize,
			   std::map<std::string,double>& modelData, 
			   std::ofstream& logFile,PetscViewer& viewerSEQ);

    // Calculation of second order derivations of the basis polynom 
    void calcPoly2ndDerivs(InputFileData* InputData,BasisPolynom* PolynomSet,
			   ShapefunctionSet* ShepardShapeSet,
			   std::vector<Particle>& ptcls,intVector& sPtcls,
			   double& x,double& y,double& z,int& supportSize,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile,PetscViewer& viewerSEQ);

};

#endif
