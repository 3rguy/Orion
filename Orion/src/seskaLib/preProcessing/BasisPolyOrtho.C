// Calculate a orthogonal basis polynom set for the requested point
// and one for suporting each particle.

#include "BasisPolyOrtho.h"

BasisPolyOrtho::BasisPolyOrtho(InputFileData* InputData,
			       BasisPolynom* PolynomSet,
			       ShapefunctionSet* ShepardShapeSet,
			       std::vector<Particle>& ptcls,
			       intVector& sPtcls,double& x,
			       double& y,double& z,
			       int& supportSize,
			       unsigned int derivationOrder,
			       std::map<std::string,double>& modelData,  
			       std::ofstream& logFile,
			       PetscViewer& viewerSEQ) {

  using namespace std;

  switch(derivationOrder) {

    // Calculation of basis polynomial only
  case 0: 

    calcPolynom(InputData,PolynomSet,ShepardShapeSet,ptcls,sPtcls,x,y,z,
		supportSize,modelData,logFile,viewerSEQ);

    break;
  
    // Calculation of the basis polynomial and its first order derivations.
  case 1:
    
    calcPoly1stDerivs(InputData,PolynomSet,ShepardShapeSet,ptcls,sPtcls,
		      x,y,z,supportSize,modelData,logFile,viewerSEQ);
    break;

    // Calculation of the basis polynomial, its first order and second 
    // order derivations.
  case 2:
 
    calcPoly2ndDerivs(InputData,PolynomSet,ShepardShapeSet,ptcls,sPtcls,x,y,z,
		      supportSize,modelData,logFile,viewerSEQ);

    break;

  default:
    logFile<<"In BasisPolyOrtho::BasisPolyOrtho polynomial derivation\n"
	   <<" order is not supported!"<<endl;
    cerr<<"In BasisPolyOrtho::BasisPolyOrtho polynomial derivation\n"
	<<" order is not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of basis polynom
void BasisPolyOrtho::calcPolynom(InputFileData* InputData,
				 BasisPolynom* PolynomSet,
				 ShapefunctionSet* ShepardShapeSet,
				 std::vector<Particle>& ptcls,
				 intVector& sPtcls,double& x,double& y,
				 double& z,int& supportSize,
				 std::map<std::string,double>& modelData,  
				 std::ofstream& logFile,
				 PetscViewer& viewerSEQ) {
  
  using namespace std;

  dbMatrix& P = PolynomSet->getBasis();

  dbVector& shepardShapes = ShepardShapeSet->getShapefunctions();

  // Determine the actual basis at particlur point and also at its 
  // supporting particles. 
  int linEQSize = P[0].size();

  basis = dbMatrix(supportSize+1);

  // Loop over all supporting particles of any point and the point itself 
  // to determine for each a basis polynomial.
  for(int i=0;i<=supportSize;i++) {
    basis[i].resize(linEQSize);

    // Loop over the basis entries.
    for(int j=1;j<linEQSize;j++) {
      basis[i][j] =  P[i][j];
      
      // sum up over all shepard shape functions at the current
      // particle by doing a loop over its neighbours.
      for(int k=0;k<supportSize;k++)
	basis[i][j] -= shepardShapes[k]*P[k][j];

    }
    
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of first order derivations of the basis polynom 
void BasisPolyOrtho::calcPoly1stDerivs(InputFileData* InputData,
				       BasisPolynom* PolynomSet,
				       ShapefunctionSet* ShepardShapeSet,
				       std::vector<Particle>& ptcls,
				       intVector& sPtcls,double& x,
				       double& y,double& z,
				       int& supportSize,
				       std::map<std::string,double>& modelData,  
				       std::ofstream& logFile,
				       PetscViewer& viewerSEQ) {

  using namespace std;

  dbMatrix& P = PolynomSet->getBasis();
  dbMatrix& dPx = PolynomSet->getXDerivBasis();
  dbMatrix& dPy = PolynomSet->getYDerivBasis();
  dbMatrix& dPz = PolynomSet->getZDerivBasis();

  dbVector& shepardShapes = ShepardShapeSet->getShapefunctions();
  dbMatrix& firstDerivShepards = ShepardShapeSet->getFirstDerivShapes();

  // Determine the actual basis at particlur point and also at its 
  // supporting particles. 
  int linEQSize = P[0].size();

  basis = dbMatrix(supportSize+1);
  xDerivBasis.resize(supportSize+1);
  yDerivBasis.resize(supportSize+1);
  zDerivBasis.resize(supportSize+1);

  // Loop over all supporting particles of any point and the point itself 
  // to determine for each the first order derivations a basis polynomial.
  for(int i=0;i<=supportSize;i++) {

    basis[i].resize(linEQSize);

    xDerivBasis[i].resize(linEQSize);
    yDerivBasis[i].resize(linEQSize);
    zDerivBasis[i].resize(linEQSize);	

    // Loop over the basis entries.
    for(int j=1;j<linEQSize;j++) {
      basis[i][j] =  P[i][j];

      if(i < supportSize) {
	xDerivBasis[i][j] = 0;
	yDerivBasis[i][j] = 0;
	zDerivBasis[i][j] = 0;
      }
      else {
	xDerivBasis[i][j] = dPx[i][j];
	yDerivBasis[i][j] = dPy[i][j];
	zDerivBasis[i][j] = dPz[i][j];
      }
      
      // sum up over all shepard shape function derivations at the current
      // particle by doing a loop over its neighbours.
      for(int k=0;k<supportSize;k++) {

	// basis polyomial
	basis[i][j] -=  shepardShapes[k]*P[k][j];

	// first derivations of the basis polynomial
	xDerivBasis[i][j] -=  firstDerivShepards[0][k]*P[k][j];
	yDerivBasis[i][j] -=  firstDerivShepards[1][k]*P[k][j];
	zDerivBasis[i][j] -=  firstDerivShepards[2][k]*P[k][j];
      }

    }
    
  }
  
}

/***********************************************************************/
/***********************************************************************/
// Calculation of second order derivations of the basis polynom 
void BasisPolyOrtho::calcPoly2ndDerivs(InputFileData* InputData,
				       BasisPolynom* PolynomSet,
				       ShapefunctionSet* ShepardShapeSet,
				       std::vector<Particle>& ptcls,
				       intVector& sPtcls,double& x,
				       double& y,double& z,
				       int& supportSize,
				       std::map<std::string,double>& modelData,  
				       std::ofstream& logFile,
				       PetscViewer& viewerSEQ) {

  using namespace std;

  dbMatrix& P = PolynomSet->getBasis();
  dbMatrix& dPx = PolynomSet->getXDerivBasis();
  dbMatrix& dPy = PolynomSet->getYDerivBasis();
  dbMatrix& dPz = PolynomSet->getZDerivBasis();
  dbMatrix& dPxx = PolynomSet->getXXDerivBasis();
  dbMatrix& dPyy = PolynomSet->getYYDerivBasis();
  dbMatrix& dPzz = PolynomSet->getZZDerivBasis();
  dbMatrix& dPxy = PolynomSet->getXYDerivBasis();
  dbMatrix& dPyz = PolynomSet->getYZDerivBasis();
  dbMatrix& dPzx = PolynomSet->getZXDerivBasis();

  dbVector& shepardShapes = ShepardShapeSet->getShapefunctions();
  dbMatrix& firstDerivShepards = ShepardShapeSet->getFirstDerivShapes();
  dbMatrix& secondDerivShepards = ShepardShapeSet->getSecondDerivShapes();
 
  // Determine the actual basis at particlur point and also at its 
  // supporting particles. 
  int linEQSize = P[0].size();

  basis.resize(supportSize+1);

  xDerivBasis.resize(supportSize+1);
  yDerivBasis.resize(supportSize+1);
  zDerivBasis.resize(supportSize+1);

  xxDerivBasis.resize(supportSize+1);
  yyDerivBasis.resize(supportSize+1);
  zzDerivBasis.resize(supportSize+1);
  xyDerivBasis.resize(supportSize+1);
  yzDerivBasis.resize(supportSize+1);
  zxDerivBasis.resize(supportSize+1);

  // Loop over all supporting particles of any point and the point itself 
  // to determine for each the second order derivations a basis polynomial.
  for(int i=0;i<=supportSize;i++) {

    basis[i].resize(linEQSize);

    xDerivBasis[i].resize(linEQSize);
    yDerivBasis[i].resize(linEQSize);
    zDerivBasis[i].resize(linEQSize);	

    xxDerivBasis[i].resize(linEQSize);
    yyDerivBasis[i].resize(linEQSize);
    zzDerivBasis[i].resize(linEQSize);	
    xyDerivBasis[i].resize(linEQSize);
    yzDerivBasis[i].resize(linEQSize);
    zxDerivBasis[i].resize(linEQSize);	
    
    // Loop over the basis entries.
    for(int j=1;j<linEQSize;j++) {
      basis[i][j] =  P[i][j];

      if(i < supportSize) {
	xDerivBasis[i][j] = 0;
	yDerivBasis[i][j] = 0;
	zDerivBasis[i][j] = 0;

	xxDerivBasis[i][j] = 0;
	yyDerivBasis[i][j] = 0;
	zzDerivBasis[i][j] = 0;
	xyDerivBasis[i][j] = 0;
	yzDerivBasis[i][j] = 0;
	zxDerivBasis[i][j] = 0;
      }
      else {
	xDerivBasis[i][j] = dPx[i][j];
	yDerivBasis[i][j] = dPy[i][j];
	zDerivBasis[i][j] = dPz[i][j];
	
	xxDerivBasis[i][j] = dPxx[i][j];
	yyDerivBasis[i][j] = dPyy[i][j];
	zzDerivBasis[i][j] = dPzz[i][j];
	xyDerivBasis[i][j] = dPxy[i][j];
	yzDerivBasis[i][j] = dPyz[i][j];
	zxDerivBasis[i][j] = dPzx[i][j];
      }
      
      // sum up over all shepard shape function derivations at the current
      // particle by doing a loop over its neighbours.
      for(int k=0;k<supportSize;k++) {

	// basis polyomial
	basis[i][j] -= shepardShapes[k]*P[k][j];

	// first derivations of the basis polynomial
	xDerivBasis[i][j] -=  firstDerivShepards[0][k]*P[k][j];
	yDerivBasis[i][j] -=  firstDerivShepards[1][k]*P[k][j];
	zDerivBasis[i][j] -=  firstDerivShepards[2][k]*P[k][j];

	// second derivations of the basis polynomial
	xxDerivBasis[i][j] -= secondDerivShepards[0][k]*P[k][j];
	yyDerivBasis[i][j] -= secondDerivShepards[1][k]*P[k][j];
	zzDerivBasis[i][j] -= secondDerivShepards[2][k]*P[k][j];

	xyDerivBasis[i][j] -= secondDerivShepards[3][k]*P[k][j];
	yzDerivBasis[i][j] -= secondDerivShepards[4][k]*P[k][j];
	zxDerivBasis[i][j] -= secondDerivShepards[5][k]*P[k][j];
      }
   
    }
 
  }

}
