// Calculate a basis polynom set for the requested 
// point - one for each particle.

#include "BasisPolyRegular.h"

BasisPolyRegular::BasisPolyRegular(InputFileData* InputData,
				   std::vector<Particle>& ptcls,
				   intVector& sPtcls,double& x,
				   double& y,double& z,
				   int& supportSize,int& linEQSize,
				   unsigned int derivationOrder,
				   std::map<std::string,double>& modelData,  
				   std::ofstream& logFile) {

  using namespace std;

  // Calculation of basis polynom 
  calcPolynom(InputData,ptcls,sPtcls,x,y,z,supportSize,linEQSize,
	      logFile);
  
  // Calculation of first order derivations of the basis polynom
  if(derivationOrder > 0)
    calcPoly1stDerivs(InputData,ptcls,sPtcls,x,y,z,supportSize,
		      linEQSize,logFile);

  // Calculation of second order derivations of the basis polynom
  if(derivationOrder > 1)
    calcPoly2ndDerivs(InputData,ptcls,sPtcls,x,y,z,supportSize,
		      linEQSize,logFile);

}

/***********************************************************************/
/***********************************************************************/
// Calculation of basis polynom
void BasisPolyRegular::calcPolynom(InputFileData* InputData,
				   std::vector<Particle>& ptcls,
				   intVector& sPtcls,double& x,double& y,
				   double& z,int& supportSize,
				   int& linEQSize,std::ofstream& logFile) {

  using namespace std;

  int order = (int)InputData->getValue("basisPolynomOrder");
  int pType = (int)InputData->getValue("basisPolynomType");
  bool dimless = (bool)InputData->getValue("dimensionlessBasisPolynom");

  double xnorm,ynorm,znorm;

  // Select choosen polynom type (Pascal,Lagrangian,Serendipity).
  switch(pType) {

    // Pascal
  case 1:

    // Select choosen polynom order.
    switch(order) {

    case 0:
      linEQSize = 1;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      basis.resize(supportSize+1);

      // Loop over all supporting particles.
      for(int i=0;i<supportSize;i++) {
	basis[i].resize(linEQSize);

	basis[i][0] = 1.0;
      }

      basis[supportSize].resize(linEQSize);
      basis[supportSize][0] = 1.0;

      break;

    case 1:
      linEQSize = 4;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {

	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  basis[i] =  PascalLinear::P(xnorm,ynorm,znorm);
	}

	basis[supportSize].resize(linEQSize);
	basis[supportSize][0] = 1.0;
      }
      else {

	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = ptcls[sPtcls[i]].getCoord(0);
	  ynorm = ptcls[sPtcls[i]].getCoord(1);
	  znorm = ptcls[sPtcls[i]].getCoord(2);
	  
	  basis[i] =  PascalLinear::P(xnorm,ynorm,znorm);
	}

	basis[supportSize] =  PascalLinear::P(x,y,z);
      }

      break;

    case 2:
      linEQSize = 10;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  basis[i] =  PascalQuadratic::P(xnorm,ynorm,znorm);
	}

	basis[supportSize].resize(linEQSize);
	basis[supportSize][0] = 1.0;
      }
      else {
	basis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {

	  basis[i] =  PascalQuadratic::P(ptcls[sPtcls[i]].getCoord(0),
					 ptcls[sPtcls[i]].getCoord(1),
					 ptcls[sPtcls[i]].getCoord(2));
	}

	basis[supportSize] = PascalQuadratic::P(x,y,z);
      }     

      break;

    case 3:
      linEQSize = 20;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  basis[i] =  PascalCubic::P(xnorm,ynorm,znorm);
	}

	basis[supportSize].resize(linEQSize);
	basis[supportSize][0] = 1.0;
      }
      else {
	basis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {

	  basis[i] =  PascalCubic::P(ptcls[sPtcls[i]].getCoord(0),
				     ptcls[sPtcls[i]].getCoord(1),
				     ptcls[sPtcls[i]].getCoord(2));
	}

	basis[supportSize] = PascalCubic::P(x,y,z);
      }

      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Lagrangian
  case 2:

    switch(order) {

    case 1:
      linEQSize = 8;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  basis[i] =  LagrangeLinear::P(xnorm,ynorm,znorm);
	}
	
	basis[supportSize].resize(linEQSize);
	basis[supportSize][0] = 1.0;
      }
      else {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {

	  basis[i] = LagrangeLinear::P(ptcls[sPtcls[i]].getCoord(0),
				       ptcls[sPtcls[i]].getCoord(1),
				       ptcls[sPtcls[i]].getCoord(2));
	}

	basis[supportSize] = LagrangeLinear::P(x,y,z);
      }

      break;

    case 2:
      linEQSize = 27;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  basis[i] =  LagrangeQuadratic::P(xnorm,ynorm,znorm);
	}

	basis[supportSize].resize(linEQSize);
	basis[supportSize][0] = 1.0;
      }
      else {
	basis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  basis[i].resize(linEQSize);

	  basis[i] =  LagrangeQuadratic::P(ptcls[sPtcls[i]].getCoord(0),
					   ptcls[sPtcls[i]].getCoord(1),
					   ptcls[sPtcls[i]].getCoord(2));
	}

	basis[supportSize] = LagrangeQuadratic::P(x,y,z);
      }

      break;
    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Serendipity
  case 3:

    switch(order) {

    case 2:
      linEQSize = 20;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  basis[i] =  SerendipityQuadratic::P(xnorm,ynorm,znorm);
	}

	basis[supportSize].resize(linEQSize);
	basis[supportSize][0] = 1.0;
      }
      else {
	basis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {

	  basis[i] =  SerendipityQuadratic::P(ptcls[sPtcls[i]].getCoord(0),
					      ptcls[sPtcls[i]].getCoord(1),
					      ptcls[sPtcls[i]].getCoord(2));
	}

	basis[supportSize] = SerendipityQuadratic::P(x,y,z);

      }
      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Bernstein
  case 4:

    switch(order) {

    case 1:
      linEQSize = 8;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	basis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  basis[i] =  BernsteinLinear::P(xnorm,ynorm,znorm);
	}

	basis[supportSize].resize(linEQSize);
	basis[supportSize][7] = 1.0;
      }
      else {
	basis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {

	  basis[i] =  BernsteinLinear::P(ptcls[sPtcls[i]].getCoord(0),
					 ptcls[sPtcls[i]].getCoord(1),
					 ptcls[sPtcls[i]].getCoord(2));
	}

	basis[supportSize] = BernsteinLinear::P(x,y,z);
      }

      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

  default:    
    logFile<<"Chosen basis polynom type isn't supported!"<<endl;      
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of first order derivations of the basis polynom 
void BasisPolyRegular::calcPoly1stDerivs(InputFileData* InputData,
					 std::vector<Particle>& ptcls,
					 intVector& sPtcls,double& x,
					 double& y,double& z,
					 int& supportSize,int& linEQSize,
					 std::ofstream& logFile) {

  using namespace std;

  int order = (int)InputData->getValue("basisPolynomOrder");
  int pType = (int)InputData->getValue("basisPolynomType");
  bool dimless = (bool)InputData->getValue("dimensionlessBasisPolynom");

  double xnorm,ynorm,znorm;

  // Select choosen polynom type (Pascal,Lagrangian).
  switch(pType) {

  case 1:

    //Select choosen polynom order.
    switch(order) {
      
      // Pascal
    case 0:
      linEQSize = 1;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    
      if(dimless) {
	xDerivBasis.resize(supportSize);
	yDerivBasis.resize(supportSize);
	zDerivBasis.resize(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  xDerivBasis[i].resize(linEQSize);
	  yDerivBasis[i].resize(linEQSize);
	  zDerivBasis[i].resize(linEQSize);	
	}

      }
      else {
	xDerivBasis.resize(supportSize+1);
	yDerivBasis.resize(supportSize+1);
	zDerivBasis.resize(supportSize+1);	

	// Loop over all supporting particles.
	for(int i=0;i<=supportSize;i++) {
	  xDerivBasis[i].resize(linEQSize);
	  yDerivBasis[i].resize(linEQSize);
	  zDerivBasis[i].resize(linEQSize);	
	}
	
      }

      break;

    case 1:
      linEQSize = 4;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xDerivBasis[i] = PascalLinear::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = PascalLinear::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = PascalLinear::dPz(xnorm,ynorm,znorm);

	}

      }
      else {
	xDerivBasis = dbMatrix(supportSize+1);
	yDerivBasis = dbMatrix(supportSize+1);
	zDerivBasis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xDerivBasis[i] = PascalLinear::dPx(ptcls[sPtcls[i]].getCoord(0),
					     ptcls[sPtcls[i]].getCoord(1),
					     ptcls[sPtcls[i]].getCoord(2));
	  yDerivBasis[i] = PascalLinear::dPy(ptcls[sPtcls[i]].getCoord(0),
					     ptcls[sPtcls[i]].getCoord(1),
					     ptcls[sPtcls[i]].getCoord(2));
	  zDerivBasis[i] = PascalLinear::dPz(ptcls[sPtcls[i]].getCoord(0),
					     ptcls[sPtcls[i]].getCoord(1),
					     ptcls[sPtcls[i]].getCoord(2));
	}
	
	xDerivBasis[supportSize] = PascalLinear::dPx(x,y,z);
	yDerivBasis[supportSize] = PascalLinear::dPy(x,y,z);
	zDerivBasis[supportSize] = PascalLinear::dPz(x,y,z);
      }

      break;

    case 2:
      linEQSize = 10;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xDerivBasis[i] = PascalQuadratic::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = PascalQuadratic::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = PascalQuadratic::dPz(xnorm,ynorm,znorm);

	}

      }
      else {
	xDerivBasis = dbMatrix(supportSize+1);
	yDerivBasis = dbMatrix(supportSize+1);
	zDerivBasis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xDerivBasis[i] = PascalQuadratic::dPx(ptcls[sPtcls[i]].getCoord(0),
						ptcls[sPtcls[i]].getCoord(1),
						ptcls[sPtcls[i]].getCoord(2));
	  yDerivBasis[i] = PascalQuadratic::dPy(ptcls[sPtcls[i]].getCoord(0),
						ptcls[sPtcls[i]].getCoord(1),
						ptcls[sPtcls[i]].getCoord(2));
	  zDerivBasis[i] = PascalQuadratic::dPz(ptcls[sPtcls[i]].getCoord(0),
						ptcls[sPtcls[i]].getCoord(1),
						ptcls[sPtcls[i]].getCoord(2));
	}
	
	xDerivBasis[supportSize] = PascalQuadratic::dPx(x,y,z);
	yDerivBasis[supportSize] = PascalQuadratic::dPy(x,y,z);
	zDerivBasis[supportSize] = PascalQuadratic::dPz(x,y,z);
      }
      
      break;
      
    case 3:
      linEQSize = 20;
      
      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
      if(dimless) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xDerivBasis[i] = PascalCubic::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = PascalCubic::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = PascalCubic::dPz(xnorm,ynorm,znorm);
	  
	}

      }
      else {
	xDerivBasis = dbMatrix(supportSize+1);
	yDerivBasis = dbMatrix(supportSize+1);
	zDerivBasis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xDerivBasis[i] = PascalCubic::dPx(ptcls[sPtcls[i]].getCoord(0),
					    ptcls[sPtcls[i]].getCoord(1),
					    ptcls[sPtcls[i]].getCoord(2));
	  yDerivBasis[i] = PascalCubic::dPy(ptcls[sPtcls[i]].getCoord(0),
					    ptcls[sPtcls[i]].getCoord(1),
					    ptcls[sPtcls[i]].getCoord(2));
	  zDerivBasis[i] = PascalCubic::dPz(ptcls[sPtcls[i]].getCoord(0),
					    ptcls[sPtcls[i]].getCoord(1),
					    ptcls[sPtcls[i]].getCoord(2));
	}

	xDerivBasis[supportSize] = PascalCubic::dPx(x,y,z);
	yDerivBasis[supportSize] = PascalCubic::dPy(x,y,z);
	zDerivBasis[supportSize] = PascalCubic::dPz(x,y,z);
      }
      
      break;
      
    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }
    
    break;

    /*******************************************************************/
    // Lagrangian
  case 2:

    switch(order) {

    case 1:
      linEQSize = 8;
      
      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xDerivBasis[i] = LagrangeLinear::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = LagrangeLinear::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = LagrangeLinear::dPz(xnorm,ynorm,znorm);
	  
	}

      }
      else {
	xDerivBasis = dbMatrix(supportSize+1);
	yDerivBasis = dbMatrix(supportSize+1);
	zDerivBasis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xDerivBasis[i] = LagrangeLinear::dPx(ptcls[sPtcls[i]].getCoord(0),
					       ptcls[sPtcls[i]].getCoord(1),
					       ptcls[sPtcls[i]].getCoord(2));
	  yDerivBasis[i] = LagrangeLinear::dPy(ptcls[sPtcls[i]].getCoord(0),
					       ptcls[sPtcls[i]].getCoord(1),
					       ptcls[sPtcls[i]].getCoord(2));
	  zDerivBasis[i] = LagrangeLinear::dPz(ptcls[sPtcls[i]].getCoord(0),
					       ptcls[sPtcls[i]].getCoord(1),
					       ptcls[sPtcls[i]].getCoord(2)); 
	}
	
	xDerivBasis[supportSize] = LagrangeLinear::dPx(x,y,z);
	yDerivBasis[supportSize] = LagrangeLinear::dPy(x,y,z);
	zDerivBasis[supportSize] = LagrangeLinear::dPz(x,y,z);
      }

      break;

    case 2:
      linEQSize = 27;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xDerivBasis[i] = LagrangeQuadratic::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = LagrangeQuadratic::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = LagrangeQuadratic::dPz(xnorm,ynorm,znorm);
	  
	}

      }
      else {
	xDerivBasis = dbMatrix(supportSize+1);
	yDerivBasis = dbMatrix(supportSize+1);
	zDerivBasis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xDerivBasis[i] = 
	    LagrangeQuadratic::dPx(ptcls[sPtcls[i]].getCoord(0),
				   ptcls[sPtcls[i]].getCoord(1),
				   ptcls[sPtcls[i]].getCoord(2));
	  yDerivBasis[i] = 
	    LagrangeQuadratic::dPy(ptcls[sPtcls[i]].getCoord(0),
				   ptcls[sPtcls[i]].getCoord(1),
				   ptcls[sPtcls[i]].getCoord(2));
	  zDerivBasis[i] = 
	    LagrangeQuadratic::dPz(ptcls[sPtcls[i]].getCoord(0),
				   ptcls[sPtcls[i]].getCoord(1),
				   ptcls[sPtcls[i]].getCoord(2));
	  
	}
	
	xDerivBasis[supportSize] = LagrangeQuadratic::dPx(x,y,z);
	yDerivBasis[supportSize] = LagrangeQuadratic::dPy(x,y,z);
	zDerivBasis[supportSize] = LagrangeQuadratic::dPz(x,y,z);
      }

      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Serendipity
  case 3:

    switch(order) {

    case 2:
      linEQSize = 20;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xDerivBasis[i] = SerendipityQuadratic::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = SerendipityQuadratic::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = SerendipityQuadratic::dPz(xnorm,ynorm,znorm);
	}

      }
      else {
	xDerivBasis = dbMatrix(supportSize+1);
	yDerivBasis = dbMatrix(supportSize+1);
	zDerivBasis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xDerivBasis[i] = 
	    SerendipityQuadratic::dPx(ptcls[sPtcls[i]].getCoord(0),
				      ptcls[sPtcls[i]].getCoord(1),
				      ptcls[sPtcls[i]].getCoord(2));
	  yDerivBasis[i] = 
	    SerendipityQuadratic::dPy(ptcls[sPtcls[i]].getCoord(0),
				      ptcls[sPtcls[i]].getCoord(1),
				      ptcls[sPtcls[i]].getCoord(2));
	  zDerivBasis[i] = 
	    SerendipityQuadratic::dPz(ptcls[sPtcls[i]].getCoord(0),
				      ptcls[sPtcls[i]].getCoord(1),
				      ptcls[sPtcls[i]].getCoord(2));
	  
	}

	xDerivBasis[supportSize] = SerendipityQuadratic::dPx(x,y,z);
	yDerivBasis[supportSize] = SerendipityQuadratic::dPy(x,y,z);
	zDerivBasis[supportSize] = SerendipityQuadratic::dPz(x,y,z);
      }
      
      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Bernstein
  case 4:

    switch(order) {

    case 1:
      linEQSize = 8;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  xDerivBasis[i] = 
	    BernsteinLinear::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = 
	    BernsteinLinear::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = 
	    BernsteinLinear::dPz(xnorm,ynorm,znorm);
	  
	}

      }
      else {
	xDerivBasis = dbMatrix(supportSize+1);
	yDerivBasis = dbMatrix(supportSize+1);
	zDerivBasis = dbMatrix(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xDerivBasis[i] = 
	    BernsteinLinear::dPx(ptcls[sPtcls[i]].getCoord(0),
				 ptcls[sPtcls[i]].getCoord(1),
				 ptcls[sPtcls[i]].getCoord(2));
	  yDerivBasis[i] = 
	    BernsteinLinear::dPy(ptcls[sPtcls[i]].getCoord(0),
				 ptcls[sPtcls[i]].getCoord(1),
				 ptcls[sPtcls[i]].getCoord(2));
	  zDerivBasis[i] = 
	    BernsteinLinear::dPz(ptcls[sPtcls[i]].getCoord(0),
				 ptcls[sPtcls[i]].getCoord(1),
				 ptcls[sPtcls[i]].getCoord(2));
	  
	}

	xDerivBasis[supportSize] = BernsteinLinear::dPx(x,y,z);
	yDerivBasis[supportSize] = BernsteinLinear::dPy(x,y,z);
	zDerivBasis[supportSize] = BernsteinLinear::dPz(x,y,z);
      }

      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
  }

    break;

  default:    
    logFile<<"Chosen basis polynom type isn't supported!"<<endl;      
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of second order derivations of the basis polynom 
void BasisPolyRegular::calcPoly2ndDerivs(InputFileData* InputData,
					 std::vector<Particle>& ptcls,
					 intVector& sPtcls,double& x,
					 double& y,double& z,
					 int& supportSize,int& linEQSize,
					 std::ofstream& logFile) {

  using namespace std;

  int order = (int)InputData->getValue("basisPolynomOrder");
  int pType = (int)InputData->getValue("basisPolynomType");
  bool dimless = (bool)InputData->getValue("dimensionlessBasisPolynom");

  double xnorm,ynorm,znorm;

  // Select choosen polynom type(Pascal,Lagrangian,Serendipity).
  switch(pType) {

    // Pascal
  case 1:

    // Select choosen polynom order.
    switch(order) {

    case 0:
      linEQSize = 1;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for chosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
}

      if(dimless) {
	xxDerivBasis.resize(supportSize);
	yyDerivBasis.resize(supportSize);
	zzDerivBasis.resize(supportSize);
	xyDerivBasis.resize(supportSize);
	yzDerivBasis.resize(supportSize);
	zxDerivBasis.resize(supportSize);
	
	for(int i=0;i<supportSize;i++) {
	  xxDerivBasis[i].resize(linEQSize);
	  yyDerivBasis[i].resize(linEQSize);
	  zzDerivBasis[i].resize(linEQSize);
	  xyDerivBasis[i].resize(linEQSize);
	  yzDerivBasis[i].resize(linEQSize);
	  zxDerivBasis[i].resize(linEQSize);
	}

      }
      else {
	xxDerivBasis.resize(supportSize+1);
	yyDerivBasis.resize(supportSize+1);
	zzDerivBasis.resize(supportSize+1);
	xyDerivBasis.resize(supportSize+1);
	yzDerivBasis.resize(supportSize+1);
	zxDerivBasis.resize(supportSize+1);
	
	for(int i=0;i<=supportSize;i++) {
	  xxDerivBasis[i].resize(linEQSize);
	  yyDerivBasis[i].resize(linEQSize);
	  zzDerivBasis[i].resize(linEQSize);
	  xyDerivBasis[i].resize(linEQSize);
	  yzDerivBasis[i].resize(linEQSize);
	  zxDerivBasis[i].resize(linEQSize);
	}

      }

      break;

    case 1:
      linEQSize = 4;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xxDerivBasis = dbMatrix(supportSize,dbVector(linEQSize));
	yyDerivBasis = dbMatrix(supportSize,dbVector(linEQSize));
	zzDerivBasis = dbMatrix(supportSize,dbVector(linEQSize));
	xyDerivBasis = dbMatrix(supportSize,dbVector(linEQSize));
	yzDerivBasis = dbMatrix(supportSize,dbVector(linEQSize));
	zxDerivBasis = dbMatrix(supportSize,dbVector(linEQSize));
      }
      else {
	xxDerivBasis = dbMatrix(supportSize+1,dbVector(linEQSize));
	yyDerivBasis = dbMatrix(supportSize+1,dbVector(linEQSize));
	zzDerivBasis = dbMatrix(supportSize+1,dbVector(linEQSize));
	xyDerivBasis = dbMatrix(supportSize+1,dbVector(linEQSize));
	yzDerivBasis = dbMatrix(supportSize+1,dbVector(linEQSize));
	zxDerivBasis = dbMatrix(supportSize+1,dbVector(linEQSize));
      }

      break;

    case 2:
      linEQSize = 10;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xxDerivBasis = dbMatrix(supportSize);
	yyDerivBasis = dbMatrix(supportSize);
	zzDerivBasis = dbMatrix(supportSize);
	xyDerivBasis = dbMatrix(supportSize);
	yzDerivBasis = dbMatrix(supportSize);
	zxDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {

	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  
	  xxDerivBasis[i] = PascalQuadratic::dPxx(xnorm,ynorm,znorm);
	  yyDerivBasis[i] = PascalQuadratic::dPyy(xnorm,ynorm,znorm);
	  zzDerivBasis[i] = PascalQuadratic::dPzz(xnorm,ynorm,znorm);
	  xyDerivBasis[i] = PascalQuadratic::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = PascalQuadratic::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = PascalQuadratic::dPzx(xnorm,ynorm,znorm);
	  
	}

      }
      else {
	xxDerivBasis.resize(supportSize+1);
	yyDerivBasis.resize(supportSize+1);
	zzDerivBasis.resize(supportSize+1);
	xyDerivBasis.resize(supportSize+1);
	yzDerivBasis.resize(supportSize+1);
	zxDerivBasis.resize(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {

	  xxDerivBasis[i] = 
	    PascalQuadratic::dPxx(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	  yyDerivBasis[i] = 
	    PascalQuadratic::dPyy(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	  zzDerivBasis[i] = 
	    PascalQuadratic::dPzz(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	  xyDerivBasis[i] = 
	    PascalQuadratic::dPxy(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	  yzDerivBasis[i] = 
	    PascalQuadratic::dPyz(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	  zxDerivBasis[i] = 
	    PascalQuadratic::dPzx(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	}

	xxDerivBasis[supportSize] = PascalQuadratic::dPxx(x,y,z);
	yyDerivBasis[supportSize] = PascalQuadratic::dPyy(x,y,z);
	zzDerivBasis[supportSize] = PascalQuadratic::dPzz(x,y,z);
	xyDerivBasis[supportSize] = PascalQuadratic::dPxy(x,y,z);
	yzDerivBasis[supportSize] = PascalQuadratic::dPyz(x,y,z);
	zxDerivBasis[supportSize] = PascalQuadratic::dPzx(x,y,z);
      }

      break;

    case 3:
      linEQSize = 20;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xxDerivBasis = dbMatrix(supportSize);
	yyDerivBasis = dbMatrix(supportSize);
	zzDerivBasis = dbMatrix(supportSize);
	xyDerivBasis = dbMatrix(supportSize);
	yzDerivBasis = dbMatrix(supportSize);
	zxDerivBasis = dbMatrix(supportSize);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  
	  xxDerivBasis[i] = PascalCubic::dPxx(xnorm,ynorm,znorm);
	  yyDerivBasis[i] = PascalCubic::dPyy(xnorm,ynorm,znorm);
	  zzDerivBasis[i] = PascalCubic::dPzz(xnorm,ynorm,znorm);
	  xyDerivBasis[i] = PascalCubic::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = PascalCubic::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = PascalCubic::dPzx(xnorm,ynorm,znorm); 
	}

      }
      else {
	xxDerivBasis.resize(supportSize+1);
	yyDerivBasis.resize(supportSize+1);
	zzDerivBasis.resize(supportSize+1);
	xyDerivBasis.resize(supportSize+1);
	yzDerivBasis.resize(supportSize+1);
	zxDerivBasis.resize(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xxDerivBasis[i] = PascalCubic::dPxx(ptcls[sPtcls[i]].getCoord(0),
					      ptcls[sPtcls[i]].getCoord(1),
					      ptcls[sPtcls[i]].getCoord(2));
	  yyDerivBasis[i] = PascalCubic::dPyy(ptcls[sPtcls[i]].getCoord(0),
					      ptcls[sPtcls[i]].getCoord(1),
					      ptcls[sPtcls[i]].getCoord(2));
	  zzDerivBasis[i] = PascalCubic::dPzz(ptcls[sPtcls[i]].getCoord(0),
					      ptcls[sPtcls[i]].getCoord(1),
					      ptcls[sPtcls[i]].getCoord(2));
	  xyDerivBasis[i] = PascalCubic::dPxy(ptcls[sPtcls[i]].getCoord(0),
					      ptcls[sPtcls[i]].getCoord(1),
					      ptcls[sPtcls[i]].getCoord(2));
	  yzDerivBasis[i] = PascalCubic::dPyz(ptcls[sPtcls[i]].getCoord(0),
					      ptcls[sPtcls[i]].getCoord(1),
					      ptcls[sPtcls[i]].getCoord(2));
	  zxDerivBasis[i] = PascalCubic::dPzx(ptcls[sPtcls[i]].getCoord(0),
					      ptcls[sPtcls[i]].getCoord(1),
					      ptcls[sPtcls[i]].getCoord(2));
	  
	}

	xxDerivBasis[supportSize] = PascalCubic::dPxx(x,y,z);
	yyDerivBasis[supportSize] = PascalCubic::dPyy(x,y,z);
	zzDerivBasis[supportSize] = PascalCubic::dPzz(x,y,z);
	xyDerivBasis[supportSize] = PascalCubic::dPxy(x,y,z);
	yzDerivBasis[supportSize] = PascalCubic::dPyz(x,y,z);
	zxDerivBasis[supportSize] = PascalCubic::dPzx(x,y,z);
      }
      
      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Lagrangian
  case 2:

    switch(order) {

    case 1:
      linEQSize = 8;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xxDerivBasis = dbMatrix(supportSize);
	yyDerivBasis = dbMatrix(supportSize);
	zzDerivBasis = dbMatrix(supportSize);
	xyDerivBasis = dbMatrix(supportSize);
	yzDerivBasis = dbMatrix(supportSize);
	zxDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  xxDerivBasis[i].resize(linEQSize);
	  yyDerivBasis[i].resize(linEQSize); 
	  zzDerivBasis[i].resize(linEQSize);

	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xyDerivBasis[i] = LagrangeLinear::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = LagrangeLinear::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = LagrangeLinear::dPzx(xnorm,ynorm,znorm);
	}

      }
      else {
	xxDerivBasis.resize(supportSize+1);
	yyDerivBasis.resize(supportSize+1);
	zzDerivBasis.resize(supportSize+1);
	xyDerivBasis.resize(supportSize+1);
	yzDerivBasis.resize(supportSize+1);
	zxDerivBasis.resize(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  xxDerivBasis[i].resize(linEQSize);
	  yyDerivBasis[i].resize(linEQSize); 
	  zzDerivBasis[i].resize(linEQSize);
	  
	  xyDerivBasis[i] = LagrangeLinear::dPxy(ptcls[sPtcls[i]].getCoord(0),
						 ptcls[sPtcls[i]].getCoord(1),
						 ptcls[sPtcls[i]].getCoord(2));
	  yzDerivBasis[i] = LagrangeLinear::dPyz(ptcls[sPtcls[i]].getCoord(0),
						 ptcls[sPtcls[i]].getCoord(1),
						 ptcls[sPtcls[i]].getCoord(2));
	  zxDerivBasis[i] = LagrangeLinear::dPzx(ptcls[sPtcls[i]].getCoord(0),
						 ptcls[sPtcls[i]].getCoord(1),
						 ptcls[sPtcls[i]].getCoord(2));
	}
	
	xxDerivBasis[supportSize].resize(linEQSize);
	yyDerivBasis[supportSize].resize(linEQSize); 
	zzDerivBasis[supportSize].resize(linEQSize);
	xyDerivBasis[supportSize] = LagrangeLinear::dPxy(x,y,z);
	yzDerivBasis[supportSize] = LagrangeLinear::dPyz(x,y,z);
	zxDerivBasis[supportSize] = LagrangeLinear::dPzx(x,y,z);
      }

      break;
    
    case 2:
      linEQSize = 27;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xxDerivBasis = dbMatrix(supportSize);
	yyDerivBasis = dbMatrix(supportSize);
	zzDerivBasis = dbMatrix(supportSize);
	xyDerivBasis = dbMatrix(supportSize);
	yzDerivBasis = dbMatrix(supportSize);
	zxDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xxDerivBasis[i] = LagrangeQuadratic::dPxy(xnorm,ynorm,znorm);
	  yyDerivBasis[i] = LagrangeQuadratic::dPyz(xnorm,ynorm,znorm);
	  zzDerivBasis[i] = LagrangeQuadratic::dPzx(xnorm,ynorm,znorm);
	  xyDerivBasis[i] = LagrangeQuadratic::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = LagrangeQuadratic::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = LagrangeQuadratic::dPzx(xnorm,ynorm,znorm);
	}

      }
      else {
	xxDerivBasis.resize(supportSize+1);
	yyDerivBasis.resize(supportSize+1);
	zzDerivBasis.resize(supportSize+1);
	xyDerivBasis.resize(supportSize+1);
	yzDerivBasis.resize(supportSize+1);
	zxDerivBasis.resize(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xxDerivBasis[i] = 
	    LagrangeQuadratic::dPxx(ptcls[sPtcls[i]].getCoord(0),
				    ptcls[sPtcls[i]].getCoord(1),
				    ptcls[sPtcls[i]].getCoord(2));
	  yyDerivBasis[i] = 
	    LagrangeQuadratic::dPyy(ptcls[sPtcls[i]].getCoord(0),
				    ptcls[sPtcls[i]].getCoord(1),
				    ptcls[sPtcls[i]].getCoord(2));
	  zzDerivBasis[i] = 
	    LagrangeQuadratic::dPzz(ptcls[sPtcls[i]].getCoord(0),
				    ptcls[sPtcls[i]].getCoord(1),
				    ptcls[sPtcls[i]].getCoord(2));
	  xyDerivBasis[i] = 
	    LagrangeQuadratic::dPxy(ptcls[sPtcls[i]].getCoord(0),
				    ptcls[sPtcls[i]].getCoord(1),
				    ptcls[sPtcls[i]].getCoord(2));
	  yzDerivBasis[i] = 
	    LagrangeQuadratic::dPyz(ptcls[sPtcls[i]].getCoord(0),
				    ptcls[sPtcls[i]].getCoord(1),
				    ptcls[sPtcls[i]].getCoord(2));
	  zxDerivBasis[i] = 
	    LagrangeQuadratic::dPzx(ptcls[sPtcls[i]].getCoord(0),
				    ptcls[sPtcls[i]].getCoord(1),
				    ptcls[sPtcls[i]].getCoord(2));
	}
	
	xxDerivBasis[supportSize] = LagrangeQuadratic::dPxx(x,y,z);
	yyDerivBasis[supportSize] = LagrangeQuadratic::dPyy(x,y,z); 
	zzDerivBasis[supportSize] = LagrangeQuadratic::dPzz(x,y,z);
	xyDerivBasis[supportSize] = LagrangeQuadratic::dPxy(x,y,z);
	yzDerivBasis[supportSize] = LagrangeQuadratic::dPyz(x,y,z);
	zxDerivBasis[supportSize] = LagrangeQuadratic::dPzx(x,y,z);
      }

      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Serendipity
  case 3:

    switch(order) {

    case 2:
      linEQSize = 20;

      // Check if amount of supporting particles is sufficient for choosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for choosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(dimless) {
	xxDerivBasis = dbMatrix(supportSize);
	yyDerivBasis = dbMatrix(supportSize);
	zzDerivBasis = dbMatrix(supportSize);
	xyDerivBasis = dbMatrix(supportSize);
	yzDerivBasis = dbMatrix(supportSize);
	zxDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);
	  
	  xxDerivBasis[i] = SerendipityQuadratic::dPxy(xnorm,ynorm,znorm);
	  yyDerivBasis[i] = SerendipityQuadratic::dPyz(xnorm,ynorm,znorm);
	  zzDerivBasis[i] = SerendipityQuadratic::dPzx(xnorm,ynorm,znorm);
	  xyDerivBasis[i] = SerendipityQuadratic::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = SerendipityQuadratic::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = SerendipityQuadratic::dPzx(xnorm,ynorm,znorm);
	}

      }
      else {

	xxDerivBasis.resize(supportSize+1);
	yyDerivBasis.resize(supportSize+1);
	zzDerivBasis.resize(supportSize+1);
	xyDerivBasis.resize(supportSize+1);
	yzDerivBasis.resize(supportSize+1);
	zxDerivBasis.resize(supportSize+1);

	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xxDerivBasis[i] = 
	    SerendipityQuadratic::dPxx(ptcls[sPtcls[i]].getCoord(0),
				       ptcls[sPtcls[i]].getCoord(1),
				       ptcls[sPtcls[i]].getCoord(2));
	  yyDerivBasis[i] = 
	    SerendipityQuadratic::dPyy(ptcls[sPtcls[i]].getCoord(0),
				       ptcls[sPtcls[i]].getCoord(1),
				       ptcls[sPtcls[i]].getCoord(2));
	  zzDerivBasis[i] = 
	    SerendipityQuadratic::dPzz(ptcls[sPtcls[i]].getCoord(0),
				       ptcls[sPtcls[i]].getCoord(1),
				       ptcls[sPtcls[i]].getCoord(2));
	  xyDerivBasis[i] = 
	    SerendipityQuadratic::dPxy(ptcls[sPtcls[i]].getCoord(0),
				       ptcls[sPtcls[i]].getCoord(1),
				       ptcls[sPtcls[i]].getCoord(2));
	  yzDerivBasis[i] = 
	    SerendipityQuadratic::dPyz(ptcls[sPtcls[i]].getCoord(0),
				       ptcls[sPtcls[i]].getCoord(1),
				       ptcls[sPtcls[i]].getCoord(2));
	  zxDerivBasis[i] = 
	    SerendipityQuadratic::dPzx(ptcls[sPtcls[i]].getCoord(0),
				       ptcls[sPtcls[i]].getCoord(1),
				       ptcls[sPtcls[i]].getCoord(2));
	  
	}
	
	xxDerivBasis[supportSize] = SerendipityQuadratic::dPxx(x,y,z);
	yyDerivBasis[supportSize] = SerendipityQuadratic::dPyy(x,y,z); 
	zzDerivBasis[supportSize] = SerendipityQuadratic::dPzz(x,y,z);
	xyDerivBasis[supportSize] = SerendipityQuadratic::dPxy(x,y,z);
	yzDerivBasis[supportSize] = SerendipityQuadratic::dPyz(x,y,z);
	zxDerivBasis[supportSize] = SerendipityQuadratic::dPzx(x,y,z);
      }
      
      break;

    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;
    
    /*******************************************************************/
    // Bernstein
  case 4:

    switch(order) {
      
    case 1:
      linEQSize = 8;
      
      // Check if amount of supporting particles is sufficient for chosen 
      // polynom order.
      if(linEQSize > supportSize) {
	logFile<<"Not enough supporting particles for chosen polynom "
	       <<"order!(existing: "<<supportSize<<" - needed: "
	       <<linEQSize<<")"<<endl;      
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
      if(dimless) {
	xxDerivBasis.resize(supportSize);
	yyDerivBasis.resize(supportSize);
	zzDerivBasis.resize(supportSize);
	xyDerivBasis.resize(supportSize);
	yzDerivBasis.resize(supportSize);
	zxDerivBasis.resize(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  xxDerivBasis[i].resize(linEQSize);
	  yyDerivBasis[i].resize(linEQSize); 
	  zzDerivBasis[i].resize(linEQSize);
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  xyDerivBasis[i] = 
	    BernsteinLinear::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = 
	    BernsteinLinear::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = 
	    BernsteinLinear::dPzx(xnorm,ynorm,znorm);
	}
	
      }
      else {
	xxDerivBasis.resize(supportSize+1);
	yyDerivBasis.resize(supportSize+1);
	zzDerivBasis.resize(supportSize+1);
	xyDerivBasis.resize(supportSize+1);
	yzDerivBasis.resize(supportSize+1);
	zxDerivBasis.resize(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  xxDerivBasis[i].resize(linEQSize);
	  yyDerivBasis[i].resize(linEQSize); 
	  zzDerivBasis[i].resize(linEQSize);
	  
	  xyDerivBasis[i] = 
	    BernsteinLinear::dPxy(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	  yzDerivBasis[i] = 
	    BernsteinLinear::dPyz(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	  zxDerivBasis[i] = 
	    BernsteinLinear::dPzx(ptcls[sPtcls[i]].getCoord(0),
				  ptcls[sPtcls[i]].getCoord(1),
				  ptcls[sPtcls[i]].getCoord(2));
	}
	
	xxDerivBasis[supportSize].resize(linEQSize);
	yyDerivBasis[supportSize].resize(linEQSize); 
	zzDerivBasis[supportSize].resize(linEQSize);
	xyDerivBasis[supportSize] = BernsteinLinear::dPxy(x,y,z);
	yzDerivBasis[supportSize] = BernsteinLinear::dPyz(x,y,z);
	zxDerivBasis[supportSize] = BernsteinLinear::dPzx(x,y,z);
	
      }

      break;
      
    default:
      logFile<<"Chosen basis polynom order isn't supported!"<<endl;      
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

  default:    
    logFile<<"Chosen basis polynom type isn't supported!"<<endl;      
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
}
