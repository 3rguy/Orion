// Calculate a basis polynom set for the requested 
// point taking into accout that that the window-function 
// is asymmetric - one for each particle.

#include "BasisPolyAsym.h"

BasisPolyAsym::BasisPolyAsym(InputFileData* InputData,
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
void BasisPolyAsym::calcPolynom(InputFileData* InputData,
				std::vector<Particle>& ptcls,
				intVector& sPtcls,double& x,double& y,
				double& z,int& supportSize,
				int& linEQSize,std::ofstream& logFile) {

  using namespace std;

  int order = (int)InputData->getValue("basisPolynomOrder");
  int pType = (int)InputData->getValue("basisPolynomType");
  int normedBasis = (int)InputData->getValue("normedBasisPolynom");


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

      if(normedBasis > 0) {

	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
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

      if(normedBasis > 0) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
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

      if(normedBasis > 0) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
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

      if(normedBasis > 0) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
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

      if(normedBasis > 0) {
	basis = dbMatrix(supportSize+1);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
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

  default:    
    logFile<<"Chosen basis polynom type isn't supported!"<<endl;      
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of first order derivations of the basis polynom 
void BasisPolyAsym::calcPoly1stDerivs(InputFileData* InputData,
				      std::vector<Particle>& ptcls,
				      intVector& sPtcls,double& x,
				      double& y,double& z,
				      int& supportSize,int& linEQSize,
				      std::ofstream& logFile) {

  using namespace std;

  int order = (int)InputData->getValue("basisPolynomOrder");
  int pType = (int)InputData->getValue("basisPolynomType");
  int normedBasis = (int)InputData->getValue("normedBasisPolynom");


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
    
      if(normedBasis > 0) {
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

      if(normedBasis > 0) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  xDerivBasis[i] = PascalLinear::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = PascalLinear::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = PascalLinear::dPz(xnorm,ynorm,znorm);

	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm/dx = 1/r considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyDXNormDerivative(xnorm,ynorm,znorm,
				     ptcls[sPtcls[i]].getRadii(),
				     xDerivBasis[i],yDerivBasis[i],
	   				   xDerivBasis[i]);

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

      if(normedBasis > 0) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  xDerivBasis[i] = PascalQuadratic::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = PascalQuadratic::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = PascalQuadratic::dPz(xnorm,ynorm,znorm);

	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm/dx = 1/r considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyDXNormDerivative(xnorm,ynorm,znorm,
				     ptcls[sPtcls[i]].getRadii(),
				     xDerivBasis[i],yDerivBasis[i],
				     xDerivBasis[i]);

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
      
      if(normedBasis > 0) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  xDerivBasis[i] = PascalCubic::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = PascalCubic::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = PascalCubic::dPz(xnorm,ynorm,znorm);
	  
	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm/dx = 1/r considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyDXNormDerivative(xnorm,ynorm,znorm,
				     ptcls[sPtcls[i]].getRadii(),
				     xDerivBasis[i],yDerivBasis[i],
				     xDerivBasis[i]);

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

      if(normedBasis > 0) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  xDerivBasis[i] = LagrangeLinear::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = LagrangeLinear::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = LagrangeLinear::dPz(xnorm,ynorm,znorm);
	  
	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm/dx = 1/r considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyDXNormDerivative(xnorm,ynorm,znorm,
				     ptcls[sPtcls[i]].getRadii(),
				     xDerivBasis[i],yDerivBasis[i],
				     xDerivBasis[i]);

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

      if(normedBasis > 0) {
	xDerivBasis = dbMatrix(supportSize);
	yDerivBasis = dbMatrix(supportSize);
	zDerivBasis = dbMatrix(supportSize);
	
	// Loop over all supporting particles.
	for(int i=0;i<supportSize;i++) {
	  
	  xnorm = x - ptcls[sPtcls[i]].getCoord(0);
	  ynorm = y - ptcls[sPtcls[i]].getCoord(1);
	  znorm = z - ptcls[sPtcls[i]].getCoord(2);

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  xDerivBasis[i] = LagrangeQuadratic::dPx(xnorm,ynorm,znorm);
	  yDerivBasis[i] = LagrangeQuadratic::dPy(xnorm,ynorm,znorm);
	  zDerivBasis[i] = LagrangeQuadratic::dPz(xnorm,ynorm,znorm);
	  
	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm/dx = 1/r considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyDXNormDerivative(xnorm,ynorm,znorm,
				     ptcls[sPtcls[i]].getRadii(),
				     xDerivBasis[i],yDerivBasis[i],
				     xDerivBasis[i]);

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

  default:    
    logFile<<"Chosen basis polynom type isn't supported!"<<endl;      
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of second order derivations of the basis polynom 
void BasisPolyAsym::calcPoly2ndDerivs(InputFileData* InputData,
				      std::vector<Particle>& ptcls,
				      intVector& sPtcls,double& x,
				      double& y,double& z,
				      int& supportSize,int& linEQSize,
				      std::ofstream& logFile) {

  using namespace std;

  int order = (int)InputData->getValue("basisPolynomOrder");
  int pType = (int)InputData->getValue("basisPolynomType");
  bool normedBasis = (bool)InputData->getValue("normedBasisPolynom");


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

      if(normedBasis > 0) {
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

      if(normedBasis > 0) {
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

      if(normedBasis > 0) {
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

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  
	  xxDerivBasis[i] = PascalQuadratic::dPxx(xnorm,ynorm,znorm);
	  yyDerivBasis[i] = PascalQuadratic::dPyy(xnorm,ynorm,znorm);
	  zzDerivBasis[i] = PascalQuadratic::dPzz(xnorm,ynorm,znorm);
	  xyDerivBasis[i] = PascalQuadratic::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = PascalQuadratic::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = PascalQuadratic::dPzx(xnorm,ynorm,znorm);
	  
	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm^2/dx^ = 1/r^ considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyD2XNormDerivative(xnorm,ynorm,znorm,
				      ptcls[sPtcls[i]].getRadii(),
				      xxDerivBasis[i],yyDerivBasis[i],
				      zzDerivBasis[i],xyDerivBasis[i],
				      yzDerivBasis[i],zxDerivBasis[i]);

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

      if(normedBasis > 0) {
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

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  
	  xxDerivBasis[i] = PascalCubic::dPxx(xnorm,ynorm,znorm);
	  yyDerivBasis[i] = PascalCubic::dPyy(xnorm,ynorm,znorm);
	  zzDerivBasis[i] = PascalCubic::dPzz(xnorm,ynorm,znorm);
	  xyDerivBasis[i] = PascalCubic::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = PascalCubic::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = PascalCubic::dPzx(xnorm,ynorm,znorm); 

	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm^2/dx^ = 1/r^ considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyD2XNormDerivative(xnorm,ynorm,znorm,
				      ptcls[sPtcls[i]].getRadii(),
				      xxDerivBasis[i],yyDerivBasis[i],
				      zzDerivBasis[i],xyDerivBasis[i],
				      yzDerivBasis[i],zxDerivBasis[i]);
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

      if(normedBasis > 0) {
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

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  xyDerivBasis[i] = LagrangeLinear::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = LagrangeLinear::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = LagrangeLinear::dPzx(xnorm,ynorm,znorm);

	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm^2/dx^ = 1/r^ considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyD2XNormDerivative(xnorm,ynorm,znorm,
				      ptcls[sPtcls[i]].getRadii(),
				      xxDerivBasis[i],yyDerivBasis[i],
				      zzDerivBasis[i],xyDerivBasis[i],
				      yzDerivBasis[i],zxDerivBasis[i]);
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

      if(normedBasis > 0) {
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

	  if(normedBasis == 2)
	    getNormedCoords(xnorm,ynorm,znorm,ptcls[sPtcls[i]].getRadii());
	  
	  xxDerivBasis[i] = LagrangeQuadratic::dPxy(xnorm,ynorm,znorm);
	  yyDerivBasis[i] = LagrangeQuadratic::dPyz(xnorm,ynorm,znorm);
	  zzDerivBasis[i] = LagrangeQuadratic::dPzx(xnorm,ynorm,znorm);
	  xyDerivBasis[i] = LagrangeQuadratic::dPxy(xnorm,ynorm,znorm);
	  yzDerivBasis[i] = LagrangeQuadratic::dPyz(xnorm,ynorm,znorm);
	  zxDerivBasis[i] = LagrangeQuadratic::dPzx(xnorm,ynorm,znorm);

	  // Loop over all supporting particles and multiply each entry 
	  // with dxnorm^2/dx^ = 1/r^ considering the chain rule 
	  // dP/dx = dP/dxnorm*dxnorm/dx

	  if(normedBasis == 2)
	    multiplyD2XNormDerivative(xnorm,ynorm,znorm,
				      ptcls[sPtcls[i]].getRadii(),
				      xxDerivBasis[i],yyDerivBasis[i],
				      zzDerivBasis[i],xyDerivBasis[i],
				      yzDerivBasis[i],zxDerivBasis[i]);
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

  default:    
    logFile<<"Chosen basis polynom type isn't supported!"<<endl;      
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
}

/************************************************************************/
// obtain the normalized coordinates
void BasisPolyAsym::getNormedCoords(double& xnorm,double& ynorm,
				    double& znorm,dbVector& radii) {

  using namespace std;

  if(xnorm > 0)
    xnorm /= radii[0];
  else
    xnorm /= radii[3];
  
  if(ynorm > 0)
    ynorm /= radii[1];
  else
    ynorm /= radii[4];
  
  if(znorm > 0)
    znorm /= radii[2];
  else
    znorm /= radii[5];

}

/************************************************************************/
// multiply 1st order polynomial derivative with dxnorm/dx, dynorm/dy
// and dznorm/dz  
void BasisPolyAsym::multiplyDXNormDerivative(double& xnorm,double& ynorm,
					     double& znorm,dbVector& radii,
					     dbVector& dPx,dbVector& dPy,
					     dbVector& dPz) {

  using namespace std;

  for(int i=0;i<dPx.size();i++) {

    if(xnorm > 0)
      dPx[i] /= radii[0];
    else
      dPx[i] /= radii[3];
  
  }

  for(int i=0;i<dPy.size();i++) {

    if(ynorm > 0)
      dPy[i] /= radii[1];
    else
      dPy[i] /= radii[4];

  }
  
  for(int i=0;i<dPz.size();i++) {

    if(znorm > 0)
      dPz[i] /= radii[2];
    else
      dPz[i] /= radii[5];

  }

}

/************************************************************************/
// multiply 2nd order polynomial derivative with dxnorm^2/dx^2,
// dynorm^2/dy^2,dznorm^2/dz^2,dxnorm/dx*dynorm/dy, dynorm/dy*dznorm/dz
// and dznorm/dz*dxnorm/dx  
void BasisPolyAsym::multiplyD2XNormDerivative(double& xnorm,double& ynorm,
					      double& znorm,dbVector& radii,
					      dbVector& dPxx,dbVector& dPyy,
					      dbVector& dPzz,dbVector& dPxy,
					      dbVector& dPyz,dbVector& dPzx) {

  using namespace std;

  for(int i=0;i<dPxx.size();i++) {

    if(xnorm > 0)
      dPxx[i] /= radii[0]*radii[0];
    else
      dPxx[i] /= radii[3]*radii[3];
  
  }

  for(int i=0;i<dPyy.size();i++) {

    if(ynorm > 0)
      dPyy[i] /= radii[1]*radii[1];
    else
      dPyy[i] /= radii[4]*radii[4];

  }
  
  for(int i=0;i<dPzz.size();i++) {

    if(znorm > 0)
      dPzz[i] /= radii[2]*radii[2];
    else
      dPzz[i] /= radii[5]*radii[5];

  }

  for(int i=0;i<dPxy.size();i++) {

    if(xnorm > 0)
      dPxy[i] /= radii[0];
    else
      dPxy[i] /= radii[3];

    if(ynorm > 0)
      dPxy[i] /= radii[1];
    else
      dPxy[i] /= radii[4];
  
  }

  for(int i=0;i<dPyz.size();i++) {

    if(ynorm > 0)
      dPyz[i] /= radii[1];
    else
      dPyz[i] /= radii[4];

    if(znorm > 0)
      dPyz[i] /= radii[2];
    else
      dPyz[i] /= radii[5];

  }
  
  for(int i=0;i<dPzx.size();i++) {

    if(znorm > 0)
      dPzx[i] /= radii[2];
    else
      dPzx[i] /= radii[5];

    if(xnorm > 0)
      dPzx[i] /= radii[0];
    else
      dPzx[i] /= radii[3];

  }

}
