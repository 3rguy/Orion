// Calculate the window function of each particle supporting the 
// requested point.

#include "WindowFunctionSet.h"


WindowFunctionSet::WindowFunctionSet(InputFileData* InputData,
				     std::vector<Particle>& ptcls,
				     intVector& sPtcls,
				     double& x,double& y,double& z,
				     int& supportSize,
				     unsigned int derivationOrder,
				     std::map<std::string,double>& modelData, 
				     std::ofstream& logFile) {

  using namespace std;


  // Calculation of the window functions 
  calcWinFunctions(InputData,ptcls,sPtcls,x,y,z,supportSize,logFile);
  
  // Calculation of first order derivation of the window functions 
  if(derivationOrder > 0)
    calcWinFunction1stDerivs(InputData,ptcls,sPtcls,x,y,z,
			     supportSize,logFile);

  // Calculation of second order derivation of the window functions 
  if(derivationOrder > 1)
    calcWinFunction2ndDerivs(InputData,ptcls,sPtcls,x,y,z,
			     supportSize,logFile);


}

// For POD Calculation
WindowFunctionSet::WindowFunctionSet(InputFileData* InputData,
				     std::vector<Particle>& ptcls,
				     intVector& sPtcls,
				     double& x,double& y,double& z,
				     int& supportSize,
				     unsigned int derivationOrder,
				     std::ofstream& logFile) {

  using namespace std;


  // Calculation of the window functions
  calcWinFunctions(InputData,ptcls,sPtcls,x,y,z,supportSize,logFile);

  // Calculation of first order derivation of the window functions
  if(derivationOrder > 0)
    calcWinFunction1stDerivs(InputData,ptcls,sPtcls,x,y,z,
			     supportSize,logFile);

  // Calculation of second order derivation of the window functions
  if(derivationOrder > 1)
    calcWinFunction2ndDerivs(InputData,ptcls,sPtcls,x,y,z,
			     supportSize,logFile);


}

/***********************************************************************/
/***********************************************************************/
// Calculation of the window functions 
void WindowFunctionSet::calcWinFunctions(InputFileData* InputData,
					 std::vector<Particle>& ptcls,
					 intVector& sPtcls,
					 double& x,double& y,double& z,
					 int& supportSize,
					 std::ofstream& logFile) {

  using namespace std;

  int winFuncNorming = (int)InputData->getValue("windowfunctionNorming");
  int wType = (int)InputData->getValue("windowFunctionType");
  int particlesNum = ptcls.size();


  int radiusMode = (int)InputData->getValue("radiusDeterminationAlgorithm");

  double dxnorm,dynorm,dznorm,normFactor,xnorm,ynorm,znorm;
  int m,n;

  double xDist,yDist,zDist,dist,maxDist;

  if(windowFuncs.size() < supportSize)
    windowFuncs.resize(supportSize);

  // Choose a window function type.
  switch(wType) {

    // constant 1:
  case 0:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;
      
      windowFuncs[i] = 1.0/normFactor;
    }

    break;

    // cubic spline
  case 1:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      CubicSplineWindowFunc::calcWinFunc(xnorm,ynorm,znorm,windowFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;



      windowFuncs[i] *= normFactor;

    }

    break;

  case 2:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuarticSplineWindowFunc::calcWinFunc(xnorm,ynorm,znorm,windowFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;



      windowFuncs[i] *= normFactor;

    }

    break;

  case 3:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      ExponentialWindowFunc::calcWinFunc(xnorm,ynorm,znorm,windowFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;



      windowFuncs[i] *= normFactor;

    }

    break;

  case 4:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      GaussWindowFunc::calcWinFunc(xnorm,ynorm,znorm,windowFuncs[i]);


      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;



      windowFuncs[i] *= normFactor;

    }

    break;

  case 5:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      TenthOrderSplineWinFunc::calcWinFunc(xnorm,ynorm,znorm,
					   windowFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;



      windowFuncs[i] *= normFactor;

    }

    break;

  case 6:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuarticSplineWindowFunc2::calcWinFunc(xnorm,ynorm,znorm,
					    windowFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;



      windowFuncs[i] *= normFactor;

    }

    break;

  case 7:

    if(radiusMode == 5) {
      cerr<<"Chosen window function type does not support influence\n "
	  <<"radius determination mode "<<radiusMode<<"!"<<endl;       
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = x - ptcls[m].getCoord(0);
      ynorm = y - ptcls[m].getCoord(1);
      znorm = z - ptcls[m].getCoord(2);

      // determine max particle distance within current particle's (m) 
      // influence zone
      maxDist = 0;

      for(int j=0;j<supportSize;j++) {
	n = sPtcls[j];

	xDist = ptcls[n].getCoord(0) - ptcls[m].getCoord(0);
	yDist = ptcls[n].getCoord(1) - ptcls[m].getCoord(1);
	zDist = ptcls[n].getCoord(2) - ptcls[m].getCoord(2);

	dist = 2.0*sqrt(pow(xDist,2.0)+pow(yDist,2.0)+pow(zDist,2.0));

	if(dist > maxDist)
	  maxDist = dist;
      }

      GaussWindowFunc2::calcWinFunc(xnorm,ynorm,znorm,maxDist,
				    ptcls[m],windowFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      windowFuncs[i] *= normFactor;

    }

    break;

  case 8:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuinticSplineWindowFunc::calcWinFunc(xnorm,ynorm,znorm,
					   windowFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      windowFuncs[i] *= normFactor;

    }

    break;

  default:
    cerr<<"Chosen window function type isn't supported!"<<endl;       
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of first order derivation of the window functions 
void WindowFunctionSet::calcWinFunction1stDerivs(InputFileData* InputData,
						 std::vector<Particle>& ptcls,
						 intVector& sPtcls,
						 double& x,double& y,
						 double& z,
						 int& supportSize,
						 std::ofstream& logFile) {

  using namespace std;

  int winFuncNorming = (int)InputData->getValue("windowfunctionNorming");
  int wType = (int)InputData->getValue("windowFunctionType");
  int particlesNum = ptcls.size();

  double dxnormFactor,dynormFactor,dznormFactor,normFactor,xnorm,ynorm,
    znorm;
  int m,n;

  double xDist,yDist,zDist,dist,maxDist;

  if(xDerivWinFuncs.size() < supportSize) {
    xDerivWinFuncs.resize(supportSize);
    yDerivWinFuncs.resize(supportSize);
    zDerivWinFuncs.resize(supportSize);
  }

  // Choose a window function type.
  switch(wType) {

    // constant 1:
  case 0:
    break;

    // cubic spline:
  case 1:

    for(int i=0;i<supportSize;i++) {

      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      CubicSplineWindowFunc::calcWinFunc1stDerivs(xnorm,ynorm,znorm,
						    xDerivWinFuncs[i],
						    yDerivWinFuncs[i],
						    zDerivWinFuncs[i]);
      
      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  case 2:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuarticSplineWindowFunc::calcWinFunc1stDerivs(xnorm,ynorm,znorm,
						    xDerivWinFuncs[i],
						    yDerivWinFuncs[i],
						    zDerivWinFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  case 3:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      ExponentialWindowFunc::calcWinFunc1stDerivs(xnorm,ynorm,znorm,
						  xDerivWinFuncs[i],
						  yDerivWinFuncs[i],
						  zDerivWinFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  case 4:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      GaussWindowFunc::calcWinFunc1stDerivs(xnorm,ynorm,znorm,
					    xDerivWinFuncs[i],
					    yDerivWinFuncs[i],
					    zDerivWinFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  case 5:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      TenthOrderSplineWinFunc::calcWinFunc1stDerivs(xnorm,ynorm,znorm,
						    xDerivWinFuncs[i],
						    yDerivWinFuncs[i],
						    zDerivWinFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  case 6:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuarticSplineWindowFunc2::calcWinFunc1stDerivs(xnorm,ynorm,znorm,
						     xDerivWinFuncs[i],
						     yDerivWinFuncs[i],
						     zDerivWinFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  case 7:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = x - ptcls[m].getCoord(0);
      ynorm = y - ptcls[m].getCoord(1);
      znorm = z - ptcls[m].getCoord(2);

      // determine max particle distance within current particle's (m) 
      // influence zone
      maxDist = 0;

      for(int j=0;j<supportSize;j++) {
	n = sPtcls[j];

	xDist = ptcls[n].getCoord(0) - ptcls[m].getCoord(0);
	yDist = ptcls[n].getCoord(1) - ptcls[m].getCoord(1);
	zDist = ptcls[n].getCoord(2) - ptcls[m].getCoord(2);

	dist = 2.0*sqrt(pow(xDist,2.0)+pow(yDist,2.0)+pow(zDist,2.0));

	if(dist > maxDist)
	  maxDist = dist;

      }

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;
      GaussWindowFunc2::calcWinFunc1stDerivs(xnorm,ynorm,znorm,maxDist,
					     ptcls[m],xDerivWinFuncs[i],
					     yDerivWinFuncs[i],
					     zDerivWinFuncs[i]);


      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;

      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  case 8:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuinticSplineWindowFunc::calcWinFunc1stDerivs(xnorm,ynorm,znorm,
						    xDerivWinFuncs[i],
						    yDerivWinFuncs[i],
						    zDerivWinFuncs[i]);

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxnormFactor = normFactor/ptcls[m].getRadius(0);
      dynormFactor = normFactor/ptcls[m].getRadius(1);
      dznormFactor = normFactor/ptcls[m].getRadius(2);

      xDerivWinFuncs[i] *= dxnormFactor;
      yDerivWinFuncs[i] *= dynormFactor;
      zDerivWinFuncs[i] *= dznormFactor;

    }

    break;

  default:
    cerr<<"Chosen window function type isn't supported!"<<endl;       
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}


/***********************************************************************/
/***********************************************************************/
// Calculation of second order derivation of the window functions
void WindowFunctionSet::calcWinFunction2ndDerivs(InputFileData* InputData,
						 std::vector<Particle>& ptcls,
						 intVector& sPtcls,
						 double& x,double& y,
						 double& z,
						 int& supportSize,
						 std::ofstream& logFile) {

  using namespace std;

  int winFuncNorming = (int)InputData->getValue("windowfunctionNorming");
  int wType = (int)InputData->getValue("windowFunctionType");
  int particlesNum = ptcls.size();


  double dxxnormFactor,dyynormFactor,dzznormFactor,dxynormFactor,
    dyznormFactor,dzxnormFactor,normFactor,xnorm,ynorm,znorm;
  int m,n;

  double xDist,yDist,zDist,dist,maxDist;

  if(xxDerivWinFuncs.size() < supportSize) {
    xxDerivWinFuncs.resize(supportSize);
    yyDerivWinFuncs.resize(supportSize);
    zzDerivWinFuncs.resize(supportSize);
    
    xyDerivWinFuncs.resize(supportSize);
    yzDerivWinFuncs.resize(supportSize);
    zxDerivWinFuncs.resize(supportSize);
  }

  // Choose a window function type.
  switch(wType) {

    // constant 1:
  case 0:
    break;

    // cubic spline:
  case 1:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      CubicSplineWindowFunc::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
						  xxDerivWinFuncs[i],
						  yyDerivWinFuncs[i],
						  zzDerivWinFuncs[i],
						  xyDerivWinFuncs[i],
						  yzDerivWinFuncs[i],
						  zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  case 2:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuarticSplineWindowFunc::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
						    xxDerivWinFuncs[i],
						    yyDerivWinFuncs[i],
						    zzDerivWinFuncs[i],
						    xyDerivWinFuncs[i],
						    yzDerivWinFuncs[i],
						    zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  case 3:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      ExponentialWindowFunc::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
						  xxDerivWinFuncs[i],
						  yyDerivWinFuncs[i],
						  zzDerivWinFuncs[i],
						  xyDerivWinFuncs[i],
						  yzDerivWinFuncs[i],
						  zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  case 4:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;
      if(winFuncNorming)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));
      else
	normFactor = 1.0;

      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      GaussWindowFunc::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
					    xxDerivWinFuncs[i],
					    yyDerivWinFuncs[i],
					    zzDerivWinFuncs[i],
					    xyDerivWinFuncs[i],
					    yzDerivWinFuncs[i],
					    zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  case 5:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      TenthOrderSplineWinFunc::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
						    xxDerivWinFuncs[i],
						    yyDerivWinFuncs[i],
						    zzDerivWinFuncs[i],
						    xyDerivWinFuncs[i],
						    yzDerivWinFuncs[i],
						    zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  case 6:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuarticSplineWindowFunc2::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
						     xxDerivWinFuncs[i],
						     yyDerivWinFuncs[i],
						     zzDerivWinFuncs[i],
						     xyDerivWinFuncs[i],
						     yzDerivWinFuncs[i],
						     zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  case 7:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = x - ptcls[m].getCoord(0);
      ynorm = y - ptcls[m].getCoord(1);
      znorm = z - ptcls[m].getCoord(2);

      // determine max particle distance within current particle's (m) 
      // influence zone
      maxDist = 0;

      for(int j=0;j<supportSize;j++) {
	n = sPtcls[j];

	xDist = ptcls[n].getCoord(0) - ptcls[m].getCoord(0);
	yDist = ptcls[n].getCoord(1) - ptcls[m].getCoord(1);
	zDist = ptcls[n].getCoord(2) - ptcls[m].getCoord(2);

	dist = 2.0*sqrt(pow(xDist,2.0)+pow(yDist,2.0)+pow(zDist,2.0));

	if(dist > maxDist)
	  maxDist = dist;

      }

      GaussWindowFunc2::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
					     maxDist,ptcls[m],
					     xxDerivWinFuncs[i],
					     yyDerivWinFuncs[i],
					     zzDerivWinFuncs[i],
					     xyDerivWinFuncs[i],
					     yzDerivWinFuncs[i],
					     zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  case 8:

    for(int i=0;i<supportSize;i++) {
      m = sPtcls[i];

      // norming by 3D domain of influence
      if(winFuncNorming == 1)
	normFactor = 1.0/(ptcls[m].getRadius(0)*ptcls[m].getRadius(1)*
			  ptcls[m].getRadius(2));

      // norming by particle weight obtained from FEM mesh
      else if(winFuncNorming == 2)
	normFactor = 1.0/ptcls[m].getWeight();

      else
	normFactor = 1.0;


      dxxnormFactor = normFactor/pow(ptcls[m].getRadius(0),2);

      dyynormFactor = normFactor/pow(ptcls[m].getRadius(1),2);

      dzznormFactor = normFactor/pow(ptcls[m].getRadius(2),2);

      dxynormFactor = normFactor/(ptcls[m].getRadius(0)*
				 ptcls[m].getRadius(1));

      dyznormFactor = normFactor/(ptcls[m].getRadius(1)*
				 ptcls[m].getRadius(2));

      dzxnormFactor = normFactor/(ptcls[m].getRadius(2)*
				 ptcls[m].getRadius(0));

      xnorm = (x - ptcls[m].getCoord(0))/ptcls[m].getRadius(0);
      ynorm = (y - ptcls[m].getCoord(1))/ptcls[m].getRadius(1);
      znorm = (z - ptcls[m].getCoord(2))/ptcls[m].getRadius(2);

      QuinticSplineWindowFunc::calcWinFunc2ndDerivs(xnorm,ynorm,znorm,
						    xxDerivWinFuncs[i],
						    yyDerivWinFuncs[i],
						    zzDerivWinFuncs[i],
						    xyDerivWinFuncs[i],
						    yzDerivWinFuncs[i],
						    zxDerivWinFuncs[i]);

      xxDerivWinFuncs[i] *= dxxnormFactor;
      yyDerivWinFuncs[i] *= dyynormFactor;
      zzDerivWinFuncs[i] *= dzznormFactor;

      xyDerivWinFuncs[i] *= dxynormFactor;
      yzDerivWinFuncs[i] *= dyznormFactor;
      zxDerivWinFuncs[i] *= dzxnormFactor;

    }

    break;

  default:
    cerr<<"Chosen window function type isn't supported!"<<endl;       
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
