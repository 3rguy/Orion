#include "ShepardShapeFunc.h"


ShepardShapeFunc::ShepardShapeFunc(InputFileData* InputData,
				   int& supportSize,
				   intVector& sPtcls,
				   std::vector<Particle>& ptcls,
				   double& x,double& y,double& z,
				   unsigned int derivationOrder,
				   std::map<std::string,double>& modelData,  
				   std::ofstream& logFile,
				   PetscViewer& viewerSEQ) {

  using namespace std;

  switch(derivationOrder) {

    // Calculate at a point for all its supporting particles their 
    // shape functions.
  case 0: 

    calcShapes(InputData,supportSize,sPtcls,ptcls,x,y,z,
	       modelData,logFile,viewerSEQ);

    break;
  

  case 1:
    
    // Calculate at a point for all its supporting particles their 
    // shape functions and their first order derivations.
    calcShape1stDerivs(InputData,supportSize,sPtcls,ptcls,x,y,z,
		       modelData,logFile,viewerSEQ);

    break;


  case 2:
 
    // Calculate at a point for all its supporting particles their 
    // shape functions and their first and second order derivations.
    calcShape2ndDerivs(InputData,supportSize,sPtcls,ptcls,x,y,z,
		       modelData,logFile,viewerSEQ);

    break;

  default:
    logFile<<"In ShepardShapeFunc::ShepardShapeFunc derivation order\n"
	   <<" is not supported!"<<endl;
    cerr<<"In ShepardShapeFunc::ShepardShapeFunc derivation order\n"
	<<" is not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions.
void ShepardShapeFunc::calcShapes(InputFileData* InputData,
				  int& supportSize,
				  intVector& sPtcls,
				  std::vector<Particle>& ptcls,
				  double& x,double& y,double& z,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile,
				  PetscViewer& viewerSEQ) {

  using namespace std;

  int particlesNum = ptcls.size();

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new WindowFunctionSet(InputData,ptcls,
						      sPtcls,x,y,z,
						      supportSize,0,
						      modelData,
						      logFile);


  dbVector& W = WFuncSet->getWindowFuncs();

#ifdef _geometryDebugMode_
  double xrad,yrad,zrad;
  logFile<<"##########################################################"<<endl;
  logFile<<"***************** shepard function set **************"<<endl;
  logFile<<"point "<<x<<" "<<y<<" "<<z<<endl;
  logFile<<"SupportSize "<<supportSize<<endl;
  for(int k=0;k<supportSize;k++) {
    logFile<<"Ptcle "<<sPtcls[k]<<": "<<ptcls[sPtcls[k]].getCoord(0)<<" "
	   <<ptcls[sPtcls[k]].getCoord(1)<<" "
	   <<ptcls[sPtcls[k]].getCoord(2)<<endl;
    logFile<<"Radii "<<ptcls[sPtcls[k]].getRadius(0)<<" "
	   <<ptcls[sPtcls[k]].getRadius(1)<<" "
	   <<ptcls[sPtcls[k]].getRadius(2)<<endl;
    xrad = x-ptcls[sPtcls[k]].getCoord(0);
    yrad = y-ptcls[sPtcls[k]].getCoord(1);
    zrad = z-ptcls[sPtcls[k]].getCoord(2);
    logFile<<"rad "<<xrad<<" "<<yrad<<" "<<zrad<<" ";
    logFile<<"distance "<<sqrt(pow(xrad,2)+pow(yrad,2)+pow(zrad,2))<<endl;
  }
  logFile<<"**************** window functions ***************"<<endl;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
#endif

  double sumW = 0;
  
  for(int i=0;i<supportSize;i++)
    sumW += W[i];


  /********************************************************************/
  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  shapefunctions.resize(supportSize);

  for(int i=0;i<supportSize;i++) {

    shapefunctions[i] = W[i]/sumW;
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************ Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapefunctions[i]<<endl;
#endif
  
  delete WFuncSet;
}

/**********************************************************************/
/**********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions.
void ShepardShapeFunc::calcShapes(InputFileData* InputData,
				  int& supportSize,
				  intVector& sPtcls,
				  std::vector<Particle>& ptcls,
				  double& x,double& y,double& z,
				  dbVector& shapes,
				  int& basisTermNum,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile,
				  PetscViewer& viewerSEQ) {

  basisTermNum = 1;

  calcShapes(InputData,supportSize,sPtcls,ptcls,x,y,z,modelData,logFile,
	     viewerSEQ);

  shapes = shapefunctions;
}

/**********************************************************************/
/**********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first order derivations.
void ShepardShapeFunc::calcShape1stDerivs(InputFileData* InputData,
					  int& supportSize,
					  intVector& sPtcls,
					  std::vector<Particle>& ptcls,
					  double& x,double& y,double& z,
					  std::map<std::string,double>& modelData,
					  std::ofstream& logFile,
					  PetscViewer& viewerSEQ) {
  
  using namespace std;

  int particlesNum = ptcls.size();

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new WindowFunctionSet(InputData,ptcls,
						      sPtcls,x,y,z,
						      supportSize,1,
						      modelData,
						      logFile);


  dbVector& W = WFuncSet->getWindowFuncs();
  dbVector& dWx = WFuncSet->getXDerivWinFuncs();
  dbVector& dWy = WFuncSet->getYDerivWinFuncs();
  dbVector& dWz = WFuncSet->getZDerivWinFuncs();

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"***************** shepard function set **************"<<endl;
  double xrad,yrad,zrad;
  logFile<<"point "<<x<<" "<<y<<" "<<z<<endl;
  logFile<<"SupportSize "<<supportSize<<endl;
  for(int k=0;k<supportSize;k++) {
    logFile<<"Ptcle "<<sPtcls[k]<<": "<<ptcls[sPtcls[k]].getCoord(0)<<" "
	   <<ptcls[sPtcls[k]].getCoord(1)<<" "
	   <<ptcls[sPtcls[k]].getCoord(2)<<endl;
    logFile<<"Radii "<<ptcls[sPtcls[k]].getRadius(0)<<" "
	   <<ptcls[sPtcls[k]].getRadius(1)<<" "
	   <<ptcls[sPtcls[k]].getRadius(2)<<endl;
    xrad = x-ptcls[sPtcls[k]].getCoord(0);
    yrad = y-ptcls[sPtcls[k]].getCoord(1);
    zrad = z-ptcls[sPtcls[k]].getCoord(2);
    logFile<<"rad "<<xrad<<" "<<yrad<<" "<<zrad<<" ";
    logFile<<"distance "<<sqrt(pow(xrad,2)+pow(yrad,2)+pow(zrad,2))<<endl;
  }
  logFile<<"**************** window functions ***************"<<endl;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
  logFile<<"*** dWx ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWx[i]<<endl;
  logFile<<"*** dWy ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWy[i]<<endl;
  logFile<<"*** dWz ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWz[i]<<endl;
#endif

  double sumW = 0;
  
  for(int i=0;i<supportSize;i++)
    sumW += W[i];


  /********************************************************************/
  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  shapefunctions.resize(supportSize);

  for(int i=0;i<supportSize;i++) {

    shapefunctions[i] = W[i]/sumW;
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************ Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapefunctions[i]<<endl;
#endif

  /********************************************************************/
  double sumW2 = pow(sumW,2);

  double sumDxW = 0;
  double sumDyW = 0;
  double sumDzW = 0;

  for(int i=0;i<supportSize;i++) {
    sumDxW += dWx[i];
    sumDyW += dWy[i];
    sumDzW += dWz[i];
  }
  
  /********************************************************************/
  // calculation of first order derivations
  firstDerivShapes.resize(3);

  for(int i=0;i<firstDerivShapes.size();i++)
    firstDerivShapes[i].resize(supportSize);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i=0;i<supportSize;i++) {
    
    firstDerivShapes[0][i] = (dWx[i]*sumW-W[i]*sumDxW)/sumW2;
    firstDerivShapes[1][i] = (dWy[i]*sumW-W[i]*sumDyW)/sumW2;
    firstDerivShapes[2][i] = (dWz[i]*sumW-W[i]*sumDzW)/sumW2;
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************** first order derivations ***************"<<endl;
  for(int i=0;i<supportSize;i++) {
    logFile<<i<<": "<<sPtcls[i]<<" "<<firstDerivShapes[0][i]<<" "
	   <<firstDerivShapes[1][i]<<" "
	   <<firstDerivShapes[2][i]<<endl;
  }
#endif


  delete WFuncSet;
}


void ShepardShapeFunc::calcShapes(InputFileData* InputData,
				  int& supportSize,
				  intVector& sPtcls,
				  std::vector<Particle>& ptcls,
				  double& x,double& y,double& z,
				  dbVector& shapes,
				  dbMatrix& dShapes,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile,
				  PetscViewer& viewerSEQ) {

  using namespace std;

  calcShape1stDerivs(InputData,supportSize,sPtcls,ptcls,x,y,z,
		     modelData,logFile,viewerSEQ);

  shapes = shapefunctions;
  dShapes = firstDerivShapes;

}

/**********************************************************************/
/**********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first and second order derivations.
void ShepardShapeFunc::calcShape2ndDerivs(InputFileData* InputData,
					  int& supportSize,
					  intVector& sPtcls,
					  std::vector<Particle>& ptcls,
					  double& x,double& y,double& z,
					  std::map<std::string,double>& modelData,
					  std::ofstream& logFile,
					  PetscViewer& viewerSEQ) {

  using namespace std;

  int particlesNum = ptcls.size();

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new WindowFunctionSet(InputData,ptcls,
						      sPtcls,x,y,z,
						      supportSize,2,
						      modelData,
						      logFile);

  dbVector& W = WFuncSet->getWindowFuncs();

  dbVector& dWx = WFuncSet->getXDerivWinFuncs();
  dbVector& dWy = WFuncSet->getYDerivWinFuncs();
  dbVector& dWz = WFuncSet->getZDerivWinFuncs();

  dbVector& dWxx = WFuncSet->getXXDerivWinFuncs();
  dbVector& dWyy = WFuncSet->getYYDerivWinFuncs();
  dbVector& dWzz = WFuncSet->getZZDerivWinFuncs();

  dbVector& dWxy = WFuncSet->getXYDerivWinFuncs();
  dbVector& dWyz = WFuncSet->getYZDerivWinFuncs();
  dbVector& dWzx = WFuncSet->getZXDerivWinFuncs();

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"***************** shepard function set **************"<<endl;
  double xrad,yrad,zrad;
  logFile<<"point "<<x<<" "<<y<<" "<<z<<endl;
  logFile<<"SupportSize "<<supportSize<<endl;
  for(int k=0;k<supportSize;k++) {
    logFile<<"Ptcle "<<sPtcls[k]<<": "<<ptcls[sPtcls[k]].getCoord(0)<<" "
	   <<ptcls[sPtcls[k]].getCoord(1)<<" "
	   <<ptcls[sPtcls[k]].getCoord(2)<<endl;
    logFile<<"Radii "<<ptcls[sPtcls[k]].getRadius(0)<<" "
	   <<ptcls[sPtcls[k]].getRadius(1)<<" "
	   <<ptcls[sPtcls[k]].getRadius(2)<<endl;
    xrad = x-ptcls[sPtcls[k]].getCoord(0);
    yrad = y-ptcls[sPtcls[k]].getCoord(1);
    zrad = z-ptcls[sPtcls[k]].getCoord(2);
    logFile<<"rad "<<xrad<<" "<<yrad<<" "<<zrad<<" ";
    logFile<<"distance "<<sqrt(pow(xrad,2)+pow(yrad,2)+pow(zrad,2))<<endl;
  }
  logFile<<"**************** window functions ***************"<<endl;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
  logFile<<"*** dWx ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWx[i]<<endl;
  logFile<<"*** dWy ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWy[i]<<endl;
  logFile<<"*** dWz ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWz[i]<<endl;
  logFile<<"*** dWxx ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWxx[i]<<endl;
  logFile<<"*** dWyy ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWyy[i]<<endl;
  logFile<<"*** dWzz ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWzz[i]<<endl;
  logFile<<"*** dWxy ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWxy[i]<<endl;
  logFile<<"*** dWyz ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWyz[i]<<endl;
  logFile<<"*** dWzx ****"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<dWzx[i]<<endl;
#endif

  double sumW = 0;
  
  for(int i=0;i<supportSize;i++)
    sumW += W[i];


  /********************************************************************/
  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  shapefunctions.resize(supportSize);

  for(int i=0;i<supportSize;i++) {

    shapefunctions[i] = W[i]/sumW;
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************ Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapefunctions[i]<<endl;
#endif

  /********************************************************************/
  // calculation of first order derivations
  double sumW2 = pow(sumW,2);

  double sumDxW = 0;
  double sumDyW = 0;
  double sumDzW = 0;

  for(int i=0;i<supportSize;i++) {
    sumDxW += dWx[i];
    sumDyW += dWy[i];
    sumDzW += dWz[i];
  }

  firstDerivShapes.resize(3);

  for(int i=0;i<firstDerivShapes.size();i++)
    firstDerivShapes[i].resize(supportSize);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i=0;i<supportSize;i++) {
    
    firstDerivShapes[0][i] = (dWx[i]*sumW-W[i]*sumDxW)/sumW2;
    firstDerivShapes[1][i] = (dWy[i]*sumW-W[i]*sumDyW)/sumW2;
    firstDerivShapes[2][i] = (dWz[i]*sumW-W[i]*sumDzW)/sumW2;
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************** first order derivations ***************"<<endl;
  for(int i=0;i<supportSize;i++) {
    logFile<<i<<" "<<sPtcls[i]<<": "<<firstDerivShapes[0][i]<<" "
	   <<firstDerivShapes[1][i]<<" "
	   <<firstDerivShapes[2][i]<<endl;
  }
#endif
  
  /********************************************************************/
  // calculation of second order derivations
  double sumW3 = pow(sumW,3);

  double sumDxW2 = pow(sumDxW,2);
  double sumDyW2 = pow(sumDyW,2);
  double sumDzW2 = pow(sumDzW,2);

  double sumDxxW = 0;
  double sumDyyW = 0;
  double sumDzzW = 0;

  double sumDxyW = 0;
  double sumDyzW = 0;
  double sumDzxW = 0;

  for(int i=0;i<supportSize;i++) {
    sumDxxW += dWxx[i];
    sumDyyW += dWyy[i];
    sumDzzW += dWzz[i];
    
    sumDxyW += dWxy[i];
    sumDyzW += dWyz[i];
    sumDzxW += dWzx[i];
  }

  secondDerivShapes.resize(6);

  for(int i=0;i<secondDerivShapes.size();i++)
    secondDerivShapes[i].resize(supportSize);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i=0;i<supportSize;i++) {

    secondDerivShapes[0][i] = (dWxx[i]*sumW2-W[i]*sumW*sumDxxW-
			       2*dWx[i]*sumW*sumDxW+2*W[i]*sumDxW2)/sumW3;

    secondDerivShapes[1][i] = (dWyy[i]*sumW2-W[i]*sumW*sumDyyW-
			       2*dWy[i]*sumW*sumDyW+2*W[i]*sumDyW2)/sumW3;

    secondDerivShapes[2][i] = (dWzz[i]*sumW2-W[i]*sumW*sumDzzW-
			       2*dWz[i]*sumW*sumDzW+2*W[i]*sumDzW2)/sumW3;

    secondDerivShapes[3][i] = 
      (dWxy[i]*sumW2-dWx[i]*sumW*sumDyW-
       dWy[i]*sumW*sumDxW+2*W[i]*sumDxW*sumDyW-W[i]*sumW*sumDxyW)/sumW3;

    secondDerivShapes[4][i] = 
      (dWyz[i]*sumW2-dWy[i]*sumW*sumDzW-
       dWz[i]*sumW*sumDyW+2*W[i]*sumDyW*sumDzW-W[i]*sumW*sumDyzW)/sumW3;

    secondDerivShapes[5][i] = 
      (dWzx[i]*sumW2-dWz[i]*sumW*sumDxW-
       dWx[i]*sumW*sumDzW+2*W[i]*sumDzW*sumDxW-W[i]*sumW*sumDzxW)/sumW3;
  }


#ifdef _geometryDebugMode_
    logFile<<"********** second order derivations **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<": "<<secondDerivShapes[0][i]<<" "
	     <<secondDerivShapes[1][i]<<" "
	     <<secondDerivShapes[2][i]<<" "
	     <<secondDerivShapes[3][i]<<" "
	     <<secondDerivShapes[4][i]<<" "
	     <<secondDerivShapes[5][i]<<endl;
  double pum = 0;
  double pumDx = 0;
  double pumDy = 0;
  double pumDz = 0;
  double pumDxx = 0;
  double pumDyy = 0;
  double pumDzz = 0;
  double pumDxy = 0;
  double pumDyz = 0;
  double pumDzx = 0;
  for(int i=0;i<supportSize;i++) {
    pum += shapefunctions[i];
    pumDx += firstDerivShapes[0][i];
    pumDy += firstDerivShapes[1][i];
    pumDz += firstDerivShapes[2][i];
    pumDxx += secondDerivShapes[0][i];
    pumDyy += secondDerivShapes[1][i];
    pumDzz += secondDerivShapes[2][i];
    pumDxy += secondDerivShapes[3][i];
    pumDyz += secondDerivShapes[4][i];
    pumDzx += secondDerivShapes[5][i];
  }
  logFile<<"pum = "<<pum<<endl;
  logFile<<"pumDx = "<<pumDx<<endl;
  logFile<<"pumDy = "<<pumDy<<endl;
  logFile<<"pumDz = "<<pumDz<<endl;
  logFile<<"pumDxx = "<<pumDxx<<endl;
  logFile<<"pumDyy = "<<pumDyy<<endl;
  logFile<<"pumDzz = "<<pumDzz<<endl;
  logFile<<"pumDxy = "<<pumDxy<<endl;
  logFile<<"pumDyz = "<<pumDyz<<endl;
  logFile<<"pumDzx = "<<pumDzx<<endl;
#endif

  delete WFuncSet;
}

void ShepardShapeFunc::calcShapes(InputFileData* InputData,
				  int& supportSize,
				  intVector& sPtcls,
				  std::vector<Particle>& ptcls,
				  double& x,double& y,double& z,
				  dbVector& shapes,
				  dbMatrix& dShapes,
				  dbMatrix& d2Shapes,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile, 
				  PetscViewer& viewerSEQ) {

  using namespace std;

  calcShape2ndDerivs(InputData,supportSize,sPtcls,ptcls,x,y,z,
		     modelData,logFile,viewerSEQ);

  shapes = shapefunctions;
  dShapes = firstDerivShapes;
  d2Shapes = secondDerivShapes;

}
