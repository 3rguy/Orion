#include "AsymShapeFunc.h"

// Calculate at a point for all its supporting particles their 
// shape functions.
void AsymShapeFunc::calcShapes(InputFileData* InputData,
			       int& supportSize,
			       intVector& sPtcls,
			       std::vector<Particle>& ptcls,
			       double& x,double& y,double& z,
			       dbVector& shapeFuncs,
			       int& basisTermNum,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerSEQ) {

  using namespace std;

  InputData->setValue("windowFunctionType",4.0);
  // InputData->getValue("windowfunctionNorming",2.0);

  int usedDims = (int)modelData["usedDimensions"];
  bool normedBasis = (bool)InputData->getValue("normedBasisPolynom");
  int wType = (int)InputData->getValue("windowFunctionType");


  if((int)InputData->getValue("radiusDeterminationAlgorithm") != 6) {
    logFile<<"In AsymShapeFunc::calcShapes input variable\n" 
	   <<"'radiusDeterminationAlgorithm' needs to be set to '6'."
	   <<endl;
    cerr<<"In AsymShapeFunc::calcShapes input variable\n" 
	   <<"'radiusDeterminationAlgorithm' needs to be set to '6'."
	<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

#ifdef _geometryDebugMode_
  double xrad,yrad,zrad;
  logFile<<"##########################################################"<<endl;
  logFile<<"point "<<x<<" "<<y<<" "<<z<<endl;
  logFile<<"SupportSize "<<supportSize<<endl;
  for(int k=0;k<supportSize;k++) {
    logFile<<"Ptcle "<<sPtcls[k]<<": "<<ptcls[sPtcls[k]].getCoord(0)<<" "
	   <<ptcls[sPtcls[k]].getCoord(1)<<" "
	   <<ptcls[sPtcls[k]].getCoord(2)<<endl;
    logFile<<"Radii "<<ptcls[sPtcls[k]].getRadius(0)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims)<<" "
	   <<ptcls[sPtcls[k]].getRadius(1)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims+1)<<" "
	   <<ptcls[sPtcls[k]].getRadius(2)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims+2)<<endl;
    xrad = x-ptcls[sPtcls[k]].getCoord(0);
    yrad = y-ptcls[sPtcls[k]].getCoord(1);
    zrad = z-ptcls[sPtcls[k]].getCoord(2);
    logFile<<"rad "<<xrad<<" "<<yrad<<" "<<zrad<<" ";
    logFile<<"distance "<<sqrt(pow(xrad,2)+pow(yrad,2)+pow(zrad,2))<<endl;
  }
#endif

  // The size of the linear equation system.
  int linEQSize;
  int particlesNum = ptcls.size();

  // Calculate a polynom base set and its derivatives for 
  // all particles.
  //InputData->setValue("normedBasisPolynom",1);

  int derivOrder = 0;

  BasisPolynom* PolynomSet = new BasisPolyAsym(InputData,ptcls,sPtcls,
					       x,y,z,supportSize,
					       linEQSize,derivOrder,
					       modelData,logFile);
  
  dbMatrix& P = PolynomSet->getBasis();
  basisTermNum = linEQSize;

  // Check the particle support
  if(linEQSize > supportSize) {
    logFile <<"In AsymShapeFunc::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.

  WindowFunctionSet* WFuncSet;

  WFuncSet = new WindowFuncAsym(InputData,ptcls,
				sPtcls,x,y,z,
				supportSize,
				derivOrder,
				modelData,
				logFile,viewerSEQ);



  dbVector& W = WFuncSet->getWindowFuncs();

#ifdef _geometryDebugMode_
  logFile<<"LinEQSize "<<linEQSize<<endl;
  logFile<<"*************** basis polynom set ****************"<<endl;
  logFile<<"Whole size "<<P.size()<<endl;
  for(int i=0;i<supportSize;i++) {
    logFile<<"Ptcle "<<sPtcls[i]<<" Size "<<P[i].size()<<endl;
    
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<P[i][j]<<endl;
  }
  logFile<<"**************** window functions ***************"<<endl;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
#endif

  /**********************************************************************/
  // make use a shifted basis and continuous approximation set-up
  if(normedBasis) {

    // Assembly the momentum matrix 'M'
    Mat M;
    intVector matIdx(linEQSize);
    dbMatrix MValues(linEQSize,dbVector(linEQSize));

    map<string,bool> mOptions;
    createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);

    for(int i=0;i<linEQSize;i++) {
      matIdx[i] = i;

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  MValues[i][j] = MValues[i][j]+P[k][i]*P[k][j]*W[k];

    }

#ifdef _geometryDebugMode_
    logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<MValues[i][j];
      logFile<<endl;
    }
#endif

    for(int i=0;i<linEQSize;i++)
      MatSetValues(M,1,&i,linEQSize,&matIdx[0],&MValues[i][0],
		   INSERT_VALUES);
	
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
    MatSetOption(M,MAT_SYMMETRIC,PETSC_TRUE);
    //  MatView(M,viewerSEQ);

    // Assemble vector P(0).
    Vec P0;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&P0);
    intVector vecIdx(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      vecIdx[i] = i;
    }

    VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
    // Create vector B.
    Vec B;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&B);
  
    // Solve the linear equation system 'MB = P(0)' and get vector B.
    KSP ksp;
    PC pc;
    Mat F;

    KSPCreate(MPI_COMM_SELF,&ksp);
    KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetUp(ksp);

  
    KSPSolve(ksp,P0,B);

    // Check whether the solving was successful or not.
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);

    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    //  VecView(B,viewerSEQ);

    destroyPETScVec(P0);

    // Calculate the shape functions.
    double* BCalc;
    VecGetArray(B,&BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      shapeFuncs[i] = 0;

      for(int j=0;j<linEQSize;j++)
	shapeFuncs[i] += P[i][j]*BCalc[j];

      shapeFuncs[i] = shapeFuncs[i]*W[i];
    }
  
#ifdef _geometryDebugMode_
    logFile<<"************ Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif
  
    // Destroy all petsc objects.
    destroyPETScSolver(ksp);

    VecRestoreArray(B,&BCalc);
    destroyPETScVec(B);
    destroyPETScMat(M);
  }

  /**********************************************************************/
  // make use of the conventional MLS approximation scheme
  else {

    // Assembly the momentum matrix 'M'.
    Mat M;
    intVector matIdx(linEQSize);
    dbMatrix MValues(linEQSize,dbVector(linEQSize));

    map<string,bool> mOptions;
    createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);

    for(int i=0;i<linEQSize;i++) {
      matIdx[i] = i;

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  MValues[i][j] = MValues[i][j]+P[k][i]*P[k][j]*W[k];

    }

#ifdef _geometryDebugMode_
    logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<MValues[i][j];
      logFile<<endl;
    }
#endif

    for(int i=0;i<linEQSize;i++)
      MatSetValues(M,1,&i,linEQSize,&matIdx[0],&MValues[i][0],
		   INSERT_VALUES);
	
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
    MatSetOption(M,MAT_SYMMETRIC,PETSC_TRUE);
    //  MatView(M,viewerSEQ);

    // Assemble vector P(0).
    Vec P0;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&P0);
    intVector vecIdx(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      vecIdx[i] = i;
    }

    VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
    // Create vector B.
    Vec B;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&B);
  
    // Solve the linear equation system 'MB = P(0)' and get vector B.
    PC pc;
    KSP ksp;
    Mat F;
    
    KSPCreate(MPI_COMM_SELF,&ksp);
    KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetUp(ksp);

  
    KSPSolve(ksp,P0,B);

    // Check whether the solving was successful or not.
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);

    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    //  VecView(B,viewerSEQ);

    destroyPETScVec(P0);

    // Calculate the shape functions.
    double* BCalc;
    VecGetArray(B,&BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      shapeFuncs[i] = 0;

      for(int j=0;j<linEQSize;j++)
	shapeFuncs[i] += P[i][j]*BCalc[j];

      shapeFuncs[i] = shapeFuncs[i]*W[i];
    }
  
#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<BCalc[i]<<endl;
    logFile<<"************ Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif

    // Destroy all petsc objects.
    destroyPETScSolver(ksp);

    VecRestoreArray(B,&BCalc);
    destroyPETScVec(B);
    destroyPETScMat(M);

  }

  delete PolynomSet,WFuncSet;
}

/************************************************************************/
/************************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first order derivations.
void AsymShapeFunc::calcShapes(InputFileData* InputData,
			       int& supportSize,
			       intVector& sPtcls,
			       std::vector<Particle>& ptcls,
			       double& x,double& y,double& z,
			       dbVector& shapeFuncs,
			       dbVector& xDerivShapes,
			       dbVector& yDerivShapes,
			       dbVector& zDerivShapes,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerSEQ) {

  using namespace std;

  InputData->setValue("windowFunctionType",4.0);
  //  InputData->getValue("windowfunctionNorming",2.0);

  int usedDims = (int)modelData["usedDimensions"];
  bool normedBasis = (bool)InputData->getValue("normedBasisPolynom");
  int wType = (int)InputData->getValue("windowFunctionType");

  if((int)InputData->getValue("radiusDeterminationAlgorithm") != 6) {
    logFile<<"In AsymShapeFunc::calcShapes input variable\n" 
	   <<"'radiusDeterminationAlgorithm' needs to be set to '6'."
	   <<endl;
    cerr<<"In AsymShapeFunc::calcShapes input variable\n" 
	   <<"'radiusDeterminationAlgorithm' needs to be set to '6'."
	<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  double xrad,yrad,zrad;
  logFile<<"point "<<x<<" "<<y<<" "<<z<<endl;
  logFile<<"SupportSize "<<supportSize<<endl;
  for(int k=0;k<supportSize;k++) {
    logFile<<"Ptcle "<<sPtcls[k]<<": "<<ptcls[sPtcls[k]].getCoord(0)<<" "
	   <<ptcls[sPtcls[k]].getCoord(1)<<" "
	   <<ptcls[sPtcls[k]].getCoord(2)<<endl;
    logFile<<"Radii "<<ptcls[sPtcls[k]].getRadius(0)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims)<<" "
	   <<ptcls[sPtcls[k]].getRadius(1)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims+1)<<" "
	   <<ptcls[sPtcls[k]].getRadius(2)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims+2)<<endl;
    xrad = x-ptcls[sPtcls[k]].getCoord(0);
    yrad = y-ptcls[sPtcls[k]].getCoord(1);
    zrad = z-ptcls[sPtcls[k]].getCoord(2);
    logFile<<"rad "<<xrad<<" "<<yrad<<" "<<zrad<<" ";
    logFile<<"distance "<<sqrt(pow(xrad,2)+pow(yrad,2)+pow(zrad,2))<<endl;
  }
#endif

  // The size of the linear equation system.
  int linEQSize;
  int particlesNum = ptcls.size();

  // Calculate a polynom base set and its derivatives for 
  // all particles.
  //InputData->setValue("normedBasisPolynom",1);

  int derivOrder = 1;

  BasisPolynom* PolynomSet = new BasisPolyAsym(InputData,ptcls,sPtcls,
					       x,y,z,supportSize,
					       linEQSize,derivOrder,
					       modelData,
					       logFile);

  dbMatrix& P = PolynomSet->getBasis();
  dbMatrix& dPx = PolynomSet->getXDerivBasis();
  dbMatrix& dPy = PolynomSet->getYDerivBasis();
  dbMatrix& dPz = PolynomSet->getZDerivBasis();

  // Check the particle support
  if(linEQSize > supportSize) {
    logFile<<"In AsymShapeFunc::calcShapes not enough particles support "
	   <<"current point! (existing: "<<supportSize<<" - needed: "
	   <<linEQSize<<")"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.

  WindowFunctionSet* WFuncSet;

  WFuncSet = new WindowFuncAsym(InputData,ptcls,
				sPtcls,x,y,z,
				supportSize,
				derivOrder,
				modelData,
				logFile,viewerSEQ);

  dbVector& W = WFuncSet->getWindowFuncs();
  dbVector& dWx = WFuncSet->getXDerivWinFuncs();
  dbVector& dWy = WFuncSet->getYDerivWinFuncs();
  dbVector& dWz = WFuncSet->getZDerivWinFuncs();

#ifdef _geometryDebugMode_
  logFile<<"LinEQSize "<<linEQSize<<endl;
  logFile<<"*************** basis polynom set ****************"<<endl;
  logFile<<"Whole size "<<P.size()<<endl;
  for(int i=0;i<supportSize;i++) {
    logFile<<"Ptcle "<<sPtcls[i]<<" Size "<<P[i].size()<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<P[i][j]<<endl;
    logFile<<"*** dPx ****"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPx[i][j]<<endl;
    logFile<<"*** dPy ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPy[i][j]<<endl;
    logFile<<"*** dPz ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPz[i][j]<<endl;
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

  /**********************************************************************/
  // make use a shifted basis and continuous approximation set-up
  if(normedBasis) {

    // Assembly the momentum matrix 'M'
    Mat M;
    intVector matIdx(linEQSize);
    dbMatrix MValues(linEQSize,dbVector(linEQSize));
    
    map<string,bool> mOptions;
    createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);

    for(int i=0;i<linEQSize;i++) {
      matIdx[i] = i;

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  MValues[i][j] = MValues[i][j]+P[k][i]*P[k][j]*W[k];
   
    }

#ifdef _geometryDebugMode_
    logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<MValues[i][j];
      logFile<<endl;
    }
#endif

    for(int i=0;i<linEQSize;i++)
      MatSetValues(M,1,&i,linEQSize,&matIdx[0],&MValues[i][0],
		   INSERT_VALUES);
	
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
    MatSetOption(M,MAT_SYMMETRIC,PETSC_TRUE);
    //  MatView(M,viewerSEQ);
  
    // Assemble vector P(0).
    Vec P0;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&P0);
    intVector vecIdx(linEQSize);

    for(int i=0;i<linEQSize;i++)
      vecIdx[i] = i;

    VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
    // Create vector B.
    Vec B;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&B);
  
    // Solve the linear equation system 'MB = P(0)' and get vector B.
    KSP ksp;
    PC pc;
    Mat F;

    KSPCreate(MPI_COMM_SELF,&ksp);
    KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetUp(ksp);
  
    KSPSolve(ksp,P0,B);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);

    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    destroyPETScVec(P0);

    // Calculate the shape functions.
    double* BCalc;
    VecGetArray(B,&BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      shapeFuncs[i] = 0;

      for(int j=0;j<linEQSize;j++)
	shapeFuncs[i] = shapeFuncs[i]+P[i][j]*BCalc[j];

      shapeFuncs[i] = shapeFuncs[i]*W[i];
    }
  
#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<BCalc[i]<<endl;
    logFile<<"************ Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the x-derivation of momentum matrix 'M'
    dbMatrix dMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	dMValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dMValues[i][j] +=
	    ((dPx[k][i]*P[k][j]+P[k][i]*dPx[k][j])*W[k]+P[k][i]*P[k][j]*dWx[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = -dM B .
    Vec dBx;
    Vec RS;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&RS);
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBx);

    dbVector rightSide(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - dMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBx);
    //  VecView(sBx,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    double* dBCalc;
    VecGetArray(dBx,&dBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	xDerivShapes[i] = xDerivShapes[i] + dPx[i][j]*BCalc[j]*W[i]
	  + P[i][j]*dBCalc[j]*W[i] + P[i][j]*BCalc[j]*dWx[i];


    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dBCalc[i]<<endl;
    logFile<<"********** x-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBx,&dBCalc);
    destroyPETScVec(dBx);

    /********************************************************************/
    // Assemble the y-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	dMValues[i][j] = 0;
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dMValues[i][j] = dMValues[i][j]
	    +((dPy[k][i]*P[k][j]+P[k][i]*dPy[k][j])*W[k]+P[k][i]*P[k][j]*dWy[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = -dM B .
    Vec dBy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - dMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBy);
    //  VecView(dBy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the y-derivation of the shape functions.
    VecGetArray(dBy,&dBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	yDerivShapes[i] = yDerivShapes[i]
	  + dPy[i][j]*BCalc[j]*W[i] + P[i][j]*dBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWy[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dBCalc[i]<<endl;
    logFile<<"************* y-Abl. d. Ansatzfunktion **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBy,&dBCalc);
    destroyPETScVec(dBy);

    /********************************************************************/
    // Assemble the z-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	dMValues[i][j] = 0;
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dMValues[i][j] = dMValues[i][j]
	    +((dPz[k][i]*P[k][j]+P[k][i]*dPz[k][j])*W[k]+P[k][i]*P[k][j]*dWz[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = -dM B .*/
    Vec dBz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - dMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBz);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the z-derivation of the shape functions.
    VecGetArray(dBz,&dBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	zDerivShapes[i] = zDerivShapes[i]
	  + dPz[i][j]*BCalc[j]*W[i] + P[i][j]*dBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWz[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dBCalc[i]<<endl;
    logFile<<"************ z-Abl. d. Ansatzfunktionen **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBz,&dBCalc);
    destroyPETScVec(dBz);

    /********************************************************************/
    // Destroy all petsc objects.
    destroyPETScSolver(ksp);

    VecRestoreArray(B,&BCalc);
    destroyPETScVec(B);
    destroyPETScVec(RS);
    destroyPETScMat(M);

  }

  /**********************************************************************/
  // make use of the conventional MLS approximation scheme
  else {

    Mat M;
    intVector matIdx(linEQSize);
    dbMatrix MValues(linEQSize,dbVector(linEQSize));

    map<string,bool> mOptions;
    createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);

    for(int i=0;i<linEQSize;i++) {
      matIdx[i] = i;

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  MValues[i][j] = MValues[i][j]+P[k][i]*P[k][j]*W[k];
   
    }

#ifdef _geometryDebugMode_
    logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<MValues[i][j];
      logFile<<endl;
    }
#endif

    for(int i=0;i<linEQSize;i++)
      MatSetValues(M,1,&i,linEQSize,&matIdx[0],&MValues[i][0],
		   INSERT_VALUES);
	
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
    MatSetOption(M,MAT_SYMMETRIC,PETSC_TRUE);
    //  MatView(M,viewerSEQ);
  
    // Assemble vector P(0).
    Vec P0;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&P0);
    intVector vecIdx(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      vecIdx[i] = i;
    }

    VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
    // Create vector B.
    Vec B;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&B);
  
    // Solve the linear equation system 'MB = P(0)' and get vector B.
    KSP ksp;
    PC pc;
    Mat F;
    
    KSPCreate(MPI_COMM_SELF,&ksp);
    KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetUp(ksp);

  
    KSPSolve(ksp,P0,B);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);

    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    destroyPETScVec(P0);

    // Calculate the shape functions.
    double* BCalc;
    VecGetArray(B,&BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      shapeFuncs[i] = 0;

      for(int j=0;j<linEQSize;j++)
	shapeFuncs[i] = shapeFuncs[i]+P[i][j]*BCalc[j];

      shapeFuncs[i] = shapeFuncs[i]*W[i];
    }
  
#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<BCalc[i]<<endl;
    logFile<<"************ Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the x-derivation of momentum matrix 'M'
    dbMatrix dMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	dMValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dMValues[i][j] += P[k][i]*P[k][j]*dWx[k];
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = dP - dM B .
    Vec dBx;
    Vec RS;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&RS);
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBx);

    dbVector rightSide(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPx[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBx);
    //  VecView(sBx,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix\n  of the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    double* dBCalc;
    VecGetArray(dBx,&dBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	xDerivShapes[i] += P[i][j]*dBCalc[j]*W[i] + P[i][j]*BCalc[j]*dWx[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dBCalc[i]<<endl;
    logFile<<"********** x-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBx,&dBCalc);
    destroyPETScVec(dBx);

    /********************************************************************/
    // Assemble the y-derivation of momentum matrix 'M'.
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	dMValues[i][j] = 0;
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dMValues[i][j] += P[k][i]*P[k][j]*dWy[k];
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = dP - dM B .
    Vec dBy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPy[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBy);
    //  VecView(dBy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix\n  of the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix\n  of the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    VecGetArray(dBy,&dBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	yDerivShapes[i] += P[i][j]*dBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWy[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dBCalc[i]<<endl;
    logFile<<"************* y-Abl. d. Ansatzfunktion **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBy,&dBCalc);
    destroyPETScVec(dBy);

    /********************************************************************/
    // Assemble the z-derivation of momentum matrix 'M'.
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	dMValues[i][j] = 0;
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dMValues[i][j] += P[k][i]*P[k][j]*dWz[k];
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = dP - dM B .
    Vec dBz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPz[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBz);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    VecGetArray(dBz,&dBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	zDerivShapes[i] += P[i][j]*dBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWz[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dBCalc[i]<<endl;
    logFile<<"************ z-Abl. d. Ansatzfunktionen **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBz,&dBCalc);
    destroyPETScVec(dBz);

    /********************************************************************/
    // Destroy all petsc objects.
    destroyPETScSolver(ksp);

    VecRestoreArray(B,&BCalc);
    destroyPETScVec(B);
    destroyPETScVec(RS);
    destroyPETScMat(M);
  }

  delete PolynomSet,WFuncSet;
}

/************************************************************************/
/************************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first and second order derivations.
void AsymShapeFunc::calcShapes(InputFileData* InputData,
			       int& supportSize,
			       intVector& sPtcls,
			       std::vector<Particle>& ptcls,
			       double& x,double& y,double& z,
			       dbVector& shapeFuncs,
			       dbVector& xDerivShapes,
			       dbVector& yDerivShapes,
			       dbVector& zDerivShapes,
			       dbVector& xxDerivShapes,
			       dbVector& yyDerivShapes,
			       dbVector& zzDerivShapes,
			       dbVector& xyDerivShapes,
			       dbVector& yzDerivShapes,
			       dbVector& zxDerivShapes,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerSEQ) {

  using namespace std;

  InputData->setValue("windowFunctionType",4.0);
  //  InputData->getValue("windowfunctionNorming",2.0);

  int usedDims = (int)modelData["usedDimensions"];
  bool normedBasis = (bool)InputData->getValue("normedBasisPolynom");
  int wType = (int)InputData->getValue("windowFunctionType");

  if((int)InputData->getValue("radiusDeterminationAlgorithm") != 6) {
    logFile<<"In AsymShapeFunc::calcShapes input variable\n" 
	   <<"'radiusDeterminationAlgorithm' needs to be set to '6'."
	   <<endl;
    cerr<<"In AsymShapeFunc::calcShapes input variable\n" 
	   <<"'radiusDeterminationAlgorithm' needs to be set to '6'."
	<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // The size of the linear equation system.
  int linEQSize;
  int particlesNum = ptcls.size();

  int derivOrder = 2;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  double xrad,yrad,zrad;
  logFile<<"point "<<x<<" "<<y<<" "<<z<<endl;
  logFile<<"SupportSize "<<supportSize<<endl;
  for(int k=0;k<supportSize;k++) {
    logFile<<"Ptcle "<<sPtcls[k]<<": "<<ptcls[sPtcls[k]].getCoord(0)<<" "
	   <<ptcls[sPtcls[k]].getCoord(1)<<" "
	   <<ptcls[sPtcls[k]].getCoord(2)<<endl;
    logFile<<"Radii "<<ptcls[sPtcls[k]].getRadius(0)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims)<<" "
	   <<ptcls[sPtcls[k]].getRadius(1)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims+1)<<" "
	   <<ptcls[sPtcls[k]].getRadius(2)<<"/ "
	   <<ptcls[sPtcls[k]].getRadius(usedDims+2)<<endl;
    xrad = x-ptcls[sPtcls[k]].getCoord(0);
    yrad = y-ptcls[sPtcls[k]].getCoord(1);
    zrad = z-ptcls[sPtcls[k]].getCoord(2);
    logFile<<"rad "<<xrad<<" "<<yrad<<" "<<zrad<<" ";
    logFile<<"distance "<<sqrt(pow(xrad,2)+pow(yrad,2)+pow(zrad,2))<<endl;
  }
#endif

  // Calculate a polynom base set and the derivatives for 
  // all particles.
  //InputData->setValue("normedBasisPolynom",1);

  BasisPolynom* PolynomSet = new BasisPolyAsym(InputData,ptcls,sPtcls,
					       x,y,z,supportSize,
					       linEQSize,derivOrder,
					       modelData,logFile);

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

  // Check the particle support
  if(linEQSize > supportSize) {
    logFile <<"In AsymShapeFunc::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.

  WindowFunctionSet* WFuncSet;

  WFuncSet = new WindowFuncAsym(InputData,ptcls,
				sPtcls,x,y,z,
				supportSize,
				derivOrder,
				modelData,
				logFile,viewerSEQ);

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
  logFile<<"LinEQSize "<<linEQSize<<endl;
  logFile<<"*************** basis polynom set ****************"<<endl;
  logFile<<"Whole size "<<P.size()<<endl;
  for(int i=0;i<supportSize;i++) {
    logFile<<"Ptcle "<<sPtcls[i]<<" Size "<<P[i].size()<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<P[i][j]<<endl;
    logFile<<"*** dPx ****"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPx[i][j]<<endl;
    logFile<<"*** dPy ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPy[i][j]<<endl;
    logFile<<"*** dPz ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPz[i][j]<<endl;
    logFile<<"*** dPxx ****"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPxx[i][j]<<endl;
    logFile<<"*** dPyy ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPyy[i][j]<<endl;
    logFile<<"*** dPzz ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPzz[i][j]<<endl;
    logFile<<"*** dPxy ****"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPxy[i][j]<<endl;
    logFile<<"*** dPyz ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPyz[i][j]<<endl;
    logFile<<"*** dPzx ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPzx[i][j]<<endl;
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

  /********************************************************************/
  // make use a shifted basis and continuous approximation set-up
  if(normedBasis) {

    // Assembly the momentum matrix 'M'
    Mat M;
    intVector matIdx(linEQSize);
    dbMatrix MValues(linEQSize,dbVector(linEQSize));
    
    map<string,bool> mOptions;
    createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);

    for(int i=0;i<linEQSize;i++) {
      matIdx[i] = i;

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  MValues[i][j] = MValues[i][j]+P[k][i]*P[k][j]*W[k];
   
    }

#ifdef _geometryDebugMode_
    logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<MValues[i][j];
      logFile<<endl;
    }
#endif

    for(int i=0;i<linEQSize;i++)
      MatSetValues(M,1,&i,linEQSize,&matIdx[0],&MValues[i][0],
		   INSERT_VALUES);
	
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
    MatSetOption(M,MAT_SYMMETRIC,PETSC_TRUE);
    //  MatView(M,viewerSEQ);
  
    // Assemble the right side vector
    Vec RS;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&RS);
    intVector vecIdx(linEQSize);

    for(int i=0;i<linEQSize;i++)
      vecIdx[i] = i;

    VecSetValues(RS,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
    // Create vector B.
    Vec B;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&B);
  
    // Solve the linear equation system 'MB = P(0)' and get vector B.
    KSP ksp;
    PC pc;
    Mat F;
    
    KSPCreate(MPI_COMM_SELF,&ksp);
    KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetUp(ksp);

  
    KSPSolve(ksp,RS,B);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);

    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the shape functions.
    double* BCalc;
    VecGetArray(B,&BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      shapeFuncs[i] = 0;

      for(int j=0;j<linEQSize;j++)
	shapeFuncs[i] = shapeFuncs[i]+P[i][j]*BCalc[j];

      shapeFuncs[i] = shapeFuncs[i]*W[i];
    }
  
#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<BCalc[i]<<endl;
    logFile<<"************ Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the x-derivation of momentum matrix 'M'
    dbMatrix dxMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++)

      for(int j=0;j<linEQSize;j++)

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dxMValues[i][j] = dxMValues[i][j]
	    +((dPx[k][i]*P[k][j]+P[k][i]*dPx[k][j])*W[k]+P[k][i]*P[k][j]*dWx[k]);
     
	}

#ifdef _geometryDebugMode_
    logFile<<"******* dx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dxMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dxMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = -dM B .
    Vec dBx;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBx);

    dbVector rightSide(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dxMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBx);
    //  VecView(sBx,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    double* dxBCalc;
    VecGetArray(dBx,&dxBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	xDerivShapes[i] = xDerivShapes[i] + dPx[i][j]*BCalc[j]*W[i]
	  + P[i][j]*dxBCalc[j]*W[i] + P[i][j]*BCalc[j]*dWx[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dxBCalc[i]<<endl;
    logFile<<"********** x-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xDerivShapes[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the y-derivation of momentum matrix 'M'
    dbMatrix dyMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++)

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dyMValues[i][j] = dyMValues[i][j]
	    +((dPy[k][i]*P[k][j]+P[k][i]*dPy[k][j])*W[k]+P[k][i]*P[k][j]*dWy[k]);
	}

#ifdef _geometryDebugMode_
    logFile<<"******* dy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dyMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dyMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = -dM B .
    Vec dBy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dyMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBy);
    //  VecView(dBy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    double* dyBCalc;
    VecGetArray(dBy,&dyBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	yDerivShapes[i] = yDerivShapes[i]
	  + dPy[i][j]*BCalc[j]*W[i] + P[i][j]*dyBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWy[i];
    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dyBCalc[i]<<endl;
    logFile<<"************* y-Abl. d. Ansatzfunktion **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yDerivShapes[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the z-derivation of momentum matrix 'M'
    dbMatrix dzMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++)

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dzMValues[i][j] = dzMValues[i][j]
	    +((dPz[k][i]*P[k][j]+P[k][i]*dPz[k][j])*W[k]+P[k][i]*P[k][j]*dWz[k]);

	}

#ifdef _geometryDebugMode_
    logFile<<"******* dz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dzMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dzMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = -dM B .
    Vec dBz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dzMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBz);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    double* dzBCalc;
    VecGetArray(dBz,&dzBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	zDerivShapes[i] = zDerivShapes[i]
	  + dPz[i][j]*BCalc[j]*W[i] + P[i][j]*dzBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWz[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dzBCalc[i]<<endl;
    logFile<<"************ z-Abl. d. Ansatzfunktionen **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zDerivShapes[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the xx-derivation of momentum matrix 'M'
    dbMatrix d2MValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  d2MValues[i][j] +=
	    ((dPxx[k][i]*P[k][j]+2.0*dPx[k][i]*dPx[k][j]+P[k][i]*dPxx[k][j])*W[k]
	     +(2.0*dPx[k][i]*P[k][j]+2.0*P[k][i]*dPx[k][j])*dWx[k]
	     +P[k][i]*P[k][j]*dWxx[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dxx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBxx = -dMxx B - 2*(dMx dBx).
    Vec dBxx;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBxx);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  2.0*dxMValues[i][j]*dxBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBxx = RS' and get vector dBxx.
    KSPSolve(ksp,RS,dBxx);
    //  VecView(sBx,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the xx-derivation of the shape functions.
    double* d2BCalc;
    VecGetArray(dBxx,&d2BCalc);

#ifndef _MLSDebugMode_
    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xxDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	xxDerivShapes[i] = xxDerivShapes[i]

	  + (d2BCalc[j]*P[i][j] + 2.0*dxBCalc[j]*dPx[i][j] + 
	     + BCalc[j]*dPxx[i][j])*W[i]

	  + 2.0*(dxBCalc[j]*P[i][j] + BCalc[j]*dPx[i][j])*dWx[i]

	  + BCalc[j]*P[i][j]*dWxx[i];

    }
#else
    for(int i=0;i<supportSize;i++) {
      xxDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	xxDerivShapes[i] = xxDerivShapes[i]

	  + (d2BCalc[j]*P[i][j] + 2.0*dxBCalc[j]*dPx[i][j] + 
	     + BCalc[j]*dPxx[i][j])*W[i]

	  + 2.0*(dxBCalc[j]*P[i][j] + BCalc[j]*dPx[i][j])*dWx[i]

	  + BCalc[j]*P[i][j]*dWxx[i];

    }
#endif

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** xx-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xxDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBxx,&d2BCalc);
    destroyPETScVec(dBxx);

    /********************************************************************/
    // Assemble the yy-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  d2MValues[i][j] +=
	    ((dPyy[k][i]*P[k][j]+2.0*dPy[k][i]*dPy[k][j]+P[k][i]*dPyy[k][j])*W[k]
	     +(2.0*dPy[k][i]*P[k][j]+2.0*P[k][i]*dPy[k][j])*dWy[k]
	     +P[k][i]*P[k][j]*dWyy[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dyy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dByy = -dMyy B - 2*(dMy dBy).
    Vec dByy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dByy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  2.0*dyMValues[i][j]*dyBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dByy = RS' and get vector dByy.
    KSPSolve(ksp,RS,dByy);
    //  VecView(sByy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the yy-derivation of the shape functions.
    VecGetArray(dByy,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yyDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	yyDerivShapes[i] = yyDerivShapes[i]

	  + (d2BCalc[j]*P[i][j] + 2.0*dyBCalc[j]*dPy[i][j] + 
	     + BCalc[j]*dPyy[i][j])*W[i]

	  + 2.0*(dyBCalc[j]*P[i][j] + BCalc[j]*dPy[i][j])*dWy[i]

	  + BCalc[j]*P[i][j]*dWyy[i];


    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** yy-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yyDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dByy,&d2BCalc);
    destroyPETScVec(dByy);

    /********************************************************************/
    // Assemble the zz-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  d2MValues[i][j] +=
	    ((dPzz[k][i]*P[k][j]+2.0*dPz[k][i]*dPz[k][j]+P[k][i]*dPzz[k][j])*W[k]
	     +(2.0*dPz[k][i]*P[k][j]+2.0*P[k][i]*dPz[k][j])*dWz[k]
	     +P[k][i]*P[k][j]*dWzz[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dzz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBzz = -dMzz B - 2*(dMz dBz).
    Vec dBzz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBzz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  2.0*dzMValues[i][j]*dzBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBzz = RS' and get vector dBzz.
    KSPSolve(ksp,RS,dBzz);
    //  VecView(sBzz,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the zz-derivation of the shape functions.
    VecGetArray(dBzz,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zzDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	zzDerivShapes[i] = zzDerivShapes[i]

	  + (d2BCalc[j]*P[i][j] + 2.0*dzBCalc[j]*dPz[i][j] 
	     + BCalc[j]*dPzz[i][j])*W[i]

	  + 2.0*(dzBCalc[j]*P[i][j] + BCalc[j]*dPz[i][j])*dWz[i]

	  + BCalc[j]*P[i][j]*dWzz[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** zz-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zzDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBzz,&d2BCalc);
    destroyPETScVec(dBzz);

    /********************************************************************/
    // Assemble the xy-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  d2MValues[i][j] +=
	    ((dPxy[k][i]*P[k][j] + dPx[k][i]*dPy[k][j]
	      +dPy[k][i]*dPx[k][j] + P[k][i]*dPxy[k][j])*W[k]

	     +(dPx[k][i]*P[k][j] + P[k][i]*dPx[k][j])*dWy[k]

	     +(dPy[k][i]*P[k][j] + P[k][i]*dPy[k][j])*dWx[k]

	     +P[k][i]*P[k][j]*dWxy[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dxy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBxy = -dMxy B - dMx dBy - dMy dBx .
    Vec dBxy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBxy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] 
	  - dxMValues[i][j]*dyBCalc[j] - dyMValues[i][j]*dxBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBxy = RS' and get vector dBxy.*/
    KSPSolve(ksp,RS,dBxy);
    //  VecView(sBxy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the xy-derivation of the shape functions.
    VecGetArray(dBxy,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xyDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	xyDerivShapes[i] = xyDerivShapes[i]

	  + (d2BCalc[j]*P[i][j] + dxBCalc[j]*dPy[i][j] + 
	     dyBCalc[j]*dPx[i][j] + BCalc[j]*dPxy[i][j])*W[i]

	  + (dxBCalc[j]*P[i][j] + BCalc[j]*dPx[i][j])*dWy[i]

	  + (dyBCalc[j]*P[i][j] + BCalc[j]*dPy[i][j])*dWx[i]

	  + BCalc[j]*P[i][j]*dWxy[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** xy-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xyDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBxy,&d2BCalc);
    destroyPETScVec(dBxy);


    /********************************************************************/
    // Assemble the yz-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  d2MValues[i][j] +=
	    ((dPyz[k][i]*P[k][j] + dPy[k][i]*dPz[k][j]
	      +dPz[k][i]*dPy[k][j] + P[k][i]*dPyz[k][j])*W[k]

	     +(dPy[k][i]*P[k][j] + P[k][i]*dPy[k][j])*dWz[k]

	     +(dPz[k][i]*P[k][j] + P[k][i]*dPz[k][j])*dWy[k]

	     +P[k][i]*P[k][j]*dWyz[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dyz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dByz = -dMyz B - dMy dBz - dMz dBy .
    Vec dByz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dByz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] 
	  - dyMValues[i][j]*dzBCalc[j] - dzMValues[i][j]*dyBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dByz = RS' and get vector dByz.
    KSPSolve(ksp,RS,dByz);
    //  VecView(sByz,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the xy-derivation of the shape functions.
    VecGetArray(dByz,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yzDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	yzDerivShapes[i] = yzDerivShapes[i]

	  + (d2BCalc[j]*P[i][j] + dyBCalc[j]*dPz[i][j] + 
	     dzBCalc[j]*dPy[i][j] + BCalc[j]*dPyz[i][j])*W[i]

	  + (dyBCalc[j]*P[i][j] + BCalc[j]*dPy[i][j])*dWz[i]

	  + (dzBCalc[j]*P[i][j] + BCalc[j]*dPz[i][j])*dWy[i]

	  + BCalc[j]*P[i][j]*dWyz[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** yz-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yzDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dByz,&d2BCalc);
    destroyPETScVec(dByz);

    /********************************************************************/
    // Assemble the zx-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  d2MValues[i][j] +=
	    ((dPzx[k][i]*P[k][j] + dPz[k][i]*dPx[k][j]
	      +dPx[k][i]*dPz[k][j] + P[k][i]*dPzx[k][j])*W[k]

	     +(dPz[k][i]*P[k][j] + P[k][i]*dPz[k][j])*dWx[k]

	     +(dPx[k][i]*P[k][j] + P[k][i]*dPx[k][j])*dWz[k]

	     +P[k][i]*P[k][j]*dWzx[k]);
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dzx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBzx = -dMzx B - dMz dBx - dMx dBz .
    Vec dBzx;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBzx);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = 0;

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] 
	  - dzMValues[i][j]*dxBCalc[j] - dxMValues[i][j]*dzBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBzx = RS' and get vector dBzx.
    KSPSolve(ksp,RS,dBzx);
    //  VecView(sBxy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In AsymShapeFunc::calcShapes calculation of moment "
	      <<"matrix's derivation of\n  the weight function failed!"
	      <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the xy-derivation of the shape functions.
    VecGetArray(dBzx,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zxDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	zxDerivShapes[i] = zxDerivShapes[i]

	  + (d2BCalc[j]*P[i][j] + dzBCalc[j]*dPx[i][j] + 
	     dxBCalc[j]*dPz[i][j] + BCalc[j]*dPzx[i][j])*W[i]

	  + (dzBCalc[j]*P[i][j] + BCalc[j]*dPz[i][j])*dWx[i]

	  + (dxBCalc[j]*P[i][j] + BCalc[j]*dPx[i][j])*dWz[i]

	  + BCalc[j]*P[i][j]*dWzx[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** zx-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zxDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBzx,&d2BCalc);
    destroyPETScVec(dBzx);

    /********************************************************************/
    // Destroy all petsc objects.
    destroyPETScSolver(ksp);

    VecRestoreArray(B,&BCalc);
    VecRestoreArray(dBx,&dxBCalc);
    VecRestoreArray(dBy,&dyBCalc);
    VecRestoreArray(dBz,&dzBCalc);

    destroyPETScVec(B);
    destroyPETScVec(dBx);
    destroyPETScVec(dBy);
    destroyPETScVec(dBz);

    destroyPETScVec(RS);
    destroyPETScMat(M);

  }
  /**********************************************************************/
  // make use of the conventional MLS approximation scheme
  else {

    Mat M;
    intVector matIdx(linEQSize);
    dbMatrix MValues(linEQSize,dbVector(linEQSize));

    map<string,bool> mOptions;
    createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);

    for(int i=0;i<linEQSize;i++) {
      matIdx[i] = i;

      for(int j=0;j<linEQSize;j++)
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  MValues[i][j] += P[k][i]*P[k][j]*W[k];
   
    }

#ifdef _geometryDebugMode_
    logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<MValues[i][j];
      logFile<<endl;
    }
#endif

    for(int i=0;i<linEQSize;i++)
      MatSetValues(M,1,&i,linEQSize,&matIdx[0],&MValues[i][0],
		   INSERT_VALUES);
	
    MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
    MatSetOption(M,MAT_SYMMETRIC,PETSC_TRUE);
    //  MatView(M,viewerSEQ);
  
    // Assemble vector P(0).
    Vec P0;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&P0);
    intVector vecIdx(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      vecIdx[i] = i;
    }

    VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
    // Create vector B.
    Vec B;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&B);
  
    // Solve the linear equation system 'MB = P(0)' and get vector B.
    KSP ksp;
    PC pc;
    Mat F;
    
    KSPCreate(MPI_COMM_SELF,&ksp);
    KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetUp(ksp);
    

    KSPSolve(ksp,P0,B);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);

    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    destroyPETScVec(P0);

    // Calculate the shape functions.
    double* BCalc;
    VecGetArray(B,&BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      shapeFuncs[i] = 0;

      for(int j=0;j<linEQSize;j++)
	shapeFuncs[i] += P[i][j]*BCalc[j];

      shapeFuncs[i] *= W[i];
    }
  
#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<BCalc[i]<<endl;
    logFile<<"************ Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the x-derivation of momentum matrix 'M'
    dbMatrix dxMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dxMValues[i][j] += P[k][i]*P[k][j]*dWx[k];
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dxMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dxMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = dP - dM B .
    Vec dBx;
    Vec RS;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&RS);
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBx);

    dbVector rightSide(linEQSize);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPx[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dxMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBx);
    //  VecView(sBx,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix\n  of the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the x-derivation of the shape functions.
    double* dxBCalc;
    VecGetArray(dBx,&dxBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	xDerivShapes[i] += P[i][j]*dxBCalc[j]*W[i] + P[i][j]*BCalc[j]*dWx[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dxBCalc[i]<<endl;
    logFile<<"********** x-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xDerivShapes[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the y-derivation of momentum matrix 'M'.
    dbMatrix dyMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dyMValues[i][j] += P[k][i]*P[k][j]*dWy[k];
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dyMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dyMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = dP - dM B .
    Vec dBy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPy[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dyMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBy);
    //  VecView(dBy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix\n  of the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix\n  of the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the y-derivation of the shape functions.
    double* dyBCalc;
    VecGetArray(dBy,&dyBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	yDerivShapes[i] += P[i][j]*dyBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWy[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dyBCalc[i]<<endl;
    logFile<<"************* y-Abl. d. Ansatzfunktion **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yDerivShapes[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the z-derivation of momentum matrix 'M'.
    dbMatrix dzMValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
      
	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++) {
	  dzMValues[i][j] += P[k][i]*P[k][j]*dWz[k];
	}
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(dzMValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<dzMValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dB = dP - dM B .
    Vec dBz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPz[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] -= dzMValues[i][j]*BCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'MdB = RS' and get vector B.
    KSPSolve(ksp,RS,dBz);
    //  VecView(B,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the z-derivation of the shape functions.
    double* dzBCalc;
    VecGetArray(dBz,&dzBCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)
	zDerivShapes[i] += P[i][j]*dzBCalc[j]*W[i]
	  + P[i][j]*BCalc[j]*dWz[i];

    }

#ifdef _geometryDebugMode_
    logFile<<"************** b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<dzBCalc[i]<<endl;
    logFile<<"************ z-Abl. d. Ansatzfunktionen **************"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zDerivShapes[i]<<endl;
#endif

    /********************************************************************/
    // Assemble the xx-derivation of momentum matrix 'M'
    dbMatrix d2MValues(linEQSize,dbVector(linEQSize));

    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  d2MValues[i][j] += P[k][i]*P[k][j]*dWxx[k];
      
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dxx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBxx = dPxx -dMxx B - 2*(dMx dBx).*/
    Vec dBxx;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBxx);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPxx[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -
	  2.0*dxMValues[i][j]*dxBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBxx = RS' and get vector dBxx.
    KSPSolve(ksp,RS,dBxx);
    //  VecView(sBx,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the xx-derivation of the shape functions.
    double* d2BCalc;
    VecGetArray(dBxx,&d2BCalc);

#ifndef _MLSDebugMode_
    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xxDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	xxDerivShapes[i] += d2BCalc[j]*P[i][j]*W[i] + 
	  2.0*dxBCalc[j]*P[i][j]*dWx[i] + BCalc[j]*P[i][j]*dWxx[i];
    }
#else
    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xxDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	xxDerivShapes[i] += d2BCalc[j]*P[i][j]*W[i] + 
	  2.0*dxBCalc[j]*P[i][j]*dWx[i] + BCalc[j]*P[i][j]*dWxx[i];
    }
#endif

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** xx-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xxDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBxx,&d2BCalc);
    destroyPETScVec(dBxx);

    /********************************************************************/
    // Assemble the yy-derivation of momentum matrix 'M'.
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  d2MValues[i][j] += P[k][i]*P[k][j]*dWyy[k];
      
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dyy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dByy = dPyy -dMyy B - 2*(dMy dBy).
    Vec dByy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dByy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPyy[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  2.0*dyMValues[i][j]*dyBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dByy = RS' and get vector dByy.
    KSPSolve(ksp,RS,dByy);
    //  VecView(sByy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the yy-derivation of the shape functions.
    VecGetArray(dByy,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yyDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	yyDerivShapes[i] += d2BCalc[j]*P[i][j]*W[i] + 
	  2.0*dyBCalc[j]*P[i][j]*dWy[i] + BCalc[j]*P[i][j]*dWyy[i];
    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** yy-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yyDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dByy,&d2BCalc);
    destroyPETScVec(dByy);

    /********************************************************************/
    // Assemble the zz-derivation of momentum matrix 'M'.
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  d2MValues[i][j] += P[k][i]*P[k][j]*dWzz[k];
      
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dzz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBzz = dPzz - dMzz B - 2*(dMz dBz).
    Vec dBzz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBzz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPzz[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  2.0*dzMValues[i][j]*dzBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBzz = RS' and get vector dBzz.
    KSPSolve(ksp,RS,dBzz);
    //  VecView(sBzz,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the zz-derivation of the shape functions.
    VecGetArray(dBzz,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zzDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	zzDerivShapes[i] += d2BCalc[j]*P[i][j]*W[i] + 
	  2.0*dzBCalc[j]*P[i][j]*dWz[i] + BCalc[j]*P[i][j]*dWzz[i];
    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** zz-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zzDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBzz,&d2BCalc);
    destroyPETScVec(dBzz);

    /********************************************************************/
    // Assemble the xy-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  d2MValues[i][j] += P[k][i]*P[k][j]*dWxy[k];
      
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dxy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBxy = dPxy -dMxy B - dMx dBy - dMy dBx .
    Vec dBxy;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBxy);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPxy[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  dxMValues[i][j]*dyBCalc[j] - dyMValues[i][j]*dxBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBxy = RS' and get vector dBxy.
    KSPSolve(ksp,RS,dBxy);
    //  VecView(sBxy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the xy-derivation of the shape functions.
    VecGetArray(dBxy,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      xyDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	xyDerivShapes[i] += d2BCalc[j]*P[i][j]*W[i] + 
	  dxBCalc[j]*P[i][j]*dWy[i] + dyBCalc[j]*P[i][j]*dWx[i] + 
	  BCalc[j]*P[i][j]*dWxy[i];
    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** xy-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<xyDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBxy,&d2BCalc);
    destroyPETScVec(dBxy);

    /********************************************************************/
    // Assemble the yz-derivation of momentum matrix 'M'
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  d2MValues[i][j] += P[k][i]*P[k][j]*dWyz[k];
      
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dyz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dByz = dPyz -dMyz B - dMy dBz - dMz dBy .
    Vec dByz;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dByz);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPyz[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  dyMValues[i][j]*dzBCalc[j] - dzMValues[i][j]*dyBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dByz = RS' and get vector dByz.
    KSPSolve(ksp,RS,dByz);
    //  VecView(sByz,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the yz-derivation of the shape functions.
    VecGetArray(dByz,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      yzDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	yzDerivShapes[i] += d2BCalc[j]*P[i][j]*W[i] + 
	  dyBCalc[j]*P[i][j]*dWz[i] + dzBCalc[j]*P[i][j]*dWy[i] + 
	  BCalc[j]*P[i][j]*dWyz[i];
    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** yz-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<yzDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dByz,&d2BCalc);
    destroyPETScVec(dByz);

    /********************************************************************/
    // Assemble the zx-derivation of momentum matrix 'M'.
    for(int i=0;i<linEQSize;i++) {

      for(int j=0;j<linEQSize;j++) {
	d2MValues[i][j] = 0;

	// Only a loop over the supported particles, since in the other 
	// case the weight function is "0".
	for(int k=0;k<supportSize;k++)
	  d2MValues[i][j] += P[k][i]*P[k][j]*dWzx[k];
      
      }
    }

#ifdef _geometryDebugMode_
    logFile<<"******* dzx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
    for(int i=0;i<linEQSize;i++) {
      for(int j=0;j<linEQSize;j++)
	if(d2MValues[i][j]<0.00000000001)
	  logFile<<" 0.000000";
	else
	  logFile<<" "<<d2MValues[i][j];
      logFile<<endl;
    }
#endif

    // Assemble the right side vector 'RS' of the linear eq. system 
    // M dBzx = dPzx -dMzx B - dMz dBx - dMx dBz .
    Vec dBzx;
    VecCreateSeq(PETSC_COMM_SELF,linEQSize,&dBzx);

    for(int i=0;i<linEQSize;i++) {
      rightSide[i] = dPzx[supportSize][i];

      for(int j=0;j<linEQSize;j++)
	rightSide[i] = rightSide[i] - d2MValues[i][j]*BCalc[j] -  
	  dzMValues[i][j]*dxBCalc[j] - dxMValues[i][j]*dzBCalc[j];
    }

    VecSetValues(RS,linEQSize,&vecIdx[0],&rightSide[0],INSERT_VALUES);

    // Solve the linear equation system 'M dBzx = RS' and get vector dBzx.
    KSPSolve(ksp,RS,dBzx);
    //  VecView(sBxy,viewerSEQ);
  
    // Check whether the solving was successful or not.
    if(reason < 1) {
      logFile <<"In EFGShapeFunc::calcShapes calculation of moment "
	      <<"matrix of\n  the weight function failed!"<<endl;
      cerr <<"In EFGShapeFunc::calcShapes calculation of moment "
	   <<"matrix of\n  the weight function failed!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Calculate the zx-derivation of the shape functions.
    VecGetArray(dBzx,&d2BCalc);

    // Only a loop over the supported particles, since in the other 
    // case the weight function is "0".
    for(int i=0;i<supportSize;i++) {
      zxDerivShapes[i] = 0;

      for(int j=0;j<linEQSize;j++)

	zxDerivShapes[i] += d2BCalc[j]*P[i][j]*W[i] + 
	  dzBCalc[j]*P[i][j]*dWx[i] + dxBCalc[j]*P[i][j]*dWz[i] + 
	  BCalc[j]*P[i][j]*dWzx[i];
    }

#ifdef _geometryDebugMode_
    logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
    for(int i=0;i<linEQSize;i++)
      logFile<<i<<" "<<d2BCalc[i]<<endl;
    logFile<<"********** zx-Abl. d. Ansatzfunktionen **********"<<endl;
    for(int i=0;i<supportSize;i++)
      logFile<<i<<" "<<sPtcls[i]<<" "<<zxDerivShapes[i]<<endl;
#endif

    VecRestoreArray(dBzx,&d2BCalc);
    destroyPETScVec(dBzx);

    /********************************************************************/
    // Destroy all petsc objects.
    destroyPETScSolver(ksp);

    VecRestoreArray(B,&BCalc);
    VecRestoreArray(dBx,&dxBCalc);
    VecRestoreArray(dBy,&dyBCalc);
    VecRestoreArray(dBz,&dzBCalc);

    destroyPETScVec(B);
    destroyPETScVec(dBx);
    destroyPETScVec(dBy);
    destroyPETScVec(dBz);

    destroyPETScVec(RS);
    destroyPETScMat(M);

  }

  delete PolynomSet,WFuncSet;
}
