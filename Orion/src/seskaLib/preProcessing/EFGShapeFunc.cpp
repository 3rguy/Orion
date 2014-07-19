#include "EFGShapeFunc.h"

// Calculate at a point for all its supporting particles their 
// shape functions.
void EFGShapeFunc::calcShapes(InputFileData* InputData,
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

  // The size of the linear equation system.
  int linEQSize;
  int particlesNum = ptcls.size();

  int derivOrder = 0;

  // Calculate a normal polynom base set and the derivation for 
  // all particles.
  InputData->setValue("dimensionlessBasisPolynom",0);

#ifdef _geometryDebugMode_
  double xrad,yrad,zrad;
  logFile<<"##########################################################"<<endl;
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
#endif

  BasisPolynom* PolynomSet = new BasisPolyRegular(InputData,ptcls,sPtcls,
						  x,y,z,supportSize,
						  linEQSize,derivOrder,
						  modelData,logFile);

  dbMatrix& P = PolynomSet->getBasis();
  basisTermNum = P[0].size();

  // Check the particle support
  if(linEQSize > supportSize) {
    logFile <<"In ShapeFunctionSets::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl;
    cerr <<"In ShapeFunctionSets::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl; 

      MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new WindowFunctionSet(InputData,ptcls,
						      sPtcls,x,y,z,
						      supportSize,0,
						      modelData,
						      logFile);
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
  logFile<<"basis at point: "<<endl;
  for(int j=0;j<linEQSize;j++)
    logFile<<supportSize<<" "<<j<<" "<<P[supportSize][j]<<endl;
  logFile<<"**************** window functions ***************"<<endl;
  double summation = 0;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++) {
    summation += W[i];
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
  }
  logFile<<"sum W = "<<summation<<endl;
#endif

  /********************************************************************/
  // Assembly the momentum matrix 'M'.
  Mat M;
  intVector matIdx(linEQSize);
  dbMatrix MValues(linEQSize,dbVector(linEQSize));

  map<string,bool> mOptions;
  createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);
  //MatCreateSeqDense(PETSC_COMM_SELF,linEQSize,linEQSize,PETSC_NULL,&M);  
  //MatSetFromOptions(M);

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
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&P0);
  intVector vecIdx(linEQSize);

  for(int i=0;i<linEQSize;i++) {
    vecIdx[i] = i;
  }

  VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
  // Create vector B.
  Vec B;
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&B);
  
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

  if(shapeFuncs.size() == 0)
    resizeArray(shapeFuncs,supportSize);

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

  delete PolynomSet,WFuncSet;
}

// Calculate at a node for all its supporting particles their
// shape functions (used in POD calculation)
void EFGShapeFunc::calcShapes(InputFileData* InputData,
			      int& supportSize,
			      intVector& sPtcls,
			      std::vector<Particle>& ptcls,
			      double& x,double& y,double& z,
			      dbVector& shapeFuncs,
			      int& basisTermNum,
			      std::ofstream& logFile) {

  using namespace std;

  // The size of the linear equation system.
  int linEQSize;
  int particlesNum = ptcls.size();

  int derivOrder = 0;

  // Calculate a normal polynom base set and the derivation for
  // all particles.
  InputData->setValue("dimensionlessBasisPolynom",0);

  BasisPolynom* PolynomSet = new BasisPolyRegular(InputData,ptcls,sPtcls,
						  x,y,z,supportSize,
						  linEQSize,derivOrder,logFile);

  dbMatrix& P = PolynomSet->getBasis();
  basisTermNum = P[0].size();

  // Check the particle support
  if(linEQSize > supportSize) {
    logFile <<"In ShapeFunctionSets::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl;
    cerr <<"In ShapeFunctionSets::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl;

      MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new WindowFunctionSet(InputData,ptcls,
						      sPtcls,x,y,z,
						      supportSize,0,
						      logFile);
  dbVector& W = WFuncSet->getWindowFuncs();

#ifdef _geometryDebugMode_
  double xrad,yrad,zrad;
  logFile<<"##########################################################"<<endl;
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
  logFile<<"LinEQSize "<<linEQSize<<endl;
  logFile<<"*************** basis polynom set ****************"<<endl;
  logFile<<"Whole size "<<P.size()<<endl;
  for(int i=0;i<supportSize;i++) {
    logFile<<"Ptcle "<<sPtcls[i]<<" Size "<<P[i].size()<<endl;

    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<P[i][j]<<endl;
  }
  logFile<<"basis at point: "<<endl;
  for(int j=0;j<linEQSize;j++)
    logFile<<supportSize<<" "<<j<<" "<<P[supportSize][j]<<endl;
  logFile<<"**************** window functions ***************"<<endl;
  double summation = 0;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++) {
    summation += W[i];
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
  }
  logFile<<"sum W = "<<summation<<endl;
#endif

  /********************************************************************/
  // Assembly the momentum matrix 'M'.
  Mat M;
  intVector matIdx(linEQSize);
  dbMatrix MValues(linEQSize,dbVector(linEQSize));

  map<string,bool> mOptions;
  createDenseSequentialPETScMat(linEQSize,linEQSize,mOptions,M,logFile);
  //MatCreateSeqDense(PETSC_COMM_SELF,linEQSize,linEQSize,PETSC_NULL,&M);
  //MatSetFromOptions(M);

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
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&P0);
  intVector vecIdx(linEQSize);

  for(int i=0;i<linEQSize;i++) {
    vecIdx[i] = i;
  }

  VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);

  // Create vector B.
  Vec B;
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&B);

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

  delete PolynomSet,WFuncSet;
}

/**********************************************************************/
/**********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first order derivations.
void EFGShapeFunc::calcShapes(InputFileData* InputData,
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

  // The size of the linear equation system.
  int linEQSize;
  int particlesNum = ptcls.size();

  int derivOrder = 1;

  // Calculate a normal polynom base set and the derivation for 
  // all particles.
  InputData->setValue("dimensionlessBasisPolynom",0);

  BasisPolynom* PolynomSet = new BasisPolyRegular(InputData,ptcls,sPtcls,
						  x,y,z,supportSize,
						  linEQSize,derivOrder,
						  modelData,logFile);

  dbMatrix& P = PolynomSet->getBasis();
  dbMatrix& dPx = PolynomSet->getXDerivBasis();
  dbMatrix& dPy = PolynomSet->getYDerivBasis();
  dbMatrix& dPz = PolynomSet->getZDerivBasis();

  // Check the particle support
  if(linEQSize > supportSize) {
    logFile <<"In ShapeFunctionSets::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl;
    cerr <<"In ShapeFunctionSets::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl; 

      MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new WindowFunctionSet(InputData,ptcls,
						      sPtcls,x,y,z,
						      supportSize,
						      derivOrder,
						      modelData,
						      logFile);
  dbVector& W = WFuncSet->getWindowFuncs();
  dbVector& dWx = WFuncSet->getXDerivWinFuncs();
  dbVector& dWy = WFuncSet->getYDerivWinFuncs();
  dbVector& dWz = WFuncSet->getZDerivWinFuncs();

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
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
  logFile<<"LinEQSize "<<linEQSize<<endl;
  logFile<<"*************** basis polynom set ****************"<<endl;
  logFile<<"Whole size "<<P.size()<<endl;
  for(int i=0;i<=supportSize;i++) {
    if(i==supportSize)
      logFile<<"INT_POINT: Size "<<P[i].size()<<endl;
    else
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

  /********************************************************************/
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
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&P0);
  intVector vecIdx(linEQSize);

  for(int i=0;i<linEQSize;i++) {
    vecIdx[i] = i;
  }

  VecSetValues(P0,linEQSize,&vecIdx[0],&P[supportSize][0],INSERT_VALUES);
  
  // Create vector B.
  Vec B;
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&B);
  
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
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&RS);
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&dBx);

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
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&dBy);

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
  VecCreateSeq(MPI_COMM_SELF,linEQSize,&dBz);

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

  delete PolynomSet,WFuncSet;
}


/**********************************************************************/
/**********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first and second order derivations.
void EFGShapeFunc::calcShapes(InputFileData* InputData,
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

  using namespace std;

  // The size of the linear equation system.
  int linEQSize;
  int particlesNum = ptcls.size();

  int derivOrder = 2;

  // Calculate a dimension-less polynom base set and the derivation for 
  // all particles.
  InputData->setValue("dimensionlessBasisPolynom",0);

  BasisPolynom* PolynomSet = new BasisPolyRegular(InputData,ptcls,sPtcls,
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
    logFile <<"In EFGShapeFunc::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl;
    cerr <<"In EFGShapeFunc::calcShapes not enough particles support "
	    <<"current point! (existing: "<<supportSize<<" - needed: "
	    <<linEQSize<<")"<<endl; 

      MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new WindowFunctionSet(InputData,ptcls,
						      sPtcls,x,y,z,
						      supportSize,
						      derivOrder,
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
  logFile<<"LinEQSize "<<linEQSize<<endl;
  logFile<<"*************** basis polynom set ****************"<<endl;
  logFile<<"Whole size "<<P.size()<<endl;
  for(int i=0;i<=supportSize;i++) {
    if(i==supportSize)
      logFile<<"INT_POINT: Size "<<P[i].size()<<endl;
    else
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
  double summation = 0;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++) {
    summation += W[i];
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
  }
  logFile<<"sum W = "<<summation<<endl;
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
  KSPCreate(PETSC_COMM_SELF,&ksp);
  KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,"lu"); // if one want to use lu-factorizing PCSetType(pc,"lu");
  //  PCLUSetUseInPlace(pc); // destroy the original matrix
  KSPSetType(ksp,"preonly"); // direct solving
  KSPSetFromOptions(ksp);
  
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

  delete PolynomSet,WFuncSet;

#ifdef _geometryDebugMode_
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
    pum += shapeFuncs[i];
    pumDx += xDerivShapes[i];
    pumDy += yDerivShapes[i];
    pumDz += zDerivShapes[i];
    pumDxx += xxDerivShapes[i];
    pumDyy += yyDerivShapes[i];
    pumDzz += zzDerivShapes[i];
    pumDxy += xyDerivShapes[i];
    pumDyz += yzDerivShapes[i];
    pumDzx += zxDerivShapes[i];
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

}
