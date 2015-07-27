#include "OrthoShapeFunc.h"

// Calculate at a point for all its supporting particles their 
// shape functions.
void
OrthoShapeFunc::calcShapes(InputFileData* InputData,int& supportSize,
                           intVector& sPtcls,std::vector<Particle>& ptcls,
                           double& x,double& y,double& z,dbVector& shapeFuncs,
                           int& basisTermNum,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile,PetscViewer& viewerSEQ) {

  using namespace std;

  // The size of the linear equation system.
  int linEQSize,modLinEQSize;
  int particlesNum = ptcls.size();

  // Calculate a regular polynom base set and the derivation for 
  // all particles.
  InputData->setValue("dimensionlessBasisPolynom",0);

  BasisPolynom* RegPolySet = new BasisPolyRegular(InputData,ptcls,sPtcls,x,y,z,
                                                  supportSize,linEQSize,0,
                                                  modelData,logFile);

  basisTermNum = linEQSize;
  modLinEQSize = linEQSize - 1;

  if(modLinEQSize < 1) {
    logFile << "In OrthoShapeFuncSet::calcShapes a constant basis\n"
        << "polynomial is not supported!" << endl;
    cerr << "In OrthoShapeFuncSet::calcShapes a constant basis\n"
        << "polynomial is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Check the particle support
  if(modLinEQSize > supportSize) {
    logFile << "In OrthoShapeFuncSet::calcShapes not enough particles support "
        << "current point! (existing: " << supportSize << " - needed: "
        << modLinEQSize << ")" << endl;
    cerr << "In OrthoShapeFuncSet::calcShapes not enough particles support "
        << "current point! (existing: " << supportSize << " - needed: "
        << modLinEQSize << ")" << endl;

    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate a Shepard shapefunction set.
  ShapefunctionSet* ShepardShapeSet = new ShepardShapeFunc(InputData,
                                                           supportSize,sPtcls,
                                                           ptcls,x,y,z,0,
                                                           modelData,logFile,
                                                           viewerSEQ);
  dbVector& shepardShapes = ShepardShapeSet->getShapefunctions();

  // Calculate the actual basis.
  BasisPolynom* PolynomSet = new BasisPolyOrtho(InputData,RegPolySet,
                                                ShepardShapeSet,ptcls,sPtcls,x,
                                                y,z,supportSize,0,modelData,
                                                logFile,viewerSEQ);

  delete RegPolySet;

  dbMatrix& P = PolynomSet->getBasis();

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new PrismaticWindowFunctionSet(InputData,ptcls,
                                                               sPtcls,x,y,z,
                                                               supportSize,0,
                                                               modelData,
                                                               logFile);

  dbVector& W = WFuncSet->getWindowFuncs();

#ifdef _geometryDebugMode_
  double xrad,yrad,zrad;
  logFile<<"######################################################"<<endl;
  logFile<<"********** ortho shapefunctions set ******************"<<endl;
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
  logFile<<"**************** window functions ***************"<<endl;
  logFile<<"Size "<<W.size()<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<W[i]<<endl;
  logFile<<"***************** shepard function set **************"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<"Ptcle "<<sPtcls[i]<<" "<<shepardShapes[i]<<endl;
#endif

  /********************************************************************/
  // Assembly the momentum matrix 'M'
  Mat M;
  intVector matIdx(modLinEQSize);
  dbMatrix MValues(modLinEQSize,dbVector(modLinEQSize));

  MatCreateSeqDense(PETSC_COMM_SELF,modLinEQSize,modLinEQSize,PETSC_NULL, &M);
  MatSetFromOptions(M);

  for(int i = 1;i < linEQSize;i++) {
    matIdx[i - 1] = i - 1;

    for(int j = 1;j < linEQSize;j++)
      
      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++)
        MValues[i - 1][j - 1] += P[k][i] * P[k][j] * W[k];

  }

#ifdef _geometryDebugMode_
  logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(MValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<MValues[i][j];
    logFile<<endl;
  }
#endif

  for(int i = 0;i < modLinEQSize;i++)
    MatSetValues(M,1, &i,modLinEQSize, &matIdx[0], &MValues[i][0],
                 INSERT_VALUES);

  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
  MatOption mOption = MAT_SYMMETRIC;
  PetscBool flag = PETSC_TRUE;
  MatSetOption(M,mOption,flag);
  //  MatView(M,viewerSEQ);

  // Assemble the right side vector.
  Vec RS;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &RS);
  intVector vecIdx(modLinEQSize);

  for(int i = 0;i < modLinEQSize;i++) {
    vecIdx[i] = i;
  }
  
  VecSetValues(RS,modLinEQSize, &vecIdx[0], &P[supportSize][1],INSERT_VALUES);

  // Create vector B.
  Vec B;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &B);
  
  // Solve the linear equation system 'MB = RS' and get vector B.
  KSP ksp;
  PC pc;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
  KSPGetPC(ksp, &pc);
  PCSetType(pc,"lu"); // if one want to use lu-factorizing PCSetType(pc,"lu");
  //PCLUSetUseInPlace(pc); // destroy the original matrix
  KSPSetType(ksp,"preonly"); // direct solving
  KSPSetFromOptions(ksp);
  
  KSPSolve(ksp,RS,B);

  // Check whether the solving was successful or not.
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  //  VecView(B,viewerSEQ);

  // Calculate the shape functions.
  double* BCalc;
  VecGetArray(B, &BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    shapeFuncs[i] = 0;

    // a_j p_j
    for(int j = 1;j < linEQSize;j++)
      shapeFuncs[i] += P[i][j] * BCalc[j - 1];

    shapeFuncs[i] *= W[i];

    // Shepard part
    shapeFuncs[i] += shepardShapes[i];
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************ Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif
  
  // Destroy all petsc objects.
  destroyPETScSolver(ksp);

  VecRestoreArray(B, &BCalc);
  destroyPETScVec(B);
  destroyPETScMat(M);
  destroyPETScVec(RS);

  delete PolynomSet,ShepardShapeSet,WFuncSet;
}

/**********************************************************************/
/**********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first order derivations.
void
OrthoShapeFunc::calcShapes(InputFileData* InputData,int& supportSize,
                           intVector& sPtcls,std::vector<Particle>& ptcls,
                           double& x,double& y,double& z,dbVector& shapeFuncs,
                           dbVector& xDerivShapes,dbVector& yDerivShapes,
                           dbVector& zDerivShapes,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile,PetscViewer& viewerSEQ) {

  using namespace std;

  // The size of the linear equation system.
  int linEQSize,modLinEQSize;
  int particlesNum = ptcls.size();

  // Calculate a regular polynom basis set and the derivation for 
  // all particles.
  InputData->setValue("dimensionlessBasisPolynom",0);

  BasisPolynom* RegPolySet = new BasisPolyRegular(InputData,ptcls,sPtcls,x,y,z,
                                                  supportSize,linEQSize,1,
                                                  modelData,logFile);

  modLinEQSize = linEQSize - 1;

  if(modLinEQSize < 1) {
    logFile << "In OrthoShapeFuncSet::calcShapes a constant basis\n"
        << "polynomial is not supported!" << endl;
    cerr << "In OrthoShapeFuncSet::calcShapes a constant basis\n"
        << "polynomial is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Check the particle support
  if(modLinEQSize > supportSize) {
    logFile << "In OrthoShapeFuncSet::calcShapes not enough particles support "
        << "current point! (existing: " << supportSize << " - needed: "
        << modLinEQSize << ")" << endl;
    cerr << "In OrthoShapeFuncSet::calcShapes not enough particles support "
        << "current point! (existing: " << supportSize << " - needed: "
        << modLinEQSize << ")" << endl;

    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate a Shepard shapefunction set.
  ShapefunctionSet* ShepardShapeSet = new ShepardShapeFunc(InputData,
                                                           supportSize,sPtcls,
                                                           ptcls,x,y,z,1,
                                                           modelData,logFile,
                                                           viewerSEQ);
  dbVector& shepardShapes = ShepardShapeSet->getShapefunctions();
  dbMatrix& firstDerivShepards = ShepardShapeSet->getFirstDerivShapes();

  // Calculate the actual basis and its first derivations.
  BasisPolynom* PolynomSet = new BasisPolyOrtho(InputData,RegPolySet,
                                                ShepardShapeSet,ptcls,sPtcls,x,
                                                y,z,supportSize,1,modelData,
                                                logFile,viewerSEQ);

  delete RegPolySet;

  dbMatrix& P = PolynomSet->getBasis();
  dbMatrix& dPx = PolynomSet->getXDerivBasis();
  dbMatrix& dPy = PolynomSet->getYDerivBasis();
  dbMatrix& dPz = PolynomSet->getZDerivBasis();

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new PrismaticWindowFunctionSet(InputData,ptcls,
                                                               sPtcls,x,y,z,
                                                               supportSize,1,
                                                               modelData,
                                                               logFile);

  dbVector& W = WFuncSet->getWindowFuncs();
  dbVector& dWx = WFuncSet->getXDerivWinFuncs();
  dbVector& dWy = WFuncSet->getYDerivWinFuncs();
  dbVector& dWz = WFuncSet->getZDerivWinFuncs();

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"********** ortho shapefunctions set ******************"<<endl;
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

  /********************************************************************/
  // Assembly the momentum matrix 'M'
  Mat M;
  intVector matIdx(modLinEQSize);
  dbMatrix MValues(modLinEQSize,dbVector(modLinEQSize));

  MatCreateSeqDense(PETSC_COMM_SELF,modLinEQSize,modLinEQSize,PETSC_NULL, &M);
  MatSetFromOptions(M);

  for(int i = 1;i < linEQSize;i++) {
    matIdx[i - 1] = i - 1;

    for(int j = 1;j < linEQSize;j++)
      
      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++)
        MValues[i - 1][j - 1] += P[k][i] * P[k][j] * W[k];

  }

#ifdef _geometryDebugMode_
  logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(MValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<MValues[i][j];
    logFile<<endl;
  }
#endif

  for(int i = 0;i < modLinEQSize;i++)
    MatSetValues(M,1, &i,modLinEQSize, &matIdx[0], &MValues[i][0],
                 INSERT_VALUES);

  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
  MatOption mOption = MAT_SYMMETRIC;
  PetscBool flag = PETSC_TRUE;
  MatSetOption(M,mOption,flag);
  //  MatView(M,viewerSEQ);

  // Assemble the right side vector.
  Vec RS;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &RS);
  intVector vecIdx(modLinEQSize);

  for(int i = 0;i < modLinEQSize;i++) {
    vecIdx[i] = i;
  }
  
  VecSetValues(RS,modLinEQSize, &vecIdx[0], &P[supportSize][1],INSERT_VALUES);

  // Create vector B.
  Vec B;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &B);
  
  // Solve the linear equation system 'MB = RS' and get vector B.
  KSP ksp;
  PC pc;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
  KSPGetPC(ksp, &pc);
  PCSetType(pc,"lu"); // if one want to use lu-factorizing PCSetType(pc,"lu");
  //PCLUSetUseInPlace(pc); // destroy the original matrix
  KSPSetType(ksp,"preonly"); // direct solving
  KSPSetFromOptions(ksp);
  
  KSPSolve(ksp,RS,B);

  // Check whether the solving was successful or not.
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  //  VecView(B,viewerSEQ);

  // Calculate the shape functions.
  double* BCalc;
  VecGetArray(B, &BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    shapeFuncs[i] = 0;

    // a_j p_j
    for(int j = 1;j < linEQSize;j++)
      shapeFuncs[i] += P[i][j] * BCalc[j - 1];

    shapeFuncs[i] *= W[i];

    // Shepard part
    shapeFuncs[i] += shepardShapes[i];
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************ Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif

  /********************************************************************/
  // x-derivation
  // Assemble the x-derivation of momentum matrix 'M'.
  dbMatrix dMValues(modLinEQSize,dbVector(modLinEQSize));

  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      dMValues[i - 1][j - 1] = 0;

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        dMValues[i - 1][j - 1] += (dPx[k][i] * P[k][j] + P[k][i] * dPx[k][j])
              * W[k] + P[k][i] * P[k][j] * dWx[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(dMValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<dMValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dB = -dM B.
  Vec dBx;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBx);

  dbVector rightSide(modLinEQSize);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPx[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] -= dMValues[i][j] * BCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'MdB = dP - dMB' and get vector B.
  KSPSolve(ksp,RS,dBx);
  //  VecView(sBx,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the x-derivation of the shape functions.
  double* dBCalc;
  VecGetArray(dBx, &dBCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    xDerivShapes[i] = 0;

    // da_j p_j + a_j dp_j
    for(int j = 1;j < linEQSize;j++)
      xDerivShapes[i] += dPx[i][j] * BCalc[j - 1] * W[i]
                                                      + P[i][j] * dBCalc[j - 1] * W[i] + P[i][j] * BCalc[j - 1] * dWx[i];

    // Shepard part
    xDerivShapes[i] += firstDerivShepards[0][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<dBCalc[i]<<endl;
  logFile<<"********** x-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<xDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dBx, &dBCalc);
  destroyPETScVec(dBx);

  /********************************************************************/
  // y-derivation
  // Assemble the y-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      dMValues[i - 1][j - 1] = 0;
      
      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        dMValues[i - 1][j - 1] += (dPy[k][i] * P[k][j] + P[k][i] * dPy[k][j])
              * W[k] + P[k][i] * P[k][j] * dWy[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(dMValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<dMValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dB = dP - dM B.
  Vec dBy;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBy);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPy[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] -= dMValues[i][j] * BCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'MdB = RS' and get vector B.
  KSPSolve(ksp,RS,dBy);
  //  VecView(dBy,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the x-derivation of the shape functions.
  VecGetArray(dBy, &dBCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    yDerivShapes[i] = 0;

    // da_j p_j + a_j dp_j
    for(int j = 1;j < linEQSize;j++)
      yDerivShapes[i] += dPy[i][j] * BCalc[j - 1] * W[i]
                                                      + P[i][j] * dBCalc[j - 1] * W[i] + P[i][j] * BCalc[j - 1] * dWy[i];

    // Shepard part
    yDerivShapes[i] += firstDerivShepards[1][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<dBCalc[i]<<endl;
  logFile<<"************* y-Abl. d. Ansatzfunktion **************"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<yDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dBy, &dBCalc);
  destroyPETScVec(dBy);

  /********************************************************************/
  // z-derivation
  // Assemble the z-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      dMValues[i - 1][j - 1] = 0;
      
      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        dMValues[i - 1][j - 1] += (dPz[k][i] * P[k][j] + P[k][i] * dPz[k][j])
              * W[k] + P[k][i] * P[k][j] * dWz[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(dMValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<dMValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dB = dP - dM B.
  Vec dBz;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBz);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPz[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] -= dMValues[i][j] * BCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'MdB = RS' and get vector B.
  KSPSolve(ksp,RS,dBz);
  //  VecView(B,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the z-derivation of the shape functions.
  VecGetArray(dBz, &dBCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    zDerivShapes[i] = 0;

    // da_j p_j + a_j dP_j
    for(int j = 1;j < linEQSize;j++)
      zDerivShapes[i] += dPz[i][j] * BCalc[j - 1] * W[i]
                                                      + P[i][j] * dBCalc[j - 1] * W[i] + P[i][j] * BCalc[j - 1] * dWz[i];

    // Shepard part
    zDerivShapes[i] += firstDerivShepards[2][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<dBCalc[i]<<endl;
  logFile<<"************ z-Abl. d. Ansatzfunktionen **************"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<zDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dBz, &dBCalc);
  destroyPETScVec(dBz);

  /********************************************************************/
  // Destroy all petsc objects.
  destroyPETScSolver(ksp);

  VecRestoreArray(B, &BCalc);
  destroyPETScVec(B);
  destroyPETScVec(RS);
  destroyPETScMat(M);

  delete PolynomSet,ShepardShapeSet,WFuncSet;
}

/**********************************************************************/
/**********************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions and their first and second order derivations.
void
OrthoShapeFunc::calcShapes(InputFileData* InputData,int& supportSize,
                           intVector& sPtcls,std::vector<Particle>& ptcls,
                           double& x,double& y,double& z,dbVector& shapeFuncs,
                           dbVector& xDerivShapes,dbVector& yDerivShapes,
                           dbVector& zDerivShapes,dbVector& xxDerivShapes,
                           dbVector& yyDerivShapes,dbVector& zzDerivShapes,
                           dbVector& xyDerivShapes,dbVector& yzDerivShapes,
                           dbVector& zxDerivShapes,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile,PetscViewer& viewerSEQ) {

  using namespace std;

  // The size of the linear equation system.
  int linEQSize,modLinEQSize;
  int particlesNum = ptcls.size();

  // Calculate a regular polynom base set and the derivation for 
  // all particles.
  InputData->setValue("dimensionlessBasisPolynom",0);

  BasisPolynom* RegPolySet = new BasisPolyRegular(InputData,ptcls,sPtcls,x,y,z,
                                                  supportSize,linEQSize,2,
                                                  modelData,logFile);

  modLinEQSize = linEQSize - 1;

  if(modLinEQSize < 1) {
    logFile << "In OrthoShapeFuncSet::calcShapes a constant basis\n"
        << "polynomial is not supported!" << endl;
    cerr << "In OrthoShapeFuncSet::calcShapes a constant basis\n"
        << "polynomial is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

#ifdef _geometryDebugMode_
  logFile<<"#######################################################"<<endl;
  logFile<<"**************** ortho shape function set *************"<<endl;
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
  logFile<<"************* regular basis polynom set **************"<<endl;
  dbMatrix& PReg = RegPolySet->getBasis();
  dbMatrix& dPRegx = RegPolySet->getXDerivBasis();
  dbMatrix& dPRegy = RegPolySet->getYDerivBasis();
  dbMatrix& dPRegz = RegPolySet->getZDerivBasis();
  dbMatrix& dPRegxx = RegPolySet->getXXDerivBasis();
  dbMatrix& dPRegyy = RegPolySet->getYYDerivBasis();
  dbMatrix& dPRegzz = RegPolySet->getZZDerivBasis();
  dbMatrix& dPRegxy = RegPolySet->getXYDerivBasis();
  dbMatrix& dPRegyz = RegPolySet->getYZDerivBasis();
  dbMatrix& dPRegzx = RegPolySet->getZXDerivBasis();
  logFile<<"Whole size "<<PReg.size()<<endl;
  for(int i=0;i<=supportSize;i++) {
    if(i != supportSize)
      logFile<<"Ptcle "<<sPtcls[i]<<" Size "<<PReg[i].size()<<endl;
    else
      logFile<<"Point "<<PReg[i].size()<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<PReg[i][j]<<endl;
    logFile<<"*** dPRegx ****"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegx[i][j]<<endl;
    logFile<<"*** dPRegy ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegy[i][j]<<endl;
    logFile<<"*** dPRegz ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegz[i][j]<<endl;
    logFile<<"*** dPRegxx ****"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegxx[i][j]<<endl;
    logFile<<"*** dPRegyy ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegyy[i][j]<<endl;
    logFile<<"*** dPRegzz ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegzz[i][j]<<endl;
    logFile<<"*** dPRegxy ****"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegxy[i][j]<<endl;
    logFile<<"*** dPRegyz ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegyz[i][j]<<endl;
    logFile<<"*** dPRegzx ***"<<endl;
    for(int j=0;j<linEQSize;j++)
      logFile<<i<<" "<<j<<" "<<dPRegzx[i][j]<<endl;
  }
#endif

  // Check the particle support
  if(modLinEQSize > supportSize) {
    logFile << "In OrthoShapeFuncSet::calcShapes not enough particles support "
        << "current point! (existing: " << supportSize << " - needed: "
        << modLinEQSize << ")" << endl;
    cerr << "In OrthoShapeFuncSet::calcShapes not enough particles support "
        << "current point! (existing: " << supportSize << " - needed: "
        << modLinEQSize << ")" << endl;

    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate a Shepard shapefunction set.
  ShapefunctionSet* ShepardShapeSet = new ShepardShapeFunc(InputData,
                                                           supportSize,sPtcls,
                                                           ptcls,x,y,z,2,
                                                           modelData,logFile,
                                                           viewerSEQ);
  dbVector& shepardShapes = ShepardShapeSet->getShapefunctions();
  dbMatrix& firstDerivShepards = ShepardShapeSet->getFirstDerivShapes();
  dbMatrix& secondDerivShepards = ShepardShapeSet->getSecondDerivShapes();

  // Calculate the actual used basis, its first and second order 
  // derivations.
  BasisPolynom* PolynomSet = new BasisPolyOrtho(InputData,RegPolySet,
                                                ShepardShapeSet,ptcls,sPtcls,x,
                                                y,z,supportSize,2,modelData,
                                                logFile,viewerSEQ);

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

  delete RegPolySet;

  // Calculate the window functions and their derivation for each
  // particle which this gauss points supports.
  WindowFunctionSet* WFuncSet = new PrismaticWindowFunctionSet(InputData,ptcls,
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
  logFile<<"#######################################################"<<endl;
  logFile<<"************ orthogonal basis polynom set *************"<<endl;
  logFile<<"Whole size "<<P.size()<<endl;
  for(int i=0;i<=supportSize;i++) {
    if(i != supportSize)
      logFile<<"Ptcle "<<sPtcls[i]<<" Size "<<P[i].size()<<endl;
    else
      logFile<<"Point "<<P[i].size()<<endl;
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
  // Assembly the momentum matrix 'M'
  Mat M;
  intVector matIdx(modLinEQSize);
  dbMatrix MValues(modLinEQSize,dbVector(modLinEQSize));

  MatCreateSeqDense(PETSC_COMM_SELF,modLinEQSize,modLinEQSize,PETSC_NULL, &M);
  MatSetFromOptions(M);

  for(int i = 1;i < linEQSize;i++) {
    matIdx[i - 1] = i - 1;

    for(int j = 1;j < linEQSize;j++)
      
      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++)
        MValues[i - 1][j - 1] += P[k][i] * P[k][j] * W[k];

  }

#ifdef _geometryDebugMode_
  logFile<<"********** Momentenmatrix der Wichtungsfunktion ********"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(MValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<MValues[i][j];
    logFile<<endl;
  }
#endif

  for(int i = 0;i < modLinEQSize;i++)
    MatSetValues(M,1, &i,modLinEQSize, &matIdx[0], &MValues[i][0],
                 INSERT_VALUES);

  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);
  
  MatOption mOption = MAT_SYMMETRIC;
  PetscBool flag = PETSC_TRUE;
  MatSetOption(M,mOption,flag);
  //MatView(M,viewerSEQ);

  // Assemble the right side vector.
  Vec RS;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &RS);
  intVector vecIdx(modLinEQSize);

  for(int i = 0;i < modLinEQSize;i++) {
    vecIdx[i] = i;
  }
  
  VecSetValues(RS,modLinEQSize, &vecIdx[0], &P[supportSize][1],INSERT_VALUES);

  // Create vector B.
  Vec B;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &B);
  
  // Solve the linear equation system 'MB = RS' and get vector B.
  KSP ksp;
  PC pc;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
  KSPGetPC(ksp, &pc);
  PCSetType(pc,"lu"); // if one want to use lu-factorizing PCSetType(pc,"lu");
  //PCLUSetUseInPlace(pc); // destroy the original matrix
  KSPSetType(ksp,"preonly"); // direct solving
  KSPSetFromOptions(ksp);
  
  KSPSolve(ksp,RS,B);

  // Check whether the solving was successful or not.
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!(" << reason << ")" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!(" << reason << ")" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  //  VecView(B,viewerSEQ);

  // Calculate the shape functions.
  double* BCalc;
  VecGetArray(B, &BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    shapeFuncs[i] = 0;

    // a_j p_j
    for(int j = 1;j < linEQSize;j++)
      shapeFuncs[i] += P[i][j] * BCalc[j - 1];

    shapeFuncs[i] *= W[i];

    // Shepard part
    shapeFuncs[i] += shepardShapes[i];
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************ Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<shapeFuncs[i]<<endl;
#endif

  /********************************************************************/
  // x-derivation
  // Assemble the x-derivation of momentum matrix 'M'.
  dbMatrix dxMValues(modLinEQSize,dbVector(modLinEQSize));

  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        dxMValues[i - 1][j - 1] += (dPx[k][i] * P[k][j] + P[k][i] * dPx[k][j])
              * W[k] + P[k][i] * P[k][j] * dWx[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(dxMValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<dxMValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dB = dP -dM B.
  Vec dBx;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBx);

  dbVector rightSide(modLinEQSize);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPx[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] -= dxMValues[i][j] * BCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'MdB = dP - dMB' and get vector B.
  KSPSolve(ksp,RS,dBx);
  //  VecView(sBx,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the x-derivation of the shape functions.
  double* dxBCalc;
  VecGetArray(dBx, &dxBCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    xDerivShapes[i] = 0;

    // da_j p_j + a_j dp_j
    for(int j = 1;j < linEQSize;j++)
      xDerivShapes[i] += dPx[i][j] * BCalc[j - 1] * W[i]
                                                      + P[i][j] * dxBCalc[j - 1] * W[i] + P[i][j] * BCalc[j - 1] * dWx[i];

    // Shepard part
    xDerivShapes[i] += firstDerivShepards[0][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<dxBCalc[i]<<endl;
  logFile<<"********** x-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<xDerivShapes[i]<<endl;
#endif

  /********************************************************************/
  // y-derivation
  dbMatrix dyMValues(modLinEQSize,dbVector(modLinEQSize));

  // Assemble the y-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      dyMValues[i - 1][j - 1] = 0;
      
      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        dyMValues[i - 1][j - 1] += (dPy[k][i] * P[k][j] + P[k][i] * dPy[k][j])
              * W[k] + P[k][i] * P[k][j] * dWy[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(dyMValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<dyMValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dB = dP - dM B.
  Vec dBy;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBy);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPy[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] -= dyMValues[i][j] * BCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'MdB = RS' and get vector B.
  KSPSolve(ksp,RS,dBy);
  //  VecView(dBy,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the x-derivation of the shape functions.
  double* dyBCalc;
  VecGetArray(dBy, &dyBCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    yDerivShapes[i] = 0;

    // da_j p_j + a_j dp_j
    for(int j = 1;j < linEQSize;j++)
      yDerivShapes[i] += dPy[i][j] * BCalc[j - 1] * W[i]
                                                      + P[i][j] * dyBCalc[j - 1] * W[i] + P[i][j] * BCalc[j - 1] * dWy[i];

    // Shepard part
    yDerivShapes[i] += firstDerivShepards[1][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<dyBCalc[i]<<endl;
  logFile<<"************* y-Abl. d. Ansatzfunktion **************"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<yDerivShapes[i]<<endl;
#endif

  /********************************************************************/
  // z-derivation
  dbMatrix dzMValues(modLinEQSize,dbVector(modLinEQSize));

  // Assemble the z-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      dzMValues[i - 1][j - 1] = 0;
      
      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        dzMValues[i - 1][j - 1] += (dPz[k][i] * P[k][j] + P[k][i] * dPz[k][j])
              * W[k] + P[k][i] * P[k][j] * dWz[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(dzMValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<dzMValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dB = dP - dM B.
  Vec dBz;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBz);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPz[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] -= dzMValues[i][j] * BCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'MdB = RS' and get vector B.
  KSPSolve(ksp,RS,dBz);
  //  VecView(B,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the z-derivation of the shape functions.
  double* dzBCalc;
  VecGetArray(dBz, &dzBCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    zDerivShapes[i] = 0;

    // da_j p_j + a_j dP_j
    for(int j = 1;j < linEQSize;j++)
      zDerivShapes[i] += dPz[i][j] * BCalc[j - 1] * W[i]
                                                      + P[i][j] * dzBCalc[j - 1] * W[i] + P[i][j] * BCalc[j - 1] * dWz[i];

    // Shepard part
    zDerivShapes[i] += firstDerivShepards[2][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<dzBCalc[i]<<endl;
  logFile<<"************ z-Abl. d. Ansatzfunktionen **************"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<zDerivShapes[i]<<endl;
#endif
  
  /********************************************************************/
  // xx-derivation
  // Assemble the xx-derivation of momentum matrix 'M'.
  dbMatrix d2MValues(modLinEQSize,dbVector(modLinEQSize));

  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        d2MValues[i - 1][j - 1] += (dPxx[k][i] * P[k][j]
                                                      + 2.0 * dPx[k][i] * dPx[k][j] + P[k][i] * dPxx[k][j]) * W[k]
                                                                                                                + (2.0 * dPx[k][i] * P[k][j] + 2.0 * P[k][i] * dPx[k][j]) * dWx[k]
                                                                                                                                                                                + P[k][i] * P[k][j] * dWxx[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dxx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
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
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBxx);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPxx[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] = rightSide[i] - d2MValues[i][j] * BCalc[j]
                                                            - 2.0 * dxMValues[i][j] * dxBCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'M dBxx = RS' and get vector dBxx.
  KSPSolve(ksp,RS,dBxx);
  //  VecView(sBx,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the xx-derivation of the shape functions.
  double* d2BCalc;
  VecGetArray(dBxx, &d2BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    xxDerivShapes[i] = 0;

    for(int j = 1;j < linEQSize;j++)

      xxDerivShapes[i] = xxDerivShapes[i]

                                       + (d2BCalc[j - 1] * P[i][j] + 2.0 * dxBCalc[j - 1] * dPx[i][j]
                                                                                                   + +BCalc[j - 1] * dPxx[i][j]) * W[i]

                                                                                                                                     + 2.0 * (dxBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPx[i][j]) * dWx[i]

                                                                                                                                                                                                         + BCalc[j - 1] * P[i][j] * dWxx[i];

    xxDerivShapes[i] += secondDerivShepards[0][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<d2BCalc[i]<<endl;
  logFile<<"********** xx-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<xxDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dBxx, &d2BCalc);
  destroyPETScVec(dBxx);

  /********************************************************************/
  // yy-derivation
  // Assemble the yy-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      d2MValues[i - 1][j - 1] = 0;

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        d2MValues[i - 1][j - 1] += (dPyy[k][i] * P[k][j]
                                                      + 2.0 * dPy[k][i] * dPy[k][j] + P[k][i] * dPyy[k][j]) * W[k]
                                                                                                                + (2.0 * dPy[k][i] * P[k][j] + 2.0 * P[k][i] * dPy[k][j]) * dWy[k]
                                                                                                                                                                                + P[k][i] * P[k][j] * dWyy[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dyy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
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
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dByy);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPyy[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] = rightSide[i] - d2MValues[i][j] * BCalc[j]
                                                            - 2.0 * dyMValues[i][j] * dyBCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'M dByy = RS' and get vector dByy.
  KSPSolve(ksp,RS,dByy);
  //  VecView(sByy,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the yy-derivation of the shape functions.
  VecGetArray(dByy, &d2BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    yyDerivShapes[i] = 0;

    for(int j = 1;j < linEQSize;j++)

      yyDerivShapes[i] = yyDerivShapes[i]

                                       + (d2BCalc[j - 1] * P[i][j] + 2.0 * dyBCalc[j - 1] * dPy[i][j]
                                                                                                   + +BCalc[j - 1] * dPyy[i][j]) * W[i]

                                                                                                                                     + 2.0 * (dyBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPy[i][j]) * dWy[i]

                                                                                                                                                                                                         + BCalc[j - 1] * P[i][j] * dWyy[i];

    yyDerivShapes[i] += secondDerivShepards[1][i];

  }

#ifdef _geometryDebugMode_
  logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<d2BCalc[i]<<endl;
  logFile<<"********** yy-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<yyDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dByy, &d2BCalc);
  destroyPETScVec(dByy);

  /********************************************************************/
  // zz-derivation
  // Assemble the zz-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      d2MValues[i - 1][j - 1] = 0;

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        d2MValues[i - 1][j - 1] += (dPzz[k][i] * P[k][j]
                                                      + 2.0 * dPz[k][i] * dPz[k][j] + P[k][i] * dPzz[k][j]) * W[k]
                                                                                                                + (2.0 * dPz[k][i] * P[k][j] + 2.0 * P[k][i] * dPz[k][j]) * dWz[k]
                                                                                                                                                                                + P[k][i] * P[k][j] * dWzz[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dzz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
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
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBzz);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPzz[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] = rightSide[i] - d2MValues[i][j] * BCalc[j]
                                                            - 2.0 * dzMValues[i][j] * dzBCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'M dBzz = RS' and get vector dBzz.
  KSPSolve(ksp,RS,dBzz);
  //  VecView(sBzz,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the zz-derivation of the shape functions.
  VecGetArray(dBzz, &d2BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    zzDerivShapes[i] = 0;

    for(int j = 1;j < linEQSize;j++)

      zzDerivShapes[i] = zzDerivShapes[i]

                                       + (d2BCalc[j - 1] * P[i][j] + 2.0 * dzBCalc[j - 1] * dPz[i][j]
                                                                                                   + BCalc[j - 1] * dPzz[i][j]) * W[i]

                                                                                                                                    + 2.0 * (dzBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPz[i][j]) * dWz[i]

                                                                                                                                                                                                        + BCalc[j - 1] * P[i][j] * dWzz[i];

    zzDerivShapes[i] += secondDerivShepards[2][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<d2BCalc[i]<<endl;
  logFile<<"********** zz-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<zzDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dBzz, &d2BCalc);
  destroyPETScVec(dBzz);

  /********************************************************************/
  // xy-derivation
  // Assemble the xy-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      d2MValues[i - 1][j - 1] = 0;

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        d2MValues[i - 1][j - 1] += (dPxy[k][i] * P[k][j] + dPx[k][i] * dPy[k][j]
                                                                              + dPy[k][i] * dPx[k][j] + P[k][i] * dPxy[k][j]) * W[k]

                                                                                                                                  + (dPx[k][i] * P[k][j] + P[k][i] * dPx[k][j]) * dWy[k]

                                                                                                                                                                                      + (dPy[k][i] * P[k][j] + P[k][i] * dPy[k][j]) * dWx[k]

                                                                                                                                                                                                                                          + P[k][i] * P[k][j] * dWxy[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dxy-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
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
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBxy);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPxy[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] = rightSide[i] - d2MValues[i][j] * BCalc[j]
                                                            - dxMValues[i][j] * dyBCalc[j] - dyMValues[i][j] * dxBCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'M dBxy = RS' and get vector dBxy.
  KSPSolve(ksp,RS,dBxy);
  //  VecView(sBxy,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the xy-derivation of the shape functions.
  VecGetArray(dBxy, &d2BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    xyDerivShapes[i] = 0;

    for(int j = 1;j < linEQSize;j++)

      xyDerivShapes[i] = xyDerivShapes[i]

                                       + (d2BCalc[j - 1] * P[i][j] + dxBCalc[j - 1] * dPy[i][j]
                                                                                             + dyBCalc[j - 1] * dPx[i][j] + BCalc[j - 1] * dPxy[i][j]) * W[i]

                                                                                                                                                           + (dxBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPx[i][j]) * dWy[i]

                                                                                                                                                                                                                         + (dyBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPy[i][j]) * dWx[i]

                                                                                                                                                                                                                                                                                       + BCalc[j - 1] * P[i][j] * dWxy[i];

    xyDerivShapes[i] += secondDerivShepards[3][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<d2BCalc[i]<<endl;
  logFile<<"********** xy-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<xyDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dBxy, &d2BCalc);
  destroyPETScVec(dBxy);

  /********************************************************************/
  // yz-derivation
  // Assemble the yz-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      d2MValues[i - 1][j - 1] = 0;

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        d2MValues[i - 1][j - 1] += (dPyz[k][i] * P[k][j] + dPy[k][i] * dPz[k][j]
                                                                              + dPz[k][i] * dPy[k][j] + P[k][i] * dPyz[k][j]) * W[k]

                                                                                                                                  + (dPy[k][i] * P[k][j] + P[k][i] * dPy[k][j]) * dWz[k]

                                                                                                                                                                                      + (dPz[k][i] * P[k][j] + P[k][i] * dPz[k][j]) * dWy[k]

                                                                                                                                                                                                                                          + P[k][i] * P[k][j] * dWyz[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dyz-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(d2MValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<d2MValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dByz = dPyz -dMyz B - dMy dBz - dMz dBy.
  Vec dByz;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dByz);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPxy[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] = rightSide[i] - d2MValues[i][j] * BCalc[j]
                                                            - dyMValues[i][j] * dzBCalc[j] - dzMValues[i][j] * dyBCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'M dByz = RS' and get vector dByz.
  KSPSolve(ksp,RS,dByz);
  //  VecView(sByz,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the xy-derivation of the shape functions.
  VecGetArray(dByz, &d2BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".
  for(int i = 0;i < supportSize;i++) {
    yzDerivShapes[i] = 0;

    for(int j = 1;j < linEQSize;j++)

      yzDerivShapes[i] = yzDerivShapes[i]

                                       + (d2BCalc[j - 1] * P[i][j] + dyBCalc[j - 1] * dPz[i][j]
                                                                                             + dzBCalc[j - 1] * dPy[i][j] + BCalc[j - 1] * dPyz[i][j]) * W[i]

                                                                                                                                                           + (dyBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPy[i][j]) * dWz[i]

                                                                                                                                                                                                                         + (dzBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPz[i][j]) * dWy[i]

                                                                                                                                                                                                                                                                                       + BCalc[j - 1] * P[i][j] * dWyz[i];

    yzDerivShapes[i] += secondDerivShepards[4][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<d2BCalc[i]<<endl;
  logFile<<"********** yz-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<yzDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dByz, &d2BCalc);
  destroyPETScVec(dByz);

  /********************************************************************/
  // zx-derivation
  // Assemble the zx-derivation of momentum matrix 'M'.
  for(int i = 1;i < linEQSize;i++) {

    for(int j = 1;j < linEQSize;j++) {
      d2MValues[i - 1][j - 1] = 0;

      // Only a loop over the supported particles, since in the other 
      // case the weight function is "0".
      for(int k = 0;k < supportSize;k++) {
        d2MValues[i - 1][j - 1] += (dPzx[k][i] * P[k][j] + dPz[k][i] * dPx[k][j]
                                                                              + dPx[k][i] * dPz[k][j] + P[k][i] * dPzx[k][j]) * W[k]

                                                                                                                                  + (dPz[k][i] * P[k][j] + P[k][i] * dPz[k][j]) * dWx[k]

                                                                                                                                                                                      + (dPx[k][i] * P[k][j] + P[k][i] * dPx[k][j]) * dWz[k]

                                                                                                                                                                                                                                          + P[k][i] * P[k][j] * dWzx[k];
      }
    }
  }

#ifdef _geometryDebugMode_
  logFile<<"******* dzx-Momentenmatrix der Wichtungsfunktion ******"<<endl;
  for(int i=0;i<modLinEQSize;i++) {
    for(int j=0;j<modLinEQSize;j++)
      if(d2MValues[i][j]<0.00000000001)
        logFile<<" 0.000000";
      else
        logFile<<" "<<d2MValues[i][j];
    logFile<<endl;
  }
#endif

  // Assemble the right side vector 'RS' of the linear eq. system 
  // M dBzx = dPzx - dMzx B - dMz dBx - dMx dBz.
  Vec dBzx;
  VecCreateSeq(PETSC_COMM_SELF,modLinEQSize, &dBzx);

  for(int i = 0;i < modLinEQSize;i++) {
    rightSide[i] = dPzx[supportSize][i + 1];

    for(int j = 0;j < modLinEQSize;j++)
      rightSide[i] = rightSide[i] - d2MValues[i][j] * BCalc[j]
                                                            - dzMValues[i][j] * dxBCalc[j] - dxMValues[i][j] * dzBCalc[j];
  }

  VecSetValues(RS,modLinEQSize, &vecIdx[0], &rightSide[0],INSERT_VALUES);

  // Solve the linear equation system 'M dBzx = RS' and get vector dBzx.
  KSPSolve(ksp,RS,dBzx);
  //  VecView(sBxy,viewerSEQ);
  
  // Check whether the solving was successful or not.
  if(reason < 1) {
    logFile << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    cerr << "In OrthoShapeFunc::calcShapes calculation of moment "
        << "matrix of\n  the weight function failed!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Calculate the zx-derivation of the shape functions.
  VecGetArray(dBzx, &d2BCalc);

  // Only a loop over the supported particles, since in the other 
  // case the weight function is "0".*/
  for(int i = 0;i < supportSize;i++) {
    zxDerivShapes[i] = 0;

    for(int j = 1;j < linEQSize;j++)

      zxDerivShapes[i] = zxDerivShapes[i]

                                       + (d2BCalc[j - 1] * P[i][j] + dzBCalc[j - 1] * dPx[i][j]
                                                                                             + dxBCalc[j - 1] * dPz[i][j] + BCalc[j - 1] * dPzx[i][j]) * W[i]

                                                                                                                                                           + (dzBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPz[i][j]) * dWx[i]

                                                                                                                                                                                                                         + (dxBCalc[j - 1] * P[i][j] + BCalc[j - 1] * dPx[i][j]) * dWz[i]

                                                                                                                                                                                                                                                                                       + BCalc[j - 1] * P[i][j] * dWzx[i];

    zxDerivShapes[i] += secondDerivShepards[5][i];
  }

#ifdef _geometryDebugMode_
  logFile<<"************** d2b-Koeffizienten ****** **********"<<endl;
  for(int i=0;i<modLinEQSize;i++)
    logFile<<i<<" "<<d2BCalc[i]<<endl;
  logFile<<"********** zx-Abl. d. Ansatzfunktionen **********"<<endl;
  for(int i=0;i<supportSize;i++)
    logFile<<i<<" "<<sPtcls[i]<<" "<<zxDerivShapes[i]<<endl;
#endif

  VecRestoreArray(dBzx, &d2BCalc);
  destroyPETScVec(dBzx);

  /********************************************************************/
  // Destroy all petsc objects.
  destroyPETScSolver(ksp);

  VecRestoreArray(B, &BCalc);
  VecRestoreArray(dBx, &dxBCalc);
  VecRestoreArray(dBy, &dyBCalc);
  VecRestoreArray(dBz, &dzBCalc);

  destroyPETScVec(B);
  destroyPETScVec(dBx);
  destroyPETScVec(dBy);
  destroyPETScVec(dBz);

  destroyPETScVec(RS);
  destroyPETScMat(M);

  delete PolynomSet,ShepardShapeSet,WFuncSet;

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
