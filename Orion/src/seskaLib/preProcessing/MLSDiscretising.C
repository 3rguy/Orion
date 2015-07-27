#include "MLSDiscretising.h"

MLSDiscretising::MLSDiscretising(InputFileData* InputData,
                                 std::ofstream& logFile) :
    MLSGaussIntegral(InputData,logFile), MLSPtcleIntegral(InputData,logFile),
    MLSShapeFuncSet(InputData,logFile), ParticleDistribution(InputData,logFile),
    BackgroundMesh(InputData,logFile), maxBoundPtcleSupport(0) {
}

/**********************************************************************/
/**********************************************************************/
// Copy the FEM background mesh, respectively nodal coordinates and
// Gauss integration scheme.
void MLSDiscretising::copyFEData(FEMGeometry* FEData,InputFileData* InputData,
                                 std::map<std::string,double>& modelData,
                                 std::ofstream& logFile) {

  using namespace std;

  int integrationMethod = (int) modelData["integrationMethod"];

  // Choose the integration method
  switch(integrationMethod) {

  // Gaussian quadrature
  case 1:

    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    BackgroundMesh::copyFEData(FEData,InputData,modelData,logFile);

    break;

    // Particle integration
  case 2:

    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    MLSPtcleIntegral::copyFEData(FEData,InputData,modelData,logFile);

    break;

  default:
    cerr << "Chosen numerical integration method isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/**********************************************************************/
/**********************************************************************/
// Determine the necessary influence radii for all particles.
void MLSDiscretising::setInfluenceRadii(InputFileData* InputData,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile,
                                        PetscViewer& viewerMPI,
                                        PetscViewer& viewerSEQ) {

  using namespace std;

  int integrationMethod = (int) modelData["integrationMethod"];
  int windowFuncShape = (int) InputData->getValue("windowFunctionShape");
  int shapeFuncType = (int) InputData->getValue("shapefunctionType");
  int radiusMethod = (int) InputData->getValue("radiusDeterminationAlgorithm");

  // Choose the integration method
  switch(integrationMethod) {

  // Calculation of the influence radii with Gaussian quadrature.
  case 1:

    MLSGaussIntegral::setInfluenceRadii(InputData,modelData,logFile,viewerMPI,
                                        viewerSEQ);

    break;

    // Calculation of the influence radii using particle integration.
  case 2:

    MLSPtcleIntegral::setInfluenceRadii(InputData,modelData,logFile,viewerMPI,
                                        viewerSEQ);

    break;

  default:
    logFile << "In MLSDiscretising::setInfluenceRadii chosen numerical\n"
        << "integration method isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/**********************************************************************/
/**********************************************************************/
// Calculation of the particle interpolants
void MLSDiscretising::calcInterpolants(InputFileData* InputData,
                                       std::map<std::string,double>& calcData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile,
                                       PetscViewer& viewerMPI,
                                       PetscViewer& viewerSEQ) {

  using namespace std;

  int integrationMethod = (int) modelData["integrationMethod"];

  // Choose the integration method
  switch(integrationMethod) {

  // Calculation of the MLS interpolanten with Gaussian quadrature.
  case 1:

    calcGaussIntShapes(InputData,calcData,modelData,logFile,viewerMPI,
                       viewerSEQ);

    break;

    // Calculation of the MLS interpolanten with particle integration.
  case 2:

    calcPtcleIntShapes(InputData,calcData,modelData,logFile,viewerMPI,
                       viewerSEQ);

    break;

  default:
    cerr << "In MLSDiscretising::calcInterpolants chosen numerical\n"
        << "integration method isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/**********************************************************************/
/**********************************************************************/
// Calculation of the MLS interpolanten with Gaussian quadrature.
void MLSDiscretising::calcGaussIntShapes(
    InputFileData* InputData,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Determine for all particles their neighbour particles.
  setPtclePtclsConn(InputData,modelData,logFile);

  // Determine the Gauss point distribution among the processors.
  setGaussDistribution(InputData,modelData,logFile);

  // Determine for all influence spheres their neighbour influence
  // spheres.
  setInflSpheresConn(InputData,logFile);

  // Determine the locally needed particles.
  setLocalPtcls(InputData,modelData,logFile,viewerMPI);

  // Distribute the degrees of freedoms to the local gauss points.
  assignGaussPtsDOF(InputData,modelData,logFile);

  // Calculate for all Gauss points and all particles a vector containing 
  // the calculated shape functions and their derivations for all 
  // supporting particles.
  MLSGaussIntegral::setAllShapeFuncs(InputData,calcData,modelData,logFile,
                                     viewerMPI,viewerSEQ);

  // Check the shape functions and their derivatives.
  if((bool) InputData->getValue("checkShapefunctions")) checkGaussShapes(
      InputData,modelData,logFile);

  // Deallocate all class arrays of particles which are not locally 
  // needed.

  deleteNonlocalPtcls(InputData,logFile);
  
}

/**********************************************************************/
/**********************************************************************/
// Calculation of the MLS interpolanten with particle integration.
void MLSDiscretising::calcPtcleIntShapes(
    InputFileData* InputData,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  // Determine for all particles their neighbour particles.
  setPtclePtclsConn(InputData,modelData,logFile);

  // Create a particle connection list needed for partitioning.
  //if((int)InputData->getValue("partitioning") == 1)
  //setPartitionGraph(InputData,modelData,logFile);

  // Determine the locally needed particles.
  setLocalPtcls(InputData,modelData,logFile,viewerMPI);

  /*********************************************************************/
  // Calculate for all particles a vector containing the calculated
  // shape functions and their derivations for all supporting particles.
  MLSPtcleIntegral::setShapeFuncsOnPtcls(InputData,calcData,modelData,logFile,
                                         viewerMPI,viewerSEQ);

}

/**********************************************************************/
/**********************************************************************/
// Rearrange particle support lists according to the new particles
// vector ordering
void MLSDiscretising::rearrangeSupportLists(
    InputFileData* InputData,intVector& newGlobalPtcls,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int integrationMethod = (int) modelData["integrationMethod"];

  // Choose the integration method
  switch(integrationMethod) {

  // Calculation of the MLS interpolanten with Gaussian quadrature.
  case 1:

    MLSGaussIntegral::rearrangeSupportLists(InputData,newGlobalPtcls,modelData,
                                            logFile);

    break;

    // Calculation of the MLS interpolanten with particle integration.
  case 2:

    MLSPtcleIntegral::rearrangeSupportLists(InputData,newGlobalPtcls,modelData,
                                            logFile);

    break;

  default:
    cerr << "Chosen numerical integration method isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
  
}

/************************************************************************/
/************************************************************************/
// Set a special matrix containing the shape function values (for r=0)
// of all particles essential and natural boundary conditions applied 
// to transform the discrete equation system to satisfy the essential 
// and natural boundary conditions.
void MLSDiscretising::setTransformingShapes(
    InputFileData* InputData,intMatrix& boundaryPtcls,
    intMatrix& localBoundPtcls,intVector& newDOFID,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {
  
  using namespace std;
  
  int usedDims = (int) modelData["usedDimensions"];
  int usedDOF = (int) modelData["usedDegreesOfFreedom"];
  map<string,double>& problemData = InputData->getProblemData();

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  intVector& localStartIdx = allExclLocalPtclsStartIdx;
  intVector& allLocalPtcls = allExclLocalPtcls;
  intVector& allGlobalLocalIdx = allGlobalExclLocalPtcleIdx;

#ifdef _forceVecModificationDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"######### set modified boundary shapefunctions #######"<<endl;
  logFile<<"******************************************************"<<endl;
  logFile<<"ptclsPortion="<<exclusiveLocalPtcls.size()<<endl;
  logFile<<"particlesNum="<<particlesNum<<endl;
  logFile<<"maxPtcleSupport="<<maxPtcleSupport<<endl;
  logFile<<"maxIndirectPtcleSupport="<<maxIndirectPtcleSupport<<endl;
  logFile<<"**************** all local particles *****************"<<endl;
  int n=0;
  for(int i=0;i<allLocalPtcls.size();i++) {
    if(i == localStartIdx[n]) {
      logFile<<"--------------------------------------------------"<<endl;
      n++;
    }
    logFile<<i<<".) ptcle "<<allLocalPtcls[i]<<endl;
  }
  logFile<<"**************** all global->local idx ***************"<<endl;
  for(int i=0;i<allGlobalLocalIdx.size();i++)
  logFile<<i<<" -> "<<allGlobalLocalIdx[i]<<endl;
  logFile<<"***************** portioning *************************"<<endl;
  for(int i=0;i<usedDOF;i++)
  logFile<<"DOF "<<i<<": bPtclePortion="
  <<localBoundPtcls[i].size()<<endl;
  logFile<<"************** local boundary particles **************"<<endl;
  for(int i=0;i<usedDOF;i++) {
    logFile<<"DOF "<<i<<endl;
    for(int j=0;j<localBoundPtcls[i].size();j++)
    logFile<<j<<".) "<<localBoundPtcls[i][j]<<endl;
  }
#endif
  
  /**********************************************************************/
  // Modify the boundary shape functions of the boundary particles by
  // inverse the matrix all shapes of all boundary particles stored.
  int ptclePortion = exclusiveLocalPtcls.size();

  intVector bIdx;
  intVector nBIdx;

  intVector bSIdx;
  intVector nBSIdx;
  
  intVector matIdx(maxPtcleSupport);
  dbVector matValues(maxPtcleSupport);
  dbVector XGlobal;

  intVector neededBoundPtcls;
  intVector dummy;
  intVector diagEntries(ptclePortion);
  intVector offDiagEntries(ptclePortion);

  int m;
  int pos;

  // Loop over all used degrees of freedom.
  for(int i = 0;i < usedDOF;i++) {
    
    // set storage scheme for the parallel coefficient matrix

    clearArray(diagEntries);
    clearArray(offDiagEntries);
    
    // loop over exclusively local particles
    for(int l = 0;l < exclusiveLocalPtcls.size();l++) {

      int& ptcle = exclusiveLocalPtcls[l];

      // boundary particle
      if(newDOFID[ptcle * usedDOF + i] < 0) {
        intVector& sPtcls = particles[ptcle].getSupportPtcls();

        for(int k = 0;k < sPtcls.size();k++) {

          // diagonal entry
          if(allGlobalLocalIdx[sPtcls[k]] >= localStartIdx[rank]
            && allGlobalLocalIdx[sPtcls[k]]
              < localStartIdx[rank] + ptclePortion)

          diagEntries[l]++;

          // off-diagonal entry
          else

          offDiagEntries[l]++;

        }

      }

      // non-boundary particle
      else

      diagEntries[l]++;

    }
    
#ifdef _forceVecModificationDebugMode_
    logFile<<"DOF "<<i<<": "<<endl;
    logFile<<"******** coefficient matrix storage scheme *********"<<endl;
    logFile<<endl;
    for(int j=0;j<diagEntries.size();j++)
    logFile<<j<<".) ptcle "<<j+localStartIdx[rank]<<": diagEntries: "
    <<diagEntries[j]<<" offDiagEntries: "<<offDiagEntries[j]
    <<endl;
    logFile<<"****************************************************"<<endl;
#endif

    // ------------------------------------------------------------------
    // Create the solution vector and the right needed for the 
    // PETSC-Solver.
    Vec X;
    Vec Y;
    createParallelPETScVec(ptclePortion,particlesNum,X);
    createParallelPETScVec(ptclePortion,particlesNum,Y);
    
    Mat A;
    map<string,bool> mOptions;
    
    createSparseParallelPETScMat(ptclePortion,ptclePortion,particlesNum,
                                 particlesNum,0,diagEntries,0,offDiagEntries,
                                 mOptions,A,logFile);

    // ------------------------------------------------------------------
    // allocate modified shape function related vectors 

    // loop over all locally needed particles where essential boundary 
    // conditions are applied

    for(int j = 0;j < localBoundPtcls[i].size();j++) {
      int& ptcle = localBoundPtcls[i][j];

      dbVector& bShps = particles[ptcle].getBShapeFuncs(i);
      dbVector& nBShps = particles[ptcle].getNBShapeFuncs(i);
      bShps.resize(maxIndirectPtcleSupport);
      nBShps.resize(maxIndirectPtcleSupport);
    }
    
    for(int j = 0;j < boundaryPtcls[i].size();j++) {

      int& ptcle = boundaryPtcls[i][j];

      intVector& sBPtcls = particles[ptcle].getSupportBPtcls(i);
      intVector& sNBPtcls = particles[ptcle].getSupportNBPtcls(i);
      sBPtcls.resize(maxIndirectPtcleSupport);
      sNBPtcls.resize(maxIndirectPtcleSupport);
    }

    // ===================================================================
    // loop over a local portion of particles and intialize Mat A row by 
    // row
    for(int j = localStartIdx[rank];j < localStartIdx[rank] + ptclePortion;
        j++) {

      int& ptcle = allLocalPtcls[j];

      // boundary particle
      if(newDOFID[ptcle * usedDOF + i] < 0) {
        intVector& sPtcls = particles[ptcle].getSupportPtcls();
        dbVector& shps = particles[ptcle].getShapeFuncs();

        // Determine the column position of matrix A entries to be stored.
        m = 0;

        for(int k = 0;k < sPtcls.size();k++)

          if(shps[k] != 0) {
            matIdx[m] = allGlobalLocalIdx[sPtcls[k]];
            matValues[m] = shps[k];
            m++;
          }

#ifdef _forceVecModificationDebugMode_
        logFile<<j<<".) Particle "<<ptcle
        <<" matIDs (size="<<sPtcls.size()<<"): "<<endl;
        for(int k=0;k<sPtcls.size();k++)
        logFile<<matIdx[k]<<" "<<matValues[k]<<endl;
        if(sPtcls.size() >
            diagEntries[j-localStartIdx[rank]]+offDiagEntries[j-localStartIdx[rank]] ||
            diagEntries[j-localStartIdx[rank]] < 1)
        logFile<<"mat-row allocation insufficient"<<endl;
#endif

        MatSetValues(A,1, &j,sPtcls.size(), &matIdx[0], &matValues[0],
                     INSERT_VALUES);

      }
      
      // non-boundary particle
      else {

#ifdef _forceVecModificationDebugMode_
        logFile<<j<<".) Particle "<<ptcle<<" matIDs: "<<j<<endl;
#endif

        MatSetValue(A,j,j,1.0,INSERT_VALUES);

      }
      
    }

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    MatSetOption(A,MAT_SYMMETRIC,PETSC_FALSE);
    MatSetOption(A,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE);
    MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE);
    
#ifdef _forceVecModificationDebugMode_
    if(particlesNum <= 512)
    MatView(A,viewerMPI);
#endif

#ifdef _checkMatrixConditioning_
    map<string,double>& problemData = InputData->getProblemData();
    checkMatrixConditioning(A,problemData,logFile);
#endif
    
    // -------------------------------------------------------------------
    // Initialize Solver.
    KSP ksp;
    PC pc;
    KSPConvergedReason reason;
    vector<string> solvingMethod;

    map<string,double> solverData = InputData->getProblemData();
    solverData["initializeKSP"] = 1;

    // iterative solving: Jacobi preconditioning and GMRES solver works best
    if( !(bool) solverData["directSolving"]) {
      //solverData["directSolver"] = 2;
      //solverData["directSolving"] = 0;
      solverData["preconditioner"] = 5; // Jacobi
      solverData["solver"] = 0; // GMRES
    }

    createParallelPETScSolver(A,solverData,solvingMethod,ksp,pc,calcData,
                              modelData,logFile);

    // ===================================================================
    // Make a loop to assemble vector Y and solve the equation system to 
    // invert A.
    intVector vecIdx(2);
    dbVector oneVec(2);
    oneVec[0] = 0.0;
    oneVec[1] = 1.0;
    
    bSIdx.resize(localBoundPtcls[i].size());
    nBSIdx.resize(localBoundPtcls[i].size());
    clearArray(bSIdx);
    clearArray(nBSIdx);

    bIdx.resize(boundaryPtcls[i].size());
    nBIdx.resize(boundaryPtcls[i].size());
    clearArray(bIdx);
    clearArray(nBIdx);

    // Loop over all particles and initialise the right hand side as 
    // identity vector -> inverse of matrix A.
    for(int j = 0;j < allLocalPtcls.size();j++) {
      int& ptcleJ = allLocalPtcls[j];
      
      if(ptcleRootList[ptcleJ] == rank) {

        if(j == 0) VecSetValues(Y,1, &j, &oneVec[1],INSERT_VALUES);

        else {
          vecIdx[0] = j - 1;
          vecIdx[1] = j;
          VecSetValues(Y,2, &vecIdx[0], &oneVec[0],INSERT_VALUES);
        }

      }
      
      VecAssemblyBegin(Y);
      VecAssemblyEnd(Y);
      
      KSPSolve(ksp,Y,X);
      KSPGetConvergedReason(ksp, &reason);

      int solverIterations;
      KSPGetIterationNumber(ksp, &solverIterations);

#ifdef _forceVecModificationDebugMode_
//       int solverIterations;
//       KSPGetIterationNumber(ksp,&solverIterations);
//       logFile<<j<<".) inverse solver iterations="<<solverIterations<<endl; 
//      KSPView(ksp,viewerMPI);
//       VecView(Y,viewerMPI);
//       VecView(X,viewerMPI);
//      MPI_Abort(MPI_COMM_WORLD,1);
#endif

      if(reason < 0) {
        logFile << "In MLSDiscretising::setTransformingShapes inversion of\n"
            << "shape function matrix failed!" << endl;
        writeKSPConvergenceReason(reason,solverIterations,logFile);
        logFile << "------------------------------------------------" << endl;
        logFile << "The particle distribution at a Dirichlet boundary "
            << "surface is bad.\n Check whether there duplicate surfaces."
            << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      getGlobalParallelPETScVec(ptclePortion,particlesNum,X,XGlobal);

      // only locally needed shape function ordinates are stored to 
      // save memory (corresponding rows of N!) 
      // note: inversion of matrix A changes the support structure

      for(int k = 0;k < localBoundPtcls[i].size();k++) {
        int& ptcleK = localBoundPtcls[i][k];
        int& idxK = allGlobalLocalIdx[ptcleK];

        if( !checkZero(XGlobal[idxK])) {

          // boundary particle (j=current column of N)
          if(newDOFID[ptcleJ * usedDOF + i] < 0) {

            dbVector& shps = particles[ptcleK].getBShapeFuncs(i);

            if(shps.size() <= bSIdx[k]) shps.resize(shps.size() * 2);

            shps[bSIdx[k]] = XGlobal[idxK];
            bSIdx[k]++;

          }

          // non-boundary particle (j=current column of N)
          else {

            dbVector& shps = particles[ptcleK].getNBShapeFuncs(i);

            if(shps.size() <= nBSIdx[k]) shps.resize(shps.size() * 2);

            shps[nBSIdx[k]] = XGlobal[idxK];
            nBSIdx[k]++;

          }

        }

      }

      // loop over all global boundary particles and store the particle 
      // support list which is needed to determine the extended particle 
      // support list later
      for(int k = 0;k < boundaryPtcls[i].size();k++) {

        int& ptcleK = boundaryPtcls[i][k];
        int& idxK = allGlobalLocalIdx[ptcleK];

        if( !checkZero(XGlobal[idxK])) {

          // boundary particle (j=current column of N)
          if(newDOFID[ptcleJ * usedDOF + i] < 0) {

            intVector& sPtcls = particles[ptcleK].getSupportBPtcls(i);

            if(sPtcls.size() <= bIdx[k]) sPtcls.resize(sPtcls.size() * 2);

            sPtcls[bIdx[k]] = ptcleJ;
            bIdx[k]++;

          }

          // non-boundary particle (j=current column of N)
          else {

            intVector& sPtcls = particles[ptcleK].getSupportNBPtcls(i);

            if(sPtcls.size() <= nBIdx[k]) sPtcls.resize(sPtcls.size() * 2);

            sPtcls[nBIdx[k]] = ptcleJ;
            nBIdx[k]++;

          }

        }

        if(maxBoundPtcleSupport < bIdx[k] + nBIdx[k]) maxBoundPtcleSupport =
          bIdx[k] + nBIdx[k];

      }

    }

    // -------------------------------------------------------------------
    // free not needed space

    for(int j = 0;j < localBoundPtcls[i].size();j++) {
      int& ptcle = localBoundPtcls[i][j];

      dbVector& bShps = particles[ptcle].getBShapeFuncs(i);
      bShps.resize(bSIdx[j]);
      
      dbVector& nBShps = particles[ptcle].getNBShapeFuncs(i);
      nBShps.resize(nBSIdx[j]);
    }
    
    for(int j = 0;j < boundaryPtcls[i].size();j++) {

      int& ptcle = boundaryPtcls[i][j];

      intVector& sBPtcls = particles[ptcle].getSupportBPtcls(i);
      sBPtcls.resize(bIdx[j]);

      intVector& sNBPtcls = particles[ptcle].getSupportNBPtcls(i);
      sNBPtcls.resize(nBIdx[j]);
    }

    // Destroy all petsc objects.
    destroyPETScSolver(ksp);
    destroyPETScMat(A);
    destroyPETScVec(X);
    destroyPETScVec(Y);
    MPI_Barrier(MPI_COMM_WORLD);
  }

#ifdef _forceVecModificationDebugMode_
  logFile<<"*********************************************************"<<endl;
  logFile<<"************ local modified boundary shapes *************"<<endl;
  logFile<<"*********************************************************"<<endl;
  logFile<<"maxBoundPtcleSupport "<<maxBoundPtcleSupport<<endl;
  for(int i=0;i<usedDOF;i++) {
    logFile<<"DOF "<<i<<" *************************************"<<endl;
    for(int j=0;j<boundaryPtcls[i].size();j++) {
      int& ptcle = boundaryPtcls[i][j];
      double PUM=0;
      logFile<<j<<".) BOUND PARTICLE "<<ptcle<<": ";
      intVector& sBPtcls = particles[ptcle].getSupportBPtcls(i);
      intVector& sNBPtcls = particles[ptcle].getSupportNBPtcls(i);
      dbVector& bShps = particles[ptcle].getBShapeFuncs(i);
      dbVector& nBShps = particles[ptcle].getNBShapeFuncs(i);
      for(int k=0;k<sBPtcls.size();k++)
      logFile<<sBPtcls[k]<<" ";
      logFile<<" || ";
      for(int k=0;k<sNBPtcls.size();k++)
      logFile<<sNBPtcls[k]<<" ";
      logFile<<endl;
      for(int k=0;k<bShps.size();k++) {
        logFile<<bShps[k]<<" ";
      }
      logFile<<" || ";
      for(int k=0;k<nBShps.size();k++) {
        logFile<<nBShps[k]<<" ";
      }
      logFile<<endl;
      //if(sNBPtcls.size() != nBShps.size()) 
      //logFile<<"Different entry number!"<<endl;
      //logFile<<"*********"<<endl;
    }
  }
#endif
#ifdef _testBoundCollocation_
  logFile<<"******************************************************"<<endl;
  logFile<<"********* test of boundary enforcement ***************"<<endl;

  // choose a particle and DOF which is specified as essential boundary 
  // particle in the input file;
  // note: best testing with curved surface

  // C++ conversion:
  // deg = deg -1
  // particle = particle -1
  int deg = 1;
  int particle = 185;

  if(deg >= usedDOF)
  deg = 0;

  if(particle >= particlesNum)
  particle = 0;

  particle = newIdx[particle];

  double one = 1.0e+03;
  double boundCond = 0.0;
  double boundDOF;

  // --------------------------------------------------------------------
  // set the modified degree of freedom vector

  for(int i=0;i<localBoundPtcls[deg].size();i++) {
    int& ptcle = localBoundPtcls[deg][i];
    intVector& sBPtcls = particles[ptcle].getSupportBPtcls(deg);
    intVector& sNBPtcls = particles[ptcle].getSupportNBPtcls(deg);

    for(int j=0;j<sBPtcls.size();j++) {
      dbVector& DOF = particles[sBPtcls[j]].getDOF();
      DOF.resize(usedDOF);
      particles[sBPtcls[j]].setDOF(deg,boundCond);
    }

    for(int j=0;j<sNBPtcls.size();j++) {
      dbVector& DOF = particles[sNBPtcls[j]].getDOF();
      DOF.resize(usedDOF);
      particles[sNBPtcls[j]].setDOF(deg,one);
    }

  }

  // --------------------------------------------------------------------
  int boundPtclsSize = boundaryPtcls[deg].size();
  dbVector newLocalDOFs(particlesNum);
  dbVector newDOFs(particlesNum);

  // compute the true boundary DOF; zero boundary function values should 
  // be achieved with modified boundary DOF = 0, since modified boundary 
  // = boundary function value 
  // note: modified non-boundary = true non-boundary DOF

  for(int i=0;i<localBoundPtcls[deg].size();i++) {

    int& ptcle = localBoundPtcls[deg][i];

    if(ptcleRootList[ptcle] != rank) continue;

    intVector& sBPtcls = particles[ptcle].getSupportBPtcls(deg);
    intVector& sNBPtcls = particles[ptcle].getSupportNBPtcls(deg);
    dbVector& bShps = particles[ptcle].getBShapeFuncs(deg);
    dbVector& nBShps = particles[ptcle].getNBShapeFuncs(deg);
    boundDOF = 0;

    for(int j=0;j<sBPtcls.size();j++)
    newLocalDOFs[ptcle] += bShps[j]*particles[sBPtcls[j]].getDOF(deg);

    for(int j=0;j<sNBPtcls.size();j++)
    newLocalDOFs[ptcle] += nBShps[j]*particles[sNBPtcls[j]].getDOF(deg);

  }

  MPI_Allreduce(&newLocalDOFs[0],&newDOFs[0],particlesNum,
      MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  for(int i=0;i<particlesNum;i++) {

    logFile<<"ptcle "<<i<<": ";

    if(newDOFs[i] != 0) {
      particles[i].setDOF(deg,newDOFs[i]);
      logFile<<newDOFs[i];
    }

    logFile<<endl;
  }

  logFile<<"-----------------------------------------------------"<<endl;
  for(int j=0;j<exclusiveLocalPtcls.size();j++)
  logFile<<"ptcle "<<exclusiveLocalPtcls[j]
  <<": "<<particles[exclusiveLocalPtcls[j]].getDOF(deg)<<endl;

  // compute the function value which should be zero, if 'particle' 
  // is a boundary particle with a boundary conditions applied

  if(findIntVecPos(particle,0,localBoundPtcls[deg].size(),
          localBoundPtcls[deg]) != -1) {

    double funcValue = 0;

    intVector& sPtcls = particles[particle].getSupportPtcls();
    dbVector& shps = particles[particle].getShapeFuncs();

    for(int j=0;j<shps.size();j++) {
      funcValue += shps[j]*particles[sPtcls[j]].getDOF(deg);
      logFile<<j<<".) "<<sPtcls[j]<<" "<<shps[j]<<" "<<particles[sPtcls[j]].getDOF(deg)<<endl;
    }

    logFile<<"funcValue = "<<funcValue<<endl;
  }

  if(InputData->getValue("maxNeumannLoadingFactor") > 0 ||
      InputData->getValue("maxDirichletLoadingFactor") > 0) {
    cerr<<"In MLSDiscretising::setTransformingShapes boundary collocation\n"
    <<"testing active, no modelling can't be done. Unset pre-processor\n"
    <<"variable '_testBoundCollocation_' !"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
#endif
  
}

/************************************************************************/
/************************************************************************/
// Assemble all degrees of freedom. 
void MLSDiscretising::getAllPtcleDOF(InputFileData* InputData,dbVector& allDOF,
                                     std::map<std::string,double>& calcData,
                                     std::map<std::string,double>& modelData,
                                     std::ofstream& logFile) {
  
  using namespace std;

  int usedDOF = (int) modelData["usedDegreesOfFreedom"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(allDOF.size() < particlesNum * usedDOF) {
    clearArray(allDOF);
    allDOF.resize(particlesNum * usedDOF);
  }
  else clearArray(allDOF);

  dbVector localAllDOF(exclusiveLocalPtcls.size() * usedDOF);

  // Loop over all local particles to gather their DOF.
  for(int i = 0;i < exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];
    dbVector& DOF = particles[ptcle].getDOF();

    // Loop over all degrees of freedom.
    for(int j = 0;j < usedDOF;j++)
      localAllDOF[i * usedDOF + j] = DOF[j];

  }
  
  // ---------------------------------------------------------------------
  // Merge the local vectors to a global one.

  mergeLocalPtcleVectors(InputData,localAllDOF,usedDOF,allDOF,calcData,
                         modelData,logFile);

#ifdef _simulationDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"******************** particle DOF ********************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"Ptcle "<<i<<" (old: "<<oldIdx[i]<<"): ";
    for(int j=0;j<usedDOF;j++)
    logFile<<allDOF[i*usedDOF+j]<<" ";
    logFile<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Assemble all degrees of freedom. 
void MLSDiscretising::printAllPtcleDOF(InputFileData* InputData,
                                       std::map<std::string,double>& calcData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {
  
  using namespace std;

  int usedDOF = (int) modelData["usedDegreesOfFreedom"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  dbVector allDOF;
  if(allDOF.size() < particlesNum * usedDOF) {
    clearArray(allDOF);
    allDOF.resize(particlesNum * usedDOF);
  }
  else clearArray(allDOF);

  dbVector localAllDOF(exclusiveLocalPtcls.size() * usedDOF);

  // Loop over all local particles to gather their DOF.
  for(int i = 0;i < exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];
    dbVector& DOF = particles[ptcle].getDOF();

    // Loop over all degrees of freedom.
    for(int j = 0;j < usedDOF;j++)
      localAllDOF[i * usedDOF + j] = DOF[j];

  }
  
  // ---------------------------------------------------------------------
  // Merge the local vectors to a global one.

  mergeLocalPtcleVectors(InputData,localAllDOF,usedDOF,allDOF,calcData,
                         modelData,logFile);

#ifdef _simulationDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"******************** particle DOF ********************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"Ptcle "<<i<<" (old: "<<oldIdx[i]<<"): ";
    for(int j=0;j<usedDOF;j++)
    logFile<<allDOF[i*usedDOF+j]<<" ";
    logFile<<endl;
  }
#endif

}

// /************************************************************************/
// /************************************************************************/
// // Assemble all step degrees of freedom. 
// void MLSDiscretising::getAllPtcleDeltaDOF(InputFileData* InputData,
// 					  dbVector& allStepDOF,
// 					  std::map<std::string,double>& calcData,
// 					  std::map<std::string,double>& modelData,
// 					  std::ofstream& logFile) {

//   using namespace std;

//   int usedDOF = (int)modelData["usedDegreesOfFreedom"];

//   int rank,size; 
//   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
//   MPI_Comm_size(MPI_COMM_WORLD,&size);

//   if(allStepDOF.size() < particlesNum*usedDOF) {
//     clearArray(allStepDOF);
//     allStepDOF.resize(particlesNum*usedDOF);
//   }
//   else
//     clearArray(allStepDOF);

//   dbVector localAllDOF(exclusiveLocalPtcls.size()*usedDOF);

//   // Loop over all local particles to gather their DOF.
//   for(int i=0;i<exclusiveLocalPtcls.size();i++) {
//     int& ptcle = exclusiveLocalPtcls[i];
//     dbVector& deltaDOF = particles[ptcle].getStepDOF();

//     if(deltaDOF.size() < usedDOF)
//       deltaDOF.resize(usedDOF);

//     // Loop over all degrees of freedom.
//     for(int j=0;j<usedDOF;j++)
//       localAllDOF[i*usedDOF+j] = deltaDOF[j];

//   }

//   // ---------------------------------------------------------------------
//   // Merge the local vectors to a global one.

//   mergeLocalPtcleVectors(InputData,localAllDOF,usedDOF,allStepDOF,calcData,
// 			 modelData,logFile);

// #ifdef _simulationDebugMode_
//   logFile<<"######################################################"<<endl;
//   logFile<<"****************** particle Step DOF ****************"<<endl;
//   for(int i=0;i<particlesNum;i++) {
//     logFile<<"Ptcle "<<i<<" (old: "<<oldIdx[i]<<"): ";
//     for(int j=0;j<usedDOF;j++)
//       logFile<<allStepDOF[i*usedDOF+j]<<" ";
//     logFile<<endl;
//   }
// #endif

// }

/************************************************************************/
/************************************************************************/
// Assemble all step degrees of freedom. 
void MLSDiscretising::getAllPtcleStepDOF(
    InputFileData* InputData,dbVector& allStepDOF,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {
  
  using namespace std;

  int usedDOF = (int) modelData["usedDegreesOfFreedom"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(allStepDOF.size() < particlesNum * usedDOF) {
    clearArray(allStepDOF);
    allStepDOF.resize(particlesNum * usedDOF);
  }
  else clearArray(allStepDOF);

  dbVector localAllDOF(exclusiveLocalPtcls.size() * usedDOF);

  // Loop over all local particles to gather their DOF.
  for(int i = 0;i < exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];
    dbVector& stepDOF = particles[ptcle].getStepDOF();

    if(stepDOF.size() < usedDOF) stepDOF.resize(usedDOF);

    // Loop over all degrees of freedom.
    for(int j = 0;j < usedDOF;j++)
      localAllDOF[i * usedDOF + j] = stepDOF[j];

  }
  
  // ---------------------------------------------------------------------
  // Merge the local vectors to a global one.

  mergeLocalPtcleVectors(InputData,localAllDOF,usedDOF,allStepDOF,calcData,
                         modelData,logFile);

#ifdef _simulationDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"****************** particle Step DOF *****************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"Ptcle "<<i<<" (old: "<<oldIdx[i]<<"): ";
    for(int j=0;j<usedDOF;j++)
    logFile<<allStepDOF[i*usedDOF+j]<<" ";
    logFile<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Assemble all deformation degrees of freedom. 
void MLSDiscretising::getAllPtcleDefDOF(InputFileData* InputData,
                                        dbVector& allDOF,
                                        std::map<std::string,double>& calcData,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile) {
  
  using namespace std;

  int defDOF = (int) modelData["deformationDegreesOfFreedom"];
  int usedDOF = (int) modelData["usedDegreesOfFreedom"];
  int stressDOF = (int) modelData["stressDegreesOfFreedom"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(allDOF.size() < particlesNum * defDOF) {
    clearArray(allDOF);
    allDOF.resize(particlesNum * defDOF);
  }
  else clearArray(allDOF);

  dbVector localAllDOF(exclusiveLocalPtcls.size() * defDOF);

  // Loop over all local particles to gather their DOF.
  for(int i = 0;i < exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];
    dbVector& DOF = particles[ptcle].getDOF();

    // Loop over all deformation degrees of freedom.
    for(int j = 0;j < defDOF;j++)
      localAllDOF[i * defDOF + j] = DOF[j];

  }
  
  // ---------------------------------------------------------------------
  // Merge the local vectors to a global one.

  mergeLocalPtcleVectors(InputData,localAllDOF,defDOF,allDOF,calcData,modelData,
                         logFile);

#ifdef _simulationDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"******************** particle DOF ********************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"Ptcle "<<i<<" (old: "<<oldIdx[i]<<"): ";
    for(int j=0;j<defDOF;j++)
    logFile<<allDOF[i*defDOF+j]<<" ";
    logFile<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Assemble all displacement degrees of freedom. 
void MLSDiscretising::getAllPtcleDispDOF(
    InputFileData* InputData,dbVector& allDOF,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {
  
  using namespace std;

  int dispDOF = (int) modelData["displacementDegreesOfFreedom"];
  int usedDOF = (int) modelData["usedDegreesOfFreedom"];
  int stressDOF = (int) modelData["stressDegreesOfFreedom"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(allDOF.size() < particlesNum * dispDOF) {
    clearArray(allDOF);
    allDOF.resize(particlesNum * dispDOF);
  }
  else clearArray(allDOF);

  dbVector localAllDOF(exclusiveLocalPtcls.size() * dispDOF);

  // Loop over all local particles to gather their DOF.
  for(int i = 0;i < exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];
    dbVector& DOF = particles[ptcle].getDOF();

    // Loop over all displacement degrees of freedom.
    for(int j = 0;j < dispDOF;j++)
      localAllDOF[i * dispDOF + j] = DOF[j];

  }

  // ---------------------------------------------------------------------
  // Merge the local vectors to a global one.

  mergeLocalPtcleVectors(InputData,localAllDOF,dispDOF,allDOF,calcData,
                         modelData,logFile);
  
#ifdef _simulationDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"******************** particle DOF ********************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"Ptcle "<<i<<" (old: "<<oldIdx[i]<<"): ";
    for(int j=0;j<dispDOF;j++)
    logFile<<allDOF[i*dispDOF+j]<<" ";
    logFile<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Return the indices of all particles any point force loading is applied.
intMatrix& MLSDiscretising::getPointForceBoundPtcleIdx(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile) {
  
  using namespace std;

  int integrationMethod = (int) modelData["integrationMethod"];

  // Choose the integration method
  switch(integrationMethod) {

  // MLS method with Gaussian quadrature.
  case 1:

    return MLSGaussIntegral::getPointForceBoundPtcleIdx();

    break;

    // MLS method with particle integration.
  case 2:

    return MLSPtcleIntegral::getPointForceBoundPtcleIdx();

    break;

  default:
    cerr << "Chosen numerical integration method isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/************************************************************************/
/************************************************************************/
// at a specified distribution of points (reference points) compute the 
// MLS approximation of a given field of matrix values known at a 
// distribution of particles
void MLSDiscretising::matrixFieldMLSApproximation(
    std::vector<Particle>& referencePoints,std::vector<Point>& points,
    std::map<std::string,double>& data,std::ofstream& logFile,
    PetscViewer& viewerSEQ) {
  
  using namespace std;

  int polynomialOrder = (int) data["basisPolynomOrder"];
  int usedDims = (int) data["usedDimensions"];

  // Re-read the input file.
  InputFileData* InputData = new InputFileData(logFile);

  InputData->setValue("radiusDeterminationAlgorithm",1);
  InputData->setValue("windowFunctionType",1);
  InputData->setValue("shapefunctionType",2);
  InputData->setValue("basisPolynomType",1);
  InputData->setValue("basisPolynomOrder",polynomialOrder);

  /*********************************************************************/
  // check the dimension of the field to be approximated is consistent
  // throughout the domain
  int numOfRows,numOfCols;
  int rows,cols;

  for(int i = 0;i < referencePoints.size();i++) {
    Particle& ptcle = referencePoints[i];
    dbMatrix& ptcleMat = ptcle.getMLSmat();

    rows = ptcleMat.size();

    if(i == 0) {
      numOfRows = rows;

      for(int j = 0;j < ptcleMat.size();j++) {

        cols = ptcleMat[j].size();

        if(j == 0) numOfCols = cols;

        else if(numOfCols != cols) {
          logFile << "In MLSDiscretising::matrixFieldMLSApproximation field\n"
              << "to be approximated is not consistent in size." << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }

    }
    else if(numOfRows != rows) {
      logFile << "In MLSDiscretising::matrixFieldMLSApproximation field\n"
          << "to be approximated is not consistent in size." << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

  if(numOfRows == 0 || numOfCols == 0) {

    logFile << "In MLSDiscretising::matrixFieldMLSApproximation no field\n"
        << "to be approximated." << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /*********************************************************************/
  // determine particle influence radii correspondingly to the distances
  // to their neighbouring particles
  //InputData->setValue("minDirectPtcleSuppReduction",0); 
  setPtcleRadsDistDepend(InputData,referencePoints,data,logFile);
  postProcInfRads(InputData,referencePoints,data,logFile);

#ifdef _anisotropyInterpolationDebugMode_
  logFile<<"#################################################"<<endl;
  logFile<<"******** approximation of matrix field **********"<<endl;
  logFile<<"*************************************************"<<endl;
  logFile<<"************** influence radii ******************"<<endl;
  for(int i=0;i<referencePoints.size();i++) {
    Particle& point = referencePoints[i];
    dbVector& coords = point.getCoords();
    logFile<<"POINT "<<i+1<<"(";
    for(int j=0;j<coords.size();j++)
    logFile<<coords[j]<<" ";
    logFile<<") radii: ";
    dbVector& radii = point.getRadii();
    for(int j=0;j<radii.size();j++)
    logFile<<radii[j]<<" ";
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // set supporting particle list for all points and compute the shape
  // function ordinates of all supporting particles
  int support,dummy;
  int maxSupport = 0;

  // loop over all points and determine their supporting particles
  for(int i = 0;i < points.size();i++) {
    Point& point = points[i];
    dbVector& coords = point.getCoords();
    intVector& suppPtcls = point.getSupportPtcls();
    allocateArray(suppPtcls,maxSupport);

    support = 0;

    for(int j = 0;j < referencePoints.size();j++) {

      // Check if current point 'i' is supported by particle 'j'.
      if(referencePoints[j].querySupported(InputData,coords,data,logFile)) {

        if(support < suppPtcls.size()) suppPtcls[support] = j;

        else suppPtcls.push_back(j);

        support++;

      }

      if(maxSupport < support) maxSupport = support;

      resizeArray(suppPtcls,support);
    }

  }

#ifdef _anisotropyInterpolationDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"***************** point support list *****************"<<endl;
  for(int i=0;i<points.size();i++) {
    Point& point = points[i];
    dbVector& coords = point.getCoords();
    intVector& suppPtcls = point.getSupportPtcls();
    dbVector& shapes = point.getShapeFuncs();
    logFile<<"POINT "<<i+1<<"(";
    for(int j=0;j<coords.size();j++)
    logFile<<coords[j]<<" ";
    logFile<<"): ";
    for(int j=0;j<suppPtcls.size();j++)
    logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
#endif

  // loop over all points and determine the particle shape function ordinates 
  for(int i = 0;i < points.size();i++) {
    Point& point = points[i];
    dbVector& coords = point.getCoords();
    intVector& suppPtcls = point.getSupportPtcls();
    dbVector& shapes = point.getShapeFuncs();
    support = suppPtcls.size();

    // compute the shape function ordinates of all supporting particles
    EFGShapeFunc::calcShapes(InputData,support,suppPtcls,referencePoints,
                             coords[0],coords[1],coords[2],shapes,dummy,data,
                             logFile,viewerSEQ);

  }

#ifdef _anisotropyInterpolationDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"***************** point shape functions **************"<<endl;
  for(int i=0;i<points.size();i++) {
    Point& point = points[i];
    dbVector& coords = point.getCoords();
    intVector& suppPtcls = point.getSupportPtcls();
    dbVector& shapes = point.getShapeFuncs();
    logFile<<"POINT "<<i+1<<"(";
    for(int j=0;j<coords.size();j++)
    logFile<<coords[j]<<" ";
    logFile<<"): ";
    for(int j=0;j<suppPtcls.size();j++)
    logFile<<"("<<suppPtcls[j]<<") "<<shapes[j]<<" ";
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // approximate all degrees of freedom at the specified field of points
  for(int i = 0;i < points.size();i++) {
    Point& point = points[i];
    intVector& suppPtcls = point.getSupportPtcls();
    dbVector& shapes = point.getShapeFuncs();
    dbMatrix& pointMat = point.getMLSmat();
    allocateArray(pointMat,numOfRows,numOfCols);
    clearArray(pointMat);

    for(int j = 0;j < suppPtcls.size();j++) {
      Particle& ptcle = referencePoints[suppPtcls[j]];
      dbMatrix& ptcleMat = ptcle.getMLSmat();

      for(int k = 0;k < ptcleMat.size();k++)
        for(int l = 0;l < ptcleMat[k].size();l++)

          pointMat[k][l] += shapes[j] * ptcleMat[k][l];

    }
    
  }
  
}
