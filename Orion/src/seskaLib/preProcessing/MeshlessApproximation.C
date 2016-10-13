#include "MeshlessApproximation.h"

MeshlessApproximation::MeshlessApproximation(
    FEMGeometry* FEData,InputFileData* InputData,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerMPI,PetscViewer& viewerSEQ) :
    BackgroundMesh(InputData,logFile), FEMDiscretising(InputData,logFile),
    MaxEntShapeFunc(InputData,logFile), MaxEntDiscretising(InputData,logFile),
    MLSGaussIntegral(InputData,logFile), MLSPtcleIntegral(InputData,logFile),
    MLSDiscretising(InputData,logFile), MLSShapeFuncSet(InputData,logFile),
    ParticleDistribution(InputData,logFile) {

  using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Copy the FEM background mesh, respectively nodal coordinates and
  // Gauss integration scheme.
  copyFEData(FEData,InputData,modelData,logFile);

#ifdef _sphericalShapeFuncDebugMode_
  vector<Particle> savedParticles = particles;
  intVector savedGlobalLocalElemIdx = globalLocalElemIdx;
  vector<FEMElement> savedNodesElements = nodesElements;
#endif

  /*********************************************************************/
  // Determine for each process a portion of particles and set a root 
  // list, which process each particle belongs to.
  setPtcleProcDistrib(InputData,logFile);

  // Determine the influence radii for all particles.
  setInfluenceRadii(InputData,modelData,logFile,viewerMPI,viewerSEQ);

  // --------------------------------------------------------------------
  // Call of test routines 
  //testSpline(logFile);
  //testShapes(InputData,particles,modelData,logFile,viewerSEQ);
  //testMLS(InputData,particles,modelData,logFile,viewerSEQ);
  //testMaxEnt(InputData,particles,modelData,logFile,viewerSEQ);
  //testUpdate(InputData,logFile);

  /*********************************************************************/
  // Calculation of the interpolanten
  calcInterpolants(InputData,calcData,modelData,logFile,viewerMPI,viewerSEQ);

#ifdef _sphericalShapeFuncDebugMode_
  // compute prismatic shape functions
  InputData->setValue("windowFunctionShape",1.0);
  vector<Particle> sphereParticles = particles;
  particles = savedParticles;
  globalLocalElemIdx = savedGlobalLocalElemIdx;
  nodesElements = savedNodesElements;
  logFile<<"###################################################"<<endl;
  logFile<<"compute prismatic shape functions"<<endl;
  logFile<<"###################################################"<<endl;
  setInfluenceRadii(InputData,modelData,logFile,viewerMPI,viewerSEQ);
  calcInterpolants(InputData,calcData,modelData,logFile,viewerMPI,viewerSEQ);
  logFile<<"###################################################"<<endl;
  logFile<<"all computed prismatic shape functions"<<endl;
  logFile<<"###################################################"<<endl;
  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    dbVector& suppShapes = particles[ptcle].getShapeFuncs();
    logFile<<i<<".) PARTICLE "<<ptcle<<" ("<<oldIdx[ptcle]<<"):: "<<endl;
    for(int j=0;j<suppPtcls.size();j++) {
      logFile<<suppPtcls[j]<<": "
      <<suppShapes[j]<<" | ";
    }
    logFile<<endl;
  }
  logFile<<"###################################################"<<endl;
  logFile<<"compare for all particles their neighbour particles"<<endl;
  logFile<<"###################################################"<<endl;
  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& sphereSuppPtcls = sphereParticles[ptcle].getSupportPtcls();
    intVector& prismaSuppPtcls = particles[ptcle].getSupportPtcls();
    string mode = "arbitrary-subvector";
    if(!compareIntVecs(mode,sphereSuppPtcls,prismaSuppPtcls,logFile)) {
      cout<<i<<".) PARTICLE "<<ptcle<<" ("<<oldIdx[ptcle]<<")"<<endl;
      cout<<"mismatch of supporting particles using spherical and prismatic window functions:";
      logFile<<i<<".) PARTICLE "<<ptcle<<" ("<<oldIdx[ptcle]<<")"<<endl;
      logFile<<"mismatch of supporting particles using spherical and prismatic window functions:";
      for(int j=0;j<sphereSuppPtcls.size();j++)
      logFile<<sphereSuppPtcls[j]<<" ";
      logFile<<endl;
      for(int j=0;j<prismaSuppPtcls.size();j++)
      logFile<<prismaSuppPtcls[j]<<" ";
      logFile<<endl;
    }
  }
  MPI_Abort(MPI_COMM_WORLD,1);
  //particles = savedParticles;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Copy the FEM background mesh, respectively nodal coordinates and
// Gauss integration scheme.
void MeshlessApproximation::copyFEData(FEMGeometry* FEData,
                                       InputFileData* InputData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {

  using namespace std;

  int shapeFuncType = (int) InputData->getValue("shapefunctionType");

  if(shapeFuncType == 1 || shapeFuncType == 2 || shapeFuncType == 3
    || shapeFuncType == 4 || shapeFuncType == 5) {

    // Copy the FEM background mesh, respectively nodal coordinates and
    // Gauss integration scheme.
    MLSDiscretising::copyFEData(FEData,InputData,modelData,logFile);
    delete FEData;
  }
  else if(shapeFuncType == 6) {

    MaxEntDiscretising::copyFEData(FEData,InputData,modelData,logFile);
    delete FEData;
  }
  else if(shapeFuncType == 7) {

    FEMDiscretising::copyFEData(FEData,InputData,modelData,logFile);
    //delete FEData; // FEData is needed as various element pointers, 
    //                  e.g. FEVolumeSet would be deleted 
  }
  else {
    logFile << "In MeshlessApproximation::copyFEData shape function\n"
        << "type '" << shapeFuncType << "' is not supported. Check input file."
        << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/**********************************************************************/
/**********************************************************************/
// Determine the necessary influence radii for all particles.
void MeshlessApproximation::setInfluenceRadii(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile,PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {
  
  using namespace std;

  int shapeFuncType = (int) InputData->getValue("shapefunctionType");

  if(shapeFuncType == 1 || shapeFuncType == 2 || shapeFuncType == 3
    || shapeFuncType == 4 || shapeFuncType == 5)

  // Determine the influence radii for all particles for MLS.
  MLSDiscretising::setInfluenceRadii(InputData,modelData,logFile,viewerMPI,
                                     viewerSEQ);

  else if(shapeFuncType == 6)

  MaxEntDiscretising::setInfluenceRadii(InputData,modelData,logFile,viewerMPI,
                                        viewerSEQ);

  else if(shapeFuncType != 7) {
    logFile << "In MeshlessApproximation::setInfluenceRadii shape function\n"
        << "type '" << shapeFuncType << "' is not supported. Check input file."
        << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/**********************************************************************/
/**********************************************************************/
// Calculation of the interpolanten
void MeshlessApproximation::calcInterpolants(
    InputFileData* InputData,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  int shapeFuncType = (int) InputData->getValue("shapefunctionType");

  if(shapeFuncType == 1 || shapeFuncType == 2 || shapeFuncType == 3
    || shapeFuncType == 4 || shapeFuncType == 5) {

    // Calculation of MLS interpolanten
    MLSDiscretising::calcInterpolants(InputData,calcData,modelData,logFile,
                                      viewerMPI,viewerSEQ);
  }
  else if(shapeFuncType == 6) {

    MaxEntDiscretising::calcInterpolants(InputData,calcData,modelData,logFile,
                                         viewerMPI,viewerSEQ);
  }
  else if(shapeFuncType == 7) {

    FEMDiscretising::calcInterpolants(InputData,calcData,modelData,logFile,
                                      viewerMPI,viewerSEQ);
  }
  else {
    logFile << "In MeshlessApproximation::calcInterpolants shape function\n"
        << "type '" << shapeFuncType << "' is not supported.\n"
        << "Check input file." << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/************************************************************************/
/************************************************************************/
// Assemble all deformation degrees of freedom. 
void MeshlessApproximation::getAllPtcleDefDOF(
    InputFileData* InputData,dbVector& allDOF,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  // Assemble all deformation degrees of freedom of an MLS approximation. 
  MLSDiscretising::getAllPtcleDefDOF(InputData,allDOF,calcData,modelData,
                                     logFile);

}

