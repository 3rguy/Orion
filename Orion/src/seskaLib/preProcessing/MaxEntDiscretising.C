#include "MaxEntDiscretising.h"

MaxEntDiscretising::MaxEntDiscretising(InputFileData* InputData,
                                       std::ofstream& logFile) :
    BackgroundMesh(InputData,logFile), MaxEntShapeFunc(InputData,logFile),
    MLSGaussIntegral(InputData,logFile), MLSPtcleIntegral(InputData,logFile),
    MLSDiscretising(InputData,logFile), MLSShapeFuncSet(InputData,logFile),
    ParticleDistribution(InputData,logFile) {
}

/**********************************************************************/
/**********************************************************************/
// Calculation of the particle interpolants
void MaxEntDiscretising::calcInterpolants(
    InputFileData* InputData,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Shape Function Test
  //    shapeFuncTest(InputData,modelData,logFile,viewerSEQ);

  //    MPI_Abort(MPI_COMM_WORLD,1);

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

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // NOTE : If approxBetaValues is de-activated, setBetaOnPtcls should
  // be enabled.

  // Calculate the MLS shape funtions and its derivatives, to be used
  // to approximate Beta and delta Beta
//   calcMlsInterpolants(InputData,calcData,modelData,logFile,viewerMPI,
// 		      viewerSEQ);

//   // Approximate the value of Beta using MLS shape functions
//   approxBetaValues(InputData,modelData,logFile);

//   // Approximate the derivates value of Beta using MLS shape functions
//   approxBetaDerivValues(InputData,modelData,logFile);

//   for(int i=0;i<particlesNum;i++) {
//     MPI_Bcast(&particles[i].getBeta(),1,
// 	      MPI_DOUBLE,ptcleRootList[i],
// 	      MPI_COMM_WORLD);

//     MPI_Bcast(&particles[i].getBetaDerivs()[0],3,
// 	      MPI_DOUBLE,ptcleRootList[i],
// 	      MPI_COMM_WORLD);
//   }

//   // Clear shape Functions and derivations of shape functions of
//   // each particle
//   clearShapeFuncs(InputData,modelData,logFile);

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  // Calculate the beta values for each particle
  setBetaOnPtcls(InputData,modelData,logFile);

  // Calculate for all Gauss points and all particles a vector containing
  // the calculated shape functions and their derivations for all
  // supporting particles.
  setAllShapeFuncs(InputData,modelData,logFile,viewerMPI,viewerSEQ);

  // Deallocate all class arrays of particles which are not locally
  // needed.
  deleteNonlocalPtcls(InputData,logFile);

}

/**********************************************************************/
/**********************************************************************/
// Calculate the beta values for each particle
void MaxEntDiscretising::setBetaOnPtcls(InputFileData* InputData,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile) {

  using namespace std;

  double gamma = (double) InputData->getValue("MaxEntGamma");
  double beta;

  for(int p = 0;p < particlesNum;p++) {

    beta = gamma / pow(particles[p].getLocalMinPtcleDistance(),2);

    // Assign Beta value to particle.
    particles[p].setBeta(beta);
  }
}

/**********************************************************************/
/**********************************************************************/
// Calculate for all Gauss points and all particles a vector containing 
// the calculated shape functions and their derivatives for all 
// supporting particles.
void MaxEntDiscretising::setAllShapeFuncs(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile,PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  // Calculate for all particles a vector containing the calculated
  // shape functions and their derivations for all supporting particles.

  setShapeFuncsOnPtcls(InputData,modelData,logFile,viewerMPI,viewerSEQ);

  // Calculate for all inner and boundary gauss points a vector
  // containing the calculated shape functions and their derivations
  // for all supporting particles.

  setShapeFuncsOnGauss(InputData,modelData,logFile,viewerMPI,viewerSEQ);

  setShapeFuncsOnBGauss(InputData,modelData,logFile,viewerMPI,viewerSEQ);

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all gauss points a vector containing the calculated
// shape functions and their derivations for all supporting particles.
void MaxEntDiscretising::setShapeFuncsOnGauss(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile,PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  unsigned int derivationOrder =
    (unsigned int) modelData["shapesDerivationOrderOnIntPoints"];

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  dbMatrix secondDerivShapes(6);

  int supportSize;

  // Loop over all local gauss points to set their shape function
  // vectors.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int maxSupportSize = 0;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"******************* gausspoint shapes ********************"<<endl;
#endif

  for(int i = 0;i < localGaussPtsNum;i++) {
    supportSize = gaussPoints[i].getSupportCounts();

    if(supportSize > maxSupportSize) {
      shapeFuncs.resize(supportSize);

      if(derivationOrder > 0)

      for(int j = 0;j < 3;j++)
        firstDerivShapes[j].resize(supportSize);

      if(derivationOrder > 1)

      for(int j = 0;j < 6;j++)
        secondDerivShapes[j].resize(supportSize);

    }

#ifdef _geometryDebugMode_
    logFile<<"Gauss point "<<gaussPoints[i].getGlobalID()<<endl;
#endif

    double& x = gaussPoints[i].getCoord(0);
    double& y = gaussPoints[i].getCoord(1);
    double& z = gaussPoints[i].getCoord(2);

    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,gaussPoints[i].getSupportPtcls(),particles,x,y,z,
                     shapeFuncs,gaussPoints[i].getLocalMinPtcleDistance(),
                     modelData,logFile,viewerSEQ);

      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);

    }
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,gaussPoints[i].getSupportPtcls(),particles,x,y,z,
                     shapeFuncs,firstDerivShapes,
                     gaussPoints[i].getLocalMinPtcleDistance(),modelData,
                     logFile,viewerSEQ);

      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      gaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else {
      logFile << "In MaxEntDiscretising::setShapeFuncsOnGauss shape function "
          << "derivations higher than second order are not supported!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

#ifdef _geometryDebugMode_
  double PUM;
  logFile<<"******************************************************"<<endl;
  logFile<<"******** calculated shape sets on gauss pts **********"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    PUM = 0;
    logFile<<"GAUSSPOINT "<<gaussPoints[i].getGlobalID()<<":: ";
    for(int j=0;j<gaussPoints[i].getSupportCounts();j++) {
      PUM += gaussPoints[i].getShapeFunc(j);
      logFile<<gaussPoints[i].getSupportPtcle(j)<<": "
      <<gaussPoints[i].getShapeFunc(j)<<" | ";
    }
    if(fabs(1.0-PUM) > 0.000000001)
    logFile<<" PUM= "<<PUM<<endl;
    else
    logFile<<endl;
    if(derivationOrder > 0) {
      dbMatrix& firstDerivs = gaussPoints[i].getFirstDerivShapes();
      logFile<<"first order derivations"<<endl;
      for(int k=0;k<3;k++) {
        for(int j=0;j<firstDerivs[k].size();j++)
        logFile<<gaussPoints[i].getSupportPtcle(j)<<": "
        <<firstDerivs[k][j]<<" |";
        logFile<<endl;
      }
    }
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all particles a vector containing the calculated
// shape functions and their derivations for all supported particles.
void MaxEntDiscretising::setShapeFuncsOnPtcls(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile,PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int derivationOrder = (int) modelData["shapesDerivationOrderOnPtcls"];

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"********************* particle shapes ********************"<<endl;
  logFile<<"**********************************************************"<<endl;
  logFile<<"derivation order="<<derivationOrder<<endl;
#endif

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  int supportSize;

  // Loop over processor's locally needed particles portion.
  for(int i = 0;i < localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    supportSize = suppPtcls.size();

    if(shapeFuncs.size() < supportSize) shapeFuncs.resize(supportSize);

    if(derivationOrder > 0)

    for(int j = 0;j < firstDerivShapes.size();j++)
      firstDerivShapes[j].resize(supportSize);

    double& x = particles[ptcle].getCoord(0);
    double& y = particles[ptcle].getCoord(1);
    double& z = particles[ptcle].getCoord(2);

    // compute the shape functions
    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,suppPtcls,particles,x,y,z,shapeFuncs,
                     particles[ptcle].getLocalMinPtcleDistance(),modelData,
                     logFile,viewerSEQ);

      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);

    }

    // compute the shape functions and their first derivatives
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,suppPtcls,particles,x,y,z,shapeFuncs,
                     firstDerivShapes,
                     particles[ptcle].getLocalMinPtcleDistance(),modelData,
                     logFile,viewerSEQ);

      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);
      particles[ptcle].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else {
      logFile << "In MaxEntDiscretising::setShapeFuncsOnPtcls shape function\n "
          << "derivations higher than second order are not supported!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

#ifdef _geometryDebugMode_
  int zeroCounter;
  double PUM;
  logFile<<"******************************************************"<<endl;
  logFile<<"*********** calculated shape sets on ptcls ***********"<<endl;
  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    dbVector& suppShapes = particles[ptcle].getShapeFuncs();
    dbMatrix& firstDerivs = particles[i].getFirstDerivShapes();
    PUM = 0;
    logFile<<"Ptcle "<<ptcle<<"("<<oldIdx[ptcle]<<"):: ";
    for(int j=0;j<suppPtcls.size();j++) {
      PUM += suppShapes[j];
      logFile<<suppPtcls[j]<<": "
      <<suppShapes[j]<<" | ";
    }
    if(fabs(1.0-PUM) > 1.0e-13)
    logFile<<" PUM= "<<PUM<<endl;
    if(derivationOrder > 0) {
      zeroCounter = 0;
      logFile<<"first order derivations"<<endl;
      for(int k=0;k<firstDerivs.size();k++) {
        for(int j=0;j<firstDerivs[k].size();j++)
        logFile<<suppPtcls[j]<<": "<<firstDerivs[k][j]<<" | ";
        logFile<<endl;
      }
      logFile<<"zeroCounter = "<<zeroCounter<<" of "<<suppPtcls.size()*3<<endl;
    }
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all boundary gauss points a vector containing the 
// calculated shape functions for all supported particles.
void MaxEntDiscretising::setShapeFuncsOnBGauss(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile,PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"*************** boundary gausspoint shapes ***************"<<endl;
#endif

  unsigned int derivationOrder =
    (unsigned int) modelData["shapesDerivationOrderOnBoundIntPoints"];

  int supportSize;

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);

  int m = 0;
  int localMaxSupportSize = 0;

  // Loop over all local boundary gauss points to set their shape
  // function vectors.
  for(int i = 0;i < localBGaussPtsNum;i++) {
    supportSize = boundGaussPoints[i].getSupportCounts();

    if(supportSize > localMaxSupportSize) {
      localMaxSupportSize = supportSize;
      shapeFuncs.resize(supportSize);

      if(derivationOrder > 0)

      for(int j = 0;j < 3;j++)
        firstDerivShapes[j].resize(supportSize);

    }

#ifdef _geometryDebugMode_
    logFile<<"Boundary gauss point "<<i<<endl;
#endif

    double& x = boundGaussPoints[i].getCoord(0);
    double& y = boundGaussPoints[i].getCoord(1);
    double& z = boundGaussPoints[i].getCoord(2);

    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,boundGaussPoints[i].getSupportPtcls(),particles,
                     x,y,z,shapeFuncs,
                     boundGaussPoints[i].getLocalMinPtcleDistance(),modelData,
                     logFile,viewerSEQ);

      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
    }
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,boundGaussPoints[i].getSupportPtcls(),particles,
                     x,y,z,shapeFuncs,firstDerivShapes,
                     boundGaussPoints[i].getLocalMinPtcleDistance(),modelData,
                     logFile,viewerSEQ);

      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      boundGaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else {
      logFile << "In MaxEntDiscretising::setShapeFuncsOnBGauss shape function "
          << "derivations higher than second order are not supported!" << endl;

      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

#ifdef _geometryDebugMode_
  double PUM;
  logFile<<"######################################################"<<endl;
  logFile<<"** calculated shapes on boundary integration points **"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    PUM = 0;
    logFile<<"Bound GAUSSPOINT "<<i<<": support = "
    <<boundGaussPoints[i].getSupportCounts()<<" ";
    for(int j=0;j<boundGaussPoints[i].getSupportCounts();j++) {
      PUM += boundGaussPoints[i].getShapeFunc(j);
      logFile<<boundGaussPoints[i].getSupportPtcle(j)<<": "<<boundGaussPoints[i].getShapeFunc(j)<<" | ";
    }
    if(fabs(1.0-PUM) > 0.000000001)
    logFile<<"PUM= "<<PUM<<endl;
    else
    logFile<<endl;
    if(derivationOrder > 0) {
      dbMatrix& firstDerivs = boundGaussPoints[i].getFirstDerivShapes();
      logFile<<"first order derivations"<<endl;
      for(int k=0;k<3;k++) {
        for(int j=0;j<firstDerivs[k].size();j++)
        logFile<<boundGaussPoints[i].getSupportPtcle(j)<<": "
        <<firstDerivs[k][j]<<" |";
        logFile<<endl;
      }
    }
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Calculate at a certain point for all its supporting particles their 
// shape functions.
void MaxEntDiscretising::calcShapeFuncs(InputFileData* InputData,
                                        intVector& sPtcls,
                                        std::vector<Particle>& particles,
                                        double& x,double& y,double& z,
                                        dbVector& shapeFuncs,
                                        double& MinPtcleDist,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile,
                                        PetscViewer& viewerSEQ) {

  using namespace std;

  int shapeFuncType = (int) InputData->getValue("shapefunctionType");

  switch(shapeFuncType) {

  case 6:

    MaxEntShapeFunc::calcShapes(InputData,sPtcls,particles,x,y,z,shapeFuncs,
                                modelData,logFile,viewerSEQ);

    break;

  default:
    logFile << "In MaxEntDiscretising::calcShapeFuncs chosen shape function\n"
        << "type isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/************************************************************************/
/************************************************************************/
// Calculate at a certain point for all its supporting particles their 
// shape functions and their first order derivations.
void MaxEntDiscretising::calcShapeFuncs(InputFileData* InputData,
                                        intVector& sPtcls,
                                        std::vector<Particle>& particles,
                                        double& x,double& y,double& z,
                                        dbVector& shapeFuncs,
                                        dbMatrix& firstDerivShapes,
                                        double& MinPtcleDist,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile,
                                        PetscViewer& viewerSEQ) {

  using namespace std;

  int shapeFuncType = (int) InputData->getValue("shapefunctionType");

  switch(shapeFuncType) {

  case 6:

    MaxEntShapeFunc::calcShapes(InputData,sPtcls,particles,x,y,z,shapeFuncs,
                                firstDerivShapes,modelData,logFile,viewerSEQ);

    break;

  default:
    logFile << "In MaxEntDiscretising::calcShapeFuncs chosen shape function\n"
        << "type isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
}

/************************************************************************/
/************************************************************************/
// Calculate the shape functions and derivatives of the shape functions
// using MLS
void MaxEntDiscretising::calcMlsInterpolants(
    InputFileData* InputData,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerMPI,PetscViewer& viewerSEQ) {

  using namespace std;

  int defaultShapeFuncType = (int) InputData->getValue("shapefunctionType");
  int defaultRadDeterAlgorithm = (int) InputData->getValue(
      "radiusDeterminationAlgorithm");

  InputData->setValue("shapefunctionType",2.0);
  InputData->setValue("radiusDeterminationAlgorithm",4.0);

  int integrationMethod = (int) modelData["integrationMethod"];

  // Choose the integration method
  switch(integrationMethod) {

  // Calculation of the MLS interpolanten with Gaussian quadrature.
  case 1:

    MLSGaussIntegral::setShapeFuncsOnPtcls(InputData,calcData,modelData,logFile,
                                           viewerMPI,viewerSEQ);

    // Check the shape functions and their derivatives.
    if((bool) InputData->getValue("checkShapefunctions")) checkGaussShapes(
        InputData,modelData,logFile);

    break;

  default:
    cerr << "In MaxEntDiscretising::calcMlsInterpolants chosen numerical\n"
        << "integration method isn't supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  InputData->setValue("shapefunctionType",defaultShapeFuncType);
  InputData->setValue("radiusDeterminationAlgorithm",defaultRadDeterAlgorithm);
}

/************************************************************************/
/************************************************************************/
// Approximate the value of Beta using MLS shape functions
void MaxEntDiscretising::approxBetaValues(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile) {

  using namespace std;

#ifdef _geometryDebugMode_
  logFile<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  logFile<<"********************* Approximating Beta ******************"<<endl;
  logFile<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
#endif

  double gamma = (double) InputData->getValue("MaxEntGamma");
  double beta;

  for(int p = 0;p < exclusiveLocalPtcls.size();p++) {
    beta = 0;
    intVector supportPtcls =
      particles[exclusiveLocalPtcls[p]].getSupportPtcls();
    for(int sp = 0;sp < supportPtcls.size();sp++) {
      beta += (gamma
        / pow(particles[supportPtcls[sp]].getLocalMinPtcleDistance(),2))
        * particles[exclusiveLocalPtcls[p]].getShapeFuncs()[sp];
    }

    // Assign Beta value to particle.
    particles[exclusiveLocalPtcls[p]].setBeta(beta);

  }

#ifdef _geometryDebugMode_
  for (int p=0;p<exclusiveLocalPtcls.size();p++) {
    beta= 0;
    intVector supportPtcls = particles[exclusiveLocalPtcls[p]].getSupportPtcls();
    logFile<<"For Particle ["<<exclusiveLocalPtcls[p]<<"]: "<<particles[exclusiveLocalPtcls[p]].getBeta()<<endl;
  }
#endif
}

/************************************************************************/
/************************************************************************/
// Approximate the derivative values of Beta using MLS shape functions
void MaxEntDiscretising::approxBetaDerivValues(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile) {

  using namespace std;

#ifdef _geometryDebugMode_
  logFile<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  logFile<<"****************** Approximating delta Beta ***************"<<endl;
  logFile<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
#endif

  double gamma = (double) InputData->getValue("MaxEntGamma");
  int usedDims = (int) modelData["usedDimensions"];
  double beta,betaDerivDim;
  dbVector ptcleBetaDeriv(usedDims);

  for(int p = 0;p < exclusiveLocalPtcls.size();p++) {

    intVector supportPtcls =
      particles[exclusiveLocalPtcls[p]].getSupportPtcls();
    dbMatrix& firstDerivs =
      particles[exclusiveLocalPtcls[p]].getFirstDerivShapes();
    for(int d = 0;d < firstDerivs.size();d++) {
      betaDerivDim = 0;
      for(int sp = 0;sp < firstDerivs[d].size();sp++) {
        beta = (gamma
          / pow(particles[supportPtcls[sp]].getLocalMinPtcleDistance(),2));
        betaDerivDim += firstDerivs[d][sp] * beta;
      }
      ptcleBetaDeriv[d] = betaDerivDim;
    }

    // assigne derivatives of beta to particle.
    particles[exclusiveLocalPtcls[p]].setBetaDerivs(ptcleBetaDeriv);

  }

#ifdef _geometryDebugMode_
  for (int p=0;p<particlesNum;p++) {
    intVector& suppPtcls = particles[p].getSupportPtcls();
    dbMatrix& firstDerivs = particles[p].getFirstDerivShapes();
    logFile<<"Particle["<<p<<"] => first order derivations: "<<endl;
    for(int k=0;k<firstDerivs.size();k++) {
      for(int j=0;j<firstDerivs[k].size();j++)
      logFile<<suppPtcls[j]<<": "<<firstDerivs[k][j]<<" | ";
      logFile<<endl;
    }
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Clear shape Functions and derivations of shape functions of
// each particle
void MaxEntDiscretising::clearShapeFuncs(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(usedDims);

  int supportSize;

  for(int i = 0;i < exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];

    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    supportSize = suppPtcls.size();

    particles[ptcle].getShapeFuncs().clear();
    shapeFuncs.resize(supportSize);
    particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);

    for(int id = 0;id < usedDims;id++) {
      particles[ptcle].getFirstDerivShapes()[id].clear();
      firstDerivShapes[id].resize(supportSize);
    }
    particles[ptcle].setFirstDerivShapes(supportSize,firstDerivShapes);
  }

}

/************************************************************************/
/************************************************************************/
// Output for debugging purposes
void MaxEntDiscretising::debugLogFile(InputFileData* InputData,
                                      std::map<std::string,double>& modelData,
                                      std::ofstream& logFile) {

  using namespace std;

  //    // Print out all the particles' distances
  //    logFile<<"=========================================================="<<endl;
  //    logFile<<"******************** Particle Distances ******************"<<endl;
  //    logFile<<"=========================================================="<<endl;
  //    for (int p=0;p<particlesNum;p++)
  //        logFile<< "Particle["<<p<<"]:\t"
  //               <<particles[p].getLocalMinPtcleDistance()<<endl;

  //    logFile<<"=========================================================="<<endl;
  //    logFile<<"**************** Particle Elements & Nodes ***************"<<endl;
  //    logFile<<"=========================================================="<<endl;
  //     for(int i=0;i<localGlobalPtcls.size();i++) {
  //         intVector elemVec = particles[localGlobalPtcls[i]].getElems();
  //         logFile<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  //         logFile<<"Considering Particle -"<<particles[localGlobalPtcls[i]].getID();
  //         logFile<<" with " <<elemVec.size()<<" alpha elements"<<endl;

  //         logFile << "elements " << elemVec[0] << " with "<<endl;
  ////         logFile <<nodesElements[elemVec[0]].getNodes().size()<< endl;

  //         for(int j=0;j<elemVec.size();j++){
  //             logFile<<"Element ["<<elemVec[j]<<"]: ";
  //             intVector betaNodes = nodesElements[elemVec[j]].getNodes();
  //             for(int k=0;k<betaNodes.size();k++){
  //                 logFile<<betaNodes[k]<<", ";
  //             }
  //             logFile<<endl;
  //         }
  //         logFile<<endl;
  //     }

  //    // Print out all the particles' distances
  //    logFile<<"==========================================================="<<endl;
  //    logFile<<"******************** Gausspoint Distances *****************"<<endl;
  //    logFile<<"==========================================================="<<endl;
  //    for (int g=0;g<localGaussPtsNum;g++)
  //        logFile<< "Particle["<<g<<"]:\t"
  //               <<gaussPoints[g].getLocalMinPtcleDistance()<<endl;

  //    // Print out all the particles' distances
  //    logFile<<"==========================================================="<<endl;
  //    logFile<<"******************* B-Gausspoint Distances ****************"<<endl;
  //    logFile<<"==========================================================="<<endl;
  //    for (int bg=0;bg<localBGaussPtsNum;bg++)
  //        logFile<< "Particle["<<bg<<"]:\t"
  //               << boundGaussPoints[bg].getLocalMinPtcleDistance()<<endl;

  //    // **************************************************************************

  //    //Outputing all shape functions
  //    logFile<<"=========================================================="<<endl;
  //    logFile<<"***************** Particle Shape Functions ***************"<<endl;
  //    logFile<<"=========================================================="<<endl;
  //    for (int p=0;p<particlesNum;p++){
  //        double total = 0;
  //        for (int p_=0;p_<particles[p].getShapeFuncs().size();p_++){
  //        logFile<< "Particle["<<p<<"]:\t"
  //               <<particles[p].getShapeFuncs()[p_]<<endl;
  //        total += particles[p].getShapeFuncs()[p_];
  //        }
  //        logFile<<"\t\t ============>Total: "<<total<<endl;
  //    }

  //    logFile<<"=========================================================="<<endl;
  //    logFile<<"*********** Particle Shape Functions Derivatives *********"<<endl;
  //    logFile<<"=========================================================="<<endl;
  //    for (int p=0;p<particlesNum;p++){
  //        intVector& suppPtcls = particles[p].getSupportPtcls();
  //        dbMatrix& firstDerivs = particles[p].getFirstDerivShapes();
  //        logFile<<"Particle["<<p<<"] => first order derivations: "<<endl;
  //        for(int k=0;k<firstDerivs.size();k++) {
  //            for(int j=0;j<firstDerivs[k].size();j++)
  //                logFile<<suppPtcls[j]<<": "<<firstDerivs[k][j]<<" | ";
  //            logFile<<endl;
  //        }
  //    }

  ////#ifdef _geometryDebugMode_
  //    int derivationOrder = 1;
  //    int zeroCounter;
  //    double PUM;
  //    logFile<<"******************************************************"<<endl;
  //    logFile<<"*********** calculated shape sets on ptcls ***********"<<endl;
  //    for(int i=0;i<localGlobalPtcls.size();i++) {
  //        int& ptcle = localGlobalPtcls[i];
  //        intVector& suppPtcls = particles[ptcle].getSupportPtcls();
  //        dbVector& suppShapes = particles[ptcle].getShapeFuncs();
  //        dbMatrix& firstDerivs = particles[i].getFirstDerivShapes();
  //        dbMatrix& secondDerivs = particles[i].getSecondDerivShapes();
  //        PUM = 0;
  //        logFile<<"Ptcle "<<ptcle<<"("<<oldIdx[ptcle]<<"):: ";
  //        for(int j=0;j<suppPtcls.size();j++) {
  //            PUM += suppShapes[j];
  //            logFile<<suppPtcls[j]<<": "
  //                  <<suppShapes[j]<<" | ";
  //        }
  //        if(fabs(1.0-PUM) > 1.0e-13)
  //            logFile<<" PUM= "<<PUM<<endl;
  //        if(derivationOrder > 0) {
  //            zeroCounter = 0;
  //            logFile<<"first order derivations"<<endl;
  //            for(int k=0;k<firstDerivs.size();k++) {
  //                for(int j=0;j<firstDerivs[k].size();j++)
  //                    logFile<<suppPtcls[j]<<": "<<firstDerivs[k][j]<<" | ";
  //                logFile<<endl;
  //            }
  //            logFile<<"zeroCounter = "<<zeroCounter<<" of "<<suppPtcls.size()*3<<endl;
  //        }
  //    }
  ////#endif

  //    // Print out all the particles' distances
  //    logFile<<"==========================================================="<<endl;
  //    logFile<<"***************** Gausspoint Shape Functions **************"<<endl;
  //    logFile<<"==========================================================="<<endl;
  //    for (int g=0;g<localGaussPtsNum;g++){
  //        for (int g_=0;g_<particles[g].getShapeFuncs().size();g_++)
  //        logFile<< "Gauss["<<g<<"]:\t"
  //               <<gaussPoints[g].getShapeFuncs()[g_]<<endl;
  //    }

  //    // Print out all the particles' distances
  //    logFile<<"==========================================================="<<endl;
  //    logFile<<"**************** B-Gausspoint Shape Functions *************"<<endl;
  //    logFile<<"==========================================================="<<endl;
  //    for (int bg=0;bg<localBGaussPtsNum;bg++){
  //        for (int bg_=0;bg_<particles[bg].getShapeFuncs().size();bg_++)
  //        logFile<< "BGauss["<<bg<<"]:\t"
  //               <<boundGaussPoints[bg].getShapeFuncs()[bg_]<<endl;
  //    }

}

/************************************************************************/
/************************************************************************/
// Shape Function Test
void MaxEntDiscretising::shapeFuncTest(InputFileData* InputData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile,
                                       PetscViewer& viewerSEQ) {

  using namespace std;

  // Clear particles array
  particles.clear();

  // Define some constants
  double blk_spacing = 1,gamma = 1.8,tolLag = 1e-6,target_zero = 1.e-5;

  //     ====================================================================
  // Mesh: Calculating own coordinates

  // Define mesh points
  dbVector meshX,meshY;
  int meshSize = 10;

  double spacing = 10 / meshSize;
  blk_spacing = spacing;

  double increment = 0;

  for(;;) {
    meshX.push_back(increment);
    meshY.push_back(increment);
    increment += spacing;

    if(increment > 10) {
      break;
    }
  }

  // Populate particles vector
  double beta = gamma / (blk_spacing * blk_spacing);
  cout << "Value of Beta: " << beta << endl;
  Particle ptcle(3);
  int counter = 1;
  for(int y = 0;y < meshY.size();y++) {
    for(int x = 0;x < meshX.size();x++) {

      ptcle.setID(counter);
      ptcle.setCoords(meshX[x],meshY[y],0);
      cout << "Particles: " << meshX[x] << "," << meshY[y] << endl;
      ptcle.setBeta(beta);
      particles.push_back(ptcle);
      counter++;

    }
  }

  //     ====================================================================
  // Mesh: Read from files

  // Read coordinate File and populate particle list
  //    double coordX, coordY;
  //    int counter=0;
  //    Particle ptcle(3);

  //    double beta = gamma/(blk_spacing*blk_spacing);
  //    cout << "Value of Beta: "<< beta << endl;

  //    ifstream readFile("coords.list");

  //    if(!readFile) {
  //        cout <<"Can't open input file maxent-input.dat!"<< std::endl;

  //    }else if (readFile){
  //        readFile >> coordX >> coordY;

  //        while (!readFile.eof( )) {

  //            cout << coordX << "," << coordY << endl;
  //            ptcle.setID(counter);
  //            ptcle.setCoords(coordX,coordY,0);
  //            ptcle.setBeta(beta);
  //            particles.push_back(ptcle);

  //            counter++;

  //            readFile >> coordX >> coordY;
  //          }
  //    }

  //    readFile.close();
  //    ====================================================================

  // Find the number of supporting particles
  double nPtcls = particles.size();
  cout << "Num of Particles: " << particles.size() << endl;

  // Calculate the influence radius
  double radius = blk_spacing * sqrt( -1 / gamma * log(target_zero));

  if(radius < (2 * blk_spacing)) {
    radius = 2 * blk_spacing;
  }

  // Find the neighbours of a particle particles
  int mainPtcle = 38;
  double distX,distY,distZ,distance;
  intVector sPtclsVec;

  for(int r = 0;r < nPtcls;r++) {
    distX = particles[mainPtcle].getCoord(0) - particles[r].getCoord(0);
    distY = particles[mainPtcle].getCoord(1) - particles[r].getCoord(1);
    distZ = particles[mainPtcle].getCoord(2) - particles[r].getCoord(2);

    distance = sqrt(pow(distX,2) + pow(distY,2) + pow(distZ,2));

    if(distance < radius) {
      sPtclsVec.push_back(r);
    }
  }

  int sPtclsNum = sPtclsVec.size();
  cout << " Number of supporting particles: " << sPtclsNum << endl;

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(2);

  if(shapeFuncs.size() < sPtclsNum) shapeFuncs.resize(sPtclsNum);

  for(int j = 0;j < firstDerivShapes.size();j++)
    firstDerivShapes[j].resize(sPtclsNum);

  double minPtcleDist = 0;

  calcShapeFuncs(InputData,sPtclsVec,particles,particles[mainPtcle].getCoord(0),
                 particles[mainPtcle].getCoord(1),
                 particles[mainPtcle].getCoord(2),shapeFuncs,firstDerivShapes,
                 minPtcleDist,modelData,logFile,viewerSEQ);

  // Store the calculated shape function set.
  particles[mainPtcle].setShapeFuncs(sPtclsNum,shapeFuncs);
  //    particles[mainPtcle].setFirstDerivShapes(sPtclsNum,firstDerivShapes);

  // Output the shape function values and its corressponding derivatives
  double shapeFuncSum = 0;

  logFile << "****************************************************" << endl;
  logFile << "****************************************************" << endl;
  logFile << "********* Test Shape Function Values ***************" << endl;
  logFile << "****************************************************" << endl;
  logFile << "****************************************************" << endl;
  for(int s = 0;s < sPtclsNum;s++) {
    logFile << "ID[" << particles[sPtclsVec[s]].getID() << "]: "
        << shapeFuncs[s] << endl;
    shapeFuncSum += shapeFuncs[s];
  }
  logFile << "------------------------------------------------" << endl;
  logFile << "Sum of Shape Function Values : " << shapeFuncSum << endl;

  logFile << "****************************************************" << endl;
  logFile << "****************************************************" << endl;
  logFile << "******* Test Deriv Shape Function Values ***********" << endl;
  logFile << "****************************************************" << endl;
  logFile << "****************************************************" << endl;
  for(int d = 0;d < sPtclsNum;d++) {
    logFile << "ID[" << particles[sPtclsVec[d]].getID() << "]: "
        << firstDerivShapes[0][d] << ",\t" << firstDerivShapes[1][d] << endl;
  }

}

