#include "MLSPtcleIntegral.h"

MLSPtcleIntegral::MLSPtcleIntegral(InputFileData* InputData,
				   std::ofstream& logFile)  :
  MLSShapeFuncSet(InputData,logFile),
  ParticleDistribution(InputData,logFile) {}

/**********************************************************************/
/**********************************************************************/
// Copy the FEM background mesh, respectively nodal coordinates and
// Gauss integration scheme.
void MLSPtcleIntegral::copyFEData(FEMGeometry* FEData,
				  InputFileData* InputData,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile) {

  using namespace std;

  // Determine the number of particles 'particlesNum'.
  particlesNum =  FEData->getNodesNum();

  // Contains the connectivity new particle identifiers to the original 
  // read ones.
  oldIdx = FEData->getOldIdx();

  // Contains the connectivity original read particle identifiers to the 
  // new sorted ones.
  newIdx = FEData->getNewIdx();

  // Set the array 'particles' containing the particle properties.
  particles = FEData->getNodesVec();  

  // Set the indices of all particles any body loading is applied.
  bodyForcePtcleIdx = FEData->getBodyForcePtcleIdx();

  // Set the indices of all particles any boundary loading is applied
  tractionBoundPtcleIdx = FEData->getTractionBoundPtcleIdx();
  lineForceBoundPtcleIdx  = FEData->getLineForceBoundPtcleIdx();
  surfacePressureBoundPtcleIdx = FEData->getSurfacePressureBoundPtcleIdx(); 
  pointForceBoundPtcleIdx = FEData->getPointForceBoundPtcleIdx();

  // Set the indices of particles deformation boundary conditions are 
  // applied.
  pointDispBoundPtcleIdx = FEData->getPointDispBoundPtcleIdx();
  lineDispBoundPtcleIdx = FEData->getLineDispBoundPtcleIdx();
  surfaceDispBoundPtcleIdx = FEData->getSurfaceDispBoundPtcleIdx();

  // Delete of all FEM data, because there is no further need.
  delete FEData;


}

/**********************************************************************/
/**********************************************************************/
// Rearrange particle support lists according to the new particles
// vector ordering
void MLSPtcleIntegral::rearrangeSupportLists(InputFileData* InputData,
					     intVector& newGlobalPtcls,
					     std::map<std::string,double>& modelData,
					     std::ofstream& logFile) {

  using namespace std;

  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];

    intVector& vec = particles[ptcle].getSupportPtcls();
    intVector tmpVec = vec;

    for(int j=0;j<vec.size();j++)
      vec[j] = newGlobalPtcls[tmpVec[j]];

    sortIntVector(vec,0,vec.size()-1);


  }

#ifdef _geometryDebugMode_
  logFile<<"********** particle - particle support list **********"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    logFile<<"Ptcle "<<i<<": ";
    for(int j=0;j<suppPtcls.size();j++)
      logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all particles a vector containing the calculated
// shape functions and their derivations for all supported particles.
void MLSPtcleIntegral::setShapeFuncsOnPtcls(InputFileData* InputData,
					    std::map<std::string,double>& calcData,
					    std::map<std::string,double>& modelData,
					    std::ofstream& logFile,
					    PetscViewer& viewerMPI,
					    PetscViewer& viewerSEQ) {

  using namespace std;

  int derivationOrder = (int)modelData["shapesDerivationOrderOnIntPoints"];

  int rank,size; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"##################### particle shapes ####################"<<endl;
#endif

  dbVector shapeFuncs(maxPtcleSupport);
  dbMatrix firstDerivShapes(3,dbVector(maxPtcleSupport));
  dbMatrix secondDerivShapes(6,dbVector(maxPtcleSupport));
  int basisTermNum,supportSize;
  
  // Loop over processor's locally needed particles portion.
  for(int i=0;i<localGlobalPtcls.size();i++) {

    int& ptcle = localGlobalPtcls[i];

    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    supportSize = suppPtcls.size();
    
#ifdef _geometryDebugMode_
    logFile<<"PARTICLE "<<i<<endl;
#endif

    double& x = particles[ptcle].getCoord(0);
    double& y = particles[ptcle].getCoord(1);
    double& z = particles[ptcle].getCoord(2);

    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,supportSize,suppPtcls,particles,x,y,z,
		     shapeFuncs,basisTermNum,modelData,logFile,viewerSEQ);
    
      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);

    }
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,supportSize,suppPtcls,particles,x,y,z,
		     shapeFuncs,firstDerivShapes,modelData,logFile,
		     viewerSEQ);
    
      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);
      particles[ptcle].setFirstDerivShapes(supportSize,firstDerivShapes);


    }
    else if(derivationOrder == 2) {

      calcShapeFuncs(InputData,supportSize,suppPtcls,particles,x,y,z,
		     shapeFuncs,firstDerivShapes,secondDerivShapes,
		     modelData,logFile,viewerSEQ);
    
      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);
      particles[ptcle].setFirstDerivShapes(supportSize,firstDerivShapes);
      particles[ptcle].setSecondDerivShapes(supportSize,secondDerivShapes);

    }
    else {
      logFile <<"In MLSPtcleIntegral::setShapeFuncsOnPtcls shape function "
	      <<"derivations higher than second order are not supported!"<<endl;
      cerr <<"In MLSPtcleIntegral::setShapeFuncsOnPtcls shape function "
	   <<"derivations higher than second order are not supported!"<<endl;

      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

//   //MPI_Barrier(MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  double PUM;
  logFile<<"######### calculated shape sets on ptcls ##########"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    dbVector& shapes = particles[i].getShapeFuncs();
    dbMatrix& firstDerivShapes = particles[i].getFirstDerivShapes();
    dbMatrix& secondDerivShapes = particles[i].getSecondDerivShapes();
    PUM = 0;
    logFile<<"PARTICLE "<<i<<": "<<endl;
    for(int j=0;j<suppPtcls.size();j++) {
      PUM += shapes[j];
      logFile<<suppPtcls[j]<<": "<<shapes[j]<<" ";
    }
    if(fabs(1.0-PUM) > 0.000000001)
      logFile<<" PUM= "<<PUM<<endl;
    else
      logFile<<endl;
    logFile<<"**** first derivations :"<<endl;
    for(int j=0;j<firstDerivShapes.size();j++)
      for(int k=0;k<suppPtcls.size();k++)
	logFile<<suppPtcls[j]<<": "<<firstDerivShapes[j][k]<<" ";
    logFile<<endl;
    logFile<<"**** second derivations :"<<endl;
    for(int j=0;j<secondDerivShapes.size();j++)
      for(int k=0;k<suppPtcls.size();k++)
	logFile<<suppPtcls[j]<<": "<<secondDerivShapes[j][k]<<" ";
    logFile<<endl;
  }

#endif

}

/************************************************************************/
/************************************************************************/
// Return all particle indices where a surface displacement boundary 
// condition is applied (point,line,surface).
intMatrix& MLSPtcleIntegral::getAllDisplacementBoundPtcleIdx(InputFileData* InputData,
							     std::map<std::string,double>& modelData,
							     std::ofstream& logFile) {


  using namespace std;

  if(allDisplacementBoundPtcleIdx.size() > 0) 
    return allDisplacementBoundPtcleIdx;

  else {

    if((bool)modelData["displacementConstraint"] || 
       modelData["rotationConstraint"]) {

      pushBackVector(allDisplacementBoundPtcleIdx,
		     surfaceDispBoundPtcleIdx);
      pushBackVector(allDisplacementBoundPtcleIdx,
		     lineDispBoundPtcleIdx);
      pushBackVector(allDisplacementBoundPtcleIdx,
		     pointDispBoundPtcleIdx);
    }

    return allDisplacementBoundPtcleIdx;
  }

}

/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where a surface force boundary 
// condition is applied (point,line,surface).
intMatrix& MLSPtcleIntegral::getAllForceBoundPtcleIdx(InputFileData* InputData,
						      std::map<std::string,double>& modelData,
						      std::ofstream& logFile) {


  using namespace std;

  if(allForceBoundPtcleIdx.size() > 0) 
    return allForceBoundPtcleIdx;

  else {

    if((bool)modelData["tractionLoad"])
      pushBackVector(allForceBoundPtcleIdx,
		     tractionBoundPtcleIdx);

    if((bool)modelData["surfacePressureLoad"])
      pushBackVector(allForceBoundPtcleIdx,
		     surfacePressureBoundPtcleIdx);
    
    if((bool)modelData["lineForceLoad"])
      pushBackVector(allForceBoundPtcleIdx,
		     lineForceBoundPtcleIdx);

    if((bool)modelData["pointForceLoad"])
      pushBackVector(allForceBoundPtcleIdx,
		     pointForceBoundPtcleIdx);

    return allForceBoundPtcleIdx;
  }

}




