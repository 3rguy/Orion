#include "FEMDiscretising.h"

FEMDiscretising::FEMDiscretising(InputFileData* InputData,
				 std::ofstream& logFile) :
  BackgroundMesh(InputData,logFile),
  ParticleDistribution(InputData,logFile) {


  using namespace std;


  if((int)InputData->getValue("shapefunctionType") == 7) {
    InputData->setValue("radiusDeterminationAlgorithm",1);
    InputData->setValue("supportComputationMode",1);
    
    InputData->setValue("model",4);              
    InputData->setValue("constitutiveLaw",1);                
    InputData->setValue("boundaryEnforcementMethod",2);              
    InputData->setValue("penaltyParameter",1.00000e+15);     
    InputData->setValue("nitschePenaltyOnly",1);    
  }

}



/**********************************************************************/
/**********************************************************************/
// Copy the FEM background mesh, respectively nodal coordinates and
// Gauss integration scheme.
void FEMDiscretising::copyFEData(FEMGeometry* FEData,
				 InputFileData* InputData,
				 std::map<std::string,double>& modelData,
				 std::ofstream& logFile) {

  using namespace std;

  nodesElements = FEData->getNodesElemsVec();
  surfaceNodesElems = FEData->getSurfaceNodesElemsVec();
  lineNodesElems = FEData->getLineNodesElemsVec();

  allLineBoundGaussPtsIdx = FEData->getAllLineGaussPtsIdx(logFile);
  allSurfaceBoundGaussPtsIdx = FEData->getAllSurfaceGaussPtsIdx(logFile);

  vector<FEMElement>& dummyElemVec = FEData->getNodesElemsVec();
  dummyElemVec.resize(0,FEMElement(0));
  dummyElemVec = FEData->getSurfaceNodesElemsVec();
  dummyElemVec.resize(0,FEMElement(0));
  dummyElemVec = FEData->getLineNodesElemsVec();
  dummyElemVec.resize(0,FEMElement(0));

  intMatrix& dummyIntMat = FEData->getAllLineGaussPtsIdx(logFile);
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getAllSurfaceGaussPtsIdx(logFile);
  resizeArray(dummyIntMat,0);

  BackgroundMesh::copyFEData(FEData,InputData,modelData,logFile);

}

/**********************************************************************/
/**********************************************************************/
// Calculation of the shape functions and their derivatives with
// respect to global coordinates.
void FEMDiscretising::calcInterpolants(InputFileData* InputData,
				       std::map<std::string,double>& calcData,
				       std::map<std::string,double>& modelData,
				       std::ofstream& logFile,
				       PetscViewer& viewerMPI,
				       PetscViewer& viewerSEQ) {

  using namespace std;


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // Determine for all nodes their connected neighbour nodes.
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

  // Calculate for all Gauss points and all nodes a vector containing
  // the calculated shape functions and their derivations for all
  // supporting nodes.
  setAllShapeFuncs(InputData,modelData,logFile,viewerMPI,viewerSEQ);

  // Deallocate all class arrays of particles which are not locally
  // needed.
  deleteNonlocalPtcls(InputData,logFile);

}

/**********************************************************************/
/**********************************************************************/
// Determine for all nodes their connected neighbour nodes.
void FEMDiscretising::setPtclePtclsConn(InputFileData* InputData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile) {
  
  using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //   int startIdx = exclusiveLocalPtcls.front();
  //   int endIdx = exclusiveLocalPtcls.back()+1;

  bool supported;
  int supportSize,minPtcleSupport;
  int localMinSupport;

  int localMaxSupport = 0;

#ifdef _geometryDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"********** calculation of ptcle-ptcle connlist **********"<<endl;
#endif

  // Loop over a local portion of particles and create for each a support list.
  for(int i=0;i<particles.size();i++) {
    Particle& ptcle = particles[i];

    intVector& suppPtcls = ptcle.getSupportPtcls();
    suppPtcls.push_back(i);
  }
  
  maxPtcleSupport = 1;
  minPtcleSupport = 1;
  
  
  logFile<<"minPtcleSupport = "<<minPtcleSupport<<endl;
  logFile<<"maxPtcleSupport = "<<maxPtcleSupport<<endl;
  
  if(rank == 0)
    cout<<"ptcleSupport: "<<minPtcleSupport<<" - "<<maxPtcleSupport<<endl;
  
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
// Determine for each Gauss point its supporting nodes.
void FEMDiscretising::setGaussPtcleConn(InputFileData* InputData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile) {

  using namespace std;

  int gaussPtcleConnect =
    (int)InputData->getValue("gaussParticleConnectivity");
  int mode = (int)InputData->getValue("supportComputationMode");

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef _geometryDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"********** calculation of gauss-ptcle connlist **********"<<endl;
  double oldTime = MPI_Wtime();
  for (int i=0;i<localGaussPtsNum;i++) {
    GaussPoint& gPoint = gaussPoints[i];
    intVector& elemInfo = gPoint.getElementInfo();
    intVector& nodes = nodesElements[elemInfo[0]].getNodes();
    logFile<<"volume-element "<<elemInfo[0]<<": ";
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<endl;
  }
#endif

  int localMinSupport;
  int localMaxSupport = maxPtcleSupport;

  // Loop over all local gauss points.
  for (int i=0;i<localGaussPtsNum;i++) {
    GaussPoint& gPoint = gaussPoints[i];
    
    intVector& suppPtcls = gPoint.getSupportPtcls();
    intVector& elemInfo = gPoint.getElementInfo();
    intVector& nodes = nodesElements[elemInfo[0]].getNodes();

    resizeArray(suppPtcls,nodes.size());

    for(int j=0;j<nodes.size();j++)
      suppPtcls[j] = newIdx[nodes[j]-1];

    if(suppPtcls.size() > localMaxSupport) {
      localMaxSupport = suppPtcls.size();
      maxGaussSupport = suppPtcls.size();
    }

  }

  /**********************************************************************/
  // Determine global maximum number of supporting particles for a
  // gauss point.
  MPI_Allreduce(&localMaxSupport,&maxGaussSupport,1,
		MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  // Determine global minimum number of supporting particles for a
  // gauss point.
  MPI_Allreduce(&localMinSupport,&minGaussSupport,1,
		MPI_INT,MPI_MIN,MPI_COMM_WORLD);


  logFile<<"minGaussSupport = "<<minGaussSupport<<endl;
  logFile<<"maxGaussSupport = "<<maxGaussSupport<<endl;

  if(rank == 0)
    cout<<"Gauss-point support: "<<minGaussSupport<<" - "<<maxGaussSupport<<endl;


#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"******** supporting particles of gauss points ********"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    logFile<<i<<".) GAUSS POINT "<<gaussPoints[i].getGlobalID()
	   <<" (supp "<<gaussPoints[i].getSupportCounts()<<"): ";
    for(int j=0;j<gaussPoints[i].getSupportCounts();j++) {
      logFile<<gaussPoints[i].getSupportPtcle(j)<<" ";
    }
    logFile<<endl;
  }
  if(rank == 0)
    cout<<"gauss finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Abort(MPI_COMM_WORLD,1);
}

/**********************************************************************/
/**********************************************************************/
// Determine for each boundary gauss point its supporting nodes.
void FEMDiscretising::setBGaussPtcleConn(InputFileData* InputData,
					 std::map<std::string,double>& modelData,
					 std::ofstream& logFile) {

  using namespace std;

  int gaussPtcleConnect =
    (int)InputData->getValue("gaussParticleConnectivity");

  intMatrix& surfaceBGaussPtsIdx = allSurfaceBoundGaussPtsIdx;
  intMatrix& lineBGaussPtsIdx = allLineBoundGaussPtsIdx;

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef _geometryDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"******* calculation of bound-gauss-ptcle connlist *******"<<endl;
  double oldTime = MPI_Wtime();
  for (int i=0;i<allSurfaceBoundGaussPtsIdx.size();i++) {
    int& idx = allSurfaceBoundGaussPtsIdx[i][0];
    logFile<<"Surface-GPt "<<idx<<": "<<endl;
    GaussPoint& gPoint = boundGaussPoints[idx];
    intVector& elemInfo = gPoint.getElementInfo();
    intVector& nodes = surfaceNodesElems[elemInfo[0]].getNodes();
    logFile<<"surface-element "<<elemInfo[0]<<": ";
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<endl;
  }
#endif

  int localMinSupport;
  int localMaxSupport = maxGaussSupport;

  // Loop over all local gauss points.
  //for (int i=0;i<localBGaussPtsNum;i++) {
  for (int i=0;i<allSurfaceBoundGaussPtsIdx.size();i++) {
    int& idx = allSurfaceBoundGaussPtsIdx[i][0];
    GaussPoint& gPoint = boundGaussPoints[idx];
    
    intVector& suppPtcls = gPoint.getSupportPtcls();
    intVector& elemInfo = gPoint.getElementInfo();
    intVector& nodes = surfaceNodesElems[elemInfo[0]].getNodes();

    resizeArray(suppPtcls,nodes.size());

    for(int j=0;j<nodes.size();j++)
      suppPtcls[j] = newIdx[nodes[j]-1];

    if(suppPtcls.size() > localMaxSupport) {
      localMaxSupport = suppPtcls.size();
      maxGaussSupport = suppPtcls.size();
    }

  }

#ifdef _geometryDebugMode_
  logFile<<"**** calculation of line bound-gauss-ptcle connlist *****"<<endl;
  for (int i=0;i<allLineBoundGaussPtsIdx.size();i++) {
    int& idx = allLineBoundGaussPtsIdx[i][0];
    logFile<<"Line-GPt "<<idx<<": "<<endl;
    GaussPoint& gPoint = boundGaussPoints[idx];
    intVector& elemInfo = gPoint.getElementInfo();
    intVector& nodes = lineNodesElems[elemInfo[0]].getNodes();
    logFile<<"line-element "<<elemInfo[0]<<": ";
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<endl;
  }
#endif

  // Loop over all local gauss points.
  //for (int i=0;i<localBGaussPtsNum;i++) {
  for (int i=0;i<allLineBoundGaussPtsIdx.size();i++) {
    int& idx = allLineBoundGaussPtsIdx[i][0];
    GaussPoint& gPoint = boundGaussPoints[idx];
    
    intVector& suppPtcls = gPoint.getSupportPtcls();
    intVector& elemInfo = gPoint.getElementInfo();
    intVector& nodes = lineNodesElems[elemInfo[0]].getNodes();

    resizeArray(suppPtcls,nodes.size());

    for(int j=0;j<nodes.size();j++)
      suppPtcls[j] = newIdx[nodes[j]-1];

    if(suppPtcls.size() > localMaxSupport) {
      localMaxSupport = suppPtcls.size();
      maxGaussSupport = suppPtcls.size();
    }

  }

  // Determine global maximum number of supporting particles for a
  // gauss point.
  MPI_Allreduce(&localMaxSupport,&maxBGaussSupport,1,
		MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  // Determine global minimum number of supporting particles for a
  // gauss point.
  MPI_Allreduce(&localMinSupport,&minBGaussSupport,1,
		MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  // Check if the number of support boundary Gauss points is decisive
  if(maxBGaussSupport > maxGaussSupport)
    maxGaussSupport = maxBGaussSupport;

  if(minBGaussSupport < minGaussSupport)
    minGaussSupport = minBGaussSupport;


#ifdef _geometryDebugMode_
  logFile<<"*********************************************************"<<endl;
  logFile<<"********* supporting particles of boundary gauss *********"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    logFile<<"Bound GAUSSPOINT "<<i<<": ";
    for(int j=0;j<boundGaussPoints[i].getSupportCounts();j++) {
      logFile<<boundGaussPoints[i].getSupportPtcle(j)<<" ";
    }
    logFile<<endl;
  }
  if(rank == 0)
    cout<<"bound gauss finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
#endif

}

/************************************************************************/
/************************************************************************/
// Determine the Gauss point distribution among the processors.
void FEMDiscretising::setGaussDistribution(InputFileData* InputData,
					   std::map<std::string,double>& modelData,
					   std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"*********** set Gauss point distribution *************"<<endl;
#endif

  intVector globalLocalIdx(globalBGaussPtsNum);

  // Partitioning of volume and boundary Gauss points using software
  // package ParMETIS.
  if((int)InputData->getValue("partitionGaussPoints") == 1 && size > 1) {

    // Modify temporarily all volume and boundary Gauss point indices
    // to get consecutive ordering locally (necessary for ParMetis!).
    int defaultPortion,startBGIdx,endBGIdx;
    intVector gaussPortions(size);
    intVector bGaussPortions(size);

    // Determine the initial boundary Gauss point distribution.
    getInitBGaussSplitting(InputData,defaultPortion,startBGIdx,endBGIdx,
			   modelData,logFile);

    localBGaussPtsNum = endBGIdx - startBGIdx;

    MPI_Allgather(&localBGaussPtsNum,1,MPI_INT,&bGaussPortions[0],1,
		  MPI_INT,MPI_COMM_WORLD);

    MPI_Allgather(&localGaussPtsNum,1,MPI_INT,&gaussPortions[0],1,
		  MPI_INT,MPI_COMM_WORLD);

    // Loop over all processors to set the global volume Gauss point start
    // and end index locally, but including the boundary Gauss points.
    int startGIdx = 0;

    int m=0;

    while(m < rank) {
      startGIdx += gaussPortions[m]+bGaussPortions[m];
      m++;
    }

    int endGIdx = startGIdx + localGaussPtsNum;

    // Loop over all local volume Gauss points to set their indices.
    for(int i=0,j=startGIdx;i<localGaussPtsNum;i++,j++)
      gaussPoints[i].setGlobalID(j);

    intVector localBGIdx(localBGaussPtsNum);

    // Loop over all local boundary Gauss points to assemble the local
    // indices.
    for(int i=0,j=endGIdx;i<localBGaussPtsNum;i++,j++)
      localBGIdx[i] = j;

    int* recvCounts = new int[size];
    int* displs = new int[size];

    int sendBuf = localBGaussPtsNum;
    MPI_Allgather(&sendBuf,1,MPI_INT,recvCounts,1,MPI_INT,MPI_COMM_WORLD);

    displs[0] = 0;

    for(int i=1;i<size;i++)
      displs[i] = displs[i-1]+recvCounts[i-1];

    intVector globalBGIdx(globalBGaussPtsNum);

    MPI_Allgatherv(&localBGIdx[0],recvCounts[rank],MPI_INT,
		   &globalBGIdx[0],recvCounts,displs,MPI_INT,
		   MPI_COMM_WORLD);

    // Loop over all global boundary Gauss points to set their indices.
    for(int i=0;i<globalBGaussPtsNum;i++)
      boundGaussPoints[i].setGlobalID(globalBGIdx[i]);

#ifdef _geometryDebugMode_
    logFile<<"*********** newly set Gauss indics *****************"<<endl;
    for(int i=0;i<localGaussPtsNum;i++)
      logFile<<i<<".) "<<gaussPoints[i].getGlobalID()<<endl;
    logFile<<"****************************************************"<<endl;
    for(int i=startBGIdx;i<endBGIdx;i++)
      logFile<<i<<".) "<<boundGaussPoints[i].getGlobalID()<<endl;
    logFile<<"****************************************************"<<endl;
    for(int i=0;i<globalBGaussPtsNum;i++)
      logFile<<i<<".) "<<boundGaussPoints[i].getGlobalID()<<endl;
#endif

    // ------------------------------------------------------------------
    // Set a vector contains the coordinates of all volume and boundary
    // Gauss points.
    intVector idx(localGaussPtsNum);
    dbVector localCoords(localGaussPtsNum*usedDims);

    // Loop over all local volume Gauss points.
    for(int i=0;i<localGaussPtsNum;i++) {
      idx[i] = gaussPoints[i].getGlobalID();

      for(int j=0;j<usedDims;j++)
	localCoords[i*usedDims+j] = gaussPoints[i].getCoord(j);

    }

#ifdef _geometryDebugMode_
    logFile<<"*********** local vol gauss coordinates ************"<<endl;
    for(int i=0;i<localGaussPtsNum;i++) {
      logFile<<i<<".) "<<idx[i]<<": ";
      for(int j=0;j<usedDims;j++)
	logFile<<localCoords[i*usedDims+j]<<" ";
      logFile<<endl;
    }
#endif

    // The index part.
    intVector allIdx(globalGaussPtsNum);

    MPI_Allgather(&localGaussPtsNum,1,MPI_INT,recvCounts,1,MPI_INT,
		  MPI_COMM_WORLD);

    displs[0] = 0;

    for(int i=1;i<size;i++)
      displs[i] = displs[i-1]+recvCounts[i-1];

    MPI_Allgatherv(&idx[0],recvCounts[rank],MPI_INT,&allIdx[0],recvCounts,
		   displs,MPI_INT,MPI_COMM_WORLD);

    // The coordinates' part.
    dbVector tmpGlobalCoords(globalGaussPtsNum*usedDims);

    for(int i=0;i<size;i++)
      recvCounts[i] *= usedDims;

    displs[0] = 0;

    for(int i=1;i<size;i++)
      displs[i] = displs[i-1]+recvCounts[i-1];

    MPI_Allgatherv(&localCoords[0],recvCounts[rank],MPI_DOUBLE,
		   &tmpGlobalCoords[0],recvCounts,displs,MPI_DOUBLE,
		   MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
    logFile<<"*********** global vol gauss coordinates ***********"<<endl;
    for(int i=0;i<globalGaussPtsNum;i++) {
      logFile<<i<<".) "<<allIdx[i]<<": ";
      for(int j=0;j<usedDims;j++)
	logFile<<tmpGlobalCoords[i*usedDims+j]<<" ";
      logFile<<endl;
    }
#endif

    // Order the coordinates according to their global indices.
    dbVector globalCoords((globalGaussPtsNum+globalBGaussPtsNum)*usedDims);

    // Loop over all global volume Gauss points.
    for(int i=0;i<globalGaussPtsNum;i++)

      for(int j=0;j<usedDims;j++)
	globalCoords[allIdx[i]*usedDims+j] = tmpGlobalCoords[i*usedDims+j];

    // Loop over all global boundary Gauss points.
    for(int i=0;i<globalBGaussPtsNum;i++)

      for(int j=0;j<usedDims;j++)
	globalCoords[boundGaussPoints[i].getGlobalID()*usedDims+j] =
	  boundGaussPoints[i].getCoord(j);


#ifdef _geometryDebugMode_
    logFile<<"************* global gauss coordinates ***************"<<endl;
    for(int i=0;i<globalGaussPtsNum+globalBGaussPtsNum;i++) {
      logFile<<"GAUSSPOINT "<<i<<".) ";
      for(int j=0;j<usedDims;j++)
	logFile<<globalCoords[i*usedDims+j]<<" ";
      logFile<<endl;
    }
#endif

    delete[] displs,recvCounts;

    // ------------------------------------------------------------------
    // Determine for each gauss point its supporting particles.
    localBGaussPtsNum = globalBGaussPtsNum;

    setGaussPtcleConn(InputData,modelData,logFile);
    setBGaussPtcleConn(InputData,modelData,logFile);

    localBGaussPtsNum = endBGIdx - startBGIdx;

    // ------------------------------------------------------------------
    intVector gaussProcList;

    // Create a gauss points connection list needed for partitioning
    // and initiate the partitioning.
    setPartitionGraph(InputData,gaussProcList,globalCoords,modelData,
		      logFile);

    // ------------------------------------------------------------------
    // Exchange the volume Gauss data between the processor that each of
    // them possesses the correct local data.
    exchangeGaussData(InputData,gaussProcList,modelData,logFile);

    // Split the boundary Gauss points between the processor that each of
    // them possesses the correct local data.
    splitBoundGaussPts(InputData,gaussProcList,globalLocalIdx,modelData,
		       logFile);

  }

  /**********************************************************************/
  // Straight forward partitioning of volume and boundary Gauss points.
  else {

    // Determine the initial boundary Gauss point distribution.
    int defaultBGaussPortion,startIdx,endIdx;
    getInitBGaussSplitting(InputData,defaultBGaussPortion,startIdx,endIdx,
			   modelData,logFile);

    localBGaussPtsNum = endIdx - startIdx;

#ifdef _geometryDebugMode_
    logFile<<"******** boundary Gauss points distribution ********"<<endl;
    logFile<<"defaultBGaussPortion = "<<defaultBGaussPortion<<endl;
    logFile<<"startIdx "<<startIdx<<" endIdx "<<endIdx<<endl;
    logFile<<"localBGaussPtsNum "<<localBGaussPtsNum
	   <<" globalBGaussPtsNum "<<globalBGaussPtsNum<<endl;
#endif


    int increment;
    vector<GaussPoint>::iterator startPos;
    vector<GaussPoint>::iterator endPos;

    // erase the leading part of the boundary Gauss point's vector.
    if(rank > 0) {

      startPos =
	boundGaussPoints.begin();
      endPos = boundGaussPoints.begin() + startIdx;

      boundGaussPoints.erase(startPos,endPos);
    }

    // erase the end part of the boundary Gauss point's vector.
    if(rank+1 < size && boundGaussPoints.size() > defaultBGaussPortion) {

      startPos = boundGaussPoints.begin() + defaultBGaussPortion;
      endPos = boundGaussPoints.end();

      boundGaussPoints.erase(startPos,endPos);
    }

    // erase the nonlocal boundary Gauss points from the corresponding
    // surface and line elements

    int count,idx;
    intVector newIntPts;

    // loop over all surface elements
    for(int i=0;i<surfaceNodesElems.size();i++) {
      
      FEMElement& elem = surfaceNodesElems[i];
      intVector& intPts = elem.getSurfaceIntegrationPts();
      
      resizeArray(newIntPts,intPts.size());
      count = 0;
      
      for(int j=0;j<intPts.size();j++) {
	idx = intPts[j];
	
	if(idx >= startIdx && idx < endIdx) {
	  newIntPts[count] = idx;
	  count++;
	}
	
      }

      resizeArray(newIntPts,count);
      intPts = newIntPts;
    }

    // loop over all line elements
    for(int i=0;i<lineNodesElems.size();i++) {

      FEMElement& elem = lineNodesElems[i];
      intVector& intPts = elem.getLineIntegrationPts();

      resizeArray(newIntPts,intPts.size());
      count = 0;

      for(int j=0;j<intPts.size();j++) {
	idx = intPts[j];

	if(idx >= startIdx && idx < endIdx) {
	  newIntPts[count] = idx;
	  count++;
	}

      }

      resizeArray(newIntPts,count);
      intPts = newIntPts;
    }

#ifdef _geometryDebugMode_
    logFile<<"******* local stored boundary Gauss points *********"<<endl;
    logFile<<"global startIdx = "<<startIdx<<" endIdx = "<<endIdx<<endl;
    logFile<<"localBGaussPtsNum = "<<localBGaussPtsNum<<" ?=  "
	   <<boundGaussPoints.size()<<endl;
    for(int i=0;i<localBGaussPtsNum;i++) {
      logFile<<i<<".) "<<boundGaussPoints[i].getGlobalID()<<"; coords: ";
      for(int j=0;j<usedDims;j++)
	logFile<<boundGaussPoints[i].getCoord(j)<<" ";
      logFile<<"; weight = "<<boundGaussPoints[i].getWeight()<<"; ";
      dbVector& normal = boundGaussPoints[i].getSurfaceNormal();
      logFile<<"normal : ";
      for(int j=0;j<normal.size();j++)
	logFile<<normal[j]<<" ";
      logFile<<endl;
    }
#endif

    // ------------------------------------------------------------------
    // Store the local portion of boundary Gauss point indices which
    // is a Dirichlet or a v. Neumann boundary condition applied to.

    globalLocalIdx.assign(globalLocalIdx.size(),-1);
    int m=0;

    for(int i=startIdx;i<endIdx;i++) {
      globalLocalIdx[i] = m;
      m++;
    }

#ifdef _geometryDebugMode_
    m=0;
    logFile<<"********* global -> local connectivity list ********"<<endl;
    for(int i=0;i<globalLocalIdx.size();i++) {
      if(globalLocalIdx[i] != -1)
	logFile<<i<<".) global: "<<boundGaussPoints[m].getGlobalID()
	       <<" -> local: "<<globalLocalIdx[i]<<endl;
      else
	logFile<<i<<".) "<<globalLocalIdx[i]
	       <<" -> not available "<<endl;
      m++;
    }
#endif

  }


  if(localGaussPtsNum == 0) {
    logFile<<"In function BackgroundMesh::setGaussDistribution volume "
	   <<"Gauss point split\n"
	   <<"amongst processes gives zero exclusively local volume "
	   <<"Gauss points."<<endl
           <<"Choose less processes."<<endl;
    cerr<<"In function BackgroundMesh::setGaussDistribution volume "
	<<"Gauss point split\n"
	<<"amongst processes gives zero exclusively local volume\n"
	<<"Gauss points at process "<<rank<<"."<<endl
        <<"Choose less processes."<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /*********************************************************************/
  // modify all vectors which contain the indices of all boundary
  // Gauss points corresponding to the new local ordering

  // surface displacement boundary conditions
  intMatrix tmpGaussPtsIdx = surfaceDispBoundGaussPtsIdx;
  int m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfaceDispBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfaceDispBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfaceDispBoundGaussPtsIdx.resize(m);

  // line displacement boundary conditions
  tmpGaussPtsIdx = lineDispBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      lineDispBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	lineDispBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  lineDispBoundGaussPtsIdx.resize(m);

  // point displacement boundary conditions
  tmpGaussPtsIdx = pointDispBoundPtcleIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      pointDispBoundPtcleIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	pointDispBoundPtcleIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  pointDispBoundPtcleIdx.resize(m);

  // --------------------------------------------------------------------
  // surface rotation boundary conditions
  tmpGaussPtsIdx = surfaceRotBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfaceRotBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfaceRotBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfaceRotBoundGaussPtsIdx.resize(m);

  // line rotation boundary conditions
  tmpGaussPtsIdx = lineRotBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      lineRotBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	lineRotBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  lineRotBoundGaussPtsIdx.resize(m);

  // point rotation boundary conditions
  tmpGaussPtsIdx = pointRotBoundPtcleIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      pointRotBoundPtcleIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	pointRotBoundPtcleIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  pointRotBoundPtcleIdx.resize(m);

  // --------------------------------------------------------------------
  // surface electric boundary conditions
  tmpGaussPtsIdx = surfaceElectricBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfaceElectricBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfaceElectricBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfaceElectricBoundGaussPtsIdx.resize(m);

  // line electric boundary conditions
  tmpGaussPtsIdx = lineElectricBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      lineElectricBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	lineElectricBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  lineElectricBoundGaussPtsIdx.resize(m);


  // point electric boundary conditions
  tmpGaussPtsIdx = pointElectricBoundPtcleIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      pointElectricBoundPtcleIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	pointElectricBoundPtcleIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  pointElectricBoundPtcleIdx.resize(m);


  // --------------------------------------------------------------------
  // surface depolarisation boundary conditions
  tmpGaussPtsIdx = surfaceDepolarisationBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfaceDepolarisationBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfaceDepolarisationBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfaceDepolarisationBoundGaussPtsIdx.resize(m);

  // line depolarisation boundary conditions
  tmpGaussPtsIdx = lineDepolarisationBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      lineDepolarisationBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	lineDepolarisationBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  lineDepolarisationBoundGaussPtsIdx.resize(m);


  // point depolarisation boundary conditions
  tmpGaussPtsIdx = pointDepolarisationBoundPtcleIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      pointDepolarisationBoundPtcleIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	pointDepolarisationBoundPtcleIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  pointDepolarisationBoundPtcleIdx.resize(m);

  // --------------------------------------------------------------------
  // surface micro boundary conditions
  tmpGaussPtsIdx = surfaceMicroBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfaceMicroBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfaceMicroBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfaceMicroBoundGaussPtsIdx.resize(m);

  // line micro boundary conditions
  tmpGaussPtsIdx = lineMicroBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      lineMicroBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	lineMicroBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  lineMicroBoundGaussPtsIdx.resize(m);

  // --------------------------------------------------------------------
  // surface pressure boundary conditions
  tmpGaussPtsIdx = surfacePressureBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfacePressureBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfacePressureBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfacePressureBoundGaussPtsIdx.resize(m);

  // traction boundary conditions
  tmpGaussPtsIdx = tractionBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      tractionBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	tractionBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  tractionBoundGaussPtsIdx.resize(m);

  // line force boundary conditions
  tmpGaussPtsIdx = lineForceBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      lineForceBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	lineForceBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  lineForceBoundGaussPtsIdx.resize(m);

  // point force boundary conditions
  //tmpGaussPtsIdx = pointForceBoundPtcleIdx;
  //m=0;

  //for(int i=0;i<tmpGaussPtsIdx.size();i++)

  // if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
  // int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
  // pointForceBoundPtcleIdx[m][0] = ID;

  //for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
  //pointForceBoundPtcleIdx[m][j] = tmpGaussPtsIdx[i][j];
  
  //m++;
  //}
  
  // pointForceBoundPtcleIdx.resize(m);
  
  // --------------------------------------------------------------------
  // surface moment boundary conditions
  tmpGaussPtsIdx = surfaceMomentBoundGaussPtsIdx;
  m=0;
  
  for(int i=0;i<tmpGaussPtsIdx.size();i++)
    
    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfaceMomentBoundGaussPtsIdx[m][0] = ID;
      
      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfaceMomentBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfaceMomentBoundGaussPtsIdx.resize(m);

  // --------------------------------------------------------------------
  // surface electric charge boundary conditions
  tmpGaussPtsIdx = surfaceElectricChargeBoundGaussPtsIdx;
  m=0;

  for(int i=0;i<tmpGaussPtsIdx.size();i++)

    if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
      int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
      surfaceElectricChargeBoundGaussPtsIdx[m][0] = ID;

      for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
	surfaceElectricChargeBoundGaussPtsIdx[m][j] = tmpGaussPtsIdx[i][j];

      m++;
    }

  surfaceElectricChargeBoundGaussPtsIdx.resize(m);


  globalLocalIdx.resize(0);
  tmpGaussPtsIdx.resize(0);

  
#ifdef _geometryDebugMode_
  logFile<<"****** local surface deformation boundary Gauss points *******"<<endl;
  for(int i=0;i<surfaceDispBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfaceDispBoundGaussPtsIdx[i].size();j++)
      logFile<<surfaceDispBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local line deformation boundary Gauss points *******"<<endl;
  for(int i=0;i<lineDispBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<lineDispBoundGaussPtsIdx[i].size();j++)
      logFile<<lineDispBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local point deformation boundary Gauss points *******"<<endl;
  for(int i=0;i<pointDispBoundPtcleIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<pointDispBoundPtcleIdx[i].size();j++)
      logFile<<pointDispBoundPtcleIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local surface rotation boundary Gauss points *******"<<endl;
  for(int i=0;i<surfaceRotBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfaceRotBoundGaussPtsIdx[i].size();j++)
      logFile<<surfaceRotBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local line rotation boundary Gauss points *******"<<endl;
  for(int i=0;i<lineRotBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<lineRotBoundGaussPtsIdx[i].size();j++)
      logFile<<lineRotBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local point rotation boundary Gauss points *******"<<endl;
  for(int i=0;i<pointRotBoundPtcleIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<pointRotBoundPtcleIdx[i].size();j++)
      logFile<<pointRotBoundPtcleIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local surface electric boundary Gauss points *******"<<endl;
  for(int i=0;i<surfaceElectricBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfaceElectricBoundGaussPtsIdx[i].size();j++)
      logFile<<surfaceElectricBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local line electric boundary Gauss points *******"<<endl;
  for(int i=0;i<lineElectricBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<lineElectricBoundGaussPtsIdx[i].size();j++)
      logFile<<lineElectricBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local point electric boundary Gauss points *******"<<endl;
  for(int i=0;i<pointElectricBoundPtcleIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<pointElectricBoundPtcleIdx[i].size();j++)
      logFile<<pointElectricBoundPtcleIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local surface depolarisation boundary Gauss points *******"<<endl;
  for(int i=0;i<surfaceDepolarisationBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfaceDepolarisationBoundGaussPtsIdx[i].size();j++)
      logFile<<surfaceDepolarisationBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local line depolarisation boundary Gauss points *******"<<endl;
  for(int i=0;i<lineDepolarisationBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<lineDepolarisationBoundGaussPtsIdx[i].size();j++)
      logFile<<lineDepolarisationBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local point depolarisation boundary Gauss points *******"<<endl;
  for(int i=0;i<pointDepolarisationBoundPtcleIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<pointDepolarisationBoundPtcleIdx[i].size();j++)
      logFile<<pointDepolarisationBoundPtcleIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local surface micro boundary Gauss points *******"<<endl;
  for(int i=0;i<surfaceMicroBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfaceMicroBoundGaussPtsIdx[i].size();j++)
      logFile<<surfaceMicroBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local line micro boundary Gauss points *******"<<endl;
  for(int i=0;i<lineMicroBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<lineMicroBoundGaussPtsIdx[i].size();j++)
      logFile<<lineMicroBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local surface pressure boundary Gauss points *******"<<endl;
  for(int i=0;i<surfacePressureBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfacePressureBoundGaussPtsIdx[i].size();j++)
      logFile<<surfacePressureBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local traction boundary Gauss points *******"<<endl;
  for(int i=0;i<tractionBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<tractionBoundGaussPtsIdx[i].size();j++)
      logFile<<tractionBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local line force boundary Gauss points *******"<<endl;
  for(int i=0;i<lineForceBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<lineForceBoundGaussPtsIdx[i].size();j++)
      logFile<<lineForceBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"****** local pointforce boundary Gauss points *******"<<endl;
  for(int i=0;i<pointForceBoundPtcleIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<pointForceBoundPtcleIdx[i].size();j++) {
      logFile<<pointForceBoundPtcleIdx[i][j]<<" "<<endl;
    }
  }
  logFile<<"****** local surface moment boundary Gauss points *******"<<endl;
  for(int i=0;i<surfaceMomentBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfaceMomentBoundGaussPtsIdx[i].size();j++)
      logFile<<surfaceMomentBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"***** local surface electric charge boundary Gauss points ****"
	 <<endl;
  for(int i=0;i<surfaceElectricChargeBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) ";
    for(int j=0;j<surfaceElectricChargeBoundGaussPtsIdx[i].size();j++)
      logFile<<surfaceElectricChargeBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
#endif

  resizeArray(allSurfaceBoundGaussPtsIdx,0);
  allSurfaceBoundGaussPtsIdx =
    getAllSurfaceGaussPtsIdx(logFile);

  resizeArray(allLineBoundGaussPtsIdx,0);
  allLineBoundGaussPtsIdx =
    getAllLineGaussPtsIdx(logFile);

  /*********************************************************************/
  // Determine for each gauss point its supporting particles.
  
  setGaussPtcleConn(InputData,modelData,logFile);
  setBGaussPtcleConn(InputData,modelData,logFile);
    
#ifdef _geometryDebugMode_
  int lowestPtcle,highestPtcle;
  lowestPtcle = particlesNum - 1;
  highestPtcle = 0;
  for(int i=0;i<localGaussPtsNum;i++) {
    intVector& sPtcls = gaussPoints[i].getSupportPtcls();
    for(int j=0;j<sPtcls.size();j++) {
      if(sPtcls[j] < lowestPtcle) lowestPtcle = sPtcls[j];
      if(sPtcls[j] > highestPtcle) highestPtcle = sPtcls[j];
    }
  }
  logFile<<"************ supported particles range ****************"<<endl;
  logFile<<lowestPtcle<<" - "<<highestPtcle<<endl;
#endif


  /*********************************************************************/
  // Create a root list for all global volume Gauss points.
  intVector localRootList(globalGaussPtsNum);
  gaussRootList = intVector(globalGaussPtsNum);

  for(int i=0;i<localGaussPtsNum;i++)
    localRootList[gaussPoints[i].getGlobalID()] = rank;

  MPI_Allreduce(&localRootList[0],&gaussRootList[0],globalGaussPtsNum,
		MPI_INT,MPI_SUM,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"localGaussPtsNum = "<<localGaussPtsNum<<endl;
  logFile<<"localBGaussPtsNum = "<<localBGaussPtsNum<<endl;
  logFile<<"globalGaussPtsNum = "<<globalGaussPtsNum<<endl;
  logFile<<"globalBGaussPtsNum = "<<globalBGaussPtsNum<<endl;
  logFile<<"************ local vol Gauss root list ***************"<<endl;
  for(int i=0;i<globalGaussPtsNum;i++)
    logFile<<"GAUSSPOINT "<<i<<" "<<localRootList[i]<<endl;
  logFile<<"************* global vol Gauss root list *************"<<endl;
  for(int i=0;i<globalGaussPtsNum;i++) {
    logFile<<"GAUSSPOINT "<<i<<" "<<gaussRootList[i]<<endl;
  }
#endif

  // Create a root list for all global boundary Gauss points.
  if(globalBGaussPtsNum > globalGaussPtsNum)
    localRootList.resize(globalBGaussPtsNum);

  clearArray(localRootList);

  bGaussRootList = intVector(globalBGaussPtsNum);

  for(int i=0;i<localBGaussPtsNum;i++)
    localRootList[boundGaussPoints[i].getGlobalID()] =
      rank;

  MPI_Allreduce(&localRootList[0],&bGaussRootList[0],globalBGaussPtsNum,
		MPI_INT,MPI_SUM,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"************ local bound Gauss root list *************"<<endl;
  for(int i=0;i<globalBGaussPtsNum;i++)
    logFile<<"Bound GAUSSPOINT "<<i<<" "<<localRootList[i]<<endl;
  logFile<<"*********** global bound Gauss root list *************"<<endl;
  for(int i=0;i<globalBGaussPtsNum;i++) {
    logFile<<"Bound GAUSSPOINT "<<i<<" "<<bGaussRootList[i]<<endl;
  }
  logFile<<"######################################################"<<endl;
  logFile<<"**************** volume Gauss points *****************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    int& matID = gaussPoints[i].getMaterialID();
    dbVector& gCoords = gaussPoints[i].getCoords();
    logFile<<i<<".) Volume GAUSSPOINT "
	   <<gaussPoints[i].getGlobalID()<<": ";
    for(int j=0;j<gCoords.size();j++)
      logFile<<gCoords[j]<<" ";
    logFile<<"matID "<<matID<<endl;
  }
  for(int i=0;i<localBGaussPtsNum;i++) {
    int& matID =  boundGaussPoints[i].getMaterialID();
    dbVector& gCoords = boundGaussPoints[i].getCoords();
    logFile<<i<<".) Bound GAUSSPOINT: "
	   <<boundGaussPoints[i].getGlobalID()<<": ";
    for(int j=0;j<gCoords.size();j++)
      logFile<<gCoords[j]<<" ";
    logFile<<" matID "<<matID<<endl;
  }
#endif

  // Finite element stuff is not needed anymore -> free space
  resizeArray(globalLocalElemIdx,0);
  //  vector<FEMElement>(0,FEMElement(0)).swap(nodesElements);


  int maxGaussPtsNum,minGaussPtsNum,maxBGaussPtsNum,minBGaussPtsNum;

  MPI_Allreduce(&localGaussPtsNum,&maxGaussPtsNum,1,MPI_INT,MPI_MAX,
		MPI_COMM_WORLD);
  MPI_Allreduce(&localGaussPtsNum,&minGaussPtsNum,1,MPI_INT,MPI_MIN,
		MPI_COMM_WORLD);
  
  if(rank == 0)
    cout<<"volume Gauss point split: "<<minGaussPtsNum<<" - "
	<<maxGaussPtsNum<<endl;
  
}

/**********************************************************************/
/**********************************************************************/
// Calculate for all Gauss points and all nodes a vector containing 
// the calculated shape functions and their derivatives for all 
// supporting nodes.
void FEMDiscretising::setAllShapeFuncs(InputFileData* InputData,
				       std::map<std::string,double>& modelData,
				       std::ofstream& logFile,
				       PetscViewer& viewerMPI,
				       PetscViewer& viewerSEQ) {

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
// Calculate for all volume Gauss points a vector containing the 
// shape functions and their derivations for their respective supporting
// nodes.
void FEMDiscretising::setShapeFuncsOnGauss(InputFileData* InputData,
					   std::map<std::string,double>& modelData,
					   std::ofstream& logFile,
					   PetscViewer& viewerMPI,
					   PetscViewer& viewerSEQ) {

  using namespace std;


  // Loop over all local gauss points to set their shape function
  // vectors.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"******************* gausspoint shapes ********************"<<endl;
#endif

  int supportSize,idx;

  // Loop over all local elements.
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];
    intVector& intPts = elem.getVolumeIntegrationPts();

    // Get the template element shape function set.
    dbMatrix& sFuncs = elem.getVolumeShapeFuncOrds();

    // Loop over all Gauss points of the current element to determine
    // the shapefunction and shapefunction derivatives ordinates of 
    // the element nodes at the Gauss point.
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      GaussPoint& gPoint = gaussPoints[idx];
      intVector& suppPtcls = gPoint.getSupportPtcls();
      dbVector& shapeFuncs = gPoint.getShapeFuncs();
      dbMatrix& firstDerivShapes = gPoint.getFirstDerivShapes();

#ifdef _geometryDebugMode_
      logFile<<"Gauss point "<<gPoint.getGlobalID()<<endl;
#endif

      supportSize = suppPtcls.size();

      // set the shapefunctions
      shapeFuncs = sFuncs[j];

      // set the shapefunction derivatives
      firstDerivShapes = gPoint.getFEShapeFuncDerivs();      
    }

  }

#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"******** calculated shape sets on gauss pts **********"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    logFile<<"GAUSSPOINT "<<gaussPoints[i].getGlobalID()<<":: ";
    for(int j=0;j<gaussPoints[i].getSupportCounts();j++) {
      logFile<<gaussPoints[i].getSupportPtcle(j)<<": "
	     <<gaussPoints[i].getShapeFunc(j)<<" | ";
    }
    dbMatrix& firstDerivs = gaussPoints[i].getFirstDerivShapes();
    logFile<<"first order derivations"<<endl;
    for(int k=0;k<3;k++) {
      for(int j=0;j<firstDerivs[k].size();j++)
	logFile<<gaussPoints[i].getSupportPtcle(j)<<": "
	       <<firstDerivs[k][j]<<" |";
      logFile<<endl;
    }
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all boundary gauss points a vector containing the 
// calculated shape functions for all supported particles.
void FEMDiscretising::setShapeFuncsOnBGauss(InputFileData* InputData,
					    std::map<std::string,double>& modelData,
					    std::ofstream& logFile,
					    PetscViewer& viewerMPI,
					    PetscViewer& viewerSEQ) {

  using namespace std;

  vector<FEMElement>& lElems = lineNodesElems;
  vector<FEMElement>& sElems = surfaceNodesElems;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"*************** boundary gausspoint shapes ***************"<<endl;
#endif

  int supportSize,idx;

  // Loop over all surface elements.
  for(int i=0;i<sElems.size();i++) {
    FEMElement& elem = sElems[i];
    intVector& intPts = elem.getSurfaceIntegrationPts();

    // Get the template element shape function set.
    dbMatrix& sFuncs = elem.getSurfaceShapeFuncOrds();

    // Loop over all Gauss points of the current element to determine
    // the shapefunction and shapefunction derivatives ordinates of 
    // the element nodes at the Gauss point.
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      GaussPoint& gPoint = boundGaussPoints[idx];
      intVector& suppPtcls = gPoint.getSupportPtcls();
      dbVector& shapeFuncs = gPoint.getShapeFuncs();
      //dbMatrix& firstDerivShapes = gPoint.getFirstDerivShapes();

#ifdef _geometryDebugMode_
      logFile<<"SGauss point "<<gPoint.getGlobalID()<<endl;
#endif

      supportSize = suppPtcls.size();

      // set the shapefunctions
      shapeFuncs = sFuncs[j];

      // set the shapefunction derivatives - not continuous at element faces

    }

  }

  /**********************************************************************/
  // Loop over all line elements.
  for(int i=0;i<lElems.size();i++) {
    FEMElement& elem = lElems[i];
    intVector& intPts = elem.getLineIntegrationPts();

    // Get the template element shape function set.
    dbMatrix& sFuncs = elem.getLineShapeFuncOrds();

    // Loop over all Gauss points of the current element to determine
    // the shapefunction and shapefunction derivatives ordinates of 
    // the element nodes at the Gauss point.
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      GaussPoint& gPoint = boundGaussPoints[idx];
      intVector& suppPtcls = gPoint.getSupportPtcls();
      dbVector& shapeFuncs = gPoint.getShapeFuncs();
      //dbMatrix& firstDerivShapes = gPoint.getFirstDerivShapes();

#ifdef _geometryDebugMode_
      logFile<<"LGauss point "<<gPoint.getGlobalID()<<endl;
#endif

      supportSize = suppPtcls.size();

      // set the shapefunctions
      shapeFuncs = sFuncs[j];

      // set the shapefunction derivatives - not continuous at element faces!


    }

  }

#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"**** calculated shape sets on boundary Gauss pts *****"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    GaussPoint& gPoint = boundGaussPoints[i];
    intVector& suppPtcls = gPoint.getSupportPtcls();
    dbVector& shapeFuncs = gPoint.getShapeFuncs();
    logFile<<"LGAUSSPOINT "<<gPoint.getGlobalID()
	   <<"(supp. "<<gPoint.getSupportCounts()<<"/"
	   <<suppPtcls.size()<<"): ";
    for(int j=0;j<shapeFuncs.size();j++)
      logFile<<shapeFuncs[j]<<" ";
    logFile<<endl;
    dbMatrix& firstDerivs = gPoint.getFirstDerivShapes();
    logFile<<"first order derivations :"<<endl;
    for(int k=0;k<3;k++) {
      for(int j=0;j<firstDerivs[k].size();j++)
	logFile<<firstDerivs[k][j]<<" ";
      logFile<<endl;
    }
  }
#endif

}
 
/**********************************************************************/
/**********************************************************************/
// Calculate for all nodes a vector containing the calculated
// shape functions and their derivations for all supported nodes.
void FEMDiscretising::setShapeFuncsOnPtcls(InputFileData* InputData,
					   std::map<std::string,double>& modelData,
					   std::ofstream& logFile,
					   PetscViewer& viewerMPI,
					   PetscViewer& viewerSEQ) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int derivationOrder =
    (int)modelData["shapesDerivationOrderOnPtcls"];

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"*********************** nodal shapes *********************"<<endl;
  logFile<<"**********************************************************"<<endl;
#endif


  int supportSize;

  // Loop over processor's locally needed particles portion.
  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    dbVector& shapeFuncs = particles[ptcle].getShapeFuncs();
    dbMatrix& firstDerivShapes = 
      particles[ptcle].getFirstDerivShapes();

#ifdef _geometryDebugMode_
    logFile<<"Node "<<ptcle<<endl;
#endif
    
    // set the shapefunctions
    shapeFuncs.push_back(1.0); 
    
    for(int j=0;j<3;j++)
      firstDerivShapes[j].resize(supportSize);
    

    // set the shapefunction derivatives - not continuous!

  }


#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"*********** calculated shape sets on ptcls ***********"<<endl;
  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcleID = localGlobalPtcls[i];
    Particle& ptcle = particles[ptcleID];
    intVector& suppPtcls = ptcle.getSupportPtcls();
    dbVector& suppShapes = ptcle.getShapeFuncs();
    dbMatrix& firstDerivs = ptcle.getFirstDerivShapes();
    logFile<<"Ptcle "<<ptcleID<<"("<<oldIdx[ptcleID]<<"):: ";
    for(int j=0;j<suppPtcls.size();j++) {
      logFile<<suppPtcls[j]<<": "
	     <<suppShapes[j]<<" | ";
    }
    logFile<<endl;
    logFile<<"first order derivations"<<endl;
    for(int k=0;k<firstDerivs.size();k++) {
      for(int j=0;j<firstDerivs[k].size();j++)
	logFile<<suppPtcls[j]<<": "<<firstDerivs[k][j]<<" | ";
      logFile<<endl;
    }
    logFile<<endl;
  }
#endif


}
