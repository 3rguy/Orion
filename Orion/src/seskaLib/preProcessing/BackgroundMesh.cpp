# include "BackgroundMesh.h"


/**********************************************************************/
/**********************************************************************/
// Copy the FEM background mesh, respectively nodal coordinates and
// Gauss integration scheme.
void BackgroundMesh::copyFEData(FEMGeometry* FEData,
                                InputFileData* InputData,
                                std::map<std::string,double>& modelData,
                                std::ofstream& logFile) {

  using namespace std;

  maxGaussPtsPerVolumeElem = getMaxGaussPtsPerVolumeElem();

  int modeRadius = (int)InputData->getValue("radiusDeterminationAlgorithm");
  int modeSupport = (int)InputData->getValue("supportComputationMode");

  bool dirichletIntegration =
    (bool)modelData["dirichletBoundaryIntegration"];

  elementRootList = FEData->getElementRootList();
  newLocalElemIdx = FEData->getNewLocalElemIdx();
  elemGaussIdx = FEData->getElemGaussIdx();

  // Element node connectivity list is only needed to calculate the
  // influence radii indepently on the particle weight.
  if(modeRadius == 4 || modeRadius == 6 || modeSupport == 2) {
    nodesElements = FEData->getNodesElemsVec();
    globalLocalElemIdx = FEData->getGlobalLocalElemIdx();
  }

  vector<FEMElement>& dummyElemVec = FEData->getNodesElemsVec();
  dummyElemVec.resize(0,FEMElement(0));

  /*********************************************************************/
  // Determine the number of particles 'particlesNum'
  particlesNum =  FEData->getNodesNum();

  // Contains the connectivity new particle identifiers to the original
  // read ones.
  oldIdx = FEData->getOldIdx();

  // Contains the connectivity original read particle identifiers to the
  // new sorted ones.
  newIdx = FEData->getNewIdx();

  //intVector& dummyIntVec = FEData->getOldIdx();
  //clearArray(dummyIntVec);
  //dummyIntVec = FEData->getNewIdx();
  //clearArray(dummyIntVec);

  // --------------------------------------------------------------------
  // Set vectors containing the indices Gauss points which have various
  // contraints applied on.

  // Indices of all boundary Gauss points any deformation boundary
  // conditions are applies
  surfaceDispBoundGaussPtsIdx = FEData->getSurfaceDispBoundGaussPtsIdx();
  lineDispBoundGaussPtsIdx = FEData->getLineDispBoundGaussPtsIdx();
  pointDispBoundPtcleIdx = FEData->getPointDispBoundPtcleIdx();

  surfaceRotBoundGaussPtsIdx = FEData->getSurfaceRotBoundGaussPtsIdx();
  lineRotBoundGaussPtsIdx = FEData->getLineRotBoundGaussPtsIdx();
  pointRotBoundPtcleIdx = FEData->getPointRotBoundPtcleIdx();

  surfaceElectricBoundGaussPtsIdx = FEData->getSurfaceElectricBoundGaussPtsIdx();
  lineElectricBoundGaussPtsIdx = FEData->getLineElectricBoundGaussPtsIdx();
  pointElectricBoundPtcleIdx = FEData->getPointElectricBoundPtcleIdx();

  surfaceDepolarisationBoundGaussPtsIdx = FEData->getSurfaceDepolarisationBoundGaussPtsIdx();
  lineDepolarisationBoundGaussPtsIdx = FEData->getLineDepolarisationBoundGaussPtsIdx();
  pointDepolarisationBoundPtcleIdx = FEData->getPointDepolarisationBoundPtcleIdx();

  surfaceMicroBoundGaussPtsIdx = FEData->getSurfaceMicroBoundGaussPtsIdx();
  lineMicroBoundGaussPtsIdx = FEData->getLineMicroBoundGaussPtsIdx();
  //pointMicroBoundPtcleIdx = FEData->getPointMicroBoundPtcleIdx();

  // Indices of all boundary Gauss points any boundary loading is
  // applied
  pointForceBoundPtcleIdx = FEData->getPointForceBoundPtcleIdx();

  lineForceBoundGaussPtsIdx =
    FEData->getLineForceBoundGaussPtsIdx();

  tractionBoundGaussPtsIdx = FEData->getTractionBoundGaussPtsIdx();

  surfacePressureBoundGaussPtsIdx =
    FEData->getSurfacePressureBoundGaussPtsIdx();

  surfaceMomentBoundGaussPtsIdx =
    FEData->getSurfaceMomentBoundGaussPtsIdx();

  surfaceElectricChargeBoundGaussPtsIdx =
    FEData->getSurfaceElectricChargeBoundGaussPtsIdx();

  // Indices of local gauss points any body forces and moments are
  // applied.
  bodyForceGaussPtsIdx = FEData->getBodyForceGaussPtsIdx();

  bodyMomentGaussPtsIdx = FEData->getBodyMomentGaussPtsIdx();

  bodyElectricChargeGaussPtsIdx = FEData->getBodyElectricChargeGaussPtsIdx();


  intMatrix& dummyIntMat = FEData->getSurfaceDispBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getLineDispBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getPointDispBoundPtcleIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getSurfaceRotBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getLineRotBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getPointRotBoundPtcleIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getSurfaceElectricBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getLineElectricBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getPointElectricBoundPtcleIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getSurfaceDepolarisationBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getLineDepolarisationBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getPointDepolarisationBoundPtcleIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getSurfaceMicroBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getLineMicroBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getPointForceBoundPtcleIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getLineForceBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getTractionBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getSurfacePressureBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getSurfaceMomentBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getSurfaceElectricChargeBoundGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getBodyForceGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getBodyMomentGaussPtsIdx();
  resizeArray(dummyIntMat,0);
  dummyIntMat = FEData->getBodyElectricChargeGaussPtsIdx();
  resizeArray(dummyIntMat,0);

  //---------------------------------------------------------------------
  // Set the array 'particles' containing the particle properties.
  particles = FEData->getNodesVec();

  vector<Particle>& dummyPtcleVec = FEData->getNodesVec();
  vector<Particle>(0,Particle(0)).swap(dummyPtcleVec);
  //dummyPtcleVec.resize(0,Particle(0));

  // Determine the global and the local number of gauss points.
  globalGaussPtsNum = FEData->getGlobalGaussPtsNum();
  localGaussPtsNum = FEData->getLocalGaussPtsNum();

  // Set the vector 'gaussPoints' containing all gauss points.
  gaussPoints = FEData->getGaussPointsVec();

  vector<GaussPoint>& dummyGPointVec = FEData->getGaussPointsVec();
  vector<GaussPoint>(0).swap(dummyGPointVec);
  //dummyGPointVec.resize(0);

  // Determine the number of boundary gauss points.
  globalBGaussPtsNum = FEData->getGlobalBGaussPtsNum();

  // Set the vector 'boundGaussPoints' containing all necessary boundary
  // gauss points.
  boundGaussPoints = FEData->getBoundGaussPtsVec();

  dummyGPointVec = FEData->getBoundGaussPtsVec();
  vector<GaussPoint>(0).swap(dummyGPointVec);
  //dummyGPointVec.resize(0);

  // copy the boundary point integration particles
  if(dirichletIntegration) {

    int oldSize = boundGaussPoints.size();
    int deltaSize = pointDispBoundPtcleIdx.size();
    boundGaussPoints.resize(oldSize + deltaSize);
    globalBGaussPtsNum += pointDispBoundPtcleIdx.size();

    for(int i=0,j=oldSize;i<pointDispBoundPtcleIdx.size();i++,j++) {

      int& ptcle = pointDispBoundPtcleIdx[i][0];
      int& weightID = pointDispBoundPtcleIdx[i][1];
      int& condID = pointDispBoundPtcleIdx[i][2];
      int& surfaceID = pointDispBoundPtcleIdx[i][3];
      double& weight = particles[ptcle].getIntWeight(weightID);
      dbVector& boundConds = particles[ptcle].getDeformationBoundConds(condID);
      blVector& affectedDOF = particles[ptcle].getDeformationBoundDOF(condID);
      dbVector& surfaceNormal = particles[ptcle].getSurfaceNormal(surfaceID);
      dbVector& coords = particles[ptcle].getCoords();
      int& matID = particles[ptcle].getMaterialID();

      int& gaussID = boundGaussPoints[j].getGlobalID();
      double& gaussWeight = boundGaussPoints[j].getIntWeight(0);
      dbVector& gaussConds = boundGaussPoints[j].getDeformationBoundConds(condID);
      blVector& gaussDOF = boundGaussPoints[j].getDeformationBoundDOF(condID);
      dbVector& gaussNormal = boundGaussPoints[j].getSurfaceNormal();
      dbVector& gaussCoords = boundGaussPoints[j].getCoords();
      int& gaussMatID = boundGaussPoints[j].getMaterialID();

      pointDispBoundPtcleIdx[i][0] = j;
      pointDispBoundPtcleIdx[i][1] = 0;

      gaussID = j;
      gaussWeight = weight;
      gaussConds = boundConds;
      gaussDOF = affectedDOF;
      gaussNormal = surfaceNormal;
      gaussCoords = coords;
      gaussMatID = matID;

    }


    // -------------------------------------------------------------------

    oldSize = boundGaussPoints.size();
    deltaSize = pointRotBoundPtcleIdx.size();
    boundGaussPoints.resize(oldSize + deltaSize);
    globalBGaussPtsNum += pointRotBoundPtcleIdx.size();

    for(int i=0,j=oldSize;i<pointRotBoundPtcleIdx.size();i++,j++) {

      int& ptcle = pointRotBoundPtcleIdx[i][0];
      int& weightID = pointRotBoundPtcleIdx[i][1];
      int& condID = pointRotBoundPtcleIdx[i][2];
      int& surfaceID = pointRotBoundPtcleIdx[i][3];
      double& weight = particles[ptcle].getIntWeight(weightID);
      dbVector& boundConds = particles[ptcle].getDeformationBoundConds(condID);
      blVector& affectedDOF = particles[ptcle].getDeformationBoundDOF(condID);
      dbVector& surfaceNormal = particles[ptcle].getSurfaceNormal(surfaceID);
      dbVector& coords = particles[ptcle].getCoords();
      int& matID = particles[ptcle].getMaterialID();

      int& gaussID = boundGaussPoints[j].getGlobalID();
      double& gaussWeight = boundGaussPoints[j].getIntWeight(0);
      dbVector& gaussConds = boundGaussPoints[j].getDeformationBoundConds(condID);
      blVector& gaussDOF = boundGaussPoints[j].getDeformationBoundDOF(condID);
      dbVector& gaussNormal = boundGaussPoints[j].getSurfaceNormal();
      dbVector& gaussCoords = boundGaussPoints[j].getCoords();
      int& gaussMatID = boundGaussPoints[j].getMaterialID();

      pointRotBoundPtcleIdx[i][0] = j;
      pointRotBoundPtcleIdx[i][1] = 0;

      gaussID = j;
      gaussWeight = weight;
      gaussConds = boundConds;
      gaussDOF = affectedDOF;
      gaussNormal = surfaceNormal;
      gaussCoords = coords;
      gaussMatID = matID;

    }

    // -------------------------------------------------------------------

    oldSize = boundGaussPoints.size();
    deltaSize = pointElectricBoundPtcleIdx.size();
    boundGaussPoints.resize(oldSize + deltaSize);
    globalBGaussPtsNum += pointElectricBoundPtcleIdx.size();

    for(int i=0,j=oldSize;i<pointElectricBoundPtcleIdx.size();i++,j++) {

      int& ptcle = pointElectricBoundPtcleIdx[i][0];
      int& weightID = pointElectricBoundPtcleIdx[i][1];
      int& condID = pointElectricBoundPtcleIdx[i][2];
      int& surfaceID = pointElectricBoundPtcleIdx[i][3];
      double& weight = particles[ptcle].getIntWeight(weightID);
      dbVector& boundConds = particles[ptcle].getDeformationBoundConds(condID);
      blVector& affectedDOF = particles[ptcle].getDeformationBoundDOF(condID);
      dbVector& surfaceNormal = particles[ptcle].getSurfaceNormal(surfaceID);
      dbVector& coords = particles[ptcle].getCoords();
      int& matID = particles[ptcle].getMaterialID();

      int& gaussID = boundGaussPoints[j].getGlobalID();
      double& gaussWeight = boundGaussPoints[j].getIntWeight(0);
      dbVector& gaussConds = boundGaussPoints[j].getDeformationBoundConds(condID);
      blVector& gaussDOF = boundGaussPoints[j].getDeformationBoundDOF(condID);
      dbVector& gaussNormal = boundGaussPoints[j].getSurfaceNormal();
      dbVector& gaussCoords = boundGaussPoints[j].getCoords();
      int& gaussMatID = boundGaussPoints[j].getMaterialID();

      pointElectricBoundPtcleIdx[i][0] = j;
      pointElectricBoundPtcleIdx[i][1] = 0;

      gaussID = j;
      gaussWeight = weight;
      gaussConds = boundConds;
      gaussDOF = affectedDOF;
      gaussNormal = surfaceNormal;
      gaussCoords = coords;
      gaussMatID = matID;

    }

    //TEMP: needed? yes, for points.
    // -------------------------------------------------------------------

    oldSize = boundGaussPoints.size();
    deltaSize = pointDepolarisationBoundPtcleIdx.size();
    boundGaussPoints.resize(oldSize + deltaSize);
    globalBGaussPtsNum += pointDepolarisationBoundPtcleIdx.size();

    for(int i=0,j=oldSize;i<pointDepolarisationBoundPtcleIdx.size();i++,j++) {

      int& ptcle = pointDepolarisationBoundPtcleIdx[i][0];
      int& weightID = pointDepolarisationBoundPtcleIdx[i][1];
      int& condID = pointDepolarisationBoundPtcleIdx[i][2];
      int& surfaceID = pointDepolarisationBoundPtcleIdx[i][3];
      double& weight = particles[ptcle].getIntWeight(weightID);
      dbVector& boundConds = particles[ptcle].getDeformationBoundConds(condID);
      blVector& affectedDOF = particles[ptcle].getDeformationBoundDOF(condID);
      dbVector& surfaceNormal = particles[ptcle].getSurfaceNormal(surfaceID);
      dbVector& coords = particles[ptcle].getCoords();
      int& matID = particles[ptcle].getMaterialID();

      int& gaussID = boundGaussPoints[j].getGlobalID();
      double& gaussWeight = boundGaussPoints[j].getIntWeight(0);
      dbVector& gaussConds = boundGaussPoints[j].getDeformationBoundConds(condID);
      blVector& gaussDOF = boundGaussPoints[j].getDeformationBoundDOF(condID);
      dbVector& gaussNormal = boundGaussPoints[j].getSurfaceNormal();
      dbVector& gaussCoords = boundGaussPoints[j].getCoords();
      int& gaussMatID = boundGaussPoints[j].getMaterialID();

      pointDepolarisationBoundPtcleIdx[i][0] = j;
      pointDepolarisationBoundPtcleIdx[i][1] = 0;

      gaussID = j;
      gaussWeight = weight;
      gaussConds = boundConds;
      gaussDOF = affectedDOF;
      gaussNormal = surfaceNormal;
      gaussCoords = coords;
      gaussMatID = matID;

    }

  }

  else {
    resizeArray(pointDispBoundPtcleIdx,0);
    resizeArray(pointRotBoundPtcleIdx,0);
    resizeArray(pointElectricBoundPtcleIdx,0);
    resizeArray(pointDepolarisationBoundPtcleIdx,0);
  }

//   // Delete of all FEM data, because there is no further need.
//   delete FEData;

}

/************************************************************************/
/************************************************************************/
// Calculate the particles' influence radii element dependent.
void BackgroundMesh::setPtcleRadsElemDepend(InputFileData* InputData,
                                            dbVector& globalInfluenceRadii,
                                            std::map<std::string,double>& modelData,
                                            std::ofstream& logFile) {


  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  bool ptcleSet;
  int neighPtcle,ptcle,support;
  int maxSupport = 0;
  double distance;

  dbVector localInfluenceRadii(particlesNum*3);
  dbMatrix distances;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  if(elementRootList.size() < globalLocalElemIdx.size()) {

    intVector localElemRootList(globalLocalElemIdx.size());
    elementRootList.resize(globalLocalElemIdx.size());

    //   intVector localElemRootList(globalLocalElemIdx.size());
    //   intVector elementRootList(globalLocalElemIdx.size());

    for(int i=0;i<globalLocalElemIdx.size();i++)

      if(globalLocalElemIdx[i] >= 0)
	localElemRootList[i] = rank;

    MPI_Allreduce(&localElemRootList[0],&elementRootList[0],
		  globalLocalElemIdx.size(),MPI_INT,MPI_SUM,
		  MPI_COMM_WORLD);

    localElemRootList.resize(0);

#ifdef _influenceRadiusDebugMode_
    logFile<<"######################################################"<<endl;
    logFile<<"****** compute element dependent influence radii *****"<<endl;
    logFile<<"******************************************************"<<endl;
    logFile<<"***************** element root list ******************"<<endl;
    for(int i=0;i<globalLocalElemIdx.size();i++)
      logFile<<"Element "<<i<<": "<<elementRootList[i]<<endl;
#endif

  }

  int sourceIdx;

  intVector recvCounter(size);
  intMatrix recvIdx(size,intVector(nodesElements.size()));

  int sendIdxSize;
  intVector sendIdx;

  int nonLocalNodesSize;
  intVector sendCounter(size);
  intMatrix sendNodes(size,intVector(exclusiveLocalPtcls.size()));

  intVector nonLocalNodes;

  intVector recvCounts(size);
  intVector displs(size);

  // Loop over processor's particles portion.
  //for(int i=0;i<endIdx;i++) {
  for(int i=0;i<particles.size();i++) {
    support = 0;
    ptcleSet = false;
    distances.resize(0);

    intVector& elements = particles[i].getElems();

    string sPrint = "Particles[i]";

    printVector(elements,sPrint.c_str(),logFile);

    // Loop over all elements the current particle belongs to.
    for(int j=0;j<elements.size();j++) {
      int& localID = globalLocalElemIdx[elements[j]];

      // check if local element
      if(localID < 0) continue;

      intVector& nodes_j = nodesElements[localID].getNodes();

      // Loop over all element's particles.
      for(int k=0;k<nodes_j.size();k++) {
	ptcle = nodes_j[k]-1;

	// Check if current neighbour particle equals current particle.
	if(ptcle == i)
	  continue;

	// Check if current neighbour particle is already set.
	else {

	  for(int m=0;m<distances.size();m++)

	    if(ptcle == distances[m][3]) {
	      ptcleSet = true;
	      break;
	    }

	}

	if(ptcleSet) {
	  ptcleSet = false;
	  continue;
	}

	if(distances.size() <= support)
	  distances.resize(support+1,dbVector(usedDims+1));

	// coordinate x
	distances[support][0] = fabs(particles[ptcle].getCoord(0)
				     -particles[i].getCoord(0));

	if(localInfluenceRadii[i*3] < distances[support][0])
	  localInfluenceRadii[i*3] = distances[support][0];

	// coordinate y
	distances[support][1]= fabs(particles[ptcle].getCoord(1)
				    -particles[i].getCoord(1));

	if(localInfluenceRadii[i*3+1] < distances[support][1])
	  localInfluenceRadii[i*3+1] = distances[support][1];

	// coordinate z
	distances[support][2] = fabs(particles[ptcle].getCoord(2)
				     -particles[i].getCoord(2));

	if(localInfluenceRadii[i*3+2] < distances[support][2])
	  localInfluenceRadii[i*3+2] = distances[support][2];

	distances[support][3] = ptcle;

	support++;
      }

      if(support > maxSupport)
	maxSupport = support;
    }

#ifdef _influenceRadiusDebugMode_
    logFile<<"PARTICLE "<<i<<": motherElems: ";
    for(int m=0;m<elements.size();m++)
      logFile<<globalLocalElemIdx[elements[m]]<<" ";
    logFile<<" own elements ptcle support = "
	   <<support<<"; "<<endl;
    for(int m=0;m<distances.size();m++)
      logFile<<"neigh ptcle "<<distances[m][3]<<" dists: "
	     <<distances[m][0]<<" "
	     <<distances[m][1]<<" "
	     <<distances[m][2]<<endl;
    logFile<<"radii: "<<localInfluenceRadii[i*usedDims]<<" "
	   <<localInfluenceRadii[i*usedDims+1]<<" "
	   <<localInfluenceRadii[i*usedDims+2]<<endl;
    logFile<<"maxSupport = "<<maxSupport<<endl;
    logFile<<"-----------------------------------"<<endl;
#endif

    /*******************************************************************/
    // incorporation of the neighbour particles mother elements and their
    // nodes

    bool motherElemsOnly =
      InputData->getValue("setPtcleRadsElemDependMotherElemsOnly");

    if(motherElemsOnly) continue;

    clearArray(recvCounter);

    // Incorporate the neighbour particles' neighbours too.
    for(int m=0;m<elements.size();m++) {
      int& localID = globalLocalElemIdx[elements[m]];

#ifdef _influenceRadiusDebugMode_
      logFile<<"motherElem "<<elements[m]<<" -> localID="<<localID
	     <<" particles:"<<endl;
#endif

      // check if local element
      if(localID >= 0) {

	intVector& nodes_n = nodesElements[localID].getNodes();

	for(int n=0;n<nodes_n.size();n++) {
	  neighPtcle = nodes_n[n]-1;

#ifdef _influenceRadiusDebugMode_
	  logFile<<neighPtcle<<" daughterElems: ";
#endif

	  if(neighPtcle == i) continue;

	  intVector& neighElems = particles[neighPtcle].getElems();

	  // Loop over all neighbour particles' elements.
	  for(int j=0;j<neighElems.size();j++) {
	    int& localID_n = globalLocalElemIdx[neighElems[j]];

#ifdef _influenceRadiusDebugMode_
	    logFile<<neighElems[j]<<"("<<localID_n<<") ";
#endif

	    // check if daughter element is local
	    if(localID_n >= 0) {

	      intVector& nodes_k = nodesElements[localID_n].getNodes();

	      // Loop over all element's particles.
	      for(int k=0;k<nodes_k.size();k++) {
		ptcle = nodes_k[k]-1;
		ptcleSet = false;

		if(ptcle == i) continue;

		// ensure that already incorporated particles are skipped.
		for(int r=0;r<distances.size();r++)

		  if(distances[r][3] == ptcle) {
		    ptcleSet = true;
		    break;
		  }

		if(ptcleSet) continue;

		if(distances.size() <= support)
		  distances.resize(support+1,dbVector(usedDims+1));

		// coordinate x
		distances[support][0] = fabs(particles[ptcle].getCoord(0)
					     -particles[i].getCoord(0));

		if(localInfluenceRadii[i*3] < distances[support][0])
		  localInfluenceRadii[i*3] = distances[support][0];

		// coordinate y
		distances[support][1] = fabs(particles[ptcle].getCoord(1)
					     -particles[i].getCoord(1));

		if(localInfluenceRadii[i*3+1] < distances[support][1])
		  localInfluenceRadii[i*3+1] = distances[support][1];

		// coordinate z
		distances[support][2] = fabs(particles[ptcle].getCoord(2)
					     -particles[i].getCoord(2));

		if(localInfluenceRadii[i*3+2] < distances[support][2])
		  localInfluenceRadii[i*3+2] = distances[support][2];

		distances[support][3] = ptcle;

		support++;
	      }

	    }

	    // daugther element not local -> nodes need to be retrieved from
	    // other processes
	    else {

	      sourceIdx = elementRootList[neighElems[j]];

	      if(recvIdx[sourceIdx].size() > recvCounter[sourceIdx])
		recvIdx[sourceIdx][recvCounter[sourceIdx]] = neighElems[j];

	      else
		recvIdx[sourceIdx].push_back(neighElems[j]);

	      recvCounter[sourceIdx]++;

	    }

	  }

#ifdef _influenceRadiusDebugMode_
	  logFile<<endl;
#endif

	}

      }

    }

    // loop over all "source" processes and remove redundent entries
    for(int j=0;j<recvIdx.size();j++)

      removeRedundantEntries(recvIdx[j],0,recvCounter[j],logFile);


#ifdef _influenceRadiusDebugMode_
    logFile<<"other elements support = "<<support<<endl;
    for(int m=0;m<distances.size();m++)
      logFile<<"neigh ptcle "<<(int)distances[m][3]
	     <<" dists: " <<distances[m][0]<<" "
	     <<distances[m][1]<<" "
	     <<distances[m][2]<<endl;
    logFile<<"radii: "<<localInfluenceRadii[i*usedDims]<<" "
	   <<localInfluenceRadii[i*usedDims+1]<<" "
	   <<localInfluenceRadii[i*usedDims+2]<<endl;
    logFile<<"non-local elements:"<<endl;
    for(int j=0;j<recvIdx.size();j++) {
      logFile<<"proc-"<<j<<": ";
      for(int k=0;k<recvCounter[j];k++)
	logFile<<recvIdx[j][k]<<" ";
      logFile<<endl;
    }
    logFile<<"------------------------------------"<<endl;
#endif


    // ==================================================================
    // incorporation of non-local elements and their nodes

    clearArray(sendCounter);

    // loop over all "source" processes
    for(int j=0;j<size;j++) {

      MPI_Allgather(&recvCounter[j],1,MPI_INT,&recvCounts[0],1,
		    MPI_INT,MPI_COMM_WORLD);

      if(j == rank) {

	sendIdxSize = 0;

	for(int k=0;k<size;k++)
	  sendIdxSize += recvCounts[k];

	if(sendIdx.size() < sendIdxSize)
	  sendIdx.resize(sendIdxSize);

      }

      displs[0] = 0;

      for(int k=1;k<size;k++)
	displs[k] = displs[k-1]+recvCounts[k-1];

      // process 'j' gets from all other processes those element indices, the
      // corresponding nodes of which it needs to send back
      MPI_Gatherv(&recvIdx[j][0],recvCounts[rank],MPI_INT,&sendIdx[0],
		  &recvCounts[0],&displs[0],MPI_INT,j,MPI_COMM_WORLD);

      // assemble the being sent nodes
      if(j == rank) {

#ifdef _influenceRadiusDebugMode_
	logFile<<"sendIdx: ";
	for(int k=0;k<size;k++) {
	  for(int l=0;l<recvCounts[k];l++)
	    logFile<<sendIdx[displs[k]+l]<<" ";
	  logFile<<"| ";
	}
	logFile<<endl;
#endif

	// loop over all "destination" processes
	for(int k=0;k<recvCounts.size();k++) {

	  int n=0;

	  //  process does not send to itself
	  if(k == rank) continue;

	  else {

	    // loop over all elements current process 'k' requests
	    for(int l=0;l<recvCounts[k];l++) {

	      int& localElemIdx = globalLocalElemIdx[sendIdx[displs[k]+l]];
	      intVector& nodes_l = nodesElements[localElemIdx].getNodes();

	      // store the nodes which need to be sent
	      for(int r=0;r<nodes_l.size();r++) {

		if(sendNodes[k].size() > sendCounter[k])
		  sendNodes[k][sendCounter[k]] = nodes_l[r];

		else
		  sendNodes[k].push_back(nodes_l[r]);

		sendCounter[k]++;
	      }

	    }

	  }

	}

      }

    }

    // loop over all "destination" processes and remove redundent entries
    for(int j=0;j<sendNodes.size();j++)

      // remove redundent entries
      removeRedundantEntries(sendNodes[j],0,sendCounter[j],logFile);


#ifdef _influenceRadiusDebugMode_
    logFile<<"------------------------------------"<<endl;
    for(int j=0;j<size;j++) {
      logFile<<"proc-"<<j<<" sendNodes: ";
      for(int k=0;k<sendCounter[j];k++)
	logFile<<sendNodes[j][k]<<" ";
      logFile<<"| ";
    }
    logFile<<endl;
#endif

    // ------------------------------------------------------------------
    // exchange the the nodal indices between the processes

    // loop over all "source" processes
    for(int j=0;j<size;j++) {

      MPI_Allgather(&sendCounter[j],1,MPI_INT,&recvCounts[0],1,
		    MPI_INT,MPI_COMM_WORLD);

      if(j == rank) {

	nonLocalNodesSize = 0;

	for(int k=0;k<size;k++)
	  nonLocalNodesSize += recvCounts[k];

	if(nonLocalNodes.size() < nonLocalNodesSize)
	  nonLocalNodes.resize(nonLocalNodesSize);

      }

      displs[0] = 0;

      for(int k=1;k<size;k++)
	displs[k] = displs[k-1]+recvCounts[k-1];

      // process 'j' gets from all other processes those element indices, the
      // corresponding nodes of which it needs to send back
      MPI_Gatherv(&sendNodes[j][0],recvCounts[rank],MPI_INT,&nonLocalNodes[0],
		  &recvCounts[0],&displs[0],MPI_INT,j,MPI_COMM_WORLD);

    }

#ifdef _influenceRadiusDebugMode_
    logFile<<"nonLocalNODES(!): ";
    for(int j=0;j<nonLocalNodesSize;j++)
      logFile<<nonLocalNodes[j]<<" ";
    logFile<<endl;
#endif

    // ------------------------------------------------------------------
    // determine the distances of particle 'i' to all "non-local" element
    // nodes

    for(int j=0;j<nonLocalNodesSize;j++) {

      ptcle = nonLocalNodes[j]-1;
      ptcleSet = false;

      if(ptcle == i) continue;

      // ensure that already incorporated particles are skipped.
      for(int r=0;r<distances.size();r++)

	if(distances[r][3] == ptcle) {
	  ptcleSet = true;
	  break;
	}

      if(ptcleSet) continue;

      if(distances.size() <= support)
	distances.resize(support+1,dbVector(usedDims+1));

      // coordinate x
      distances[support][0] = fabs(particles[ptcle].getCoord(0)
				   -particles[i].getCoord(0));

      if(localInfluenceRadii[i*3] < distances[support][0])
	localInfluenceRadii[i*3] = distances[support][0];

      // coordinate y
      distances[support][1] = fabs(particles[ptcle].getCoord(1)
				   -particles[i].getCoord(1));

      if(localInfluenceRadii[i*3+1] < distances[support][1])
	localInfluenceRadii[i*3+1] = distances[support][1];

      // coordinate z
      distances[support][2] = fabs(particles[ptcle].getCoord(2)
				   -particles[i].getCoord(2));

      if(localInfluenceRadii[i*3+2] < distances[support][2])
	localInfluenceRadii[i*3+2] = distances[support][2];

      distances[support][3] = ptcle;

      support++;
    }


#ifdef _influenceRadiusDebugMode_
    logFile<<"non-local elements support = "<<support<<endl;
    for(int m=0;m<distances.size();m++)
      logFile<<"neigh ptcle "<<(int)distances[m][3]
	     <<" dists: " <<distances[m][0]<<" "
	     <<distances[m][1]<<" "
	     <<distances[m][2]<<endl;
    logFile<<"radii: "<<localInfluenceRadii[i*usedDims]<<" "
	   <<localInfluenceRadii[i*usedDims+1]<<" "
	   <<localInfluenceRadii[i*usedDims+2]<<endl;
    logFile<<"------------------------------------"<<endl;
#endif

  }


  // Finite element stuff is not needed anymore -> free space
  //globalLocalElemIdx.resize(0);
  //nodesElements = vector<FEMElement>(0,FEMElement(usedDOF));

  // ==================================================================
  // merge locally determined influence radii

  globalInfluenceRadii.resize(particlesNum*3);

  for(int i=0;i<particlesNum*3;i+=3)

    MPI_Allreduce(&localInfluenceRadii[i],&globalInfluenceRadii[i],3,
		  MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);


  // ==================================================================
  // Apply global influence radius increasements.
  double x1Addend = InputData->getValue("x1InfluenceRadiusAddend");
  double x2Addend = InputData->getValue("x2InfluenceRadiusAddend");
  double x3Addend = InputData->getValue("x3InfluenceRadiusAddend");


  // Apply increasements additive.
  if(x1Addend != 0 || x2Addend != 0 || x3Addend != 0) {

    for(int i=0;i<particlesNum*3;i+=3) {

      // 1-direction
      if(x1Addend != 0)
	globalInfluenceRadii[i] += x1Addend;

      // 2-direction
      if(x2Addend != 0)
	globalInfluenceRadii[i+1] += x2Addend;

      // 3-direction
      if(x3Addend != 0)
	globalInfluenceRadii[i+2] += x3Addend;

    }

  }

}

/************************************************************************/
/************************************************************************/
// Calculate the particles' influence radii element dependent with
// different values in negative and positive coordinate direction
void BackgroundMesh::setPtcleRadsElemDepend2(InputFileData* InputData,
                                             std::map<std::string,double>& modelData,
                                             std::ofstream& logFile) {


  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int shapeFuncType = (int)InputData->getValue("shapefunctionType");

  if(shapeFuncType != 5) {
    logFile<<"In BackgroundMesh::setPtcleRadsElemDepend2 "
	   <<"radiusDeterminationAlgorithm '6' is not supported for\n"
	   <<"shapefunctionType '"<<shapeFuncType<<"'. Choose "
	   <<"radiusDeterminationAlgorithm '4' instead."<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  //InputData->setValue("shapefunctionType",5.0);
  InputData->setValue("plusMinusDirectionDependentRadius",1.0);

  int neighPtcle,ptcle,support;
  int maxSupport = 0;
  double distance;

  int vecSize = 2*usedDims;
  dbVector localInfluenceRadii(particlesNum*vecSize);
  dbMatrix distances;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  if(elementRootList.size() < globalLocalElemIdx.size()) {

    intVector localElemRootList(globalLocalElemIdx.size());
    elementRootList.resize(globalLocalElemIdx.size());

    //   intVector localElemRootList(globalLocalElemIdx.size());
    //   intVector elementRootList(globalLocalElemIdx.size());

    for(int i=0;i<globalLocalElemIdx.size();i++)

      if(globalLocalElemIdx[i] >= 0)
	localElemRootList[i] = rank;

    MPI_Allreduce(&localElemRootList[0],&elementRootList[0],
		  globalLocalElemIdx.size(),MPI_INT,MPI_SUM,
		  MPI_COMM_WORLD);

    localElemRootList.resize(0);

#ifdef _influenceRadiusDebugMode_
    logFile<<"######################################################"<<endl;
    logFile<<"****** compute element dependent influence radii *****"<<endl;
    logFile<<"******************************************************"<<endl;
    logFile<<"***************** element root list ******************"<<endl;
    for(int i=0;i<globalLocalElemIdx.size();i++)
      logFile<<"Element "<<i<<": "<<elementRootList[i]<<endl;
#endif

  }

  int sourceIdx;
  int posIdx,negIdx;

  intVector recvCounter(size);
  intMatrix recvIdx(size,intVector(nodesElements.size()));

  int sendIdxSize;
  intVector sendIdx;

  int nonLocalNodesSize;
  intVector sendCounter(size);
  intMatrix sendNodes(size,intVector(exclusiveLocalPtcls.size()));

  intVector nonLocalNodes;

  intVector recvCounts(size);
  intVector displs(size);

  blVector includedPtcls(particles.size());

  // Loop over processor's particles portion.

  for(int i=0;i<particles.size();i++) {

    clearArray(distances);
    resizeArray(distances,maxSupport,usedDims+1);

    support = 0;

    clearArray(includedPtcls);
    includedPtcls[i] = true;

    intVector& elements = particles[i].getElems();

    // Loop over all elements the current particle belongs to.
    for(int j=0;j<elements.size();j++) {
      int& localID = globalLocalElemIdx[elements[j]];

      // check if local element
      if(localID < 0) continue;

      intVector& nodes_j = nodesElements[localID].getNodes();

      // Loop over all element's particles.
      for(int k=0;k<nodes_j.size();k++) {
	ptcle = nodes_j[k]-1;

	// Check if current neighbour particle is already set.
	if (includedPtcls[ptcle])
	  continue;

	// not set yet
	else {

	  if(distances.size() <= support)
	    distances.resize(support+1,dbVector(usedDims+1));


	  includedPtcls[ptcle] = true;

	  // --------------------------------------------------------------
	  // coordinate x
	  distances[support][0]= particles[ptcle].getCoord(0)
	    -particles[i].getCoord(0);

	  posIdx = i*vecSize;
	  negIdx = i*vecSize + usedDims;

	  // positive direction
	  if(distances[support][0] >= 0.0 &&
	     localInfluenceRadii[posIdx] < distances[support][0])

	    localInfluenceRadii[posIdx] = distances[support][0];

	  // negative direction
	  else if(distances[support][0] < 0.0 &&
		  localInfluenceRadii[negIdx] < fabs(distances[support][0]))

	    localInfluenceRadii[negIdx] = fabs(distances[support][0]);

	  // --------------------------------------------------------------
	  // coordinate y
	  distances[support][1]= particles[ptcle].getCoord(1)
	    -particles[i].getCoord(1);

	  posIdx += 1;
	  negIdx += 1;

	  // positive direction
	  if(distances[support][1] >= 0.0 &&
	     localInfluenceRadii[posIdx] < distances[support][1])

	    localInfluenceRadii[posIdx] = distances[support][1];

	  // negative direction
	  else if(distances[support][1] < 0.0 &&
		  localInfluenceRadii[negIdx] < fabs(distances[support][1]))

	    localInfluenceRadii[negIdx] = fabs(distances[support][1]);


	  // --------------------------------------------------------------
	  // coordinate z
	  distances[support][2] = particles[ptcle].getCoord(2)
	    -particles[i].getCoord(2);

	  posIdx += 1;
	  negIdx += 1;


	  // positive direction
	  if(distances[support][2] >= 0.0 &&
	     localInfluenceRadii[posIdx] < distances[support][2])

	    localInfluenceRadii[posIdx] = distances[support][2];

	  // negative direction
	  else if(distances[support][2] < 0.0 &&
		  localInfluenceRadii[negIdx] < fabs(distances[support][2]))

	    localInfluenceRadii[negIdx] = fabs(distances[support][2]);



	  distances[support][3] = ptcle;

	  support++;
	}

	if(support > maxSupport)
	  maxSupport = support;
      }

    }

#ifdef _influenceRadiusDebugMode_
    logFile<<"PARTICLE "<<i<<": motherElems: ";
    for(int m=0;m<elements.size();m++)
      logFile<<globalLocalElemIdx[elements[m]]<<" ";
    logFile<<" own elements ptcle support = "
	   <<support<<"; "<<endl;
    for(int m=0;m<distances.size();m++)
      logFile<<"neigh ptcle "<<distances[m][3]<<" dists: "
	     <<distances[m][0]<<" "
	     <<distances[m][1]<<" "
	     <<distances[m][2]<<endl;
    logFile<<"radii: "<<localInfluenceRadii[i*usedDims]<<" "
	   <<localInfluenceRadii[i*usedDims+1]<<" "
	   <<localInfluenceRadii[i*usedDims+2]<<endl;
    logFile<<"maxSupport = "<<maxSupport<<endl;
    logFile<<"-----------------------------------"<<endl;
#endif

    /*******************************************************************/
    // incorporation of the neighbour particles mother elements and their
    // nodes

    bool motherElemsOnly =
      InputData->getValue("setPtcleRadsElemDependMotherElemsOnly");

    if(motherElemsOnly) continue;

    clearArray(recvCounter);

    // Incorporate the neighbour particles' neighbours too.
    for(int m=0;m<elements.size();m++) {
      int& localID = globalLocalElemIdx[elements[m]];

#ifdef _influenceRadiusDebugMode_
      logFile<<"motherElem "<<elements[m]<<" -> localID="<<localID
	     <<" particles:"<<endl;
#endif

      // check if local element
      if(localID >= 0) {

	intVector& nodes_n = nodesElements[localID].getNodes();

	for(int n=0;n<nodes_n.size();n++) {
	  neighPtcle = nodes_n[n]-1;

#ifdef _influenceRadiusDebugMode_
	  logFile<<neighPtcle<<" daughterElems: ";
#endif

	  if(neighPtcle == i) continue;

	  intVector& neighElems = particles[neighPtcle].getElems();

	  // Loop over all neighbour particles' elements.
	  for(int j=0;j<neighElems.size();j++) {
	    int& localID_n = globalLocalElemIdx[neighElems[j]];

#ifdef _influenceRadiusDebugMode_
	    logFile<<neighElems[j]<<"("<<localID_n<<") ";
#endif

	    // check if daughter element is local
	    if(localID_n >= 0) {

	      intVector& nodes_k = nodesElements[localID_n].getNodes();

	      // Loop over all element's particles.
	      for(int k=0;k<nodes_k.size();k++) {
		ptcle = nodes_k[k]-1;

		// Check if current neighbour particle is already set.
		if (includedPtcls[ptcle])
		  continue;

		// not set yet
		else {


		  if(distances.size() <= support)
		    distances.resize(support+1,dbVector(usedDims+1));


		  includedPtcls[ptcle] = true;

		  // -------------------------------------------------------
		  // coordinate x
		  distances[support][0] = particles[ptcle].getCoord(0)
		    -particles[i].getCoord(0);

		  posIdx = i*vecSize;
		  negIdx = i*vecSize + usedDims;

		  // positive direction
		  if(distances[support][0] >= 0.0 &&
		     localInfluenceRadii[posIdx] < distances[support][0])

		    localInfluenceRadii[posIdx] = distances[support][0];

		  // negative direction
		  else if(distances[support][0] < 0.0 &&
			  localInfluenceRadii[negIdx] < fabs(distances[support][0]))

		    localInfluenceRadii[negIdx] = fabs(distances[support][0]);

		  // -------------------------------------------------------
		  // coordinate y
		  distances[support][1]= particles[ptcle].getCoord(1)
		    -particles[i].getCoord(1);

		  posIdx += 1;
		  negIdx += 1;

		  // positive direction
		  if(distances[support][1] >= 0.0 &&
		     localInfluenceRadii[posIdx] < distances[support][1])

		    localInfluenceRadii[posIdx] = distances[support][1];

		  // negative direction
		  else if(distances[support][1] < 0.0 &&
			  localInfluenceRadii[negIdx] < fabs(distances[support][1]))

		    localInfluenceRadii[negIdx] = fabs(distances[support][1]);


		  // -------------------------------------------------------
		  // coordinate z
		  distances[support][2] = particles[ptcle].getCoord(2)
		    -particles[i].getCoord(2);

		  posIdx += 1;
		  negIdx += 1;


		  // positive direction
		  if(distances[support][2] >= 0.0 &&
		     localInfluenceRadii[posIdx] < distances[support][2])

		    localInfluenceRadii[posIdx] = distances[support][2];

		  // negative direction
		  else if(distances[support][2] < 0.0 &&
			  localInfluenceRadii[negIdx] < fabs(distances[support][2]))

		    localInfluenceRadii[negIdx] = fabs(distances[support][2]);


		  distances[support][3] = ptcle;

		  support++;
		}

	      }

	    }

	    // daugther element not local -> nodes need to be retrieved from
	    // other processes
	    else {

	      sourceIdx = elementRootList[neighElems[j]];

	      if(recvIdx[sourceIdx].size() > recvCounter[sourceIdx])
		recvIdx[sourceIdx][recvCounter[sourceIdx]] = neighElems[j];

	      else
		recvIdx[sourceIdx].push_back(neighElems[j]);

	      recvCounter[sourceIdx]++;

	    }

	  }

#ifdef _influenceRadiusDebugMode_
	  logFile<<endl;
#endif

	}

      }

    }

    // loop over all "source" processes and remove redundent entries
    for(int j=0;j<recvIdx.size();j++)

      removeRedundantEntries(recvIdx[j],0,recvCounter[j],logFile);


#ifdef _influenceRadiusDebugMode_
    logFile<<"other elements support = "<<support<<endl;
    for(int m=0;m<distances.size();m++)
      logFile<<"neigh ptcle "<<(int)distances[m][3]
	     <<" dists: " <<distances[m][0]<<" "
	     <<distances[m][1]<<" "
	     <<distances[m][2]<<endl;
    logFile<<"radii: "<<localInfluenceRadii[i*usedDims]<<" "
	   <<localInfluenceRadii[i*usedDims+1]<<" "
	   <<localInfluenceRadii[i*usedDims+2]<<endl;
    logFile<<"non-local elements:"<<endl;
    for(int j=0;j<recvIdx.size();j++) {
      logFile<<"proc-"<<j<<": ";
      for(int k=0;k<recvCounter[j];k++)
	logFile<<recvIdx[j][k]<<" ";
      logFile<<endl;
    }
    logFile<<"------------------------------------"<<endl;
#endif


    // ==================================================================
    // incorporation of non-local elements and their nodes

    clearArray(sendCounter);

    // loop over all "source" processes
    for(int j=0;j<size;j++) {

      MPI_Allgather(&recvCounter[j],1,MPI_INT,&recvCounts[0],1,
		    MPI_INT,MPI_COMM_WORLD);

      if(j == rank) {

	sendIdxSize = 0;

	for(int k=0;k<size;k++)
	  sendIdxSize += recvCounts[k];

	if(sendIdx.size() < sendIdxSize)
	  sendIdx.resize(sendIdxSize);

      }

      displs[0] = 0;

      for(int k=1;k<size;k++)
	displs[k] = displs[k-1]+recvCounts[k-1];

      // process 'j' gets from all other processes those element indices, the
      // corresponding nodes of which it needs to send back
      MPI_Gatherv(&recvIdx[j][0],recvCounts[rank],MPI_INT,&sendIdx[0],
		  &recvCounts[0],&displs[0],MPI_INT,j,MPI_COMM_WORLD);

      // assemble the being sent nodes
      if(j == rank) {

#ifdef _influenceRadiusDebugMode_
	logFile<<"sendIdx: ";
	for(int k=0;k<size;k++) {
	  for(int l=0;l<recvCounts[k];l++)
	    logFile<<sendIdx[displs[k]+l]<<" ";
	  logFile<<"| ";
	}
	logFile<<endl;
#endif

	// loop over all "destination" processes
	for(int k=0;k<recvCounts.size();k++) {

	  int n=0;

	  //  process does not send to itself
	  if(k == rank) continue;

	  else {

	    // loop over all elements current process 'k' requests
	    for(int l=0;l<recvCounts[k];l++) {

	      int& localElemIdx = globalLocalElemIdx[sendIdx[displs[k]+l]];
	      intVector& nodes_l = nodesElements[localElemIdx].getNodes();

	      // store the nodes which need to be sent
	      for(int r=0;r<nodes_l.size();r++) {

		if(sendNodes[k].size() > sendCounter[k])
		  sendNodes[k][sendCounter[k]] = nodes_l[r];

		else
		  sendNodes[k].push_back(nodes_l[r]);

		sendCounter[k]++;
	      }

	    }

	  }

	}

      }

    }

    // loop over all "destination" processes and remove redundent entries
    for(int j=0;j<sendNodes.size();j++)

      // remove redundent entries
      removeRedundantEntries(sendNodes[j],0,sendCounter[j],logFile);


#ifdef _influenceRadiusDebugMode_
    logFile<<"------------------------------------"<<endl;
    for(int j=0;j<size;j++) {
      logFile<<"proc-"<<j<<" sendNodes: ";
      for(int k=0;k<sendCounter[j];k++)
	logFile<<sendNodes[j][k]<<" ";
      logFile<<"| ";
    }
    logFile<<endl;
#endif

    // ------------------------------------------------------------------
    // exchange the the nodal indices between the processes

    // loop over all "source" processes
    for(int j=0;j<size;j++) {

      MPI_Allgather(&sendCounter[j],1,MPI_INT,&recvCounts[0],1,
		    MPI_INT,MPI_COMM_WORLD);

      if(j == rank) {

	nonLocalNodesSize = 0;

	for(int k=0;k<size;k++)
	  nonLocalNodesSize += recvCounts[k];

	if(nonLocalNodes.size() < nonLocalNodesSize)
	  nonLocalNodes.resize(nonLocalNodesSize);

      }

      displs[0] = 0;

      for(int k=1;k<size;k++)
	displs[k] = displs[k-1]+recvCounts[k-1];

      // process 'j' gets from all other processes those element indices, the
      // corresponding nodes of which it needs to send back
      MPI_Gatherv(&sendNodes[j][0],recvCounts[rank],MPI_INT,&nonLocalNodes[0],
		  &recvCounts[0],&displs[0],MPI_INT,j,MPI_COMM_WORLD);

    }

#ifdef _influenceRadiusDebugMode_
    logFile<<"nonLocalNODES(!): ";
    for(int j=0;j<nonLocalNodesSize;j++)
      logFile<<nonLocalNodes[j]<<" ";
    logFile<<endl;
#endif

    // ------------------------------------------------------------------
    // determine the distances of particle 'i' to all "non-local" element
    // nodes

    for(int j=0;j<nonLocalNodesSize;j++) {

      ptcle = nonLocalNodes[j]-1;

      // Check if current neighbour particle is already set.
      if (includedPtcls[ptcle])
	continue;

      // not set yet
      else {


	if(distances.size() <= support)
	  distances.resize(support+1,dbVector(usedDims+1));


	includedPtcls[ptcle] = true;


	// -------------------------------------------------------
	// coordinate x
	distances[support][0]= particles[ptcle].getCoord(0)
	  -particles[i].getCoord(0);

	posIdx = i*vecSize;
	negIdx = i*vecSize + usedDims;

	// positive direction
	if(distances[support][0] >= 0.0 &&
	   localInfluenceRadii[posIdx] < distances[support][0])

	  localInfluenceRadii[posIdx] = distances[support][0];

	// negative direction
	else if(distances[support][0] < 0.0 &&
		localInfluenceRadii[negIdx] < fabs(distances[support][0]))

	  localInfluenceRadii[negIdx] = fabs(distances[support][0]);

	// -------------------------------------------------------
	// coordinate y
	distances[support][1]= particles[ptcle].getCoord(1)
	  -particles[i].getCoord(1);

	posIdx += 1;
	negIdx += 1;

	// positive direction
	if(distances[support][1] >= 0.0 &&
	   localInfluenceRadii[posIdx] < distances[support][1])

	  localInfluenceRadii[posIdx] = distances[support][1];

	// negative direction
	else if(distances[support][1] < 0.0 &&
		localInfluenceRadii[negIdx] < fabs(distances[support][1]))

	  localInfluenceRadii[negIdx] = fabs(distances[support][1]);


	// -------------------------------------------------------
	// coordinate z
	distances[support][2] = particles[ptcle].getCoord(2)
	  -particles[i].getCoord(2);

	posIdx += 1;
	negIdx += 1;


	// positive direction
	if(distances[support][2] >= 0.0 &&
	   localInfluenceRadii[posIdx] < distances[support][2])

	  localInfluenceRadii[posIdx] = distances[support][2];

	// negative direction
	else if(distances[support][2] < 0.0 &&
		localInfluenceRadii[negIdx] < fabs(distances[support][2]))

	  localInfluenceRadii[negIdx] = fabs(distances[support][2]);


	distances[support][3] = ptcle;

	support++;
      }

    }


#ifdef _influenceRadiusDebugMode_
    logFile<<"non-local elements support = "<<support<<endl;
    for(int m=0;m<distances.size();m++)
      logFile<<"neigh ptcle "<<(int)distances[m][3]
	     <<" dists: " <<distances[m][0]<<" "
	     <<distances[m][1]<<" "
	     <<distances[m][2]<<endl;
    logFile<<"radii: "<<localInfluenceRadii[i*usedDims*2]<<" - "
	   <<localInfluenceRadii[i*usedDims*2+3]<<" || "
	   <<localInfluenceRadii[i*usedDims*2+1]<<" - "
	   <<localInfluenceRadii[i*usedDims*2+4]<<" || "
           <<localInfluenceRadii[i*usedDims*2+2]<<" - "
	   <<localInfluenceRadii[i*usedDims*2+5]<<endl;
    logFile<<"------------------------------------"<<endl;
#endif

  }

  // ==================================================================
  // Apply global influence radius increasements.
  double x1Addend = InputData->getValue("x1InfluenceRadiusAddend");
  double x2Addend = InputData->getValue("x2InfluenceRadiusAddend");
  double x3Addend = InputData->getValue("x3InfluenceRadiusAddend");


  // Apply increasements additive.
  if(x1Addend != 0 || x2Addend != 0 || x3Addend != 0) {

    for(int i=0;i<particlesNum*vecSize;i+=vecSize) {

      // 1-direction
      if(x1Addend != 0) {
	localInfluenceRadii[i] += x1Addend;
	localInfluenceRadii[i+3] += x1Addend;
      }

      // 2-direction
      if(x2Addend != 0) {
	localInfluenceRadii[i+1] += x2Addend;
	localInfluenceRadii[i+4] += x2Addend;
      }

      // 3-direction
      if(x3Addend != 0) {
	localInfluenceRadii[i+2] += x3Addend;
	localInfluenceRadii[i+5] += x3Addend;
      }

    }

  }

  /*********************************************************************/
  // merge locally determined influence radii

  // loop over all particles
  for(int i=0,j=0;i<particlesNum*vecSize;i+=vecSize,j++) {

    dbVector& globalInfRadii = particles[j].getRadii();
    globalInfRadii.resize(vecSize);

    MPI_Allreduce(&localInfluenceRadii[i],&globalInfRadii[0],vecSize,
		  MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  }

}

/************************************************************************/
/************************************************************************/
// Return the maximal number of Gauss points per volume element. 
int BackgroundMesh::getMaxGaussPtsPerVolumeElem() {

  using namespace std;

  if(maxGaussPtsPerVolumeElem > 0)
    return maxGaussPtsPerVolumeElem;

  else {

    int gaussPtsPerElem;

    maxGaussPtsPerVolumeElem = 0;

    for(int i=0;i<nodesElements.size();i++) {
      intVector& intPts = nodesElements[i].getVolumeIntegrationPts();
      gaussPtsPerElem = intPts.size();

      if(gaussPtsPerElem > maxGaussPtsPerVolumeElem)
	maxGaussPtsPerVolumeElem = gaussPtsPerElem;
    }
  }
}

/**********************************************************************/
/**********************************************************************/
// Determine the necessary influence radii for all particles dependent
// on the FEM background mesh.
void BackgroundMesh::setInfluenceRadii(InputFileData* InputData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile,
                                       PetscViewer& viewerMPI,
                                       PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  int mode = (int)InputData->getValue("radiusDeterminationAlgorithm");
  int minPtcleSupport = (int)InputData->getValue("minParticleSupport");
  double multiplier = InputData->getValue("influenceRadiusMultiplier");
    

  if(minPtcleSupport > particlesNum) {
    logFile<<"Too less particles in the system ("<<particlesNum<<") to "
	   <<"satisfy the minimum particle support ("
	   <<minPtcleSupport<<"!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }


  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef _influenceRadiusDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"******** Calculation of influence radii **************"<<endl;
  logFile<<"rank: "<<rank<<endl;
  logFile<<"******************** rootlist ************************"<<endl;
  for(int i=0;i<particlesNum;i++)
    logFile<<"PARTICLE "<<i<<": "<<ptcleRootList[i]<<endl;
  logFile<<"******************************************************"<<endl;
#endif

  dbVector localInfRadii;
  dbVector localInfluenceRadii;
  dbVector globalInfluenceRadii;

  switch(mode) {

    // determine particle influence radii correspondingly to the distances
    // to their neighbouring particles
  case 1:

    setPtcleRadsDistDepend(InputData,modelData,logFile);

    break;

    // ******************************************************************
    // Calculated the particle's influence radius  support dependent.
  case 2:

    setPtcleRadsSupportDepend(InputData,gaussPoints,boundGaussPoints,
			      modelData,logFile);

    for(int i=0;i<particlesNum;i++) {

      localInfRadii = particles[i].getRadii();
      dbVector& globalInfRadii = particles[i].getRadii();

      MPI_Allreduce(&localInfRadii[0],&globalInfRadii[0],3,
		    MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    }

    break;

    // ******************************************************************
    // Calculate the particle's influence radius  weight dependent.
  case 3:

    for(int i=0;i<particlesNum;i++) {

      dbVector& globalInfRadii = particles[i].getRadii();

      globalInfRadii[0] = particles[i].getWeight();
      globalInfRadii[1] = particles[i].getWeight();
      globalInfRadii[2] = particles[i].getWeight();
    }

    break;

    // ******************************************************************
    // Calculate the particles' influence radii element dependent.
  case 4:

    setPtcleRadsElemDepend(InputData,globalInfluenceRadii,modelData,logFile);

    for(int i=0,j=0;i<particlesNum*3;i+=3,j++) {

      dbVector& globalInfRadii = particles[j].getRadii();

      globalInfRadii[0] = globalInfluenceRadii[i];
      globalInfRadii[1] = globalInfluenceRadii[i+1];
      globalInfRadii[2] = globalInfluenceRadii[i+2];

    }

    break;

    // *******************************************************************
    // determine particle influence radii correspondingly to the distances
    // to their neighbouring particles, but with different values for
    // negative and positive coordinate direction
  case 5:

    setPtcleRadsDistDepend2(InputData,modelData,logFile);

    break;

    // ******************************************************************
    // Calculate the particles' influence radii element dependent with
    // different values in negative and positive coordinate direction
  case 6:

    setPtcleRadsElemDepend2(InputData,modelData,logFile);

    break;

  default:

    logFile<<"In ParticleDistribution::setInfluenceRadii Particle "
	   <<"influence radii have not been determined.\n"
	   <<"Set input datum 'radiusDeterminationAlgorithm'<1,2,3,4,5,6>"
	   <<"is not set."<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /**********************************************************************/
  // Post-process all particles' influence radii.

  postProcInfRads(InputData,modelData,logFile);
    
}

/************************************************************************/
/************************************************************************/
// Determine for all influence spheres their neighbour influence
// spheres.
void BackgroundMesh::setInflSpheresConn(InputFileData* InputData,
                                        std::ofstream& logFile) {

  using namespace std;

  intVector suppCounts(globalGaussPtsNum);
  intMatrix supportingPtcls(globalGaussPtsNum);
  int localMaxIndirectPtcleSupport = 0; //local max indirect ptcle support on each proc

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef _influenceSphereDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"********** influence spheres assembling ************"<<endl;
#endif


  for(int i=0;i<particles.size();i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    intVector& inflSpheres = particles[i].getInflSpheres();

    inflSpheres = suppPtcls;
  }

  // Determine for all local Gausspoints the counts of supporting
  // particles, the processor that possesses these and the supporting
  // particles.

  // Loop over all local Gauss points.
  for(int i=0;i<localGaussPtsNum;i++) {
    suppCounts[gaussPoints[i].getGlobalID()] =
      gaussPoints[i].getSupportCounts();
    supportingPtcls[gaussPoints[i].getGlobalID()] =
      gaussPoints[i].getSupportPtcls();
  }

  // Determine for all global Gausspoints the counts of supporting
  // particles.
  for(int i=0;i<globalGaussPtsNum;i++)
    MPI_Bcast(&suppCounts[i],1,MPI_INT,gaussRootList[i],MPI_COMM_WORLD);

#ifdef _influenceSphereDebugMode_
  logFile<<"*********** global supported counts by gauss *********"<<endl;
  for(int i=0;i<globalGaussPtsNum;i++)
    logFile<<i<<" "<<suppCounts[i]<<endl;
#endif

  // Determine for all global Gausspoints the supporting particles.
  for(int i=0;i<globalGaussPtsNum;i++) {

    supportingPtcls[i].resize(suppCounts[i]);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&supportingPtcls[i][0],suppCounts[i],
	      MPI_INT,gaussRootList[i],MPI_COMM_WORLD);
  }

#ifdef _influenceSphereDebugMode_
  logFile<<"******* supporting particles by global gauss *********"<<endl;
  for(int i=0;i<globalGaussPtsNum;i++) {
    logFile<<"GAUSSPOINT "<<i<<": ";
    for(int j=0;j<supportingPtcls[i].size();j++)
      logFile<<supportingPtcls[i][j]<<" ";
    logFile<<endl;
  }
#endif

  // Loop over all local volume Gausspoints.
  for(int i=0;i<localGaussPtsNum;i++) {
    int globalGPIdx = gaussPoints[i].getGlobalID();

    // Loop over all supporting particles of current Gausspoint.
    for(int j=0;j<supportingPtcls[globalGPIdx].size();j++) {
      intVector& inflSpheres =
        particles[supportingPtcls[globalGPIdx][j]].getInflSpheres();

#ifdef _influenceSphereDebugMode_
      logFile<<"*************************************"<<endl;
      logFile<<globalGPIdx<<" particle "<<supportingPtcls[globalGPIdx][j]<<": ";
      for(int k=0;k<inflSpheres.size();k++)
        logFile<<inflSpheres[k]<<" ";
      logFile<<endl;
#endif

      // Loop over all supporting particles of current Gausspoint.
      for(int k=0;k<supportingPtcls[globalGPIdx].size();k++) {

#ifdef _influenceSphereDebugMode_
        logFile<<"check ptcle "<<supportingPtcls[globalGPIdx][k]<<": ";
#endif

        if(inflSpheres.size() == 0)
          inflSpheres.push_back(supportingPtcls[globalGPIdx][k]);

        // Loop over all influencing spheres of particle
        // 'supportingPtcls[i][j]'.
        for(int l=0;l<inflSpheres.size();l++) {

#ifdef _influenceSphereDebugMode_
          logFile<<inflSpheres[l]<<" ";
#endif

          if(supportingPtcls[globalGPIdx][k] == inflSpheres[l])
            break;
          else if(supportingPtcls[globalGPIdx][k] < inflSpheres[l]) {
            inflSpheres.insert(inflSpheres.begin()+l,supportingPtcls[globalGPIdx][k]);
            break;
          }
          else if(l+1 == inflSpheres.size()) {
            inflSpheres.insert(inflSpheres.end(),supportingPtcls[globalGPIdx][k]);
            break;
          }

        }

#ifdef _influenceSphereDebugMode_
        logFile<<endl;
#endif

      }

      if(inflSpheres.size() > localMaxIndirectPtcleSupport)
        localMaxIndirectPtcleSupport = inflSpheres.size();

    }


  }

#ifdef _influenceSphereDebugMode_
  logFile<<"******** influencing spheres incl vol Gauss ********"<<endl;
  logFile<<"maxIndirectPtcleSupport "<<maxIndirectPtcleSupport<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& inflParticles = particles[i].getInflSpheres();
    logFile<<"PARTICLE "<<i<<": ";
    for(int j=0;j<inflParticles.size();j++)
      logFile<<inflParticles[j]<<" ";
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // Determine for all local Gausspoints the counts of supporting
  // particles, the processor that possesses these and the supporting 
  // particles.

  if(globalBGaussPtsNum > globalGaussPtsNum) {
    suppCounts.resize(globalBGaussPtsNum);
    supportingPtcls.resize(globalBGaussPtsNum);
  }

  // Loop over all local Gauss points.
  for(int i=0;i<localBGaussPtsNum;i++) {
    suppCounts[boundGaussPoints[i].getGlobalID()] = 
      boundGaussPoints[i].getSupportCounts();
    supportingPtcls[boundGaussPoints[i].getGlobalID()] = 
      boundGaussPoints[i].getSupportPtcls();
  }

  // Determine for all global Gausspoints the counts of supporting
  // particles.
  for(int i=0;i<globalBGaussPtsNum;i++)
    MPI_Bcast(&suppCounts[i],1,MPI_INT,bGaussRootList[i],MPI_COMM_WORLD);

#ifdef _influenceSphereDebugMode_
  logFile<<"******** global supported counts by bound gauss ******"<<endl;
  for(int i=0;i<globalBGaussPtsNum;i++)
    logFile<<i<<" "<<suppCounts[i]<<endl;
#endif

  // Determine for all global boundary Gausspoints the supporting 
  // particles.
  for(int i=0;i<globalBGaussPtsNum;i++) {
    supportingPtcls[i].resize(suppCounts[i]);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&supportingPtcls[i][0],suppCounts[i],
	      MPI_INT,bGaussRootList[i],MPI_COMM_WORLD);

  }

#ifdef _influenceSphereDebugMode_
  logFile<<"****** supporting particles by global bGauss *******"<<endl;
  for(int i=0;i<supportingPtcls.size();i++) {
    logFile<<"Bound GAUSSPOINT "<<i<<": ";
    for(int j=0;j<supportingPtcls[i].size();j++)
      logFile<<supportingPtcls[i][j]<<" ";
    logFile<<endl;
  }
#endif

  // Loop over local boundary Gauss points.
  for(int i=0;i<localBGaussPtsNum;i++) {
    int globalBGPIdx = boundGaussPoints[i].getGlobalID();


    // Loop over all supporting particles of current boundary Gausspoint.
    for(int j=0;j<supportingPtcls[globalBGPIdx].size();j++) {
      intVector& inflSpheres = 
        particles[supportingPtcls[globalBGPIdx][j]].getInflSpheres();

#ifdef _influenceSphereDebugMode_
      logFile<<"****************"<<endl;
      logFile<<globalBGPIdx<<" particle "<<supportingPtcls[globalBGPIdx][j]<<": ";
      for(int k=0;k<inflSpheres.size();k++)
	logFile<<inflSpheres[k]<<" ";
      logFile<<endl;
#endif

      // Loop over all supporting particles of current boundary Gausspoint.
      for(int k=0;k<supportingPtcls[globalBGPIdx].size();k++) {

#ifdef _influenceSphereDebugMode_
        logFile<<"check ptcle "<<supportingPtcls[globalBGPIdx][k]<<": ";
#endif

	if(inflSpheres.size() == 0)
          inflSpheres.push_back(supportingPtcls[globalBGPIdx][k]);
	
	// Loop over all influencing spheres of particle 
	// 'suppPtcls[i][j]'.
	for(int l=0;l<inflSpheres.size();l++) {

#ifdef _influenceSphereDebugMode_
	  logFile<<inflSpheres[l]<<" ";
#endif

          if(supportingPtcls[globalBGPIdx][k] == inflSpheres[l])
	    break;
          else if(supportingPtcls[globalBGPIdx][k] < inflSpheres[l]) {
            inflSpheres.insert(inflSpheres.begin()+l,supportingPtcls[globalBGPIdx][k]);
	    break;
	  }
	  else if(l+1 == inflSpheres.size()) {
            inflSpheres.insert(inflSpheres.end(),supportingPtcls[globalBGPIdx][k]);
	    break;
	  }

	}

#ifdef _influenceSphereDebugMode_
	logFile<<endl;
#endif

      }

      if(inflSpheres.size() > localMaxIndirectPtcleSupport)
        localMaxIndirectPtcleSupport = inflSpheres.size();

    }
    
  }

  /***********************************************************************/
  // broadcast and merge locally determined influencing spheres with all 
  // other processes

  int localSpheresSize,globalSpheresSize,allGlobalSpheresSize;
  intVector recvCounts(size);
  intVector displs(size);
  intVector allGlobalSpheres;


  for(int i=0;i<particlesNum;i++) {
    intVector& inflSpheres = particles[i].getInflSpheres();

    localSpheresSize = inflSpheres.size();
    MPI_Allreduce(&localSpheresSize,&globalSpheresSize,1,MPI_INT,MPI_MAX,
		  MPI_COMM_WORLD);


    MPI_Allgather(&localSpheresSize,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		  MPI_COMM_WORLD);

    allGlobalSpheresSize = 0;

    for(int j=0;j<size;j++) 
      allGlobalSpheresSize += recvCounts[j];

    resizeArray(allGlobalSpheres,allGlobalSpheresSize);

    displs[0] = 0;
      
    for(int j=1;j<size;j++) 
      displs[j] = displs[j-1]+recvCounts[j-1];


    MPI_Allgatherv(&inflSpheres[0],recvCounts[rank],MPI_INT,
		   &allGlobalSpheres[0],&recvCounts[0],
		   &displs[0],MPI_INT,MPI_COMM_WORLD);

    // extract the new influencing spheres
    removeRedundantEntries(allGlobalSpheres,0,allGlobalSpheresSize,
			   logFile);
    resizeArray(allGlobalSpheres,allGlobalSpheresSize);
    
    // ----
    // store the influencing spheres
    inflSpheres = allGlobalSpheres;


    if(maxIndirectPtcleSupport < allGlobalSpheresSize)
      maxIndirectPtcleSupport = allGlobalSpheresSize;
    
  }



//   // get the maximum local indirect support
//   MPI_Allreduce(&localMaxIndirectPtcleSupport,&maxIndirectPtcleSupport,1,
// 		MPI_INT,MPI_MAX,MPI_COMM_WORLD);

//   //create inflspheres 2d vector for inflspheres on each proc
//   intMatrix allInflSpheres;

//   for(int i=0;i<particles.size();i++){
//     allocateArray(allInflSpheres,size,maxIndirectPtcleSupport);

//     // get influencing spheres on local processor
//     intVector& inflSpheres = particles[i].getInflSpheres();

//     //should have memory reallocated at end of iteration
//     int inflSpheresArray[maxIndirectPtcleSupport];

//     //convert to array
//     for(int j=0;j<inflSpheres.size();j++)
//       inflSpheresArray[j] = inflSpheres[j];
//     for(int j=inflSpheres.size();j<maxIndirectPtcleSupport;j++)
//       inflSpheresArray[j] = 0;

//     // create inflspheres 2d array for inflspheres on each proc
//     // should have memory reallocated at end of iteration
//     int allInflSpheresArray[size][maxIndirectPtcleSupport];

//     //first buffer the local inflSpheres objects,
//     // vec may be large but consistent size needed
//     inflSpheres.resize(maxIndirectPtcleSupport);

//     //gather all local inflSphere objects into global inflSphere matrix
//     //each row is a proc
//     MPI_Allgather(&inflSpheresArray[0],maxIndirectPtcleSupport,MPI_INT,
// 		  &allInflSpheresArray[0],maxIndirectPtcleSupport,MPI_INT,MPI_COMM_WORLD);

//     MPI_Barrier(MPI_COMM_WORLD); //may not be necessary


//     //convert arrays back to vectors for resizing simplicity
//     for(int k=0;k<size;k++){
//       for(int j=0;j<maxIndirectPtcleSupport;j++){
// 	if(j>0 && allInflSpheresArray[k][j] == 0){
// 	  allInflSpheres[k].resize(j);
// 	  break;
// 	}
// 	else
// 	  allInflSpheres[k][j] = allInflSpheresArray[k][j];
//       }
//     }



//     //combine local infl spheres from all procs
//     combineSortedVecs(allInflSpheres);
//     //convert to vector
//     inflSpheres = allInflSpheres[0];
//     //set maximum indirect ptcle support
//     if(inflSpheres.size() > maxIndirectPtcleSupport)
//       maxIndirectPtcleSupport = inflSpheres.size();

//   }
//   resizeArray(allInflSpheres,0);


  /***********************************************************************/
  // Check whether any particles support no Gauss point. 
  for(int i=0;i<particlesNum;i++) {
    intVector& inflSpheres = particles[i].getInflSpheres();

    if(inflSpheres.size() == 0) {
      logFile<<"BackgroundMesh::setInflSpheresConn some particles "
	     <<"don't support any Gauss point!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

#ifdef _geometryDebugMode_
  logFile<<"****************************************************"<<endl;
  logFile<<"******** influencing spheres incl bound Gauss ******"<<endl;
  int minSupp,minPtcle;
  logFile<<"maxIndirectPtcleSupport "<<maxIndirectPtcleSupport<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& inflParticles = particles[i].getInflSpheres();
    if(i == 0) {
      minSupp = inflParticles.size();
      minPtcle = i;
    }
    else if(inflParticles.size() < minSupp) {
      minSupp = inflParticles.size();
      minPtcle = i;
    }
    logFile<<"PARTICLE "<<i<<": ";
    for(int j=0;j<inflParticles.size();j++)
      logFile<<inflParticles[j]<<" ";
    logFile<<endl;
  }
  logFile<<"minSupp="<<minSupp<<" -> particle "<<minPtcle<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine for all particles their neighbour particles.
void BackgroundMesh::setPtclePtclsConn(InputFileData* InputData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int supportSize,minPtcleSupport,alphaPtcle,betaPtcle;
  int localMinSupport, minPtcleDistMode=1;
  double localMinPtcleDist;

  int localMaxSupport = 0;

#ifdef _supportDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"********** calculation of ptcle-ptcle connlist **********"<<endl;
  intMatrix allSupportLists(particlesNum);
  double oldTime = MPI_Wtime();
#endif



  int mode = (int)InputData->getValue("supportComputationMode");

  // Straight forward support list computation by checking for each particle
  // the distance to all other particles in the system
  if(mode == 1) {

    // Loop over a local portion of particles and create for each a
    // support list.
    for(int i=0;i<exclusiveLocalPtcls.size();i++) {
      int& ptcle = exclusiveLocalPtcls[i];

      dbVector& coords = particles[ptcle].getCoords();
      intVector& suppPtcls = particles[ptcle].getSupportPtcls();
      suppPtcls.resize(localMaxSupport);

      supportSize=0;

#ifdef _supportDebugMode_
      logFile<<"PARTICLE "<<i<<"**********"<<endl;
#endif

      setPointPtclsConn(InputData,coords,supportSize,suppPtcls,localMinPtcleDist,
			modelData,logFile);


      particles[ptcle].setLocalMinPtcleDistance(localMinPtcleDist);

      // -----------------------------------------------------------------

      if(supportSize > localMaxSupport)
	localMaxSupport = supportSize;

      if(i == 0)
	localMinSupport = supportSize;

      else if(supportSize < localMinSupport)
	localMinSupport = supportSize;

#ifdef _supportDebugMode_
      logFile<<"localMinSupport = "<<localMinSupport<<endl;
      logFile<<"localMaxSupport = "<<localMaxSupport<<endl;
#endif

    }

    // -------------------------------------------------------------------
    // broadcast the locally determined particle support lists to the other
    // processes

    int globalSize;

    for(int i=0;i<particlesNum;i++) {
      intVector& suppPtcls = particles[i].getSupportPtcls();

      globalSize = suppPtcls.size();

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&globalSize,1,MPI_INT,ptcleRootList[i],
		MPI_COMM_WORLD);

      if(rank != ptcleRootList[i])
	suppPtcls.resize(globalSize);

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&suppPtcls[0],globalSize,MPI_INT,ptcleRootList[i],
		MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&particles[i].getLocalMinPtcleDistance(),1,
		MPI_DOUBLE,ptcleRootList[i],
		MPI_COMM_WORLD);

    }

  }

  /**********************************************************************/
  // Determine for all particles their neighbour particles making use of
  // the finite element mesh.

  else if(mode == 2) {

    blVector includedPtcls(particlesNum);
    blVector includedElems(globalLocalElemIdx.size());
    intVector allSupportSize(particlesNum);
    intMatrix neededElemCounter(particlesNum,intVector(size));
    intMatrix3 neededElemIdx(particlesNum,intMatrix(size));
    dbMatrix localMinPtcleDist(particlesNum);
    blMatrix visitedBetaPtcls_(exclusiveLocalPtcls.size(),blVector(particlesNum));
    dbMatrix tempVisitedBetaPtcls(particlesNum);
    bool visited;

    // Define Variables for finding localMinDistance
    double ptcleDist; dbVector minPtcleDistVec(particlesNum);

    int sourceIdx;

    if(elementRootList.size() < globalLocalElemIdx.size()) {

      intVector localElemRootList(globalLocalElemIdx.size());
      elementRootList.resize(globalLocalElemIdx.size());

      for(int i=0;i<globalLocalElemIdx.size();i++)

	if(globalLocalElemIdx[i] >= 0)
	  localElemRootList[i] = rank;

      MPI_Allreduce(&localElemRootList[0],&elementRootList[0],
		    globalLocalElemIdx.size(),MPI_INT,MPI_SUM,
		    MPI_COMM_WORLD);

      resizeArray(localElemRootList,0);
    }

#ifdef _supportDebugMode_
    logFile<<"***************** element root list ***************"<<endl;
    for(int i=0;i<globalLocalElemIdx.size();i++){
      logFile<<"Element "<<i<<": "<<elementRootList[i]<<endl;
      logFile<<"GlobalLocalElemIdx : "<< globalLocalElemIdx[i]<<endl;
    }
#endif


    // -------------------------------------------------------------------
    // Loop over all particles.
    for(int i=0;i<exclusiveLocalPtcls.size();i++) {

      alphaPtcle = exclusiveLocalPtcls[i];

#ifdef _supportDebugMode_
      logFile<<"ALPHA-ptcle "<<alphaPtcle<<endl;
#endif

      intVector& suppPtcls = particles[alphaPtcle].getSupportPtcls();
      dbVector& coords = particles[alphaPtcle].getCoords();
      intVector& alphaElems = particles[alphaPtcle].getElems();
      minPtcleDistVec[alphaPtcle] = numeric_limits<double>::max();


      suppPtcls.resize(localMaxSupport);

      if(allSupportSize[alphaPtcle] < suppPtcls.size())
	suppPtcls[allSupportSize[alphaPtcle]] = alphaPtcle;

      else
	suppPtcls.push_back(alphaPtcle);

      includedPtcls[alphaPtcle] = true;
      allSupportSize[alphaPtcle] += 1;


      if(neededElemCounter[alphaPtcle].size() == 0)
	neededElemCounter[alphaPtcle].resize(size);

      if(neededElemIdx[alphaPtcle].size() == 0)
	neededElemIdx[alphaPtcle].resize(size);


      if(neededElemIdx[alphaPtcle][0].size() == 0) {

	for(int i=0;i<size;i++)
	  neededElemIdx[alphaPtcle][i].resize(nodesElements.size());

      }


      // Loop over all alpha-elements the current alpha-particle
      // belongs to.
      for(int j=0;j<alphaElems.size();j++) {
	int& localID = globalLocalElemIdx[alphaElems[j]];

#ifdef _supportDebugMode_
	logFile<<"alpha-elem "<<alphaElems[j]<<endl;
#endif

	// check if alpha-element is a local one
	if(localID >= 0) {
	  includedElems[alphaElems[j]] = true;

	  // Retreive beta Nodes for a particular local alpha element
	  intVector& betaNodes = nodesElements[localID].getNodes();

	  // Loop over the alpha element's portion of particles.
	  for(int k=0;k<betaNodes.size();k++) {
	    betaPtcle = betaNodes[k]-1;

	    if (includedElems[alphaElems[j]] == true){
	      setElemPtcleSupport(InputData,betaPtcle,
				  includedPtcls,includedElems,
				  allSupportSize[alphaPtcle],suppPtcls,
				  coords,neededElemCounter[alphaPtcle],
				  neededElemIdx[alphaPtcle],modelData,
				  logFile);
	    }

	    // Verify that the Beta particle has not been
	    // considered before or is not the alpha
	    // particle itself !
	    if (visitedBetaPtcls_[i][betaPtcle] == false && betaPtcle != alphaPtcle){
	      visitedBetaPtcls_[i][betaPtcle] = true ;
	      tempVisitedBetaPtcls[alphaPtcle].push_back(betaPtcle);
	      if(particles[betaPtcle].querySupported(InputData,
						     coords,ptcleDist,
						     modelData,logFile)) {

		// Add ptcle distance to list
		localMinPtcleDist[alphaPtcle].push_back(ptcleDist);

		// Check if particle distance is
		// smaller than the minimum one
		if (minPtcleDistVec[alphaPtcle]>ptcleDist){
		  minPtcleDistVec[alphaPtcle] = ptcleDist;
		}
	      }
	    }
	  }
	}

	// -------------------------------------------------------------------
	// element already checked.
	else if(localID >= 0 && includedElems[alphaElems[j]]) continue;

	// element not checked yet and nonlocal
	else if(localID < 0 && !includedElems[alphaElems[j]]) {

	  sourceIdx = elementRootList[alphaElems[j]];

	  if(neededElemIdx[alphaPtcle][sourceIdx].size() >
	     neededElemCounter[alphaPtcle][sourceIdx])

	    neededElemIdx[alphaPtcle][sourceIdx]
	      [neededElemCounter[alphaPtcle][sourceIdx]] = alphaElems[j];

	  else
	    neededElemIdx[alphaPtcle][sourceIdx].push_back(alphaElems[j]);

	  neededElemCounter[alphaPtcle][sourceIdx]++;

#ifdef _supportDebugMode_
	  for(int j=0;j<size;j++)
	    for(int k=0;k<neededElemCounter[alphaPtcle][j];k++)
	      logFile<<neededElemIdx[alphaPtcle][j][k]<<" ";
	  logFile<<endl;
#endif

	}


      }

      //resizeArray(suppPtcls,supportSize);

      // clear arrays includedPtcls and includedElems
      for(int j=0;j<allSupportSize[alphaPtcle];j++)
	includedPtcls[suppPtcls[j]] = false;

      clearArray(includedElems);

      if(allSupportSize[alphaPtcle] > localMaxSupport) {
	localMaxSupport = allSupportSize[alphaPtcle];
	maxPtcleSupport = allSupportSize[alphaPtcle];
      }

      // free not needed memory
      for(int j=0;j<size;j++)
	resizeArray(neededElemIdx[alphaPtcle][j],
		    neededElemCounter[alphaPtcle][j]);

    }


#ifdef _supportDebugMode_
    logFile<<"********** temp particle support list *********"<<endl;
    for(int i=0;i<particlesNum;i++) {
      intVector& suppPtcls = particles[i].getSupportPtcls();
      logFile<<"Ptcle "<<i<<": ";
      for(int j=0;j<allSupportSize[i];j++)
	logFile<<suppPtcls[j]<<" ";
      logFile<<endl;
    }
    logFile<<"************** nonlocal elements **************"<<endl;
    for(int i=0;i<particlesNum;i++) {
      logFile<<"Ptcle "<<i<<" - nonlocal elems: ";
      for(int j=0;j<size;j++) {
	logFile<<"proc "<<j<<": ";
	for(int k=0;k<neededElemCounter[i][j];k++)
	  logFile<<neededElemIdx[i][j][k]<<" ";
	logFile<<" | ";
      }
      logFile<<endl;
    }
    if(rank == 0)
      cout<<"ptcle support list computing - first round finished in "
	  <<MPI_Wtime()-oldTime<<" secs"<<endl;
    int round = 2;
#endif

    // -------------------------------------------------------------------
    // incorporate the nonlocal nodes

    bool flag;
    int localNeededElems,globalNeededElems,exclusiveCounter =0;
    intVector nonLocalNodes,nonLocalBetaNodes;
    globalNeededElems = 1;

    // loop over all particles until all processes have their particle
    // support list assembled
    while(globalNeededElems != 0) {

      localNeededElems = 0;

      for(int i=0;i<particlesNum;i++) {

	alphaPtcle = i;
	intVector& suppPtcls = particles[alphaPtcle].getSupportPtcls();

#ifdef _supportDebugMode_
	logFile<<"ALPHA-ptcle "<<alphaPtcle<<endl;
#endif

	flag = getNonLocalSupport(InputData,alphaPtcle,
				  includedPtcls,
				  neededElemCounter[alphaPtcle],
				  neededElemIdx[alphaPtcle],
				  nonLocalNodes,nonLocalBetaNodes,
				  modelData,logFile);


	// set array includedElems
	for(int j=0;j<size;j++)

	  for(int k=0;k<neededElemCounter[alphaPtcle][j];k++)
	    includedElems[neededElemIdx[alphaPtcle][j][k]] = true;

	clearArray(neededElemCounter[i]);

	// loop over the non local Beta particles to find
	// the additional supporting particles' distances
	// with respect to the alpha particle
	dbVector& coords = particles[alphaPtcle].getCoords();

	if (ptcleRootList[alphaPtcle] == rank){

	  for(int l=0;l<nonLocalBetaNodes.size();l++) {
	    betaPtcle = nonLocalBetaNodes[l]-1;

	    if (visitedBetaPtcls_[exclusiveCounter][betaPtcle] == false && betaPtcle != alphaPtcle){
	      visitedBetaPtcls_[exclusiveCounter][betaPtcle] = true ;

	      // Find distance between alpha and beta particle
	      if(particles[betaPtcle].querySupported(InputData,
						     coords,ptcleDist,
						     modelData,logFile)) {

		// Add ptcle distance to list
		localMinPtcleDist[alphaPtcle].push_back(ptcleDist);

		// Check if particle distance is
		// smaller than the minimum one
		if (minPtcleDistVec[alphaPtcle]>ptcleDist){
		  minPtcleDistVec[alphaPtcle] = ptcleDist;
		}

	      }
	    }
	  }
	  exclusiveCounter++;
	}
	nonLocalBetaNodes.clear();

	if(flag) {

	  // set the currently supporting particles
	  for(int j=0;j<allSupportSize[alphaPtcle];j++)
	    includedPtcls[suppPtcls[j]] = true;

	  // Loop over the alpha particle's nonlocal nodes
	  for(int k=0;k<nonLocalNodes.size();k++) {
	    betaPtcle = nonLocalNodes[k]-1;

	    setElemPtcleSupport(InputData,betaPtcle,
				includedPtcls,includedElems,
				allSupportSize[alphaPtcle],suppPtcls,
				coords,neededElemCounter[alphaPtcle],
				neededElemIdx[alphaPtcle],
				modelData,logFile);
	  }

	  // clear arrays includedPtcls and includedElems
	  for(int j=0;j<allSupportSize[alphaPtcle];j++)
	    includedPtcls[suppPtcls[j]] = false;

	  clearArray(includedElems);


	  if(allSupportSize[alphaPtcle] > localMaxSupport) {
	    localMaxSupport = allSupportSize[alphaPtcle];
	    maxPtcleSupport = allSupportSize[alphaPtcle];
	  }

	  if(allSupportSize[alphaPtcle] < localMinSupport)
	    localMinSupport = allSupportSize[alphaPtcle];

	}

	resizeArray(suppPtcls,allSupportSize[alphaPtcle]);
	sortIntVector(suppPtcls,0,suppPtcls.size()-1);

	// ---------------
	// loop over all processes and determine the total number of still
	// needed non-local elements
	for(int j=0;j<size;j++)

	  localNeededElems += neededElemCounter[alphaPtcle][j];

	// free not needed memory
	for(int j=0;j<size;j++)
	  resizeArray(neededElemIdx[alphaPtcle][j],
		      neededElemCounter[alphaPtcle][j]);

      }

      MPI_Allreduce(&localNeededElems,&globalNeededElems,1,
		    MPI_INT,MPI_SUM,MPI_COMM_WORLD);

#ifdef _supportDebugMode_
      if(rank == 0)
	cout<<"ptcle support list computing - "<<round<<". round finished in "
	    <<MPI_Wtime()-oldTime<<" secs"<<endl;
      round++;
#endif

    }


    double sum;
    // Compute the average distance between particles
    for (int m=0;m<particlesNum;m++){
      sum=0;

      if(!localMinPtcleDist[m].empty()){
	for (int q=0;q<localMinPtcleDist[m].size();q++){
	  sum += localMinPtcleDist[m][q];
	}

	if(minPtcleDistMode==1)
	  // Use the average suppirting particle distances
	  particles[m].
	    setLocalMinPtcleDistance(sum/localMinPtcleDist[m].size());

	else if (minPtcleDistMode==2)
	  // Use the minimum supporting particle distances
	  particles[m].
	    setLocalMinPtcleDistance(minPtcleDistVec[m]);
      }
    }

    // -------------------------------------------------------------------
    // Merge the particle support list from all processes

#ifdef _supportDebugMode_
    logFile<<"******* temp merged particle support list ******"<<endl;
#endif

    // Loop over all particles.
    for(int i=0;i<particlesNum;i++) {

      intVector& suppPtcls = particles[i].getSupportPtcls();

      supportSize = suppPtcls.size();

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&supportSize,1,MPI_INT,ptcleRootList[i],
		MPI_COMM_WORLD);

      if(rank != ptcleRootList[i])
	suppPtcls.resize(supportSize);

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&suppPtcls[0],supportSize,MPI_INT,ptcleRootList[i],
		MPI_COMM_WORLD);

#ifdef _supportDebugMode_
      logFile<<"Ptcle "<<i<<": ";
      for(int j=0;j<suppPtcls.size();j++)
	logFile<<suppPtcls[j]<<" ";
      logFile<<endl;
#endif

      if(supportSize > localMaxSupport) {
	localMaxSupport = supportSize;
	maxPtcleSupport = supportSize;
      }

      if(i == 0)
	localMinSupport = supportSize;

      else if(supportSize < localMinSupport)
	localMinSupport = supportSize;

      // broadcast the locally determined 'localMinPtcleDistance'
      // (needed for MAXENT)
      MPI_Bcast(&particles[i].getLocalMinPtcleDistance(),1,
		MPI_DOUBLE,ptcleRootList[i],
		MPI_COMM_WORLD);

#ifdef _supportDebugMode_
      logFile<<"localMinSupport = "<<localMinSupport<<endl;
      logFile<<"localMaxSupport = "<<localMaxSupport<<endl;
#endif

    }
  }

  else {
    logFile<<"In BackgroundMesh::setPtclePtclsConn particle - "
	   <<"particle\n computation mode "<<mode<<" is not supported!"
	   <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /**********************************************************************/

  // Determine global maximum number of supporting particles for a
  // particle.
  MPI_Allreduce(&localMaxSupport,&maxPtcleSupport,1,
		MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  // Determine global minimum number of supporting particles for a
  // particle.
  MPI_Allreduce(&localMinSupport,&minPtcleSupport,1,
		MPI_INT,MPI_MIN,MPI_COMM_WORLD);


  logFile<<"minPtcleSupport = "<<minPtcleSupport<<endl;
  logFile<<"maxPtcleSupport = "<<maxPtcleSupport<<endl;

  if(rank == 0)
    cout<<"particle support:    "<<minPtcleSupport<<" - "<<maxPtcleSupport<<endl;

#ifdef _geometryDebugMode_
  logFile<<"********** supporting particles of particles *********"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    logFile<<"Ptcle "<<i<<": ";
    for(int j=0;j<suppPtcls.size();j++)
      logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
  logFile<<"********** local min particle distances *********"<<endl;
  for(int i=0;i<particlesNum;i++)
    logFile<<"Ptcle "<<i<<": "
	   <<particles[i].getLocalMinPtcleDistance()<<endl;
#endif
#ifdef _supportDebugMode_
  if(rank == 0)
    cout<<"ptcle support list computing finished in "
	<<MPI_Wtime()-oldTime<<" secs"<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine for a particle its neighbour particles making use of the
// finite element mesh.
void BackgroundMesh::setElemPtcleSupport(InputFileData* InputData,
                                         int& betaPtcle,
                                         blVector& includedPtcls,
                                         blVector& includedElems,
                                         int& supportSize,
                                         intVector& suppPtcls,
                                         dbVector& coords,
                                         intVector& neededElemCounter,
                                         intMatrix& neededElemIdx,
                                         std::map<std::string,double>& modelData,
                                         std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef _supportDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"****** compute connectivity list of a particle *******"<<endl;
  logFile<<"******************************************************"<<endl;
#endif

  int sourceIdx;

  if(neededElemCounter.size() == 0)
    neededElemCounter.resize(size);

  if(neededElemIdx.size() == 0)
    neededElemIdx.resize(size);


  if(neededElemIdx[0].size() == 0) {

    for(int i=0;i<size;i++)
      neededElemIdx[i].resize(nodesElements.size());

  }

  intVector recvCounts(size);
  intVector displs(size);


  int gammaPtcle;

  // Check if alpha-particle is supported by current beta-particle
  // beta-particle is not included yet.
  if(!includedPtcls[betaPtcle] && particles[betaPtcle].
     querySupported(InputData,coords,modelData,logFile)) {


    if(supportSize < suppPtcls.size())
      suppPtcls[supportSize] = betaPtcle;

    else
      suppPtcls.push_back(betaPtcle);

    includedPtcls[betaPtcle] = true;
    supportSize += 1;

#ifdef _supportDebugMode_
    logFile<<"beta-ptcle "<<betaPtcle<<"(included)"<<endl;
#endif

    intVector& betaElems = particles[betaPtcle].getElems();

    // Loop over all beta-elements the current beta-particle belongs to.
    for(int i=0;i<betaElems.size();i++) {
      int& localID = globalLocalElemIdx[betaElems[i]];

#ifdef _supportDebugMode_
      logFile<<"beta-elem "<<betaElems[i]<<endl;
#endif

      // check if beta-element is a local one and not checked yet
      if(localID >= 0 && !includedElems[betaElems[i]]) {

	includedElems[betaElems[i]] = true;
	intVector& gammaNodes = nodesElements[localID].getNodes();

	// Loop over the gamma element's portion of particles.
	for(int j=0;j<gammaNodes.size();j++) {
	  gammaPtcle = gammaNodes[j]-1;


	  setElemPtcleSupport(InputData,gammaPtcle,includedPtcls,
			      includedElems,supportSize,suppPtcls,coords,
			      neededElemCounter,neededElemIdx,modelData,
			      logFile);

	}

      }

      // -------------------------------------------------------------------
      // element already checked.
      else if(localID >= 0 && includedElems[betaElems[i]]) continue;

      // element not checked yet and nonlocal
      else if(localID < 0 && !includedElems[betaElems[i]]) {

	sourceIdx = elementRootList[betaElems[i]];

	if(neededElemIdx[sourceIdx].size() > neededElemCounter[sourceIdx])
	  neededElemIdx[sourceIdx][neededElemCounter[sourceIdx]] = betaElems[i];

	else
	  neededElemIdx[sourceIdx].push_back(betaElems[i]);

	neededElemCounter[sourceIdx]++;

#ifdef _supportDebugMode_
	for(int j=0;j<size;j++)
	  for(int k=0;k<neededElemCounter[j];k++)
	    logFile<<neededElemIdx[j][k]<<" ";
	logFile<<endl;
#endif

      }

    }

  }

  /**********************************************************************/
  // assemble the list of elements this process needs from the other
  // processes

  // loop over all "source" processes and remove redundent entries
  for(int i=0;i<size;i++)

    removeRedundantEntries(neededElemIdx[i],0,neededElemCounter[i],logFile);

}

/**********************************************************************/
/**********************************************************************/
// Determine for a point a portion of nonlocal supporting particles.
bool BackgroundMesh::getNonLocalSupport(InputFileData* InputData,
                                        intVector& neededElemCounter,
                                        intMatrix& neededElemIdx,
                                        intVector& nonLocalNodes,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile) {


  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  intVector recvCounts(size);
  intVector displs(size);

#ifdef _supportDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"************ determine nonlocal particles ************"<<endl;
  logFile<<"******************************************************"<<endl;
#endif

  int sendIdxSize;
  intVector sendIdx;

  int nonLocalNodesSize;
  intVector sendCounter(size);
  intMatrix sendNodes(size,intVector(exclusiveLocalPtcls.size()));


  // loop over all "source" processes to get the non-local nodes
  for(int i=0;i<size;i++) {

    int localCounts = neededElemCounter[i];

    // get the number of elements needed from process "j"
    MPI_Allgather(&localCounts,1,MPI_INT,&recvCounts[0],1,MPI_INT,MPI_COMM_WORLD);

    // allocate the sendIdx array
    if(i == rank) {

      sendIdxSize = 0;

      for(int k=0;k<size;k++)
	sendIdxSize += recvCounts[k];

      if(sendIdx.size() < sendIdxSize)
	sendIdx.resize(sendIdxSize);

    }

    displs[0] = 0;

    for(int k=1;k<size;k++)
      displs[k] = displs[k-1]+recvCounts[k-1];

    // process 'i' gets from all other processes those element indices, the
    // corresponding nodes of which it needs to send back
    MPI_Gatherv(&neededElemIdx[i][0],recvCounts[rank],MPI_INT,&sendIdx[0],
		&recvCounts[0],&displs[0],MPI_INT,i,MPI_COMM_WORLD);

    // assemble the nodes to be sent
    if(i == rank) {

#ifdef _supportDebugMode_
      logFile<<"requested elems: ";
      for(int j=0;j<sendIdxSize;j++)
	logFile<<sendIdx[j]<<" ";
      logFile<<endl;
#endif

      // loop over all "destination" processes
      for(int k=0;k<size;k++) {

	int n=0;

	//  process does not send to itself
	if(k == rank) continue;

	else {

	  // loop over all elements current process 'k' requests
	  for(int l=0;l<recvCounts[k];l++) {

	    int& localElemIdx = globalLocalElemIdx[sendIdx[displs[k]+l]];
	    intVector& nodes_l = nodesElements[localElemIdx].getNodes();

	    // store the nodes which need to be sent
	    for(int r=0;r<nodes_l.size();r++) {

	      if(sendNodes[k].size() > sendCounter[k])
		sendNodes[k][sendCounter[k]] = nodes_l[r];

	      else
		sendNodes[k].push_back(nodes_l[r]);

	      sendCounter[k]++;
	    }

	  }

	}

      }

    }

  }

  // loop over all "destination" processes and remove redundent entries
  for(int i=0;i<sendNodes.size();i++)

    // remove redundent entries
    removeRedundantEntries(sendNodes[i],0,sendCounter[i],logFile);


#ifdef _supportDebugMode_
  logFile<<"requested nodes: ";
  for(int j=0;j<size;j++) {
    for(int k=0;k<sendCounter[j];k++)
      logFile<<sendNodes[j][k]<<" ";
    logFile<<" | ";
  }
  logFile<<endl;
#endif

  // ------------------------------------------------------------------
  // exchange the the nodal indices between the processes

  // loop over all "source" processes
  for(int i=0;i<size;i++) {

    MPI_Allgather(&sendCounter[i],1,MPI_INT,&recvCounts[0],1,
		  MPI_INT,MPI_COMM_WORLD);

    if(i == rank) {

      nonLocalNodesSize = 0;

      for(int k=0;k<size;k++)
	nonLocalNodesSize += recvCounts[k];

      if(nonLocalNodes.size() < nonLocalNodesSize)
	nonLocalNodes.resize(nonLocalNodesSize);

    }

    displs[0] = 0;

    for(int k=1;k<size;k++)
      displs[k] = displs[k-1]+recvCounts[k-1];

    // process 'i' gets from all other processes those element indices, the
    // corresponding nodes of which it needs to send back
    MPI_Gatherv(&sendNodes[i][0],recvCounts[rank],MPI_INT,&nonLocalNodes[0],
		&recvCounts[0],&displs[0],MPI_INT,i,MPI_COMM_WORLD);

  }

#ifdef _supportDebugMode_
  logFile<<"received nonlocal nodes(!): ";
  for(int i=0;i<nonLocalNodesSize;i++)
    logFile<<nonLocalNodes[i]<<" ";
  logFile<<endl;
#endif


  if(nonLocalNodesSize > 0)

    return true;

  else

    return false;

}

/**********************************************************************/
/**********************************************************************/
// Determine for a point a portion of nonlocal supporting particles.
bool BackgroundMesh::getNonLocalSupport(InputFileData* InputData,
                                        int& alphaPtcle,
                                        blVector& includedPtcls,
                                        intVector& neededElemCounter,
                                        intMatrix& neededElemIdx,
                                        intVector& nonLocalNodes,
                                        intVector& nonLocalBetaNodes,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile) {


  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  intVector recvCounts(size);
  intVector displs(size);
  int sendAlpha = alphaPtcle + 1;

#ifdef _supportDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"************ determine nonlocal particles ************"<<endl;
  logFile<<"******************************************************"<<endl;
#endif

  int sendIdxSize;
  intVector sendIdx;

  int nonLocalNodesSize,nonLocalBetaNodesSize;
  intVector sendCounter(size);
  intMatrix sendNodes(size,intVector(exclusiveLocalPtcls.size()));
  intVector sendBetaCounter(size);
  intMatrix sendBetaNodes(size,intVector(exclusiveLocalPtcls.size()));
  intMatrix3 sendElemNodes(size,intMatrix(globalLocalElemIdx.size()));

  //globalLocalElemIdx
  // loop over all "source" processes to get the non-local nodes
  for(int i=0;i<size;i++) {

    int localCounts = neededElemCounter[i];

    // get the number of elements needed from process "j"
    MPI_Allgather(&localCounts,1,MPI_INT,&recvCounts[0],1,MPI_INT,MPI_COMM_WORLD);

    // allocate the sendIdx array
    if(i == rank) {

      sendIdxSize = 0;

      for(int k=0;k<size;k++)
	sendIdxSize += recvCounts[k];

      if(sendIdx.size() < sendIdxSize)
	sendIdx.resize(sendIdxSize);

    }

    displs[0] = 0;

    for(int k=1;k<size;k++)
      displs[k] = displs[k-1]+recvCounts[k-1];

    // process 'i' gets from all other processes those element indices, the
    // corresponding nodes of which it needs to send back
    MPI_Gatherv(&neededElemIdx[i][0],recvCounts[rank],MPI_INT,&sendIdx[0],
		&recvCounts[0],&displs[0],MPI_INT,i,MPI_COMM_WORLD);

    // assemble the nodes to be sent
    if(i == rank) {

#ifdef _supportDebugMode_
      logFile<<"requested elems: ";
      for(int j=0;j<sendIdxSize;j++)
	logFile<<sendIdx[j]<<" ";
      logFile<<endl;
#endif

      // loop over all "destination" processes
      for(int k=0;k<size;k++) {

	int n=0;

	//  process does not send to itself
	if(k == rank) continue;

	else {

	  // loop over all elements current process 'k' requests
	  for(int l=0;l<recvCounts[k];l++) {

	    int& localElemIdx = globalLocalElemIdx[sendIdx[displs[k]+l]];
	    intVector& nodes_l = nodesElements[localElemIdx].getNodes();

	    // store the nodes which need to be sent
	    for(int r=0;r<nodes_l.size();r++) {

	      sendElemNodes[k][sendIdx[displs[k]+l]].push_back(nodes_l[r]);

	      if(sendNodes[k].size() > sendCounter[k])
		sendNodes[k][sendCounter[k]] = nodes_l[r];

	      else
		sendNodes[k].push_back(nodes_l[r]);

	      sendCounter[k]++;
	    }

	    // -------------------------
	    // assemble beta particles only
	    int pos = findIntVecPos(sendAlpha,0,nodes_l.size(),nodes_l);

	    // non-local alpha element with beta particles
	    if(pos != -1){

	      // store the beta nodes which need to be sent
	      for(int r=0;r<nodes_l.size();r++) {

		if(sendBetaNodes[k].size() > sendBetaCounter[k])
		  sendBetaNodes[k][sendBetaCounter[k]] = nodes_l[r];

		else
		  sendBetaNodes[k].push_back(nodes_l[r]);

		sendBetaCounter[k]++;
	      }
	    }
	  }

	}

      }

    }

  }

  // loop over all "destination" processes and remove redundent entries
  for(int i=0;i<sendNodes.size();i++)

    // remove redundent entries
    removeRedundantEntries(sendNodes[i],0,sendCounter[i],logFile);

  // loop over all "destination" processes and remove redundent entries
  for(int i=0;i<sendBetaNodes.size();i++)

    // remove redundent entries
    removeRedundantEntries(sendBetaNodes[i],0,sendBetaCounter[i],logFile);


#ifdef _supportDebugMode_
  logFile<<"requested nodes: ";
  for(int j=0;j<size;j++) {
    for(int k=0;k<sendCounter[j];k++)
      logFile<<sendNodes[j][k]<<" ";
    logFile<<" | ";
  }
  logFile<<endl;
  logFile<<"requested beta nodes: ";
  for(int j=0;j<size;j++) {
    for(int k=0;k<sendBetaCounter[j];k++)
      logFile<<sendBetaNodes[j][k]<<" ";
    logFile<<" | ";
  }
  logFile<<endl;
#endif

  /********************************************************************/
  // exchange the the nodal indices between the processes

  // loop over all "source" processes
  for(int i=0;i<size;i++) {

    MPI_Allgather(&sendCounter[i],1,MPI_INT,&recvCounts[0],1,
		  MPI_INT,MPI_COMM_WORLD);

    if(i == rank) {

      nonLocalNodesSize = 0;

      for(int k=0;k<size;k++)
	nonLocalNodesSize += recvCounts[k];

      if(nonLocalNodes.size() < nonLocalNodesSize)
	nonLocalNodes.resize(nonLocalNodesSize);

    }

    displs[0] = 0;

    for(int k=1;k<size;k++)
      displs[k] = displs[k-1]+recvCounts[k-1];

    // process 'i' gets from all other processes those element indices, the
    // corresponding nodes of which it needs to send back
    MPI_Gatherv(&sendNodes[i][0],recvCounts[rank],MPI_INT,&nonLocalNodes[0],
		&recvCounts[0],&displs[0],MPI_INT,i,MPI_COMM_WORLD);

  }

#ifdef _supportDebugMode_
  logFile<<"received nonlocal nodes(!): ";
  for(int i=0;i<nonLocalNodesSize;i++)
    logFile<<nonLocalNodes[i]<<" ";
  logFile<<endl;
#endif

  // -------------------------------------------------------------------
  // loop over all "source" processes
  for(int i=0;i<size;i++) {

    MPI_Allgather(&sendBetaCounter[i],1,MPI_INT,&recvCounts[0],1,
		  MPI_INT,MPI_COMM_WORLD);

    if(i == rank) {

      nonLocalBetaNodesSize = 0;

      for(int k=0;k<size;k++)
	nonLocalBetaNodesSize += recvCounts[k];

      if(nonLocalBetaNodes.size() < nonLocalBetaNodesSize)
	nonLocalBetaNodes.resize(nonLocalBetaNodesSize);

    }

    displs[0] = 0;

    for(int k=1;k<size;k++)
      displs[k] = displs[k-1]+recvCounts[k-1];

    // process 'i' gets from all other processes those element indices, the
    // corresponding beta nodes of which it needs to send back
    MPI_Gatherv(&sendBetaNodes[i][0],recvCounts[rank],MPI_INT,&nonLocalBetaNodes[0],
		&recvCounts[0],&displs[0],MPI_INT,i,MPI_COMM_WORLD);

  }

#ifdef _supportDebugMode_
  logFile<<"received nonlocal beta nodes(!): ";
  for(int i=0;i<nonLocalBetaNodesSize;i++)
    logFile<<nonLocalBetaNodes[i]<<" ";
  logFile<<endl;
#endif

  if(nonLocalNodesSize > 0)

    return true;

  else

    return false;

}

/**********************************************************************/
/**********************************************************************/
// Determine for each particle its included gauss points.
void BackgroundMesh::setPtcleGaussConn(InputFileData* InputData,
                                       dbVector& allGCoords,
                                       intVector& inclPtsStartPos,
                                       intVector& inclPtsCounts,
                                       intVector& inclPts,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];

  dbVector radius(usedDims);

  int ptclePortion = exclusiveLocalPtcls.size();
  int gaussPtsPerElem = getMaxGaussPtsPerVolumeElem();
  int idxSize = maxGaussSupport*gaussPtsPerElem*ptclePortion;

  intVector iCounts(ptclePortion);
  intVector iPoints(idxSize);

  int oldIdxSize = idxSize;
  idxSize = 0;

  //   int ptcleStartIdx = exclusiveLocalPtcls.front();
  //   int ptcleEndIdx = exclusiveLocalPtcls.back()+1;

#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"********* particle -> Gauss connectivity list ********"<<endl;
  logFile<<"ptclePortion "<<ptclePortion<<endl;
  logFile<<"globalGaussPtsNum = "<<globalGaussPtsNum
	 <<"; globalBGaussPtsNum = "<<globalBGaussPtsNum<<endl;
#endif

  int idx;
  dbVector gCoords(usedDims);

  // Loop over particle portion.
  for(int i=0;i<exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];

    idx = 0;

    // Loop over all volume and boundary Gauss points.
    for(int j=0;j<globalGaussPtsNum+globalBGaussPtsNum;j++) {

      for(int k=0;k<usedDims;k++)
	gCoords[k] = allGCoords[j*usedDims+k];

      // Check if current point is supported.
      if(particles[ptcle].querySupported(InputData,gCoords,modelData,logFile)) {

	if(idxSize < oldIdxSize) {
	  iPoints[idxSize] = j;
	  idxSize++;
	  idx++;
	}
	else {
	  logFile <<"Too less data buffer has been allocated for vector "
		  <<"iPoints in \n function BackgroundMesh::"
		  <<"setPtcleGaussConn(existing: "<<oldIdxSize<<" - "
		  <<"needed: "<<idxSize+1<<")!"<< endl;
	  iPoints.push_back(j);
	  idxSize++;
	  idx++;
	}

      }

    }

    iCounts[i] = idx;
  }

  /*********************************************************************/
  // Create the vectors inclPts, inclStartPos and inclPtsCounts using
  // the local ones.
  int rank;
  int size;

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // The inclPtsCounts' part.
  int* recvCounts = new int[size];
  int* displs = new int[size];

  int sendBuf = ptclePortion;
  MPI_Allgather(&sendBuf,1,MPI_INT,recvCounts,1,MPI_INT,MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  inclPtsCounts = intVector(particlesNum);

  MPI_Allgatherv(&iCounts[0],recvCounts[rank],MPI_INT,
		 &inclPtsCounts[0],recvCounts,displs,MPI_INT,
		 MPI_COMM_WORLD);

  // The inclPtsStartPos' part.

  // Loop over all particles.
  inclPtsStartPos = intVector(particlesNum);

  for(int j=0;j<particlesNum;j++) {

    if(j != 0)
      inclPtsStartPos[j] = inclPtsStartPos[j-1] + inclPtsCounts[j-1];
    else
      inclPtsStartPos[j] = 0;
  }

  // The inclPoints' part.
  sendBuf = idxSize;

  MPI_Allgather(&sendBuf,1,MPI_INT,recvCounts,1,MPI_INT,MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  int globalIdxSize;
  MPI_Allreduce(&idxSize,&globalIdxSize,1,MPI_INT,MPI_SUM,
		MPI_COMM_WORLD);
  inclPts = intVector(globalIdxSize);

  MPI_Allgatherv(&iPoints[0],recvCounts[rank],MPI_INT,&inclPts[0],
		 recvCounts,displs,MPI_INT,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"******************** included points **********************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"PARTICLE "<<i<<") ";
    for(int j=inclPtsStartPos[i];j<inclPtsStartPos[i]+inclPtsCounts[i];j++) {
      logFile<<" "<<inclPts[j];
    }
    logFile<<endl;
  }
#endif

  delete[] displs,recvCounts;
  //  iCounts.~intVector();
  //  iPoints.~intVector();

}

/**********************************************************************/
/**********************************************************************/
// Determine for each Gauss point its supporting particles.
void BackgroundMesh::setGaussPtcleConn(InputFileData* InputData,
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

#ifdef _supportDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"********** calculation of gauss-ptcle connlist **********"<<endl;
  intMatrix allSupportLists(localGaussPtsNum);
  double oldTime = MPI_Wtime();
#endif

  dbVector radius(3);
  int supportCounts,localMinSupport;

  int localMaxSupport = maxPtcleSupport;


  // straight forward support list computation
  if(mode == 1) {

    // Loop over all local gauss points.
    for (int i=0;i<localGaussPtsNum;i++) {

      dbVector& coords = gaussPoints[i].getCoords();
      intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
      intVector& suppGhostPtcls =
	gaussPoints[i].getSupportBoundGhostPtcls();

#ifdef _supportDebugMode_
      logFile<<i<<".) GAUSS POINT "<<gaussPoints[i].getGlobalID()<<endl;
#endif

      setPointPtclsConn(InputData,coords,supportCounts,suppPtcls,modelData,
			logFile);

      // -----------------------------------------------------------------

      if(supportCounts + suppGhostPtcls.size() <
	 InputData->getValue("gaussParticleConnectivity")) {

	logFile<<"In BackgroundMesh::setGaussPtcleConn, too less particles\n"
	       <<"contain local gauss point "<<i
	       <<" at processor \n"<<rank<<"( necessary: "
	       <<gaussPtcleConnect<<" existing: "<<supportCounts
	       <<")!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(supportCounts + suppGhostPtcls.size() > localMaxSupport)
	localMaxSupport = supportCounts;


      if(i == 0)
	localMinSupport = supportCounts;

      else if(supportCounts  + suppGhostPtcls.size() < localMinSupport)
	localMinSupport = supportCounts;

    }

  }

  /**********************************************************************/
  // element dependent support list computation
  else if(mode == 2) {

    int betaPtcle;

    printVector(globalLocalElemIdx,"globalLocalElemIdx",logFile);
    printVector(elementRootList,"elementRootList",logFile);
    logFile << "localGaussPtsNum: " << localGaussPtsNum << endl;

    blVector includedPtcls(particlesNum);
    blVector includedElems(globalLocalElemIdx.size());
    intVector allSupportSize(localGaussPtsNum);
    intMatrix neededElemCounter(localGaussPtsNum,intVector(size));
    intMatrix3 neededElemIdx(localGaussPtsNum,intMatrix(size));

    if(elementRootList.size() < globalLocalElemIdx.size()) {

      intVector localElemRootList(globalLocalElemIdx.size());
      elementRootList.resize(globalLocalElemIdx.size());

      for(int i=0;i<globalLocalElemIdx.size();i++)

	if(globalLocalElemIdx[i] >= 0)
	  localElemRootList[i] = rank;

      MPI_Allreduce(&localElemRootList[0],&elementRootList[0],
		    globalLocalElemIdx.size(),MPI_INT,MPI_SUM,
		    MPI_COMM_WORLD);

      localElemRootList.resize(0);

#ifdef _supportDebugMode_
      logFile<<"***************** element root list ***************"<<endl;
      for(int i=0;i<globalLocalElemIdx.size();i++)
	logFile<<"Element "<<i<<": "<<elementRootList[i]<<endl;
#endif

    }

    // Loop over all local gauss points.
    for (int i=0;i<localGaussPtsNum;i++) {

      dbVector& coords = gaussPoints[i].getCoords();
      intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
      intVector& suppGhostPtcls =
    		  	  	  	  gaussPoints[i].getSupportBoundGhostPtcls();

      suppPtcls.resize(localMaxSupport);

      intVector& elemInfo = gaussPoints[i].getElementInfo();
      intVector& betaNodes = nodesElements[elemInfo[0]].getNodes();
      includedElems[elemInfo[0]] = true;

#ifdef _supportDebugMode_
      logFile<<i<<".) GAUSS POINT "<<gaussPoints[i].getGlobalID()
	     <<": localElemID "<<elemInfo[0]<<endl;
#endif

      // Loop over the beta element's portion of particles.
      for(int j=0;j<betaNodes.size();j++) {
	betaPtcle = betaNodes[j]-1;

#ifdef _supportDebugMode_
	logFile<<"BETA-ptcle "<<betaPtcle<<endl;
#endif


	setElemPtcleSupport(InputData,betaPtcle,
			    includedPtcls,includedElems,
			    allSupportSize[i],suppPtcls,coords,
			    neededElemCounter[i],
			    neededElemIdx[i],modelData,logFile);


      }

      // clear arrays includedPtcls and includedElems
      for(int j=0;j<allSupportSize[i];j++)
	includedPtcls[suppPtcls[j]] = false;

      clearArray(includedElems);

      if(allSupportSize[i] > localMaxSupport) {
	localMaxSupport = allSupportSize[i];
	maxGaussSupport = allSupportSize[i];
      }

      // free not needed memory
      for(int j=0;j<size;j++)
	resizeArray(neededElemIdx[i][j],
		    neededElemCounter[i][j]);

    }


#ifdef _supportDebugMode_
    logFile<<"********* temp gauss-particle support list *********"<<endl;
    for(int i=0;i<localGaussPtsNum;i++) {
      intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
      logFile<<"GPoint "<<i<<"(supp "<<allSupportSize[i]<<"): ";
      for(int j=0;j<allSupportSize[i];j++)
	logFile<<suppPtcls[j]<<" ";
      logFile<<endl;
    }
    logFile<<"***************** nonlocal elements ****************"<<endl;
    for(int i=0;i<localGaussPtsNum;i++) {
      logFile<<"GPoint "<<i<<" - nonlocal elems: ";
      for(int j=0;j<size;j++) {
	logFile<<"proc "<<j<<": ";
	for(int k=0;k<neededElemCounter[i][j];k++)
	  logFile<<neededElemIdx[i][j][k]<<" ";
	logFile<<" | ";
      }
      logFile<<endl;
    }
    if(rank == 0)
      cout<<"gauss point support list computing - first round\n"
	  <<"finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
    int round = 2;
#endif

    // -------------------------------------------------------------------
    // incorporate the nonlocal nodes

    bool flag;
    int localNeededElems,globalNeededElems;
    intVector nonLocalNodes;
    intVector dummyVec(size);
    intMatrix dummyMat(size);

    int maxLocalGaussPtsNum,localGlobalDifference;
    MPI_Allreduce(&localGaussPtsNum,&maxLocalGaussPtsNum,1,
		  MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    localGlobalDifference = maxLocalGaussPtsNum - localGaussPtsNum;


    globalNeededElems = 1;

    // loop over all particles until all processes have their particle
    // support list assembled
    while(globalNeededElems != 0) {

      localNeededElems = 0;

      for(int i=0;i<maxLocalGaussPtsNum;i++) {

	if(i < localGaussPtsNum) {

#ifdef _supportDebugMode_
	  logFile<<"Gauss Point "<<i<<endl;
#endif

	  flag = getNonLocalSupport(InputData,neededElemCounter[i],
				    neededElemIdx[i],nonLocalNodes,
				    modelData,logFile);


	  // set array includedElems
	  for(int j=0;j<size;j++)

	    for(int k=0;k<neededElemCounter[i][j];k++)
	      includedElems[neededElemIdx[i][j][k]] = true;

	  clearArray(neededElemCounter[i]);

	  if(flag) {

	    dbVector& coords = gaussPoints[i].getCoords();
	    intVector& suppPtcls = gaussPoints[i].getSupportPtcls();

	    // set the currently supporting particles
	    for(int j=0;j<allSupportSize[i];j++)
	      includedPtcls[suppPtcls[j]] = true;

	    // Loop over the alpha particle's nonlocal nodes
	    for(int k=0;k<nonLocalNodes.size();k++) {
	      betaPtcle = nonLocalNodes[k]-1;

	      setElemPtcleSupport(InputData,betaPtcle,
				  includedPtcls,includedElems,
				  allSupportSize[i],suppPtcls,
				  coords,neededElemCounter[i],
				  neededElemIdx[i],modelData,logFile);

	    }


	    // clear arrays includedPtcls and includedElems
	    for(int j=0;j<allSupportSize[i];j++)
	      includedPtcls[suppPtcls[j]] = false;

	    clearArray(includedElems);

	  }

	  // ---------------
	  // loop over all processes and determine the total number of still
	  // needed non-local elements
	  for(int j=0;j<size;j++)

	    localNeededElems += neededElemCounter[i][j];

	  // free not needed memory
	  for(int j=0;j<size;j++)
	    resizeArray(neededElemIdx[i][j],
			neededElemCounter[i][j]);


	}

	// ----------------
	// local process has a smaller number of Gauss points than some of
	// the other processes
	else {

#ifdef _supportDebugMode_
	  logFile<<"Dummy Gauss Point "<<i<<endl;
#endif

	  flag = getNonLocalSupport(InputData,dummyVec,dummyMat,
				    nonLocalNodes,modelData,logFile);

	}

      }

      MPI_Allreduce(&localNeededElems,&globalNeededElems,1,
		    MPI_INT,MPI_SUM,MPI_COMM_WORLD);

#ifdef _supportDebugMode_
      if(rank == 0)
	cout<<"gauss point support list computing - "<<round<<". round "
	    <<"finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
      round++;
#endif

    }


    // ----------------------------------------------------------------
    // set the final size of Gauss point support lists
    for(int i=0;i<localGaussPtsNum;i++) {

      intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
      intVector& suppGhostPtcls =
	gaussPoints[i].getSupportBoundGhostPtcls();

      resizeArray(suppPtcls,allSupportSize[i]);
      sortIntVector(suppPtcls,0,allSupportSize[i]-1);

      if(allSupportSize[i] + suppGhostPtcls.size() <
	 InputData->getValue("gaussParticleConnectivity")) {

	logFile<<"In BackgroundMesh::setGaussPtcleConn, too less particles\n"
	       <<"contain local gauss point "<<i
	       <<" at processor \n"<<rank<<"( necessary: "
	       <<gaussPtcleConnect<<" existing: "
	       <<allSupportSize[i] + suppGhostPtcls.size()
	       <<")!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(allSupportSize[i] + suppGhostPtcls.size() > localMaxSupport)
	localMaxSupport = allSupportSize[i];


      if(i == 0)
	localMinSupport = allSupportSize[i];

      else if(allSupportSize[i]  + suppGhostPtcls.size() < localMinSupport)
	localMinSupport = allSupportSize[i];

    }

#ifdef _supportDebugMode_
    string compareMode("arbitrary-subvector");
    intMatrix allSupportLists(localGaussPtsNum);
    for(int i=0;i<localGaussPtsNum;i++) {
      intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
      allSupportLists[i] = suppPtcls;
      suppPtcls.resize(0);
    }
    InputData->setValue("supportComputationMode",1.0);
    setGaussPtcleConn(InputData,modelData,logFile);
    for(int i=0;i<localGaussPtsNum;i++) {
      intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
      if(!compareIntVecs(compareMode,allSupportLists[i],
			 suppPtcls,logFile)) {
	logFile<<"different:"<<endl;
	for(int j=0;j<allSupportLists[i].size();j++)
	  logFile<<allSupportLists[i][j]<<" ";
	logFile<<endl;
	for(int j=0;j<suppPtcls.size();j++)
	  logFile<<suppPtcls[j]<<" ";
	logFile<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    InputData->setValue("supportComputationMode",2.0);
#endif

  }
  else {

    logFile<<"In BackgroundMesh::setGaussPtcleConn Gauss point - "
	   <<"particle\n computation mode "<<mode<<" is not supported!"
	   <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);

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
#endif
#ifdef _supportDebugMode_
  if(rank == 0)
    cout<<"gauss finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Abort(MPI_COMM_WORLD,1);
}

/**********************************************************************/
/**********************************************************************/
// Determine for each boundary gauss point its supporting particles.
void BackgroundMesh::setBGaussPtcleConn(InputFileData* InputData,
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

  dbVector radius(3);
  int supportCounts,localMinSupport;

#ifdef _geometryDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"******* calculation of bound-gauss-ptcle connlist *******"<<endl;
  //intMatrix allSupportLists(localBGaussPtsNum);
  double oldTime = MPI_Wtime();
#elif defined _supportDebugMode_
  double oldTime = MPI_Wtime();
#endif

  int localMaxSupport = maxGaussSupport;


  // straight forward support list computation
  if(mode == 1) {

    // Loop over all boundary gauss points.
    for (int i=0;i<localBGaussPtsNum;i++) {

      dbVector& coords = boundGaussPoints[i].getCoords();
      intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();
      intVector& suppGhostPtcls =
	boundGaussPoints[i].getSupportBoundGhostPtcls();


      setPointPtclsConn(InputData,coords,supportCounts,suppPtcls,modelData,
			logFile);


      if(supportCounts + suppGhostPtcls.size() <
	 InputData->getValue("gaussParticleConnectivity")) {

	logFile<<"In BackgroundMesh::setBGaussPtcleConn too less\n" 
	       <<"particles contain boundary gauss point "<<i
	       <<" at processor \n"<<rank<<"( necessary: "
	       <<gaussPtcleConnect<<" existing: "<<supportCounts
	       <<")!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(supportCounts  + suppGhostPtcls.size() > localMaxSupport)
	localMaxSupport = supportCounts;


      if(i == 0)
	localMinSupport = supportCounts;

      else if(supportCounts  + suppGhostPtcls.size() < localMinSupport)
	localMinSupport = supportCounts;

    }

  }

  /**********************************************************************/
  // element dependent support list computation
  else if(mode == 2) {

    int betaPtcle,sourceIdx;

    blVector includedPtcls(particlesNum);
    blVector includedElems(globalLocalElemIdx.size());
    intVector allSupportSize(localBGaussPtsNum);
    intMatrix neededElemCounter(localBGaussPtsNum,intVector(size));
    intMatrix3 neededElemIdx(localBGaussPtsNum,intMatrix(size));

    if(elementRootList.size() < globalLocalElemIdx.size()) {

      intVector localElemRootList(globalLocalElemIdx.size());
      elementRootList.resize(globalLocalElemIdx.size());

      for(int i=0;i<globalLocalElemIdx.size();i++)

	if(globalLocalElemIdx[i] >= 0)
	  localElemRootList[i] = rank;

      MPI_Allreduce(&localElemRootList[0],&elementRootList[0],
		    globalLocalElemIdx.size(),MPI_INT,MPI_SUM,
		    MPI_COMM_WORLD);

      localElemRootList.resize(0);

#ifdef _supportDebugMode_
      logFile<<"***************** element root list ***************"<<endl;
      for(int i=0;i<globalLocalElemIdx.size();i++)
	logFile<<"Element "<<i<<": "<<elementRootList[i]<<endl;
#endif

    }

    // Loop over all local boundary Gauss points.
    for (int i=0;i<localBGaussPtsNum;i++) {

      dbVector& coords = boundGaussPoints[i].getCoords();
      intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();
      intVector& suppGhostPtcls =
	boundGaussPoints[i].getSupportBoundGhostPtcls();

      suppPtcls.resize(localMaxSupport);

      intVector& elemInfo = boundGaussPoints[i].getElementInfo();
      int& globalID = elemInfo[2];
      int& localID = globalLocalElemIdx[globalID];

#ifdef _supportDebugMode_
      logFile<<i<<".) BOUND-GAUSS POINT "<<boundGaussPoints[i].getGlobalID()
	     <<": globalElemID "<<globalID<<" localElemID "<<localID<<endl;
#endif

      // check if current element is a local one
      if(localID >= 0 ) {

	intVector& betaNodes = nodesElements[localID].getNodes();
	includedElems[globalID] = true;

	// Loop over the beta element's portion of particles.
	for(int j=0;j<betaNodes.size();j++) {
	  betaPtcle = betaNodes[j]-1;

#ifdef _supportDebugMode_
	  logFile<<"BETA-ptcle "<<betaPtcle<<endl;
#endif


	  setElemPtcleSupport(InputData,betaPtcle,
			      includedPtcls,includedElems,
			      allSupportSize[i],suppPtcls,coords,
			      neededElemCounter[i],
			      neededElemIdx[i],modelData,logFile);


	}

	// clear arrays includedPtcls and includedElems
	for(int j=0;j<allSupportSize[i];j++)
	  includedPtcls[suppPtcls[j]] = false;

	clearArray(includedElems);

	if(allSupportSize[i] > localMaxSupport) {
	  localMaxSupport = allSupportSize[i];
	  maxGaussSupport = allSupportSize[i];
	}

	// free not needed memory
	for(int j=0;j<size;j++)
	  resizeArray(neededElemIdx[i][j],
		      neededElemCounter[i][j]);

      }


      // -----------------
      // element non-local
      else  {

	sourceIdx = elementRootList[globalID];

	neededElemIdx[i][sourceIdx].push_back(globalID);
	neededElemCounter[i][sourceIdx]++;

#ifdef _supportDebugMode_
	for(int j=0;j<size;j++)
	  for(int k=0;k<neededElemCounter[i][j];k++)
	    logFile<<neededElemIdx[i][j][k]<<" ";
	logFile<<endl;
#endif

      }

    }


#ifdef _supportDebugMode_
    logFile<<"********* temp gauss-particle support list *********"<<endl;
    for(int i=0;i<localBGaussPtsNum;i++) {
      intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();
      logFile<<"GPoint "<<i<<"(supp "<<allSupportSize[i]<<"): ";
      for(int j=0;j<allSupportSize[i];j++)
	logFile<<suppPtcls[j]<<" ";
      logFile<<endl;
    }
    logFile<<"***************** nonlocal elements ****************"<<endl;
    for(int i=0;i<localBGaussPtsNum;i++) {
      logFile<<"GPoint "<<i<<" - nonlocal elems: ";
      for(int j=0;j<size;j++) {
	logFile<<"proc "<<j<<": ";
	for(int k=0;k<neededElemCounter[i][j];k++)
	  logFile<<neededElemIdx[i][j][k]<<" ";
	logFile<<" | ";
      }
      logFile<<endl;
    }
    if(rank == 0)
      cout<<"gauss point support list computing - first round\n"
	  <<"finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
    int round = 2;
#endif

    // -------------------------------------------------------------------
    // incorporate the nonlocal nodes

    bool flag;
    int localNeededElems,globalNeededElems;
    intVector nonLocalNodes;
    intVector dummyVec(size);
    intMatrix dummyMat(size);

    int maxLocalGaussPtsNum,localGlobalDifference;
    MPI_Allreduce(&localBGaussPtsNum,&maxLocalGaussPtsNum,1,
		  MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    localGlobalDifference = maxLocalGaussPtsNum - localBGaussPtsNum;


    globalNeededElems = 1;

    // loop over all particles until all processes have their particle
    // support list assembled
    while(globalNeededElems != 0) {

      localNeededElems = 0;

      for(int i=0;i<maxLocalGaussPtsNum;i++) {

	if(i < localBGaussPtsNum) {

#ifdef _supportDebugMode_
	  logFile<<"Bound-Gauss Point "<<i<<endl;
#endif

	  flag = getNonLocalSupport(InputData,neededElemCounter[i],
				    neededElemIdx[i],nonLocalNodes,
				    modelData,logFile);


	  // set array includedElems
	  for(int j=0;j<size;j++)

	    for(int k=0;k<neededElemCounter[i][j];k++)
	      includedElems[neededElemIdx[i][j][k]] = true;

	  clearArray(neededElemCounter[i]);

	  if(flag) {

	    dbVector& coords = boundGaussPoints[i].getCoords();
	    intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();

	    // set the currently supporting particles
	    for(int j=0;j<allSupportSize[i];j++)
	      includedPtcls[suppPtcls[j]] = true;

	    // Loop over the alpha particle's nonlocal nodes
	    for(int k=0;k<nonLocalNodes.size();k++) {
	      betaPtcle = nonLocalNodes[k]-1;

	      setElemPtcleSupport(InputData,betaPtcle,
				  includedPtcls,includedElems,
				  allSupportSize[i],suppPtcls,
				  coords,neededElemCounter[i],
				  neededElemIdx[i],modelData,logFile);

	    }


	    // clear arrays includedPtcls and includedElems
	    for(int j=0;j<allSupportSize[i];j++)
	      includedPtcls[suppPtcls[j]] = false;

	    clearArray(includedElems);

	  }

	  // ---------------
	  // loop over all processes and determine the total number of still
	  // needed non-local elements
	  for(int j=0;j<size;j++)

	    localNeededElems += neededElemCounter[i][j];


	  // free not needed memory
	  for(int j=0;j<size;j++)
	    resizeArray(neededElemIdx[i][j],
			neededElemCounter[i][j]);

	}

	// ----------------
	// local process has a smaller number of Gauss points than some of
	else {

#ifdef _supportDebugMode_
	  logFile<<"Dummy Gauss Point "<<i<<endl;
#endif

	  flag = getNonLocalSupport(InputData,dummyVec,dummyMat,
				    nonLocalNodes,modelData,logFile);

	}

      }


      MPI_Allreduce(&localNeededElems,&globalNeededElems,1,
		    MPI_INT,MPI_SUM,MPI_COMM_WORLD);

#ifdef _supportDebugMode_
      if(rank == 0)
	cout<<"bound-gauss support list computing - "<<round<<". round "
	    <<"finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
      round++;
#endif

    }


    // ----------------------------------------------------------------
    // set the final size of Gauss point support lists
    for(int i=0;i<localBGaussPtsNum;i++) {

      intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();
      intVector& suppGhostPtcls =
	boundGaussPoints[i].getSupportBoundGhostPtcls();

      resizeArray(suppPtcls,allSupportSize[i]);
      sortIntVector(suppPtcls,0,allSupportSize[i]-1);

      if(allSupportSize[i] + suppGhostPtcls.size() <
	 InputData->getValue("gaussParticleConnectivity")) {

	logFile<<"In BackgroundMesh::setBGaussPtcleConn, too less particles\n"
	       <<"contain local gauss point "<<i
	       <<" at processor \n"<<rank<<"( necessary: "
	       <<gaussPtcleConnect<<" existing: "<<supportCounts
	       <<")!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(allSupportSize[i] + suppGhostPtcls.size() > localMaxSupport)
	localMaxSupport = allSupportSize[i];


      if(i == 0)
	localMinSupport = minGaussSupport;

      else if(allSupportSize[i]  + suppGhostPtcls.size() < localMinSupport)
	localMinSupport = allSupportSize[i];

    }

#ifdef _supportDebugMode_
    string compareMode("arbitrary-subvector");
    intMatrix allSupportLists(localBGaussPtsNum);
    for(int i=0;i<localBGaussPtsNum;i++) {
      intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();
      allSupportLists[i] = suppPtcls;
      suppPtcls.resize(0);
    }
    InputData->setValue("supportComputationMode",1.0);
    setBGaussPtcleConn(InputData,modelData,logFile);
    for(int i=0;i<localBGaussPtsNum;i++) {
      intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();
      if(!compareIntVecs(compareMode,allSupportLists[i],
			 suppPtcls,logFile)) {
	logFile<<"different:"<<endl;
	for(int j=0;j<allSupportLists[i].size();j++)
	  logFile<<allSupportLists[i][j]<<" ";
	logFile<<endl;
	for(int j=0;j<suppPtcls.size();j++)
	  logFile<<suppPtcls[j]<<" ";
	logFile<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    InputData->setValue("supportComputationMode",2.0);
#endif

  }
  else {

    logFile<<"In BackgroundMesh::setBGaussPtcleConn Gauss point - "
	   <<"particle\n computation mode "<<mode<<" is not supported!"
	   <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);

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
#endif
#ifdef _supportDebugMode_
  if(rank == 0)
    cout<<"bound gauss finished in "<<MPI_Wtime()-oldTime<<" secs"<<endl;
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Abort(MPI_COMM_WORLD,1);
}

/************************************************************************/
/************************************************************************/
// Determine the locally needed particles.
void BackgroundMesh::setLocalPtcls(InputFileData* InputData,
                                   std::map<std::string,double>& modelData,
                                   std::ofstream& logFile,
                                   PetscViewer& viewerMPI) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int usedDims = (int)modelData["usedDimensions"];
  int integrationMethod = (int)modelData["integrationMethod"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  // Determine which rows have to be stored local respectively which rows
  // get entries during the local integration (ghost entries!).


  intVector localPtcls(particlesNum);

  // Gauss quadrature
  if(integrationMethod == 1) {

    // Loop over all local volume gauss points.
    for(int i=0;i<gaussPoints.size();i++) {
      intVector& sPtcls = gaussPoints[i].getSupportPtcls();

      // Loop over all supporting particles of current gauss point.
      for(int j=0;j<sPtcls.size();j++) {
	localPtcls[sPtcls[j]] = 1;
      }

    }


    // Loop over all boundary gauss points.
    for(int i=0;i<boundGaussPoints.size();i++) {
      intVector& sPtcls = boundGaussPoints[i].getSupportPtcls();

      // Loop over all supporting particles of current gauss point.
      for(int j=0;j<sPtcls.size();j++) {
	localPtcls[sPtcls[j]] = 1;
      }

    }

  }

  // Particle integration.
  else if(integrationMethod == 2) {

    // Loop over the local portion of particle integration points.
    for(int i=0;i<ptcleRootList.size();i++) {

      if(ptcleRootList[i] == rank) {

	intVector& sPtcls = particles[i].getSupportPtcls();

	// Loop over all supporting particles of current particle point.
	for(int j=0;j<sPtcls.size();j++)
	  localPtcls[sPtcls[j]] = 1;

      }

    }

  }

#ifdef _geometryDebugMode_
  logFile<<"*********** preliminary local Particles  *************"<<endl;
  for(int i=0;i<localPtcls.size();i++)
    logFile<<"Ptcle "<<i<<": "<<localPtcls[i]<<endl;
#endif

  /*********************************************************************/
  // Mapping local particles --> global ones
  // and global particles -> local ones.
  int m=0;
  localGlobalPtcls = intVector((int)particlesNum/size);
  globalLocalPtcleIdx.resize(particlesNum);
  globalLocalPtcleIdx.assign(globalLocalPtcleIdx.size(),-1);

  // Loop over all particles.
  for(int i=0;i<particlesNum;i++) {

    if(localPtcls[i] != 0) {

      if(m >= localGlobalPtcls.size()) {
	localGlobalPtcls.push_back(i);
	globalLocalPtcleIdx[i] = i;
	m++;
      }
      else {
	localGlobalPtcls[m] = i;
	globalLocalPtcleIdx[i] = i;
	m++;
      }

    }

  }

  if(localGlobalPtcls.size() > m)
    localGlobalPtcls.resize(m);


#ifdef _geometryDebugMode_
  logFile<<"************** local->global ptcle indices **********"<<endl;
  for(int i=0;i<localGlobalPtcls.size();i++)
    logFile<<"local ptcle number: "<<i<<" global one: "
	   <<localGlobalPtcls[i]<<endl;
  logFile<<"************** global->local ptcle indices **********"<<endl;
  for(int i=0;i<globalLocalPtcleIdx.size();i++)
    logFile<<"global ptcle number: "<<i<<" local one: "
	   <<globalLocalPtcleIdx[i]<<endl;
#endif


  /*********************************************************************/
  // determine an EXCLUSIVELY local portion of particles

  intMatrix allLocalPtcls(size,intVector(particlesNum));

  allLocalPtcls[rank] = localPtcls;
  resizeArray(localPtcls,0);

  for(int i=0;i<size;i++)

    MPI_Bcast(&allLocalPtcls[i][0],particlesNum,MPI_INT,i,MPI_COMM_WORLD);


#ifdef _geometryDebugMode_
  logFile<<"*************** all local particles  ****************"<<endl;
  for(int i=0;i<allLocalPtcls[0].size();i++) {
    dbVector& coords = particles[i].getCoords();
    logFile<<"Ptcle "<<i<<": procs: ";
    for(int j=0;j<size;j++)
      if(allLocalPtcls[j][i] == 1)
	logFile<<j<<" ";
    logFile<<"coords: ";
    for(int j=0;j<coords.size();j++)
      logFile<<coords[j]<<" ";
    logFile<<endl;
  }
#endif

  int currentMinPortion,proc,procCount;
  double diff;
  intVector maxPtcle(size),minPtcle(size),localPortionCounts(size);
  intVector procsInNeed(size);

  if(exclusiveLocalPtcls.size() < particlesNum)
    exclusiveLocalPtcls.resize(particlesNum);

  for(int i=0;i<size;i++)
    minPtcle[i] = particlesNum;

  int defaultPortion = (int)ceil((double)particlesNum/size);

  // -------------------------------------------------------------------
  // loop over all particles and store those which are
  // required at only one processes

  // maintain a continuous sequence of particles on each process

  int algorithm = 1;

  if(algorithm == 1) {

    for(int i=0;i<particlesNum;i++) {

      proc = -1;
      m=0;

      // check how many processors need particle 'i'
      for(int j=0;j<size;j++) {

	if(allLocalPtcls[j][i] == 1) {
	  proc = j;
	  m++;
	}

      }

      // particle is needed by only one processor
      if(m == 1) {

	for(int j=0;j<size;j++)
	  allLocalPtcls[j][i] = 0;

	// store the current particle
	if(proc == rank )
	  exclusiveLocalPtcls[localPortionCounts[proc]] = i;

	// determine the range
	if(i < minPtcle[proc])
	  minPtcle[proc] = i;

	if(i > maxPtcle[proc])
	  maxPtcle[proc] = i;

	localPortionCounts[proc]++;

      }

    }

    // ---------------------------------------------------------------------
    // then fill up the default portion of particles for each processor

    for(int i=0;i<particlesNum;i++) {

      proc = -1;

      // check whether a processor which needs particle 'i'
      // has got its default portion yet
      for(int j=0;j<size;j++) {

	if(allLocalPtcls[j][i] == 1 &&
	   localPortionCounts[j] < defaultPortion) {
	  proc = j;
	  break;
	}

      }

      // processor 'proc' needs particle 'i' and hasn't got its default
      // portion yet
      if(proc != -1) {

	for(int j=0;j<size;j++)
	  allLocalPtcls[j][i] = 0;

	// store the current particle
	if(proc == rank)
	  exclusiveLocalPtcls[localPortionCounts[proc]] = i;

	// determine the range
	if(i < minPtcle[proc])
	  minPtcle[proc] = i;

	if(i > maxPtcle[proc])
	  maxPtcle[proc] = i;

	localPortionCounts[proc]++;

      }

    }

  }

  // strong emphasis to equalize portion but mixing of particles over
  // processes
  else {

    int procSupport = 1;

    while(procSupport <= size) {

      for(int i=0;i<particlesNum;i++) {

	procCount = 0;

	// check how many processors need particle 'i'
	for(int j=0;j<size;j++) {

	  if(allLocalPtcls[j][i] == 1) {
	    procsInNeed[procCount] = j;
	    procCount++;
	  }

	}

	// treat those particles 'i' which are needed only by 'procSupport'
	// processes
	if(procCount == procSupport) {

	  proc = -1;
	  currentMinPortion = 0;

	  // loop over processes which need current particle 'i'
	  for(int j=0;j<procCount;j++) {

	    if(localPortionCounts[procsInNeed[j]] < currentMinPortion ||
	       proc == -1) {

	      proc = procsInNeed[j];
	      currentMinPortion = localPortionCounts[procsInNeed[j]];
	    }

	  }


	  if(proc != -1) {

	    for(int j=0;j<size;j++)
	      allLocalPtcls[j][i] = 0;

	    // store the current particle
	    if(proc == rank)
	      exclusiveLocalPtcls[localPortionCounts[proc]] = i;

	    // determine the range
	    if(i < minPtcle[proc])
	      minPtcle[proc] = i;

	    if(i > maxPtcle[proc])
	      maxPtcle[proc] = i;

	    localPortionCounts[proc]++;

	  }

	}

      }

      procSupport++;

    }

  }

#ifdef _geometryDebugMode_
  logFile<<"*** intermediate exclusive local particle portion  ***"<<endl;
  for(int j=0;j<size;j++)
    logFile<<"proc "<<j<<": minPtcle="<<minPtcle[j]
	   <<" maxPtcle="<<maxPtcle[j]<<endl;
  logFile<<"------------------------------------------------------"<<endl;
  for(int i=0;i<localPortionCounts[rank];i++)
    logFile<<i<<".) "<<exclusiveLocalPtcls[i]<<endl;
  logFile<<"------------------------------------------------------"<<endl;
  logFile<<"---------------- remaining particles -----------------"<<endl;
  for(int i=0;i<particlesNum;i++) {
    for(int j=size-1;j>=0;j--) {
      if(allLocalPtcls[j][i] == 1)
	logFile<<"ptcle "<<i<<": rank="<<j<<endl;
    }
  }
#endif

  // ---------------------------------------------------------------------
  // loop over all remaining particles and store them where they fit in
  // the best

  int distributingMode =
    (int)InputData->getValue("parallelPtcleDistributingMode");

  intVector potentialHost(size);

  for(int i=0;i<particlesNum;i++) {

    proc = -1;
    diff = particlesNum;

    for(int j=0;j<size;j++)
      potentialHost[j] = particlesNum;

    // (1) location optimized (default) (2) computing load balanced
    if(distributingMode == 1) {

      // loop over all processes and determine which of those gets the
      // current particle
      for(int j=size-1;j>=0;j--) {

#ifdef _geometryDebugMode_
	if(allLocalPtcls[j][i] == 1)
	  logFile<<"ptcle "<<i<<": rank="
		 <<j<<" centre="<<(minPtcle[j] + maxPtcle[j])/2.0
		 <<" diff="<<fabs(i - (minPtcle[j] + maxPtcle[j])/2.0)<<endl;
#endif

	// check whether particle 'i' needs to be locally stored and
	// determine its distance to the centre of processor's 'j'
	// range of already stored particles
	if(allLocalPtcls[j][i] == 1 &&
	   fabs(i - (minPtcle[j] + maxPtcle[j])/2.0) < diff) {

	  proc = j;
	  diff = fabs(i - (minPtcle[j] + maxPtcle[j])/2.0);

	}

      }


      // particle 'i' needs to be stored
      if(proc > -1) {

	for(int j=0;j<size;j++)
	  allLocalPtcls[j][i] = 0;

	if(proc == rank)
	  exclusiveLocalPtcls[localPortionCounts[proc]] = i;

	// determine the range
	if(i < minPtcle[proc])
	  minPtcle[proc] = i;

	if(i > maxPtcle[proc])
	  maxPtcle[proc] = i;

	localPortionCounts[proc]++;

#ifdef _geometryDebugMode_
	logFile<<"=> assigned to proc "<<proc<<endl;
#endif

      }

    }

    // distribution more computing load balanced
    else {

      // loop over all processes and determine which of those gets the
      // current particle
      for(int j=size-1;j>=0;j--) {

#ifdef _geometryDebugMode_
	if(allLocalPtcls[j][i] == 1)
	  logFile<<"ptcle "<<i<<": rank="
		 <<j<<" centre="<<(minPtcle[j] + maxPtcle[j])/2.0
		 <<" diff="<<fabs(i - (minPtcle[j] + maxPtcle[j])/2.0)<<endl;
#endif

	// check whether particle 'i' needs to be locally stored
	if(allLocalPtcls[j][i] == 1) {

	  potentialHost[j] = localPortionCounts[j];
	  proc = size; // not the actual processor!
	}

      }


      // particle 'i' needs to be stored
      if(proc > -1) {

	proc = findMinValue(potentialHost);

	for(int j=0;j<size;j++)
	  allLocalPtcls[j][i] = 0;

	if(proc == rank)
	  exclusiveLocalPtcls[localPortionCounts[proc]] = i;

	// determine the range
	if(i < minPtcle[proc])
	  minPtcle[proc] = i;

	if(i > maxPtcle[proc])
	  maxPtcle[proc] = i;

	localPortionCounts[proc]++;

#ifdef _geometryDebugMode_
	logFile<<"=> assigned to proc "<<proc<<endl;
#endif

      }

    }

  }

  resizeArray(allLocalPtcls,0);
  resizeArray(exclusiveLocalPtcls,localPortionCounts[rank]);

  sortIntVector(exclusiveLocalPtcls,0,exclusiveLocalPtcls.size()-1);

  if(exclusiveLocalPtcls.size() == 0) {
    logFile<<"In function ParticleDistribution::setLocalPtcls particle split\n"
	   <<"amongst processes gives zero exclusively local particles.\n"
	   <<"Choose less processes."<<endl;
    if(rank == 0)
      cerr<<"In function ParticleDistribution::setLocalPtcls particle split\n"
	  <<"amongst processes gives zero exclusively local particles\n"
	  <<"at process "<<rank<<". Choose less processes."<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }



#ifdef _geometryDebugMode_
  logFile<<"******* final exclusive local particle portion  *******"<<endl;
  for(int j=0;j<size;j++)
    logFile<<"proc "<<j<<": minPtcle="<<minPtcle[j]
	   <<" maxPtcle="<<maxPtcle[j]<<endl;
  logFile<<"------------------------------------------------------"<<endl;
  for(int i=0;i<exclusiveLocalPtcls.size();i++)
    logFile<<i<<".) "<<exclusiveLocalPtcls[i]<<endl;
#endif

  // set consecutive ordering of the exclusive local particles from all
  // processors
  setAllExclLocalPtcleOrder(InputData,logFile);


  /**********************************************************************/
  // set the particle root list
  intVector localRootList(particlesNum);

  for(int i=0;i<exclusiveLocalPtcls.size();i++)
    localRootList[exclusiveLocalPtcls[i]] = rank;

  clearArray(ptcleRootList);
  MPI_Allreduce(&localRootList[0],&ptcleRootList[0],particlesNum,
		MPI_INT,MPI_SUM,MPI_COMM_WORLD);



  int maxPtcleNum,minPtcleNum,ptcleSum;
  int localPtcleNum = exclusiveLocalPtcls.size();

  MPI_Allreduce(&localPtcleNum,&maxPtcleNum,1,MPI_INT,MPI_MAX,
		MPI_COMM_WORLD);
  MPI_Allreduce(&localPtcleNum,&minPtcleNum,1,MPI_INT,MPI_MIN,
		MPI_COMM_WORLD);

  MPI_Allreduce(&localPtcleNum,&ptcleSum,1,MPI_INT,MPI_SUM,
		MPI_COMM_WORLD);

  if(ptcleSum != particlesNum) {
    logFile<<"In function BackgroundMesh::setLocalPtcls incomplete "
	   <<"particle split."<<endl;
    if(rank == 0)
      cerr<<"In function BackgroundMesh::setLocalPtcls incomplete "
	  <<"particle split."<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  if(rank == 0)
    cout<<"particle split: "<<minPtcleNum<<" - "<<maxPtcleNum<<endl;


#ifdef _geometryDebugMode_
  logFile<<"************* exclusive local Particles  *************"<<endl;
  for(int i=0;i<exclusiveLocalPtcls.size();i++) {
    logFile<<i<<".) PARTICLE "<<exclusiveLocalPtcls[i]
	   <<"; member of localGlobalPtcls: ";
    vector<int>::iterator pos;
    pos = find(localGlobalPtcls.begin(),localGlobalPtcls.end(),
	       exclusiveLocalPtcls[i]);

    if(pos == localGlobalPtcls.end()) {
      logFile<<globalLocalPtcleIdx[exclusiveLocalPtcls[i]]<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else
      logFile<<globalLocalPtcleIdx[exclusiveLocalPtcls[i]]<<" ok"<<endl;
  }
  logFile<<"************* particle root list  *******************"<<endl;
  for(int i=0;i<ptcleRootList.size();i++)
    logFile<<"Ptcle "<<i<<": "<<ptcleRootList[i]<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Rearrange Gauss point support lists according to the new particles
// vector ordering
void BackgroundMesh::rearrangeSupportLists(InputFileData* InputData,
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

  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& vec = particles[ptcle].getInflSpheres();
    intVector tmpVec = vec;

    for(int j=0;j<vec.size();j++)
      vec[j] = newGlobalPtcls[tmpVec[j]];

    sortIntVector(vec,0,vec.size()-1);
  }

  for(int i=0;i<gaussPoints.size();i++) {
    GaussPoint& gPoint = gaussPoints[i];

    intVector& vec = gPoint.getSupportPtcls();
    intVector tmpVec = vec;

    for(int j=0;j<vec.size();j++)
      vec[j] = newGlobalPtcls[tmpVec[j]];

    sortIntVector(vec,0,vec.size()-1);

  }

  for(int i=0;i<boundGaussPoints.size();i++) {
    GaussPoint& gPoint = boundGaussPoints[i];

    intVector& vec = gPoint.getSupportPtcls();
    intVector tmpVec = vec;

    for(int j=0;j<vec.size();j++)
      vec[j] = newGlobalPtcls[tmpVec[j]];

    sortIntVector(vec,0,vec.size()-1);

  }

#ifdef _geometryDebugMode_

  logFile<<"*****************************************************"<<endl;
  logFile<<"********** particle - particle support list **********"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    logFile<<"Ptcle "<<i<<": ";
    for(int j=0;j<suppPtcls.size();j++)
      logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
  logFile<<"****************** influencing spheres **************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& inflParticles = particles[i].getInflSpheres();
    logFile<<"PARTICLE "<<i<<": ";
    for(int j=0;j<inflParticles.size();j++)
      logFile<<inflParticles[j]<<" ";
    logFile<<endl;
  }
  logFile<<"*********** gauss - particle support list ************"<<endl;
  for(int i=0;i<gaussPoints.size();i++) {
    intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
    logFile<<"GPoint "<<gaussPoints[i].getGlobalID()<<": ";
    for(int j=0;j<suppPtcls.size();j++)
      logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
  logFile<<"********* bound gauss - particle support list ********"<<endl;
  for(int i=0;i<boundGaussPoints.size();i++) {
    intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();
    logFile<<"bGPoint "<<boundGaussPoints[i].getGlobalID()<<": ";
    for(int j=0;j<suppPtcls.size();j++)
      logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Create a vector containing all global gauss coordinates.
dbVector BackgroundMesh::getAllGaussCoords() {

  using namespace std;

  dbVector globalCoords(globalGaussPtsNum*3);
  dbVector tmpGlobalCoords(globalGaussPtsNum*3);
  dbVector localCoords(localGaussPtsNum*3);
  intVector idx(localGaussPtsNum);
  intVector allIdx(globalGaussPtsNum);

  for(int i=0;i<localGaussPtsNum;i++) {
    localCoords[i*3] = gaussPoints[i].getCoord(0);
    localCoords[i*3+1] = gaussPoints[i].getCoord(1);
    localCoords[i*3+2] = gaussPoints[i].getCoord(2);
    idx[i] = gaussPoints[i].getGlobalID();
  }

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  intVector recvCounts(size);
  intVector displs(size);

  // The coordinates' part.
  int sendBuf = localGaussPtsNum*3;

  MPI_Allgather(&sendBuf,1,MPI_INT,&recvCounts[0],1,MPI_INT,MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  MPI_Allgatherv(&localCoords[0],recvCounts[rank],MPI_DOUBLE,
		 &tmpGlobalCoords[0],&recvCounts[0],&displs[0],MPI_DOUBLE,
		 MPI_COMM_WORLD);

  // The index part.
  MPI_Allgather(&localGaussPtsNum,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  MPI_Allgatherv(&idx[0],recvCounts[rank],MPI_INT,&allIdx[0],&recvCounts[0],
		 &displs[0],MPI_INT,MPI_COMM_WORLD);

  // Order the coordinates according to their global indices.
  for(int i=0;i<globalGaussPtsNum;i++) {
    globalCoords[allIdx[i]*3] = tmpGlobalCoords[i*3];
    globalCoords[allIdx[i]*3+1] = tmpGlobalCoords[i*3+1];
    globalCoords[allIdx[i]*3+2] = tmpGlobalCoords[i*3+2];
  }


  return globalCoords;
}

/**********************************************************************/
/**********************************************************************/
// Distribute the degrees of freedoms to the local gauss points.
void BackgroundMesh::assignGaussPtsDOF(InputFileData* InputData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  // Loop over all local gauss points to distribute the degrees of
  // freedom to them.
  int supportSize;

  for(int i=0;i<localGaussPtsNum;i++) {

    // Set the size of vector globalDOF of each gauss point.
    gaussPoints[i].setGlobalDOFSize(usedDOF);

    supportSize = gaussPoints[i].getSupportCounts();

    // Loop over all supported particles of each gauss point and set
    // its vector globalDOF.
    for(int j=0;j<supportSize;j++) {

      for(int k=0;k<usedDOF;k++) {
	int dof = gaussPoints[i].getSupportPtcle(j)*usedDOF+k;
	gaussPoints[i].setGlobalDOF(j*usedDOF+k,dof);
      }

    }

  }

#ifdef _geometryDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"GAUSS POINTS DOF"<<endl;
  int lowestDOF,highestDOF;
  lowestDOF = particlesNum - 1;
  highestDOF = 0;
  for(int i=0;i<localGaussPtsNum;i++) {
    intVector& globalDof = gaussPoints[i].getAllGlobalDOF();
    intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
    logFile<<i<<".) GAUSSPOINT "<<gaussPoints[i].getGlobalID()
	   <<":size "<<globalDof.size()<<" ";
    for(int j=0,k=0;j<globalDof.size();j+=usedDOF,k++) {
      if(globalDof[j] < lowestDOF)
	lowestDOF = globalDof[j];
      if(globalDof[j] > highestDOF)
	highestDOF = globalDof[j];
      logFile<<suppPtcls[k]<<": ";
      for(int k=0;k<usedDOF;k++)
	logFile<<globalDof[j+k]<<" ";
    }
    logFile<<endl;
  }
  logFile<<"***************** DOF range ********************"<<endl;
  logFile<<lowestDOF<<" - "<<highestDOF<<endl;
#endif

  /**********************************************************************/
  // Loop over all local boundary gauss points to distribute the degrees
  // of freedom to them.
  for(int i=0;i<localBGaussPtsNum;i++) {
    // Set the size of vector globalDOF of each gauss point.
    boundGaussPoints[i].setGlobalDOFSize(usedDOF);

    supportSize = boundGaussPoints[i].getSupportCounts();

    // Loop over all supported particles of each gauss point and set
    // its vector globalDOF.
    for(int j=0;j<supportSize;j++) {

      for(int k=0;k<usedDOF;k++) {
	int dof = boundGaussPoints[i].getSupportPtcle(j)*usedDOF+k;
	boundGaussPoints[i].setGlobalDOF(j*usedDOF+k,dof);
      }

    }

  }

#ifdef _geometryDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"BOUNDARY GAUSS POINTS DOF"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    intVector& globalDof = boundGaussPoints[i].getAllGlobalDOF();
    logFile<<i<<".) Bound GAUSSPOINT "<<boundGaussPoints[i].getGlobalID()
	   <<".) size "<<globalDof.size()<<" ";
    for(int j=0;j<globalDof.size();j++)
      logFile<<globalDof[j]<<" ";
    logFile<<endl;
  }
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
}

/************************************************************************/
/************************************************************************/
// Determine the Gauss point distribution among the processors.
void BackgroundMesh::setGaussDistribution(InputFileData* InputData,
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

    // erase the leading part of the boundary Gauss point's vector.
    int increment;
    vector<GaussPoint>::iterator startPos;
    vector<GaussPoint>::iterator endPos;

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
  //   if(localBGaussPtsNum == 0 && modelData["boundaryEnforcementMethod"] != 1) {
  //     logFile<<"In function BackgroundMesh::setGaussDistribution "
  // 	   <<"boundary Gauss point split\n"
  // 	   <<"amongst processes gives zero exclusively local volume "
  // 	   <<"Gauss points."<<endl
  // 	   <<"Choose less processes."<<endl;
  //     cerr<<"In function BackgroundMesh::setGaussDistribution "
  // 	<<"boundary Gauss point split\n"
  // 	<<"amongst processes gives zero exclusively local boundary \n"
  // 	<<"Gauss points at process "<<rank<<"."<<endl
  // 	<<"Choose less processes."<<endl;
  //     MPI_Abort(MPI_COMM_WORLD,1);
  //   }

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

  ///if(globalLocalIdx[tmpGaussPtsIdx[i][0]] != -1) {
    // int& ID = globalLocalIdx[tmpGaussPtsIdx[i][0]];
    //pointForceBoundPtcleIdx[m][0] = ID;

    //for(int j=1;j<tmpGaussPtsIdx[i].size();j++)
    //pointForceBoundPtcleIdx[m][j] = tmpGaussPtsIdx[i][j];

    //m++;
    //}

    //pointForceBoundPtcleIdx.resize(m);

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
    vector<FEMElement>(0,FEMElement(0)).swap(nodesElements);


    int maxGaussPtsNum,minGaussPtsNum,maxBGaussPtsNum,minBGaussPtsNum;

    MPI_Allreduce(&localGaussPtsNum,&maxGaussPtsNum,1,MPI_INT,MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&localGaussPtsNum,&minGaussPtsNum,1,MPI_INT,MPI_MIN,
                  MPI_COMM_WORLD);

    if(rank == 0)
      cout<<"volume Gauss point split: "<<minGaussPtsNum<<" - "
	  <<maxGaussPtsNum<<endl;

}

/************************************************************************/
/************************************************************************/
// Determine for each gauss point its supporting particles and set up
// a gauss point connectivity list used for partitioning.
void BackgroundMesh::setPartitionGraph(InputFileData* InputData,
                                       intVector& gaussProcList,
                                       dbVector& allGCoords,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int gaussPtcleConnect =
    (int)InputData->getValue("gaussParticleConnectivity");
  double multiplier = InputData->getValue("influenceRadiusMultiplier");

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"*********** partitioning graph set up ****************"<<endl;
  logFile<<"******************************************************"<<endl;
#endif

  // Determine for each particle its included gauss points.
  intVector inclPtsStartPos,inclPtsCounts,inclPts;

  setPtcleGaussConn(InputData,allGCoords,inclPtsStartPos,inclPtsCounts,
		    inclPts,modelData,logFile);

  /*********************************************************************/
  // Set up the gauss point particle connectivity list.
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // Determine the initial boundary Gauss point distribution.
  int defaultPortion,startIdx,endIdx;
  getInitBGaussSplitting(InputData,defaultPortion,startIdx,endIdx,
			 modelData,logFile);

  // Calculate the reserved buffer space for the unreduced graph list.
  int gaussPtsPerElem = getMaxGaussPtsPerVolumeElem();

  int idxSize = (int)((localGaussPtsNum+localBGaussPtsNum)
		      *maxPtcleSupport*maxGaussSupport*gaussPtsPerElem
		      /pow(maxPtcleSupport+maxGaussSupport,1.0/2.0));

#ifdef _geometryDebugMode_
  logFile<<"*********************************************************"<<endl;
  logFile<<"******************** CONN LIST set up *******************"<<endl;
  logFile<<"localBGaussPtsNum = "<<localBGaussPtsNum<<endl;
  logFile<<"startIdx: "<<startIdx<<"; endIdx: "<<endIdx<<endl;
  logFile<<"connlist size: "<<idxSize<<endl;
#endif

  bool indicator;
  int idx,particle;
  intVector cCounts(localGaussPtsNum+localBGaussPtsNum);
  intVector tmpConnList(idxSize);

  int oldIdxSize = idxSize;
  idxSize = 0;

  // Loop over all local volume Gauss points.
  for(int i=0;i<localGaussPtsNum;i++) {

#ifdef _geometryDebugMode_
    logFile<<"****************"<<endl;
    logFile<<"GAUSS POINT "<<i<<": "<<endl;
#endif

    idx = 0;

    // Loop over its supported particles.
    for(int j=0;j<gaussPoints[i].getSupportCounts();j++) {

      particle = gaussPoints[i].getSupportPtcle(j);

#ifdef _geometryDebugMode_
      logFile<<"Ptcle "<<particle<<" incl start "
	     <<inclPtsStartPos[particle]<<" end "
	     <<inclPtsStartPos[particle]+inclPtsCounts[particle]<<endl;
#endif

      // Loop over all gauss points containing a single particle.
      for(int k=inclPtsStartPos[particle];k<inclPtsStartPos[particle]+
	    inclPtsCounts[particle];k++) {

#ifdef _geometryDebugMode_
	logFile<<"check "<<" "<<inclPts[k]<<"; idxSizes "<<idx<<" "
	       <<idxSize<<"; already set: "<<endl;
#endif

	if(gaussPoints[i].getGlobalID() == inclPts[k]) {
	  indicator = false;
#ifdef _geometryDebugMode_
	  logFile<<endl;
#endif
	}
	else {
	  indicator = true;

	  // Check, if current gauss point is already contained in list.
	  for(int l=idxSize-idx;l<idxSize;l++) {

#ifdef _geometryDebugMode_
	    logFile<<tmpConnList[l]<<" ";
#endif

	    if(tmpConnList[l] == inclPts[k]) {
	      indicator = false;

#ifdef _geometryDebugMode_
	      logFile<<endl;
#endif
	      break;
	    }

	  }

	}
	// Write current gauss point in the list.
	if(indicator) {

#ifdef _geometryDebugMode_
	  logFile<<" <- write this point"<<endl;
#endif

	  if(idxSize < oldIdxSize)
	    tmpConnList[idxSize] = inclPts[k];
	  else {
	    tmpConnList.push_back(inclPts[k]);

	    logFile <<"Too less data buffer has been allocated for vector "
		    <<"tmpConnList in \n function BackgroundMesh::"
		    <<"setPartitionGraph(existing: "<<oldIdxSize<<" - "
		    <<"needed: "<<idxSize+1<<")!"<< endl;
	  }

	  idxSize++;
	  idx++;
	}

      }

    }

    if(idx < 4) {
      logFile <<"Too less neighbour gauss points in partitioning graph "
	      <<"list for gauss point "<<gaussPoints[i].getGlobalID()
	      <<"(existing "<<idx<<" - needed 4 minimum)!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else
      cCounts[i] = idx;
  }

#ifdef _geometryDebugMode_
  int startPos=0;
  logFile<<"***************** local tmpConnList ******************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    logFile<<"VERTEX "<<i<<" -> "<<gaussPoints[i].getGlobalID()<<".) "
	   <<cCounts[i]<<": ";
    for(int j=startPos;j<startPos+cCounts[i];j++)
      logFile<<tmpConnList[j]<<" ";
    logFile<<endl;
    startPos += cCounts[i];
  }
  logFile<<"************************************"<<endl;
  for(int i=localGaussPtsNum,k=startIdx;i<localGaussPtsNum+
        localBGaussPtsNum;i++,k++) {
    logFile<<"VERTEX "<<i<<" -> "<<boundGaussPoints[k].getGlobalID()
	   <<".) "<<cCounts[i]<<": ";
    for(int j=startPos;j<startPos+cCounts[i];j++)
      logFile<<tmpConnList[j]<<" ";
    logFile<<endl;
    startPos += cCounts[i];
  }
#endif

  // Loop over a portion of boundary gauss points.
  int m=0;

  for(int i=startIdx;i<endIdx;i++) {

#ifdef _geometryDebugMode_
    logFile<<"**********************"<<endl;
    logFile<<"BOUND GAUSS POINT"<<i<<":"<<endl;
#endif

    idx = 0;

    // Loop over its supported particles.
    for(int j=0;j<boundGaussPoints[i].getSupportCounts();j++) {

      particle = boundGaussPoints[i].getSupportPtcle(j);

#ifdef _geometryDebugMode_
      logFile<<"Ptcle "<<particle<<" incl start "
	     <<inclPtsStartPos[particle]<<" end "
	     <<inclPtsStartPos[particle]+inclPtsCounts[particle]<<endl;
#endif

      // Loop over all gauss points containing a single particle.
      for(int k=inclPtsStartPos[particle];k<inclPtsStartPos[particle]+
	    inclPtsCounts[particle];k++) {

#ifdef _geometryDebugMode_
	logFile<<"check "<<idx<<" "<<idxSize<<" "<<inclPts[k]<<endl;
#endif

	if(boundGaussPoints[i].getGlobalID() == inclPts[k]) {
	  indicator = false;
#ifdef _geometryDebugMode_
	  logFile<<endl;
#endif
	}
	else {
	  indicator = true;

	  // Check, if current gauss point is already contained in list.
	  for(int l=idxSize-idx;l<idxSize;l++) {

#ifdef _geometryDebugMode_
	    logFile<<tmpConnList[l]<<" ";
#endif

	    if(tmpConnList[l] == inclPts[k]) {
	      indicator = false;
#ifdef _geometryDebugMode_
	      logFile<<endl;
#endif
	      break;
	    }

	  }

	}
	// Write current gauss point in the list.
	if(indicator) {
#ifdef _geometryDebugMode_
	  logFile<<" <- write this point"<<endl;
#endif

	  if(idxSize < oldIdxSize)
	    tmpConnList[idxSize] = inclPts[k];
	  else {
	    tmpConnList.push_back(inclPts[k]);

	    logFile <<"Too less data buffer has been allocated for vector "
		    <<"tmpConnList in \n function BackgroundMesh::"
		    <<"setPartitionGraph(existing: "<<oldIdxSize<<" - "
		    <<"needed: "<<idxSize+1<<")!"<< endl;
	  }

	  idxSize++;
	  idx++;
	}

      }

    }

    if(idx < 4) {
      logFile <<"Too less neighbour gauss points in partitioning graph "
	      <<"list for gauss point "<<boundGaussPoints[i].getGlobalID()
	      <<"(existing "<<idx<<" - needed 4 minimum)!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else
      cCounts[localGaussPtsNum+m] = idx;

    m++;
  }

#ifdef _geometryDebugMode_
  startPos=0;
  logFile<<"***************** local tmpConnList ******************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    logFile<<"VERTEX "<<i<<" -> "<<gaussPoints[i].getGlobalID()<<".) "
	   <<cCounts[i]<<": ";
    for(int j=startPos;j<startPos+cCounts[i];j++)
      logFile<<tmpConnList[j]<<" ";
    logFile<<endl;
    startPos += cCounts[i];
  }
  logFile<<"************************************"<<endl;
  for(int i=localGaussPtsNum,k=startIdx;i<localGaussPtsNum+
        localBGaussPtsNum;i++,k++) {
    logFile<<"VERTEX "<<i<<" -> "<<boundGaussPoints[k].getGlobalID()
	   <<".) "<<cCounts[i]<<": ";
    for(int j=startPos;j<startPos+cCounts[i];j++)
      logFile<<tmpConnList[j]<<" ";
    logFile<<endl;
    startPos += cCounts[i];
  }
  logFile<<"reserved buffer for tmpConnlist: "<<oldIdxSize<<" -> needed: "<<idxSize<<endl;
#endif

  //inclPts.~intVector();
  //inclPtsStartPos.~intVector();
  //inclPtsCounts.~intVector();

  /*********************************************************************/
  // Set the vertex weights
  intVector weights;
  weights = cCounts;

  /*********************************************************************/
  // Reduce the local gauss particle connectivity list.
  int startPosOld = 0;
  int startPosNew = 0;

  double connListMultiplier = InputData->getValue("graphListReducing");

  if(connListMultiplier > 1) {
    logFile <<"Input datum 'graphListReducing' mustn't be greater than "
	    <<" 1!"<<endl;
    connListMultiplier = 1.0;
  }

  idx = 0;

  // Determine the max entry number for one graph line.
  for(int i=0;i<localGaussPtsNum+localBGaussPtsNum;i++)

    if(idx < cCounts[i])
      idx = cCounts[i];


  intVector pointNum(idx);
  dbVector distance(idx);

  // Calculate the reserved buffer space for the reduced graph list.
  idxSize = (int)ceil(connListMultiplier*idxSize*1.1);

  intVector connList(idxSize);

#ifdef _geometryDebugMode_
  logFile<<"reserved buffer for connlist: "<<idxSize<<endl;
#endif

  oldIdxSize = idxSize;
  idxSize = 0;
  int globIdx;

  // Loop over all local gauss points.
  for(int i=0;i<localGaussPtsNum;i++){

    globIdx = gaussPoints[i].getGlobalID();
    int k = 0;

    // Loop over the connectivity entries of each gauss point to
    // determine the distance towards.
    while(k < cCounts[i]) {
      int& point = tmpConnList[startPosOld+k];
      distance[k] =
	sqrt(pow(allGCoords[globIdx*3]-allGCoords[point*3],2)
	     +pow(allGCoords[globIdx*3+1]-allGCoords[point*3+1],2)
	     +pow(allGCoords[globIdx*3+2]-allGCoords[point*3+2],2));
      pointNum[k] = point;
      k++;
    }

#ifdef _geometryDebugMode_
    logFile<<i<<" *************** unsorted graphline *****************"<<endl;
    for(int n=0;n<cCounts[i];n++)
      logFile<<pointNum[n]<<" "<<distance[n]<<endl;
#endif

    // Sort for each gauss points its connected gauss points according
    // to their distance towards.
    sortValuesIdx(distance,pointNum,0,cCounts[i]-1);

#ifdef _geometryDebugMode_
    logFile<<i<<" *************** sorted graphline *******************"<<endl;
    for(int n=0;n<cCounts[i];n++)
      logFile<<pointNum[n]<<" ";
    logFile<<endl;
#endif

    startPosNew += cCounts[i];

    // Determine for each gauss points the necessary number of gauss
    // points neighbours.
    cCounts[i] = (int)(ceil(connListMultiplier
			    *cCounts[i]));
    if(cCounts[i] < 4) cCounts[i] = 4;

#ifdef _geometryDebugMode_
    logFile<<i<<" ************* write points in localConn **********"<<endl;
#endif

    k = 0;

    // Loop over the nearest placed nessary number of Gauss points
    // neighbours.
    while(k < cCounts[i]) {


      if(idxSize < oldIdxSize) {
	connList[idxSize] = pointNum[k];
	idxSize++;
	k++;

#ifdef _geometryDebugMode_
	logFile<<pointNum[k]<<" ";
#endif
      }
      else {
	connList.push_back(pointNum[k]);
	idxSize++;
	k++;

#ifdef _geometryDebugMode_
	logFile<<pointNum[k]<<" ";
	logFile<<"Too less data buffer has been allocated for vector "
	       <<"connList in \n function BackgroundMesh::"
	       <<"setPartitionGraph(existing: "<<oldIdxSize<<" - "
	       <<"needed: "<<idxSize+1<<")!"<< endl;
#endif
      }

    }
#ifdef _geometryDebugMode_
    logFile<<endl;
#endif

    startPosOld = startPosNew;
  }

  m = localGaussPtsNum;

  // Loop over a portion of boundary Gauss points.
  for(int i=startIdx;i<endIdx;i++){

    globIdx = boundGaussPoints[i].getGlobalID();
    int k = 0;

    // Loop over the connectivity entries of each gauss point to
    // determine the distance towards.
    while(k < cCounts[m]) {
      int& point = tmpConnList[startPosOld+k];
      distance[k] =
	sqrt(pow(allGCoords[globIdx*3]-allGCoords[point*3],2)
	     +pow(allGCoords[globIdx*3+1]-allGCoords[point*3+1],2)
	     +pow(allGCoords[globIdx*3+2]-allGCoords[point*3+2],2));
      pointNum[k] = point;
      k++;
    }

#ifdef _geometryDebugMode_
    logFile<<i<<" *************** unsorted graphline *****************"<<endl;
    for(int n=0;n<cCounts[m];n++)
      logFile<<pointNum[n]<<" "<<distance[n]<<endl;
#endif

    // Sort for each gauss points its connected gauss points according
    // to their distance towards.
    sortValuesIdx(distance,pointNum,0,cCounts[m]-1);

#ifdef _geometryDebugMode_
    logFile<<i<<" *************** sorted graphline *******************"<<endl;
    for(int n=0;n<cCounts[m];n++)
      logFile<<pointNum[n]<<" ";
    logFile<<endl;
#endif

    startPosNew += cCounts[m];

    // Determine for each gauss points the necessary number of gauss
    // points neighbours.
    cCounts[m] = (int)(ceil(connListMultiplier
			    *cCounts[m]));
    if(cCounts[m] < 4) cCounts[m] = 4;
#ifdef _geometryDebugMode
    logFile<<i<<" ************* write points in localConn **********"<<endl;
#endif

    k = 0;

    // Loop over the nearest placed nessary number of Gauss points
    // neighbours.
    while(k < cCounts[m]) {


      if(idxSize < oldIdxSize) {
	connList[idxSize] = pointNum[k];
	idxSize++;
	k++;

#ifdef _geometryDebugMode_
	logFile<<pointNum[k]<<" ";
#endif
      }
      else {
	connList.push_back(pointNum[k]);
	idxSize++;
	k++;

#ifdef _geometryDebugMode_
	logFile<<pointNum[k]<<" ";
	logFile<<"Too less data buffer has been allocated for vector "
	       <<"connList in \n function BackgroundMesh::"
	       <<"setPartitionGraph(existing: "<<oldIdxSize<<" - "
	       <<"needed: "<<idxSize+1<<")!"<< endl;
#endif
      }

    }

#ifdef _geometryDebugMode_
    logFile<<endl;
#endif

    startPosOld = startPosNew;
    m++;
  }

#ifdef _geometryDebugMode_
  logFile<<"reserved buffer for connlist: "<<oldIdxSize<<" -> needed: "<<idxSize<<endl;
#endif

  //  allgCoords.~dbVector();
  //distance.~dbVector();
  //pointNum.~intVector();
  //tmpConnList.~intVector();

  /********************************************************************/
  // The connListStartPos' part.
  //  int* recvCounts = new int[size];
  //  int* displs = new int[size];

  // Loop over all local gauss points.
  intVector connListStartPos(localGaussPtsNum+localBGaussPtsNum+1);

  for(int i=0;i<=localGaussPtsNum+localBGaussPtsNum;i++) {

    if(i != 0)
      connListStartPos[i] = connListStartPos[i-1]
	+ cCounts[i-1];
    else
      connListStartPos[i] = 0;
  }

  // The globalGPtsRanges' part.
  intVector globalGPtsRanges(size+1);
  intVector tmpRanges(size);

  int sendBuf = localGaussPtsNum+localBGaussPtsNum;
  MPI_Allgather(&sendBuf,1,MPI_INT,&tmpRanges[0],1,
		MPI_INT,MPI_COMM_WORLD);

  globalGPtsRanges[0] = 0;

  // Loop over all processors.
  for(int i=1;i<=size;i++)
    globalGPtsRanges[i]= globalGPtsRanges[i-1] + tmpRanges[i-1];

#ifdef _geometryDebugMode_
  int sumElems = 0;
  logFile<<"******************************************************"<<endl;
  logFile<<"*************** localConnStartPos ********************"<<endl;
  for(int i=0;i<=localGaussPtsNum+localBGaussPtsNum;i++)
    logFile<<i<<".) "<<connListStartPos[i]<<endl;

  logFile<<"******************* local connList *******************"<<endl;
  logFile<<"connList size = "<<connList.size()<<endl;
  for(int i=0;i<localGaussPtsNum+localBGaussPtsNum;i++) {
    logFile<<"VERTEX "<<i<<".) "<<cCounts[i]<<": ";
    sumElems += cCounts[i];
    for(int j=connListStartPos[i];j<connListStartPos[i]+cCounts[i];j++)
      logFile<<connList[j]<<" ";
    logFile<<endl;
  }
  logFile<<"**************** localVertexWeights ******************"<<endl;
  for(int i=0;i<localGaussPtsNum+localBGaussPtsNum;i++)
    logFile<<i<<" "<<weights[i]<<endl;
  logFile<<"**************** globalGPtsRanges ********************"<<endl;
  for(int i=0;i<=size;i++)
    logFile<<globalGPtsRanges[i]<<" ";
  logFile<<endl;
  logFile<<"sumElems = "<<sumElems<<" = "<<"connListEndPos = "
	 <<connListStartPos[localGaussPtsNum+localBGaussPtsNum]<<endl;
#endif

  // Determine the array containing the global indices of the local
  // gauss points - after partitioning it contains for each gauss point
  // the processor identifier to which it newly belongs.
  gaussProcList = intVector(localGaussPtsNum+localBGaussPtsNum);

  for(int i=0;i<localGaussPtsNum;i++)
    gaussProcList[i] = gaussPoints[i].getGlobalID();

  m=localGaussPtsNum;

  for(int i=startIdx;i<endIdx;i++) {
    gaussProcList[m] = boundGaussPoints[i].getGlobalID();
    m++;
  }

  /*******************************************************************/
  // Partition the gauss points up on the processors.

#ifdef _geometryDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"************** ParMETIS Partitioning ****************"<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"gaussProcList size = "<<gaussProcList.size()<<endl;
  logFile<<"**************** old gauss proc list  ***************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++)
    logFile<<i<<".) "<<gaussPoints[i].getGlobalID()<<": "
	   <<gaussProcList[i]<<endl;
  logFile<<"*****************************************************"<<endl;
  for(int i=localGaussPtsNum,j=startIdx;i<gaussProcList.size();i++,j++)
    logFile<<i<<".) "<<boundGaussPoints[j].getGlobalID()<<": "
	   <<gaussProcList[i]<<endl;
#endif


  partitionGaussPoints(InputData,connList,weights,connListStartPos,
		       globalGPtsRanges,gaussProcList,logFile);

#ifdef _geometryDebugMode_
  logFile<<"************** new gauss proc list *******************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++)
    logFile<<i<<".) "<<gaussPoints[i].getGlobalID()<<": "
	   <<gaussProcList[i]<<endl;
  logFile<<"******************************************************"<<endl;
  for(int i=localGaussPtsNum,j=startIdx;i<gaussProcList.size();i++,j++)
    logFile<<i<<".) "<<boundGaussPoints[j].getGlobalID()<<": "
	   <<gaussProcList[i]<<endl;
#endif

  //  cCounts.~intVector();
  //  delete[] displs,recvCounts;

}

/**********************************************************************/
/**********************************************************************/
// Partition the gauss points to the processors.
void BackgroundMesh::partitionGaussPoints(InputFileData* InputData,
                                          intVector& connList,
                                          intVector& weights,
                                          intVector& connListStartPos,
                                          intVector& globalGPtsRanges,
                                          intVector& gaussProcList,
                                          std::ofstream& logFile) {

  using namespace std;

  int edgeCut,nparts,size;
  //int options[4];
  intVector options(4);
  int wgtFlag = 2;
  int numFlag = 0;
  MPI_Comm comm;

  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  nparts = size;
  options[0] = 0;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  logFile<<"In BackgroundMesh::partitionGaussPoints ParMETIS is "
	 <<"deactivated.!"<<endl;
  //MPI_Abort(MPI_COMM_WORLD,1);

  //  ParMETIS_PartKway(&globalGPtsRanges[0],&connListStartPos[0],
  //		    &connList[0],&weights[0],NULL,&wgtFlag,&numFlag,
  //		    &nparts,&options[0],&edgeCut,&gaussProcList[0],&comm);

  MPI_Comm_free(&comm);

#ifdef _geometryDebugMode_
  logFile<<"edgeCut "<<edgeCut<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Exchange the volume gauss data between the processor that each of
// them possesses the correct local data.
void BackgroundMesh::exchangeGaussData(InputFileData* InputData,
                                       intVector& gaussProcList,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {

  using namespace std;

  int rank,size;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  /**********************************************************************/
  // Get from the other processors the array containing the processor
  // identifiers to which each gauss point newly belongs.
  intVector allLocalGPtsPortions(size);
  intVector allLocalDispls(size);
  intVector allOtherProcLists(globalGaussPtsNum+globalBGaussPtsNum);

  int allLocalGPtsNum = localGaussPtsNum+localBGaussPtsNum;

  MPI_Allgather(&allLocalGPtsNum,1,MPI_INT,&allLocalGPtsPortions[0],1,
		MPI_INT,MPI_COMM_WORLD);

  allLocalDispls[0] = 0;

  for(int i=1;i<size;i++)
    allLocalDispls[i] = allLocalDispls[i-1]+allLocalGPtsPortions[i-1];

  MPI_Allgatherv(&gaussProcList[0],allLocalGPtsPortions[rank],MPI_INT,
		 &allOtherProcLists[0],&allLocalGPtsPortions[0],
		 &allLocalDispls[0],MPI_INT,MPI_COMM_WORLD);

  /**********************************************************************/
  // Assemble vector which contains the local portions of volume Gauss
  // points of each processor.
  intVector localGPtsPortions(size);
  intVector localDispls(size);

  MPI_Allgather(&localGaussPtsNum,1,MPI_INT,&localGPtsPortions[0],1,
		MPI_INT,MPI_COMM_WORLD);

  localDispls[0] = 0;

  for(int i=1;i<size;i++)
    localDispls[i] = localDispls[i-1]+localGPtsPortions[i-1];

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"********* volume gauss point data exchanging *********"<<endl;
  logFile<<"******************************************************"<<endl;
  logFile<<"************* local volume Gauss portions ************"<<endl;
  for(int i=0;i<size;i++)
    logFile<<"rank "<<i<<": "<<localGPtsPortions[i]<<endl;
  logFile<<"***************** local Gauss portions ***************"<<endl;
  for(int i=0;i<size;i++)
    logFile<<"rank "<<i<<": "<<allLocalGPtsPortions[i]<<endl;
  logFile<<"*************** global processor list ****************"<<endl;
  for(int i=0;i<allOtherProcLists.size();i++)
    logFile<<i<<": "<<allOtherProcLists[i]<<endl;
#endif

  /**********************************************************************/
  // Assemble the local vector containing the global gauss indices.
  intVector globalVolGaussIdx(globalGaussPtsNum);
  intVector localGlobalIdx(localGaussPtsNum);

  for(int i=0;i<localGaussPtsNum;i++) {
    localGlobalIdx[i] = gaussPoints[i].getGlobalID();
  }

  // Get from the other processors the array containing the global gauss
  // indices of their local gauss vector - needed as tags for
  // communication!
  MPI_Allgatherv(&localGlobalIdx[0],localGPtsPortions[rank],MPI_INT,
		 &globalVolGaussIdx[0],&localGPtsPortions[0],&localDispls[0],
		 MPI_INT,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"****************** local globalIdx *******************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    logFile<<i<<" "<<localGlobalIdx[i]<<endl;
  }
  logFile<<"***************** global globalIdx *******************"<<endl;
  for(int i=0;i<globalGaussPtsNum;i++)
    logFile<<i<<".): "<<globalVolGaussIdx[i]<<endl;
  logFile<<"********** global gauss -> processor list ************"<<endl;
  for(int i=0;i<size;i++)
    for(int j=0;j<localGPtsPortions[i];j++)
      logFile<<"GAUSSPOINT "<<globalVolGaussIdx[localDispls[i]+j]<<" -> "
	     <<allOtherProcLists[allLocalDispls[i]+j]<<endl;
#endif

  localGlobalIdx.resize(0);

  /**********************************************************************/
  // Determine the sending and receiving Gauss points and so the new size
  // of vector 'gaussPoints' too.
  int gPoint,idx;
  int newLocalGaussPtsNum = 0;
  intVector recvGaussIdx;  // contains the receiving global indices
  recvGaussIdx.reserve(localGaussPtsNum);
  intVector sendGaussIdx; // contains the sending global indices
  sendGaussIdx.reserve(localGaussPtsNum);

  int recv=0;
  int send=0;

  // Loop over all volume Gauss points.

  // Loop over all processors.
  for(int i=0;i<size;i++) {

    // Loop over a local portions of volume Gauss points.
    for(int j=0;j<localGPtsPortions[i];j++) {
      gPoint = allLocalDispls[i]+j;
      idx = localDispls[i]+j;

      // Current Gauss point will be a local one.
      if(allOtherProcLists[gPoint] == rank) {

	// Current Gauss point is still not a local one.
	if(gPoint < allLocalDispls[rank] ||
	   gPoint >= allLocalDispls[rank]+allLocalGPtsPortions[rank]) {

	  if(recv < recvGaussIdx.size())
	    recvGaussIdx[recv] = globalVolGaussIdx[idx];
	  else
	    recvGaussIdx.push_back(globalVolGaussIdx[idx]);

	  recv++;
	}

	newLocalGaussPtsNum++;
      }

      // Current Gauss point will not be a local one.
      else {

	// Current Gauss point will not be a local one any more.
	if(gPoint >= allLocalDispls[rank] &&
	   gPoint < allLocalDispls[rank]+allLocalGPtsPortions[rank]) {

	  send++;
	}

      }

    }

  }

  recvGaussIdx.resize(recv);
  sendGaussIdx.resize(send);

  // Determine the communication size.
  int commSize = recv + send;

#ifdef _geometryDebugMode_
  logFile<<"************* receiving Gauss points *****************"<<endl;
  for(int i=0;i<recvGaussIdx.size();i++)
    logFile<<i<<".): "<<recvGaussIdx[i]<<endl;
  logFile<<"Communication amount "<<commSize<<" particles "
	 <<sendGaussIdx.size()<<" sending "<<recvGaussIdx.size()
         <<" receiving"<<endl;
  logFile<<"************* allOtherProcLists **********************"<<endl;
  for(int i=0;i<allOtherProcLists.size();i++)
    logFile<<i<<" "<<allOtherProcLists[i]<<endl;
#endif

  // --------------------------------------------------------------------
  // Determine the global integration point indices body forces are
  // applied.
  int sendIdx = bodyForceGaussPtsIdx.size();
  int tmpGlobalSize;
  intVector idxDispls(size);
  intVector recvIdx(size);

  MPI_Allreduce(&sendIdx,&tmpGlobalSize,1,
		MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  MPI_Allgather(&sendIdx,1,MPI_INT,&recvIdx[0],1,MPI_INT,MPI_COMM_WORLD);

  idxDispls[0] = 0;

  for(int i=1;i<size;i++)
    idxDispls[i] = idxDispls[i-1]+recvIdx[i-1];

  intVector localGlobalBodyForceIdx(bodyForceGaussPtsIdx.size());

  for(int i=0;i<bodyForceGaussPtsIdx.size();i++)
    localGlobalBodyForceIdx[i] =
      gaussPoints[bodyForceGaussPtsIdx[i][0]].getGlobalID();

  intVector tmpGlobalBodyForceIdx(tmpGlobalSize);

  MPI_Allgatherv(&localGlobalBodyForceIdx[0],recvIdx[rank],MPI_INT,
		 &tmpGlobalBodyForceIdx[0],&recvIdx[0],&idxDispls[0],
		 MPI_INT,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"*********** old bodyForceGaussPtsIdx ****************"<<endl;
  for(int i=0;i<bodyForceGaussPtsIdx.size();i++)
    logFile<<i<<" local "<<bodyForceGaussPtsIdx[i][0]
	   <<" global int point "<<localGlobalBodyForceIdx[i]
	   <<" with vol load"<<endl;
#endif

  // Prepare vector 'bodyForceGaussPtsIdx' for a new initialisation.
  intMatrix oldBodyForceGaussPtsIdx = bodyForceGaussPtsIdx;
  bodyForceGaussPtsIdx.resize(newLocalGaussPtsNum);
  clearArray(bodyForceGaussPtsIdx);

  intVector globalBodyForceIdx(globalGaussPtsNum+globalBGaussPtsNum);

  for(int i=0;i<tmpGlobalBodyForceIdx.size();i++)
    globalBodyForceIdx[tmpGlobalBodyForceIdx[i]] = 1;

#ifdef _geometryDebugMode_
  logFile<<"************* tmpGlobalBodyForceIdx **************"<<endl;
  for(int i=0;i<tmpGlobalBodyForceIdx.size();i++)
    logFile<<"global GAUSS POINT "<<tmpGlobalBodyForceIdx[i]<<endl;
  logFile<<"************* globalBodyForceIdx *****************"<<endl;
  logFile<<"tmpGlobalSize "<<tmpGlobalSize<<" <= globalGaussPtsNum "
	 <<globalGaussPtsNum<<endl;
  for(int i=0;i<globalBodyForceIdx.size();i++)
    logFile<<"global int point "<<i
	   <<" bool applied load "<<globalBodyForceIdx[i]<<endl;
#endif

  localGlobalBodyForceIdx.resize(0);
  tmpGlobalBodyForceIdx.resize(0);


  // --------------------------------------------------------------------
  // Determine the global integration point indices body moments are
  // applied.

  if(bodyMomentGaussPtsIdx.size() > 0) {
    logFile<<"In BackgroundMesh::exchangeGaussData body moment loads "
	   <<"are\n not supported yet!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // --------------------------------------------------------------------
  // Determine the global integration point indices body electric
  // charge loads are applied.
  sendIdx = bodyElectricChargeGaussPtsIdx.size();

  MPI_Allreduce(&sendIdx,&tmpGlobalSize,1,
		MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  MPI_Allgather(&sendIdx,1,MPI_INT,&recvIdx[0],1,MPI_INT,MPI_COMM_WORLD);

  idxDispls[0] = 0;

  for(int i=1;i<size;i++)
    idxDispls[i] = idxDispls[i-1]+recvIdx[i-1];

  intVector localGlobalBodyChargeIdx(bodyElectricChargeGaussPtsIdx.size());

  for(int i=0;i<bodyElectricChargeGaussPtsIdx.size();i++)
    localGlobalBodyChargeIdx[i] =
      gaussPoints[bodyElectricChargeGaussPtsIdx[i][0]].getGlobalID();

  intVector tmpGlobalBodyChargeIdx(tmpGlobalSize);

  MPI_Allgatherv(&localGlobalBodyChargeIdx[0],recvIdx[rank],MPI_INT,
		 &tmpGlobalBodyChargeIdx[0],&recvIdx[0],&idxDispls[0],
		 MPI_INT,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"*********** old bodyElectricChargeGaussPtsIdx ****************"<<endl;
  for(int i=0;i<bodyElectricChargeGaussPtsIdx.size();i++)
    logFile<<i<<" local "<<bodyElectricChargeGaussPtsIdx[i][0]
	   <<" global int point "<<localGlobalBodyChargeIdx[i]
	   <<" with vol load"<<endl;
#endif

  // Prepare vector 'bodyElectricChargeGaussPtsIdx' for a new initialisation.
  intMatrix oldBodyChargeGaussPtsIdx = bodyElectricChargeGaussPtsIdx;
  bodyElectricChargeGaussPtsIdx.resize(newLocalGaussPtsNum);
  clearArray(bodyElectricChargeGaussPtsIdx);

  intVector globalBodyChargeIdx(globalGaussPtsNum+globalBGaussPtsNum);

  for(int i=0;i<tmpGlobalBodyChargeIdx.size();i++)
    globalBodyChargeIdx[tmpGlobalBodyChargeIdx[i]] = 1;

#ifdef _geometryDebugMode_
  logFile<<"************* tmpGlobalBodyChargeIdx **************"<<endl;
  for(int i=0;i<tmpGlobalBodyChargeIdx.size();i++)
    logFile<<"global GAUSS POINT "<<tmpGlobalBodyChargeIdx[i]<<endl;
  logFile<<"************* globalBodyChargeIdx *****************"<<endl;
  logFile<<"tmpGlobalSize "<<tmpGlobalSize<<" <= globalGaussPtsNum "
	 <<globalGaussPtsNum<<endl;
  for(int i=0;i<globalBodyChargeIdx.size();i++)
    logFile<<"global int point "<<i
	   <<" bool applied load "<<globalBodyChargeIdx[i]<<endl;
#endif

  globalVolGaussIdx.resize(0);
  localGlobalBodyChargeIdx.resize(0);
  tmpGlobalBodyChargeIdx.resize(0);


  /*********************************************************************/
  // Restructure the Gauss point vector and its staying Gauss points.
  int usedDims = (int)modelData["usedDimensions"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];

  dbMatrix sendCoords(sendGaussIdx.size());
  intMatrix sendBodyForceDOF(sendGaussIdx.size());
  dbMatrix sendBodyForces(sendGaussIdx.size());
  intMatrix sendBodyChargeDOF(sendGaussIdx.size());
  dbMatrix sendBodyCharges(sendGaussIdx.size());
  dbVector sendWeights(sendGaussIdx.size());
  int gaussIdx=0;
  int vecEndIdx = localGaussPtsNum - 1;
  int bodyForceGPtIdx = 0;
  int bodyChargeGPtIdx = 0;
  int m=0;

  // Loop over the assemble of local Gauss points.
  for(int i=0;i<localGaussPtsNum;i++) {

    // Current Gauss point 'i' will not stay.
    if(gaussProcList[i] != rank) {

      // Store Gauss point data temporarily to get sent.
      sendGaussIdx[m] = gaussPoints[i].getGlobalID();
      sendCoords[m] = gaussPoints[i].getCoords();
      sendWeights[m] = gaussPoints[i].getWeight();

      if(globalBodyForceIdx[sendGaussIdx[m]] == 1) {

	blVector& affectedDOF =  gaussPoints[i].getBodyForceDOF();

	for(int k=0;k<affectedDOF.size();k++)
	  sendBodyForceDOF[m][k] = (int)affectedDOF[k];

	sendBodyForces[m] = gaussPoints[i].getBodyForce();
      }

      if(globalBodyChargeIdx[sendGaussIdx[m]] == 1) {

	blVector& affectedDOF =  gaussPoints[i].getBodyElectricChargeDOF();

	for(int k=0;k<affectedDOF.size();k++)
	  sendBodyChargeDOF[m][k] = (int)affectedDOF[k];

	sendBodyCharges[m] = gaussPoints[i].getBodyElectricCharge();
      }

      m++;

      while(vecEndIdx > i) {

	// Copy this Gauss point at the end of the Gauss point vector to
	// position 'i'
	if(gaussProcList[vecEndIdx] == rank) {

	  gaussPoints[i].setGlobalID(gaussPoints[vecEndIdx].getGlobalID());
	  gaussPoints[i].setCoords(gaussPoints[vecEndIdx].getCoords());
	  gaussPoints[i].setWeight(gaussPoints[vecEndIdx].getWeight());

	  if(globalBodyForceIdx[gaussPoints[vecEndIdx].getGlobalID()] == 1) {
	    blVector& affectedDOF = gaussPoints[i].getBodyForceDOF();
	    dbVector& conditions = gaussPoints[i].getBodyForce();

	    // default intialisation if necessary
	    if(conditions.size() < usedDims) {
	      conditions.resize(usedDims);
	      affectedDOF.resize(usedDims);
	    }

	    affectedDOF = gaussPoints[vecEndIdx].getBodyForceDOF();
	    conditions = gaussPoints[vecEndIdx].getBodyForce();
	    bodyForceGaussPtsIdx[bodyForceGPtIdx][0] = i; // point ID
	    bodyForceGaussPtsIdx[bodyForceGPtIdx][1] = 0; // weight ID
	    bodyForceGaussPtsIdx[bodyForceGPtIdx][2] = 0; // load ID
	    bodyForceGPtIdx++;
	  }


	  if(globalBodyChargeIdx[gaussPoints[vecEndIdx].getGlobalID()] == 1) {
	    blVector& affectedDOF = gaussPoints[i].getBodyElectricChargeDOF();
	    dbVector& conditions = gaussPoints[i].getBodyElectricCharge();

	    // default intialisation if necessary
	    if(conditions.size() < usedDims) {
	      conditions.resize(usedDims);
	      affectedDOF.resize(usedDims);
	    }

	    affectedDOF = gaussPoints[vecEndIdx].getBodyElectricChargeDOF();
	    conditions = gaussPoints[vecEndIdx].getBodyElectricCharge();
	    bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][0] = i; // point ID
	    bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][1] = 0; // weight ID
	    bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][2] = 0; // load ID
	    bodyChargeGPtIdx++;
	  }

	  gaussIdx++;
	  vecEndIdx--;
	  break;
	}

	vecEndIdx--;
      }

    }

    // Current Gauss point will stay.
    else if(gaussProcList[i] == rank && i <= vecEndIdx) {

      if(globalBodyForceIdx[gaussPoints[i].getGlobalID()] == 1) {
	bodyForceGaussPtsIdx[bodyForceGPtIdx][0] = i; // point ID
	bodyForceGaussPtsIdx[bodyForceGPtIdx][1] = 0; // weight ID
	bodyForceGaussPtsIdx[bodyForceGPtIdx][2] = 0; // load ID
	bodyForceGPtIdx++;
      }

      if(globalBodyChargeIdx[gaussPoints[i].getGlobalID()] == 1) {
	bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][0] = i; // point ID
	bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][1] = 0; // weight ID
	bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][2] = 0; // load ID
	bodyChargeGPtIdx++;
      }

      gaussIdx++;
    }

  }

#ifdef _geometryDebugMode_
  logFile<<"********** restructured gauss point vector ***********"<<endl;
  logFile<<"newLoalGaussPtsNum = "<<newLocalGaussPtsNum<<endl;
  logFile<<"gaussIdx = "<<gaussIdx<<endl;
  logFile<<"bodyForceGPtIdx = "<<bodyForceGPtIdx<<endl;
  logFile<<"bodyChargeGPtIdx = "<<bodyChargeGPtIdx<<endl;
  for(int i=0;i<gaussIdx;i++)
    logFile<<i<<".) GAUSS POINT "<<gaussPoints[i].getGlobalID()
	   <<" coords: "<<gaussPoints[i].getCoord(0)<<" "
	   <<gaussPoints[i].getCoord(1)<<" "<<gaussPoints[i].getCoord(2)
	   <<" weight = "<<gaussPoints[i].getWeight()<<endl;
  logFile<<"*************** sending Gauss points *****************"<<endl;
  for(int i=0;i<sendGaussIdx.size();i++)
    logFile<<i<<".): "<<sendGaussIdx[i]<<endl;
  logFile<<"temp stored data set: "<<m<<endl;
  for(int i=0;i<m;i++) {
    logFile<<i<<".) GAUSS POINT "<<sendGaussIdx[i]<<" coords: "
	   <<sendCoords[i][0]<<" "<<sendCoords[i][1]<<" "<<sendCoords[i][2]
	   <<" weight "<<sendWeights[i]<<" loads ";
    for(int j=0;j<sendBodyForces[i].size();j++)
      logFile<<"DOF "<<sendBodyForceDOF[i][j]
	     <<" cond "<<sendBodyForces[i][j]<<" ";
    logFile<<endl;
    for(int j=0;j<sendBodyCharges[i].size();j++)
      logFile<<"DOF "<<sendBodyChargeDOF[i][j]
	     <<" cond "<<sendBodyCharges[i][j]<<" ";
    logFile<<endl;
  }
  logFile<<"**************** bodyForceGaussPtsIdx ****************"<<endl;
  for(int i=0;i<bodyForceGaussPtsIdx.size();i++)
    logFile<<i<<" local "<<bodyForceGaussPtsIdx[i][0]
	   <<" global int point with body force"<<endl;
  logFile<<"**************** bodyElectricChargeGaussPtsIdx *******"<<endl;
  for(int i=0;i<bodyElectricChargeGaussPtsIdx.size();i++)
    logFile<<i<<" local "<<bodyElectricChargeGaussPtsIdx[i][0]
	   <<" global int point with body charge"<<endl;
#endif

  // Resize Gauss point vector.
  if(newLocalGaussPtsNum != localGaussPtsNum)
    gaussPoints.resize(newLocalGaussPtsNum);

  localGaussPtsNum = newLocalGaussPtsNum;

  /**********************************************************************/
  // Loop to exchange the gauss points between the processor according
  // to the new calculated partitioning.
  int gEE = 6;      // number of exchanging elements of class 'GaussPoint'
  int recvd = 0;
  int sent = 0;
  int i=0;

  MPI_Request* reqs = new MPI_Request[commSize*gEE];
  MPI_Status* stats = new MPI_Status[commSize*gEE];

  intVector dummyDOF(usedDims);

  m=0;

#ifdef _geometryDebugMode_
  logFile<<"****************** data exchange *********************"<<endl;
#endif

  // Loop over all Gauss point has to be sent and received.
  while (m < sendGaussIdx.size() || m < recvGaussIdx.size()) {

    // Set a receive command.
    if(m < recvGaussIdx.size()) {

      // Receive data from another processor.

      // global index
      MPI_Irecv(&gaussPoints[gaussIdx].getGlobalID(),1,MPI_INT,
		MPI_ANY_SOURCE,recvGaussIdx[m]*gEE,
		MPI_COMM_WORLD,&reqs[sent+recvd]);
      recvd++;

      // coordinates
      MPI_Irecv(&gaussPoints[gaussIdx].getCoord(0),usedDims,
		MPI_DOUBLE,MPI_ANY_SOURCE,recvGaussIdx[m]*gEE+1,
		MPI_COMM_WORLD,&reqs[sent+recvd]);
      recvd++;

      // weight
      MPI_Irecv(&gaussPoints[gaussIdx].getWeight(),1,MPI_DOUBLE,
		MPI_ANY_SOURCE,recvGaussIdx[m]*gEE+2,
		MPI_COMM_WORLD,&reqs[sent+recvd]);
      recvd++;

      // body force
      if(globalBodyForceIdx[recvGaussIdx[m]] == 1) {

	blVector& affectedDOF = gaussPoints[gaussIdx].getBodyForceDOF();
	dbVector& conditions = gaussPoints[gaussIdx].getBodyForce();

	// default intialisation if necessary
	if(conditions.size() < usedDims) {
	  affectedDOF.resize(usedDims);
	  conditions.resize(usedDims);
	}

	MPI_Irecv(&dummyDOF[0],usedDims,MPI_INT,
		  MPI_ANY_SOURCE,recvGaussIdx[m]*gEE+3,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	for(int k=0;k<dummyDOF.size();k++)
	  affectedDOF[k] = (bool)dummyDOF[k];

	MPI_Irecv(&conditions[0],usedDims,MPI_DOUBLE,
		  MPI_ANY_SOURCE,recvGaussIdx[m]*gEE+3,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	bodyForceGaussPtsIdx[bodyForceGPtIdx][0] = gaussIdx; // point ID
	bodyForceGaussPtsIdx[bodyForceGPtIdx][1] = 0; // weight ID
	bodyForceGaussPtsIdx[bodyForceGPtIdx][2] = 0; // load ID
	bodyForceGPtIdx++;
	recvd++;
      }

      // material ID
      MPI_Irecv(&gaussPoints[gaussIdx].getMaterialID(),1,MPI_INT,
		MPI_ANY_SOURCE,recvGaussIdx[m]*gEE+4,
		MPI_COMM_WORLD,&reqs[sent+recvd]);
      recvd++;

      // body electric charge
      if(globalBodyChargeIdx[recvGaussIdx[m]] == 1) {

	blVector& affectedDOF = gaussPoints[gaussIdx].getBodyElectricChargeDOF();
	dbVector& conditions = gaussPoints[gaussIdx].getBodyElectricCharge();

	// default intialisation if necessary
	if(conditions.size() < usedDims) {
	  affectedDOF.resize(usedDims);
	  conditions.resize(usedDims);
	}

	MPI_Irecv(&dummyDOF[0],usedDims,MPI_INT,
		  MPI_ANY_SOURCE,recvGaussIdx[m]*gEE+5,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	for(int k=0;k<dummyDOF.size();k++)
	  affectedDOF[k] = (bool)dummyDOF[k];

	MPI_Irecv(&conditions[0],usedDims,MPI_DOUBLE,
		  MPI_ANY_SOURCE,recvGaussIdx[m]*gEE+5,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][0] = gaussIdx; // point ID
	bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][1] = 0; // weight ID
	bodyElectricChargeGaussPtsIdx[bodyChargeGPtIdx][2] = 0; // load ID
	bodyChargeGPtIdx++;
	recvd++;
      }

#ifdef _geometryDebugMode_
      logFile<<m+1<<".) received point "<<recvGaussIdx[m]
	     <<" -> local point "<<gaussIdx<<endl;
      logFile<<"message req "<<sent+recvd<<endl;
#endif

      gaussIdx++;
    }

    // Set a send command.
    if(m < sendGaussIdx.size()) {

      // Send data to another processor and replace it by data from
      // another processor.

      // global index
      MPI_Isend(&sendGaussIdx[m],1,MPI_INT,
		allOtherProcLists[sendGaussIdx[m]],
		sendGaussIdx[m]*gEE,MPI_COMM_WORLD,
		&reqs[sent+recvd]);
      sent++;

      // coordinates
      MPI_Isend(&sendCoords[m][0],usedDims,MPI_DOUBLE,
		allOtherProcLists[sendGaussIdx[m]],
		sendGaussIdx[m]*gEE+1,
		MPI_COMM_WORLD,&reqs[sent+recvd]);
      sent++;

      // weight
      MPI_Isend(&sendWeights[m],1,MPI_DOUBLE,
		allOtherProcLists[sendGaussIdx[m]],
		sendGaussIdx[m]*gEE+2,
		MPI_COMM_WORLD,&reqs[sent+recvd]);
      sent++;

      // body force
      if(globalBodyForceIdx[sendGaussIdx[m]] == 1) {

	MPI_Isend(&sendBodyForceDOF[m][0],usedDims,MPI_INT,
		  allOtherProcLists[sendGaussIdx[m]],
		  sendGaussIdx[m]*gEE+3,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	MPI_Isend(&sendBodyForces[m][0],usedDims,MPI_DOUBLE,
		  allOtherProcLists[sendGaussIdx[m]],
		  sendGaussIdx[m]*gEE+3,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	sent++;
      }

      // material ID
      MPI_Isend(&sendGaussIdx[m],1,MPI_INT,
		allOtherProcLists[sendGaussIdx[m]],
		sendGaussIdx[m]*gEE+4,MPI_COMM_WORLD,
		&reqs[sent+recvd]);
      sent++;

      // body electric charge
      if(globalBodyChargeIdx[sendGaussIdx[m]] == 1) {

	MPI_Isend(&sendBodyChargeDOF[m][0],usedDims,MPI_INT,
		  allOtherProcLists[sendGaussIdx[m]],
		  sendGaussIdx[m]*gEE+5,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	MPI_Isend(&sendBodyCharges[m][0],usedDims,MPI_DOUBLE,
		  allOtherProcLists[sendGaussIdx[m]],
		  sendGaussIdx[m]*gEE+5,
		  MPI_COMM_WORLD,&reqs[sent+recvd]);

	sent++;
      }

#ifdef _geometryDebugMode_
      logFile<<m+1<<".) sent point "<<sendGaussIdx[m]<<endl;
      logFile<<"message req "<<sent+recvd<<endl;
#endif

    }

#ifdef _geometryDebugMode_
    logFile<<"commIndex "<<m<<" maxCommIndex "
	   <<(sendGaussIdx.size() < recvGaussIdx.size() ?
	      recvGaussIdx.size() : sendGaussIdx.size())
	   <<" newGaussIdx "<<gaussIdx<<" sentCount "<<sent
	   <<" receiveCount "<<recvd<<endl;
#endif

    m++;
  }

  MPI_Waitall(sent+recvd,reqs,stats);


  // Resize vectors.
  bodyForceGaussPtsIdx.resize(bodyForceGPtIdx);
  bodyElectricChargeGaussPtsIdx.resize(bodyChargeGPtIdx);

  /*********************************************************************/
  // Compute for the local volume Gauss points a globally consecutive
  // ordering
  MPI_Allgather(&localGaussPtsNum,1,MPI_INT,&localGPtsPortions[0],1,
		MPI_INT,MPI_COMM_WORLD);

  // Loop over all processors to set the global volume Gauss point start
  // and end index locally, but excluding the boundary Gauss points.
  int startGIdx = 0;

  m=0;

  while(m < rank) {
    startGIdx += localGPtsPortions[m];
    m++;
  }

  int endGIdx = startGIdx + localGaussPtsNum;

  // Loop over the local Gauss points of current processor.
  for(int i=0,j=startGIdx;j<endGIdx;i++,j++)
    gaussPoints[i].setGlobalID(j);



#ifdef _geometryDebugMode_
  logFile<<rank<<" sent "<<sent<<" received "<<recvd<<endl;
  logFile<<"new localGaussPtsNum "<<localGaussPtsNum<<endl;
#endif
#ifdef _partitioningDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"******** new assembled gauss point data **************"<<endl;
  logFile<<"********************* coords *************************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    logFile<<i<<".) global index "<<gaussPoints[i].getGlobalID()<<" ";
    for(int j=0;j<usedDims;j++)
      logFile<<gaussPoints[i].getCoord(j)<<" ";
    logFile<<endl;
  }
#endif
#ifdef _geometryDebugMode_
  logFile<<"********************* weights ************************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++)
    logFile<<i<<".) global index "<<gaussPoints[i].getGlobalID()<<" "
	   <<gaussPoints[i].getWeight()<<endl;
  logFile<<"*********** local body force integration points *****"<<endl;
  for(int i=0;i<bodyForceGaussPtsIdx.size();i++) {
    int& ID = bodyForceGaussPtsIdx[i][0];
    blVector& affectedDOF = gaussPoints[ID].getBodyForceDOF();
    dbVector& conditions = gaussPoints[ID].getBodyForce();
    logFile<<i<<".) global index "<<gaussPoints[ID].getGlobalID()<<" ";
    for(int j=0;j<conditions.size();j++)
      logFile<<"DOF "<<affectedDOF[j]<<" cond "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"*********** local body elec charge integration points *****"<<endl;
  for(int i=0;i<bodyElectricChargeGaussPtsIdx.size();i++) {
    int& ID = bodyElectricChargeGaussPtsIdx[i][0];
    blVector& affectedDOF = gaussPoints[ID].getBodyElectricChargeDOF();
    dbVector& conditions = gaussPoints[ID].getBodyElectricCharge();
    logFile<<i<<".) global index "<<gaussPoints[ID].getGlobalID()<<" ";
    for(int j=0;j<conditions.size();j++)
      logFile<<"DOF "<<affectedDOF[j]<<" cond "<<conditions[j]<<" ";
    logFile<<endl;
  }
#endif

  delete[] reqs,stats;
}

/************************************************************************/
/************************************************************************/
// Determine the initial boundary Gauss point distribution.
void BackgroundMesh::getInitBGaussSplitting(InputFileData* InputData,
                                            int& defaultPortion,
                                            int& startIdx,int& endIdx,
                                            std::map<std::string,double>& modelData,
                                            std::ofstream& logFile) {

  using namespace std;

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  defaultPortion = (int)ceil((double)globalBGaussPtsNum/size);

  if(defaultPortion*(rank+1) <= globalBGaussPtsNum) {
    startIdx = defaultPortion*rank;
    endIdx = defaultPortion*(rank+1);
  }

  else if(defaultPortion*rank <= globalBGaussPtsNum) {
    startIdx = defaultPortion*rank;
    endIdx = globalBGaussPtsNum;
  }

  else {
    startIdx = globalBGaussPtsNum;
    endIdx = globalBGaussPtsNum;
  }

}

/**********************************************************************/
/**********************************************************************/
// Split the boundary Gauss points between the processor that each of
// them possesses the correct local data.
void BackgroundMesh::splitBoundGaussPts(InputFileData* InputData,
                                        intVector& gaussProcList,
                                        intVector& globalLocalIdx,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  intVector allOtherProcLists(globalBGaussPtsNum);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int* displs = new int[size];
  int* recvCounts = new int[size];

  /**********************************************************************/
  // Get from the other processors the array containing the processor
  // identifiers to which each gauss point newly belongs.

  // Boundary Gauss points portions.
  MPI_Allgather(&localBGaussPtsNum,1,MPI_INT,recvCounts,1,MPI_INT,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  int idx = gaussProcList.size()-localBGaussPtsNum;

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"******** boundary gauss point data splitting *********"<<endl;
  logFile<<"******************************************************"<<endl;
  logFile<<"*********** local gauss -> processor list ************"<<endl;
  logFile<<"localBGaussPtsNum = "<<localBGaussPtsNum<<endl;
  for(int i=idx,j=0;i<gaussProcList.size();i++,j++)
    logFile<<"bound GAUSSPOINT "<<j<<": "<<gaussProcList[i]<<endl;
#endif

  MPI_Allgatherv(&gaussProcList[idx],recvCounts[rank],MPI_INT,
		 &allOtherProcLists[0],recvCounts,displs,MPI_INT,
		 MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"********** global gauss -> processor list ************"<<endl;
  for(int i=0;i<allOtherProcLists.size();i++)
    logFile<<i<<".) "<<boundGaussPoints[i].getGlobalID()<<": "
	   <<allOtherProcLists[i]<<endl;
#endif

  /**********************************************************************/
  // Erase the non locally stored boundary Gauss points and set up a
  // connectivity list of global locally stored Gauss points.
  globalLocalIdx.assign(globalLocalIdx.size(),-1);

  int m=0;

  // Loop over all global Gauss points.
  for(int i=0;i<allOtherProcLists.size();i++) {

    if(allOtherProcLists[i] != rank)
      boundGaussPoints.erase(boundGaussPoints.begin()+m);

    else {
      globalLocalIdx[i] = m;
      m++;
    }

  }

  localBGaussPtsNum = boundGaussPoints.size();

  /*********************************************************************/
  // Compute for the local boundary Gauss points a global consecutive
  // ordering.
  MPI_Allgather(&localBGaussPtsNum,1,MPI_INT,recvCounts,1,MPI_INT,
		MPI_COMM_WORLD);

  // Loop over all processors to set the global volume Gauss point start
  // and end index locally, but excluding the boundary Gauss points.
  int startGIdx = 0;

  m=0;

  while(m < rank) {
    startGIdx += recvCounts[m];
    m++;
  }

  int endGIdx = startGIdx + localBGaussPtsNum;

  // Loop over the local Gauss points of current processor.
  for(int i=0,j=startGIdx;j<endGIdx;i++,j++)
    boundGaussPoints[i].setGlobalID(j);

#ifdef _partitioningDebugMode_
  int usedDims = (int)modelData["usedDimensions"];
  m=0;
  logFile<<"localBGaussPtsNum "<<localBGaussPtsNum
	 <<" globalBGaussPtsNum "<<globalBGaussPtsNum<<endl;
  logFile<<"********* local stored boundary Gauss points *********"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    logFile<<i<<".) "<<boundGaussPoints[i].getGlobalID()<<" coords ";
    for(int j=0;j<usedDims;j++)
      logFile<<boundGaussPoints[i].getCoord(j)<<" ";
    logFile<<endl;
  }
#endif
#ifdef _geometryDebugMode_
  logFile<<"********* global -> local connectivity list **********"<<endl;
  for(int i=0;i<globalLocalIdx.size();i++) {
    logFile<<i<<".) global: ";
    if(globalLocalIdx[i] != -1) {
      logFile<<boundGaussPoints[m].getGlobalID()<<" -> local: "
	     <<globalLocalIdx[i]<<endl;
      m++;
    }
    else
      logFile<<"-1 -> local: "<<globalLocalIdx[i]<<endl;
  }
#endif

  delete[] displs,recvCounts;
}

/************************************************************************/
/************************************************************************/
// Merge the local portions of element-Gauss point entries
void BackgroundMesh::mergeLocalElemGaussVectors(InputFileData* InputData,
                                                dbVector& localVec,
                                                int gaussEntries,
                                                dbMatrix3& globalMatrix,
                                                std::map<std::string,double>& calcData,
                                                std::map<std::string,double>& modelData,
                                                std::ofstream& logFile) {

  using namespace std;

  int globalElemNum = elementRootList.size();
  int localElemNum = elemGaussIdx.size();

  std::map<std::string,double>& backGroundMeshInfo =
    InputData->getBackGroundMeshInfo();
  int gaussPtsPerVolElem =
    backGroundMeshInfo["gaussPointsPerVolumeElement"];


  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  // establish a matrix containing for each process a vector with only
  // those global element indices the process locally deals with

  int localSize,globalSize;
  localSize = globalUnsortedElemIdx.size();
  MPI_Allreduce(&localSize,&globalSize,1,
		MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  intVector recvCounts(size);
  intVector displs(size);

  int sendCount;

  int n;


#ifdef _geometryDebugMode_
  logFile<<"globalSize="<<globalSize<<endl;
  n=0;
  logFile<<"*************** local entries vector ***************"<<endl;
  for(int i=0;i<localElemNum;i++) {
    logFile<<i<<".) element "<<endl;
    for(int j=0;j<gaussPtsPerVolElem;j++) {
      logFile<<j<<".) GPoint: ";
      for(int k=0;k<gaussEntries;k++) {
	logFile<<localVec[n]<<" ";
	n++;
      }
      logFile<<endl;
    }
  }
#endif

  // not determined yet
  if(globalSize == 0) {

    intVector localIdx(localElemNum);
    intVector globalIdx(globalElemNum);

    n=0;

    // loop over all elements: oldGlobalIdx 'i' = newLocalElemIdx[i]
    for(int i=0;i<newLocalElemIdx.size();i++) {

      if(newLocalElemIdx[i] != -1) {

	localIdx[n] = i;
	n++;
      }

    }

    sendCount = localIdx.size();

    MPI_Allgather(&sendCount,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		  MPI_COMM_WORLD);

    displs[0] = 0;

    for(int j=1;j<size;j++)
      displs[j] = displs[j-1]+recvCounts[j-1];

    MPI_Gatherv(&localIdx[0],localIdx.size(),MPI_INT,
		&globalIdx[0],&recvCounts[0],&displs[0],
		MPI_INT,0,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
    n=0;
    logFile<<"**************** tmp global elem idx ****************"<<endl;
    if(rank == 0) {
      for(int i=0;i<size;i++) {
	logFile<<"process "<<i<<":"<<endl;
	for(int j=displs[i];j<displs[i]+recvCounts[i];j++)
	  logFile<<j<<".) elem: "<<globalIdx[j]<<endl;
      }
    }
#endif


    if(rank == 0) {
      globalUnsortedElemIdx.resize(size);

      for(int i=0;i<size;i++) {

	globalUnsortedElemIdx[i].resize(recvCounts[i]);
	globalUnsortedElemIdx[i].assign(globalIdx.begin()+displs[i],
					globalIdx.begin()+displs[i]
					+recvCounts[i]);

      }

    }

#ifdef _geometryDebugMode_
    logFile<<"****************** global elem idx ****************"<<endl;
    if(rank == 0) {
      for(int i=0;i<globalUnsortedElemIdx.size();i++) {
	logFile<<"process "<<i<<":"<<endl;
	for(int j=0;j<globalUnsortedElemIdx[i].size();j++)
	  logFile<<j<<".) elem: "<<globalUnsortedElemIdx[i][j]<<endl;
      }
    }
#endif

  }

  /***********************************************************************/
  // get the local entries from processes to processor '0'

  dbVector globalVec;

  if(rank == 0)
    globalVec.resize(globalElemNum*gaussPtsPerVolElem*gaussEntries);

  sendCount = localVec.size();


  MPI_Allgather(&sendCount,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int j=1;j<size;j++)
    displs[j] = displs[j-1]+recvCounts[j-1];

  MPI_Gatherv(&localVec[0],recvCounts[rank],MPI_DOUBLE,
	      &globalVec[0],&recvCounts[0],&displs[0],
	      MPI_DOUBLE,0,MPI_COMM_WORLD);

  resizeArray(localVec,0);

#ifdef _geometryDebugMode_
  n=0;
  logFile<<"************* tmp global entries vector *************"<<endl;
  if(rank == 0) {
    for(int i=0;i<globalElemNum;i++) {
      logFile<<i<<".) element "<<endl;
      for(int j=0;j<gaussPtsPerVolElem;j++) {
	logFile<<j<<".) GPoint: ";
	for(int k=0;k<gaussEntries;k++) {
	  logFile<<globalVec[n]<<" ";
	  n++;
	}
	logFile<<endl;
      }
    }
  }
#endif


  /***********************************************************************/
  // assemble the matrix containing for each global element and its Gauss
  // points the entries stored the latter

  if(rank == 0) {

    if(globalMatrix.size()  == 0)
      allocateArray(globalMatrix,globalElemNum,gaussPtsPerVolElem,
		    gaussEntries);

    //clearArray(globalMatrix);

    int sIdx;
    int eIdx=0;

    // loop over all processes
    for(int i=0;i<globalUnsortedElemIdx.size();i++) {

      // loop over all processes global elem indices
      for(int j=0;j<globalUnsortedElemIdx[i].size();j++) {

	int& globalElem = globalUnsortedElemIdx[i][j];

	for(int k=0;k<gaussPtsPerVolElem;k++) {

	  sIdx = eIdx;
	  eIdx += gaussEntries;

	  globalMatrix[globalElem][k].assign(globalVec.begin()+sIdx,
					     globalVec.begin()+eIdx);

	}

      }


    }

    resizeArray(globalVec,0);

#ifdef _geometryDebugMode_
    n=0;
    logFile<<"**************** global entries vector **************"<<endl;
    for(int i=0;i<globalMatrix.size();i++) {
      logFile<<i<<".) element "<<endl;
      for(int j=0;j<globalMatrix[i].size();j++) {
	logFile<<j<<".) GPoint: ";
	for(int k=0;k<globalMatrix[i][j].size();k++) {
	  logFile<<globalMatrix[i][j][k]<<" ";
	}
	logFile<<endl;
      }
    }
#endif

  }

}

/**********************************************************************/
/**********************************************************************/
// Create a matrix containing the penalty parameters of all boundary
// Gauss points - local and otherwise.
void BackgroundMesh::getAllBoundPenaltyParameters(InputFileData* InputData,
                                                  dbVector& allPenaltyParams,
                                                  std::map<std::string,double>& modelData,
                                                  std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];

  intMatrix& deformationBGaussPtsIdx =
    getAllDefBoundGaussPtsIdx(InputData,modelData,logFile);
  intMatrix& electricBGaussPtsIdx =
    getAllElectricBoundGaussPtsIdx(InputData,modelData,logFile);
  intMatrix& depolarisationBGaussPtsIdx =
    getAllDepolarisationBoundGaussPtsIdx(InputData,modelData,logFile);

  dbVector localPenaltyParams((deformationBGaussPtsIdx.size()
			       +electricBGaussPtsIdx.size()
			       +depolarisationBGaussPtsIdx.size())*usedDOF);

  // Loop over portion of deformation boundary Gauss points.
  for(int p=0;p<deformationBGaussPtsIdx.size();p++) {
    GaussPoint& bGPoint = boundGaussPoints[deformationBGaussPtsIdx[p][0]];
    dbVector& beta = bGPoint.getPenaltyParameters();

    for(int k=0;k<defDOF;k++)
      localPenaltyParams[p*usedDOF+k] = beta[k];

  }

  int electricStartDOF = (int)modelData["electricPotentialDOFID"];
  int electricEndDOF =
    electricStartDOF + (int)modelData["electricDegreesOfFreedom"];

  // Loop over portion of electric boundary Gauss points.
  for(int p=0;p<electricBGaussPtsIdx.size();p++) {
    GaussPoint& bGPoint = boundGaussPoints[electricBGaussPtsIdx[p][0]];
    dbVector& beta = bGPoint.getPenaltyParameters();

    for(int k=electricStartDOF;k<electricEndDOF;k++)
      localPenaltyParams[p*usedDOF+k] = beta[k];

  }

  int depolarisationStartDOF = (int)modelData["depolarisationDOFID"];
  int depolarisationEndDOF =
    depolarisationStartDOF + (int)modelData["depolarisationDegreesOfFreedom"];

  // Loop over portion of depolarisation boundary Gauss points.
  for(int p=0;p<depolarisationBGaussPtsIdx.size();p++) {
    GaussPoint& bGPoint = boundGaussPoints[depolarisationBGaussPtsIdx[p][0]];
    dbVector& beta = bGPoint.getPenaltyParameters();

    for(int k=depolarisationStartDOF;k<depolarisationEndDOF;k++)
      localPenaltyParams[p*usedDOF+k] = beta[k];

  }

#ifdef _geometryDebugMode_
  int m=0;
  logFile<<"****************************************************"<<endl;
  logFile<<"********** local bound penalty parameters **********"<<endl;
  for(int i=0;i<localPenaltyParams.size();i+=usedDOF) {
    logFile<<m<<".) ";
    for(int k=0;k<usedDOF;k++)
      logFile<<localPenaltyParams[i+k]<<" ";
    logFile<<endl;
    m++;
  }
#endif

  intVector recvCounts(size);
  intVector displs(size);

  int sendBuf = localPenaltyParams.size();
  MPI_Allgather(&sendBuf,1,MPI_INT,&recvCounts[0],1,MPI_INT,MPI_COMM_WORLD);

  int numOfGlobalParams;
  MPI_Allreduce(&sendBuf,&numOfGlobalParams,1,MPI_INT,MPI_SUM,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1] + recvCounts[i-1];

  allPenaltyParams.resize(numOfGlobalParams);

  MPI_Allgatherv(&localPenaltyParams[0],recvCounts[rank],MPI_DOUBLE,
		 &allPenaltyParams[0],&recvCounts[0],&displs[0],MPI_DOUBLE,
		 MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  int n=0;
  logFile<<"****************************************************"<<endl;
  logFile<<"*********** all bound penalty parameters ***********"<<endl;
  for(int i=0;i<allPenaltyParams.size();i+=usedDOF) {
    logFile<<n<<".) ";
    for(int k=0;k<usedDOF;k++)
      logFile<<allPenaltyParams[i+k]<<" ";
    logFile<<endl;
    n++;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Create a vector containing the history variables of all Gauss points.
void BackgroundMesh::getAllGPointHistoryVariables(InputFileData* InputData,
                                                  dbVector& allHistoryVariables,
						  int& historySize,
                                                  std::map<std::string,double>& modelData,
                                                  std::ofstream& logFile) {

  using namespace std;

  int n=0;

  dbMatrix dummy = gaussPoints[0].getPlasticityHistory();
  historySize = dummy.size()*dummy[0].size();

  // volume Gauss point history assembling

  dbVector localHistoryVariables(gaussPoints.size()*historySize);

  // loop over all local volume Gauss points
  for(int i=0;i<gaussPoints.size();i++) {

    dbMatrix& history = gaussPoints[i].getPlasticityHistory();

    for(int r=0;r<history.size();r++) {

      for(int s=0;s<history[r].size();s++) {

	localHistoryVariables[n] = history[r][s];
	n++;

      }

    }

  }

#ifdef _calculationDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"******** local volume history variables **************"<<endl;
  for(int i=0;i<gaussPoints.size();i++) {
    logFile<<i<<".) GAUSSPOINT "<<gaussPoints[i].getGlobalID()<<" : "
	   <<endl;
    dbMatrix& history = gaussPoints[i].getPlasticityHistory();
    for(int r=0;r<history.size();r++) {
      for(int s=0;s<history[r].size();s++)
	logFile<<history[r][s]<<" ";
      logFile<<endl;
    }
  }
#endif

  // ====================================================================
  // Assemble the volume Gauss point history variables from all
  // processors.

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  intVector recvCounts(size);
  intVector displs(size);

  int sendBuf = localHistoryVariables.size();
  MPI_Allgather(&sendBuf,1,MPI_INT,&recvCounts[0],1,MPI_INT,MPI_COMM_WORLD);

  int numOfGlobalData;
  MPI_Allreduce(&sendBuf,&numOfGlobalData,1,MPI_INT,MPI_SUM,
		MPI_COMM_WORLD);
  allHistoryVariables.resize(numOfGlobalData);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1] + recvCounts[i-1];

  // process '0' acts as root
  MPI_Gatherv(&localHistoryVariables[0],recvCounts[rank],MPI_DOUBLE,
	      &allHistoryVariables[0],&recvCounts[0],&displs[0],MPI_DOUBLE,
	      0,MPI_COMM_WORLD);


#ifdef _calculationDebugMode_
  n=0;
  logFile<<"******************************************************"<<endl;
  logFile<<"***** all volume Gauss history variables *************"<<endl;
  for(int i=0;i<allHistoryVariables.size();i+=historySize) {
    logFile<<n<<".) ";
    for(int k=0;k<historySize;k++)
      logFile<<allHistoryVariables[i+k]<<" ";
    logFile<<endl;
    n++;
  }
#endif

  /**********************************************************************/
  // boundary Gauss point history assembling

  localHistoryVariables.resize(boundGaussPoints.size()*historySize);
  n=0;

  // loop over all local boundary Gauss points
  for(int i=0;i<boundGaussPoints.size();i++) {

    dbMatrix& history = boundGaussPoints[i].getPlasticityHistory();

    for(int r=0;r<history.size();r++) {

      for(int s=0;s<history[r].size();s++) {

	localHistoryVariables[n] = history[r][s];
	n++;

      }

    }

  }

#ifdef  _calculationDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"******* local boundary history variables *************"<<endl;
  for(int i=0;i<boundGaussPoints.size();i++) {
    logFile<<i<<".) BOUND-GAUSSPOINT "<<boundGaussPoints[i].getGlobalID()
	   <<" : "<<endl;
    dbMatrix& history = boundGaussPoints[i].getPlasticityHistory();
    for(int r=0;r<history.size();r++) {
      for(int s=0;s<history[r].size();s++)
	logFile<<history[r][s]<<" ";
      logFile<<endl;
    }
  }
#endif

  // ====================================================================
  // Assemble the boundary Gauss point history variables from all
  // processors

  sendBuf = localHistoryVariables.size();
  MPI_Allgather(&sendBuf,1,MPI_INT,&recvCounts[0],1,MPI_INT,MPI_COMM_WORLD);

  int numOfBGlobalData;
  MPI_Allreduce(&sendBuf,&numOfBGlobalData,1,MPI_INT,MPI_SUM,
		MPI_COMM_WORLD);
  allHistoryVariables.resize(numOfGlobalData+numOfBGlobalData);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1] + recvCounts[i-1];

  // process '0' acts as root
  MPI_Gatherv(&localHistoryVariables[0],recvCounts[rank],MPI_DOUBLE,
	      &allHistoryVariables[numOfGlobalData],&recvCounts[0],
	      &displs[0],MPI_DOUBLE,0,MPI_COMM_WORLD);


#ifdef _calculationDebugMode_
  n=0;
  logFile<<"****************************************************"<<endl;
  logFile<<"*********** all Gauss history variables ************"<<endl;
  for(int i=0;i<numOfGlobalData;i+=historySize) {
    logFile<<n<<".) GAUSSPOINT :";
    for(int k=0;k<historySize;k++)
      logFile<<allHistoryVariables[i+k]<<" ";
    logFile<<endl;
    n++;
  }
  n=0;
  logFile<<"***************************************************"<<endl;
  for(int i=numOfGlobalData;i<allHistoryVariables.size();i+=historySize) {
    logFile<<n<<".) BOUND-GAUSSPOINT :";
    for(int k=0;k<historySize;k++)
      logFile<<allHistoryVariables[i+k]<<" ";
    logFile<<endl;
    n++;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where a displacement boundary
// condition is applied (point,line,surface).
intMatrix& BackgroundMesh::getAllDispBoundGaussPtsIdx(InputFileData* InputData,
                                                      std::map<std::string,double>& modelData,
                                                      std::ofstream& logFile) {


  using namespace std;

  if(allDisplacementBoundGaussPtsIdx.size() > 0)
    return allDisplacementBoundGaussPtsIdx;

  else {

    if((bool)modelData["displacementConstraint"]) {

      pushBackVector(allDisplacementBoundGaussPtsIdx,
		     surfaceDispBoundGaussPtsIdx);
      pushBackVector(allDisplacementBoundGaussPtsIdx,
		     lineDispBoundGaussPtsIdx);
      pushBackVector(allDisplacementBoundGaussPtsIdx,
		     pointDispBoundPtcleIdx);
    }

#ifdef _geometryDebugMode_
    logFile<<"***************************************************"<<endl;
    logFile<<"****** all displacement gauss point indices *******"<<endl;
    for(int i=0;i<allDisplacementBoundGaussPtsIdx.size();i++) {
      logFile<<i<<".) local ID "<<allDisplacementBoundGaussPtsIdx[i][0]
	     <<": ";
      for(int j=1;j<allDisplacementBoundGaussPtsIdx[i].size();j++)
	logFile<<allDisplacementBoundGaussPtsIdx[i][j]<<" ";
      logFile<<endl;
    }
#endif

    return allDisplacementBoundGaussPtsIdx;
  }

}

/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where a rotation boundary
// condition is applied (point,line,surface).
intMatrix& BackgroundMesh::getAllRotBoundGaussPtsIdx(InputFileData* InputData,
                                                     std::map<std::string,double>& modelData,
                                                     std::ofstream& logFile) {


  using namespace std;

  if(allRotationBoundGaussPtsIdx.size() > 0)
    return allRotationBoundGaussPtsIdx;

  else {

    if((bool)modelData["rotationConstraint"]) {

      pushBackVector(allRotationBoundGaussPtsIdx,
		     surfaceRotBoundGaussPtsIdx);
      pushBackVector(allRotationBoundGaussPtsIdx,
		     lineRotBoundGaussPtsIdx);
      pushBackVector(allRotationBoundGaussPtsIdx,
		     pointRotBoundPtcleIdx);
    }

#ifdef _geometryDebugMode_
    logFile<<"***************************************************"<<endl;
    logFile<<"******** all rotation gauss point indices *********"<<endl;
    for(int i=0;i<allRotationBoundGaussPtsIdx.size();i++) {
      logFile<<i<<".) local ID "<<allRotationBoundGaussPtsIdx[i][0]
	     <<": ";
      for(int j=1;j<allRotationBoundGaussPtsIdx[i].size();j++)
	logFile<<allRotationBoundGaussPtsIdx[i][j]<<" ";
      logFile<<endl;
    }
#endif

    return allRotationBoundGaussPtsIdx;
  }

}
/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where an electric Dirichlet boundary
// condition is applied.
intMatrix& BackgroundMesh::
getAllElectricBoundGaussPtsIdx(InputFileData* InputData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile) {


  using namespace std;

  if(allElectricBoundGaussPtsIdx.size() > 0)
    return allElectricBoundGaussPtsIdx;

  else {

    if((bool)modelData["electricConstraint"]) {

      pushBackVector(allElectricBoundGaussPtsIdx,
		     surfaceElectricBoundGaussPtsIdx);
      pushBackVector(allElectricBoundGaussPtsIdx,
		     lineElectricBoundGaussPtsIdx);
      pushBackVector(allElectricBoundGaussPtsIdx,
		     pointElectricBoundPtcleIdx);
    }

#ifdef _geometryDebugMode_
    logFile<<"***************************************************"<<endl;
    logFile<<"******** all electric gauss point indices *********"<<endl;
    for(int i=0;i<allElectricBoundGaussPtsIdx.size();i++) {
      logFile<<i<<".) local ID "<<allElectricBoundGaussPtsIdx[i][0]
	     <<": ";
      for(int j=1;j<allElectricBoundGaussPtsIdx[i].size();j++)
	logFile<<allElectricBoundGaussPtsIdx[i][j]<<" ";
      logFile<<endl;
    }
#endif

    return allElectricBoundGaussPtsIdx;
  }

}

/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where an depolarisation Dirichlet boundary
// condition is applied.
intMatrix& BackgroundMesh::
getAllDepolarisationBoundGaussPtsIdx(InputFileData* InputData,
                                     std::map<std::string,double>& modelData,
                                     std::ofstream& logFile) {


  using namespace std;

  if(allDepolarisationBoundGaussPtsIdx.size() > 0)
    return allDepolarisationBoundGaussPtsIdx;

  else {

    if((bool)modelData["depolarisationConstraint"]) {

      pushBackVector(allDepolarisationBoundGaussPtsIdx,
		     surfaceDepolarisationBoundGaussPtsIdx);
      pushBackVector(allDepolarisationBoundGaussPtsIdx,
		     lineDepolarisationBoundGaussPtsIdx);
      pushBackVector(allDepolarisationBoundGaussPtsIdx,
		     pointDepolarisationBoundPtcleIdx);
    }

#ifdef _geometryDebugMode_
    logFile<<"***************************************************"<<endl;
    logFile<<"******** all depolarisation gauss point indices *********"<<endl;
    for(int i=0;i<allDepolarisationBoundGaussPtsIdx.size();i++) {
      logFile<<i<<".) local ID "<<allDepolarisationBoundGaussPtsIdx[i][0]
	     <<": ";
      for(int j=1;j<allDepolarisationBoundGaussPtsIdx[i].size();j++)
	logFile<<allDepolarisationBoundGaussPtsIdx[i][j]<<" ";
      logFile<<endl;
    }
#endif

    return allDepolarisationBoundGaussPtsIdx;
  }

}

/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where a micro boundary
// condition is applied (point,line,surface).
intMatrix& BackgroundMesh::getAllMicroBoundGaussPtsIdx(InputFileData* InputData,
                                                       std::map<std::string,double>& modelData,
                                                       std::ofstream& logFile) {


  using namespace std;

  if(allMicroBoundGaussPtsIdx.size() > 0)
    return allMicroBoundGaussPtsIdx;

  else {

    if((bool)modelData["microConstraint"]) {

      pushBackVector(allMicroBoundGaussPtsIdx,
		     surfaceMicroBoundGaussPtsIdx);
      pushBackVector(allMicroBoundGaussPtsIdx,
		     lineMicroBoundGaussPtsIdx);
      //pushBackVector(allMicroBoundGaussPtsIdx,
      //		     pointDispBoundPtcleIdx);
    }

#ifdef _geometryDebugMode_
    logFile<<"***************************************************"<<endl;
    logFile<<"****** all micro gauss point indices *******"<<endl;
    for(int i=0;i<allMicroBoundGaussPtsIdx.size();i++) {
      logFile<<i<<".) local ID "<<allMicroBoundGaussPtsIdx[i][0]
	     <<": ";
      for(int j=1;j<allMicroBoundGaussPtsIdx[i].size();j++)
	logFile<<allMicroBoundGaussPtsIdx[i][j]<<" ";
      logFile<<endl;
    }
#endif

    return allMicroBoundGaussPtsIdx;
  }

}

/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where any deformation boundary
// condition is applied.
intMatrix& BackgroundMesh::getAllDefBoundGaussPtsIdx(InputFileData* InputData,
                                                     std::map<std::string,double>& modelData,
                                                     std::ofstream& logFile) {


  using namespace std;

  if(allDeformationBoundGaussPtsIdx.size() > 0)
    return allDeformationBoundGaussPtsIdx;

  else {

    if((bool)modelData["displacementConstraint"]) {

      pushBackVector(allDeformationBoundGaussPtsIdx,
		     surfaceDispBoundGaussPtsIdx);
      pushBackVector(allDeformationBoundGaussPtsIdx,
		     lineDispBoundGaussPtsIdx);
      pushBackVector(allDeformationBoundGaussPtsIdx,
		     pointDispBoundPtcleIdx);
    }

    for(int i=0;i<allDeformationBoundGaussPtsIdx.size();i++)

      allocateArray(allDeformationBoundGaussPtsIdx[i],4);


    int position;
    intVector indexVec(4);

    if((bool)modelData["rotationConstraint"]) {

      // surface rotation constraints
      for(int i=0;i<surfaceRotBoundGaussPtsIdx.size();i++) {

	position = findIntMatPos(surfaceRotBoundGaussPtsIdx[i][0],0,
				 allDeformationBoundGaussPtsIdx.size(),0,
				 allDeformationBoundGaussPtsIdx);

	// displacement and rotation constraints are set (indexVec[3] = 2)
	if(position != -1)
	  allDeformationBoundGaussPtsIdx[i][3] = 2;

	// only rotation constraints are set (indexVec[3] = 1)
	else {

	  for(int j=0;j<3;j++)
	    indexVec[j] = surfaceRotBoundGaussPtsIdx[i][j];

	  indexVec[3] = 1;
	  pushBackVector(allDeformationBoundGaussPtsIdx,indexVec);

	}

      }

      // line rotation constraints
      for(int i=0;i<lineRotBoundGaussPtsIdx.size();i++) {

	position = findIntMatPos(lineRotBoundGaussPtsIdx[i][0],0,
				 allDeformationBoundGaussPtsIdx.size(),0,
				 allDeformationBoundGaussPtsIdx);

	if(position != -1)
	  allDeformationBoundGaussPtsIdx[i][3] = 2;

	else {

	  for(int j=0;j<3;j++)
	    indexVec[j] = lineRotBoundGaussPtsIdx[i][j];

	  indexVec[3] = 1;
	  pushBackVector(allDeformationBoundGaussPtsIdx,indexVec);

	}

      }

      // point rotation constraints
      for(int i=0;i<pointRotBoundPtcleIdx.size();i++) {

	position = findIntMatPos(pointRotBoundPtcleIdx[i][0],0,
				 allDeformationBoundGaussPtsIdx.size(),0,
				 allDeformationBoundGaussPtsIdx);

	if(position != -1)
	  allDeformationBoundGaussPtsIdx[i][3] = 2;

	else {

	  for(int j=0;j<3;j++)
	    indexVec[j] = pointRotBoundPtcleIdx[i][j];

	  indexVec[3] = 1;
	  pushBackVector(allDeformationBoundGaussPtsIdx,indexVec);

	}

      }

    }

#ifdef _geometryDebugMode_
    logFile<<"***************************************************"<<endl;
    logFile<<"****** all deformation gauss point indices ********"<<endl;
    for(int i=0;i<allDeformationBoundGaussPtsIdx.size();i++) {
      logFile<<i<<".) "<<allDeformationBoundGaussPtsIdx[i][0]<<": ";
      for(int j=1;j<allDeformationBoundGaussPtsIdx[i].size();j++)
	logFile<<allDeformationBoundGaussPtsIdx[i][j]<<" ";
      logFile<<endl;
    }
#endif

    return allDeformationBoundGaussPtsIdx;
  }

}

/************************************************************************/
/************************************************************************/
// Return all boundary Gauss indices where any force boundary condition
// is applied (point,line,surface).
intMatrix& BackgroundMesh::getAllForceBoundGaussPtsIdx(InputFileData* InputData,
                                                       std::map<std::string,double>& modelData,
                                                       std::ofstream& logFile) {


  using namespace std;

  if(allForceBoundGaussPtsIdx.size() > 0)
    return allForceBoundGaussPtsIdx;

  else {

    if((bool)modelData["tractionLoad"])
      pushBackVector(allForceBoundGaussPtsIdx,
		     tractionBoundGaussPtsIdx);

    if((bool)modelData["surfacePressureLoad"])
      pushBackVector(allForceBoundGaussPtsIdx,
		     surfacePressureBoundGaussPtsIdx);

    if((bool)modelData["lineForceLoad"])
      pushBackVector(allForceBoundGaussPtsIdx,
		     lineForceBoundGaussPtsIdx);

    if((bool)modelData["pointForceLoad"])
      pushBackVector(allForceBoundGaussPtsIdx,
		     pointForceBoundPtcleIdx);

    if((bool)modelData["surfaceMomentLoad"])
      pushBackVector(allForceBoundGaussPtsIdx,
		     surfaceMomentBoundGaussPtsIdx);

    if((bool)modelData["surfaceElectricChargeLoad"])
      pushBackVector(allForceBoundGaussPtsIdx,
		     surfaceElectricChargeBoundGaussPtsIdx);

  }

#ifdef _geometryDebugMode_
  logFile<<"***************************************************"<<endl;
  logFile<<"****** all force boundary gauss point indices *****"<<endl;
  for(int i=0;i<allForceBoundGaussPtsIdx.size();i++) {
    logFile<<i<<".) "<<allForceBoundGaussPtsIdx[i][0]<<": ";
    for(int j=1;j<allForceBoundGaussPtsIdx[i].size();j++)
      logFile<<allForceBoundGaussPtsIdx[i][j]<<" ";
    logFile<<endl;
  }
#endif

  return allForceBoundGaussPtsIdx;
}

/**********************************************************************/
/**********************************************************************/
// Return all surface integration Gauss point indices.
intMatrix& BackgroundMesh::getAllSurfaceGaussPtsIdx(std::ofstream& logFile) {

  using namespace std;

  if(allSurfaceBoundGaussPtsIdx.size() > 0)
    return allSurfaceBoundGaussPtsIdx;

  else {

    intMatrix& mat = allSurfaceBoundGaussPtsIdx;

    mat.insert(mat.begin(),surfaceDispBoundGaussPtsIdx.begin(),
	       surfaceDispBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),surfaceRotBoundGaussPtsIdx.begin(),
	       surfaceRotBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),surfaceElectricBoundGaussPtsIdx.begin(),
	       surfaceElectricBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),surfaceDepolarisationBoundGaussPtsIdx.begin(),
	       surfaceDepolarisationBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),tractionBoundGaussPtsIdx.begin(),
	       tractionBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),surfacePressureBoundGaussPtsIdx.begin(),
	       surfacePressureBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),surfaceMomentBoundGaussPtsIdx.begin(),
	       surfaceMomentBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),surfaceElectricChargeBoundGaussPtsIdx.begin(),
	       surfaceElectricChargeBoundGaussPtsIdx.end());


#ifdef _FEdebugMode_
    logFile<<"------------------------------------------------------"<<endl;
    logFile<<"------------- tmp allSurfaceGaussPtsIdx --------------"<<endl;
    for(int i=0;i<mat.size();i++)
      logFile<<i<<".) Gauss Point "<<mat[i][0]<<endl;
#endif

    removeRedundantEntries(mat,0);

#ifdef _FEdebugMode_
    logFile<<"------------------------------------------------------"<<endl;
    logFile<<"------------ final allSurfaceGaussPtsIdx -------------"<<endl;
    for(int i=0;i<mat.size();i++)
      logFile<<i<<".) Gauss Point "<<mat[i][0]<<endl;
#endif

    return mat;
  }

}

/************************************************************************/
/************************************************************************/
// Return all line integration Gauss point indices.
intMatrix& BackgroundMesh::getAllLineGaussPtsIdx(std::ofstream& logFile) {

  using namespace std;

  if(allLineBoundGaussPtsIdx.size() > 0)
    return allLineBoundGaussPtsIdx;

  else {

    intMatrix& mat = allLineBoundGaussPtsIdx;

    mat.insert(mat.begin(),lineDispBoundGaussPtsIdx.begin(),
	       lineDispBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),lineRotBoundGaussPtsIdx.begin(),
	       lineRotBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),lineElectricBoundGaussPtsIdx.begin(),
	       lineElectricBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),lineDepolarisationBoundGaussPtsIdx.begin(),
	       lineDepolarisationBoundGaussPtsIdx.end());

    mat.insert(mat.begin(),lineForceBoundGaussPtsIdx.begin(),
	       lineForceBoundGaussPtsIdx.end());

#ifdef _FEdebugMode_
    logFile<<"------------------------------------------------------"<<endl;
    logFile<<"------------- tmp allLineGaussPtsIdx -----------------"<<endl;
    for(int i=0;i<mat.size();i++)
      logFile<<i<<".) Gauss Point "<<mat[i][0]<<endl;
#endif

    removeRedundantEntries(mat,0);

#ifdef _FEdebugMode_
    logFile<<"------------------------------------------------------"<<endl;
    logFile<<"------------- final allLineGaussPtsIdx ---------------"<<endl;
    for(int i=0;i<mat.size();i++)
      logFile<<i<<".) Gauss Point "<<mat[i][0]<<endl;
#endif

    return mat;
  }

}

/************************************************************************/
/************************************************************************/
// Calculate displacements and and rotations on all gauss points
// from all processors.
void BackgroundMesh::getAllGaussDeforms(InputFileData* InputData,
                                        dbVector& deformations,
                                        std::map<std::string,double>& modelData,
                                        std::ofstream& logFile) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  intVector idx(localGaussPtsNum);
  dbVector localDefs(localGaussPtsNum*usedDOF);

  /*********************************************************************/
  // Loop over all local gauss points to determine their deformations.
  for(int i=0;i<localGaussPtsNum;i++) {
    intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
    idx[i] = gaussPoints[i].getGlobalID();

    // Loop over its supporting particles.
    for(int j=0;j<suppPtcls.size();j++)

      // Loop over all degrees of freedom.
      for(int k=0;k<usedDOF;k++)
	localDefs[i*usedDOF+k] += gaussPoints[i].getShapeFunc(j)
	  *particles[suppPtcls[j]].getDOF(k);

  }

#ifdef _postProcDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"########### local gauss pts displacements ##########"<<endl;
  double PUM;
  for(int i=0;i<localGaussPtsNum;i++) {
    intVector& suppPtcls = gaussPoints[i].getSupportPtcls();
    PUM=0;
    for(int j=0;j<suppPtcls.size();j++)
      PUM+=gaussPoints[i].getShapeFunc(j);
    logFile<<"locIdx "<<i<<" globIdx "<<idx[i]<<" ";
    for(int j=0;j<usedDOF;j++)
      logFile<<localDefs[i*usedDOF+j]<<" ";
    if(fabs(1.0-PUM) > 0.000000001)
      logFile<<"PUM = "<<PUM;
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // Assemble all local gauss points displacements from all processors.
  dbVector tmpGlobalDefs(globalGaussPtsNum*usedDOF);
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int* recvCounts = new int[size];
  int* displs = new int[size];

  // The coordinates' part.
  int sendBuf = localDefs.size();

  MPI_Allgather(&sendBuf,1,MPI_INT,recvCounts,1,MPI_INT,MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  MPI_Gatherv(&localDefs[0],localDefs.size(),MPI_DOUBLE,
	      &tmpGlobalDefs[0],recvCounts,displs,MPI_DOUBLE,
	      0,MPI_COMM_WORLD);

  // The index part.
  intVector allIdx(globalGaussPtsNum);

  MPI_Allgather(&localGaussPtsNum,1,MPI_INT,recvCounts,1,MPI_INT,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  MPI_Gatherv(&idx[0],recvCounts[rank],MPI_INT,&allIdx[0],recvCounts,
	      displs,MPI_INT,0,MPI_COMM_WORLD);

#ifdef _postProcDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"#### unsorted global gauss pts displacements #######"<<endl;
  for(int i=0;i<globalGaussPtsNum;i++) {
    logFile<<i<<" globIdx "<<allIdx[i]<<" ";
    for(int j=0;j<usedDOF;j++)
      logFile<<tmpGlobalDefs[i*usedDOF+j]<<" ";
    logFile<<endl;
  }
#endif

  // Order the coordinates according to their global indices.
  if(deformations.size() < globalGaussPtsNum*usedDOF) {
    clearArray(deformations);
    deformations.resize(globalGaussPtsNum*usedDOF);
  }
  else
    clearArray(deformations);

  for(int i=0;i<globalGaussPtsNum;i++)

    for(int j=0;j<usedDOF;j++)
      deformations[allIdx[i]*usedDOF+j] = tmpGlobalDefs[i*usedDOF+j];


#ifdef _postProcDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"###### sorted global gauss pts displacements #######"<<endl;
  for(int i=0;i<globalGaussPtsNum;i++) {
    logFile<<i<<" ";
    for(int j=0;j<usedDOF;j++)
      logFile<<deformations[i*usedDOF+j]<<" ";
    logFile<<endl;
  }
#endif

  delete[] recvCounts,displs;

}

/************************************************************************/
/************************************************************************/
// Calculate rotations of one integration point.
void BackgroundMesh::getIntPointRotation(InputFileData* InputData,
                                         GaussPoint& gPoint,
                                         dbVector& allDOF,
                                         dbVector& rotation,
                                         std::map<std::string,double>& modelData,
                                         std::ofstream& logFile) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int usedDim = (int)modelData["usedDimensions"];

  if(rotation.size() < usedDim) {
    clearArray(rotation);
    rotation.resize(usedDim);
  }
  else
    clearArray(rotation);

  if(allDOF.size() < particlesNum*usedDOF) {
    logFile <<"In BackgroundMesh::getIntPointRotation can not "
	    <<"calculate the rotation of integration point "
	    <<gPoint.getGlobalID()
	    <<",\n because the vector containing the particle DOF "
            <<"has too less entries!"<<endl;
    cerr <<"In BackgroundMesh::getIntPointRotation can not "
	 <<"calculate the rotation of integration point "
	 <<gPoint.getGlobalID()
	 <<",\n because the vector containing the particle DOF "
         <<"has too less entries!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Get for this particle all its supporting neighbours in their
  // shape function ordinates.
  intVector& suppPtcls = gPoint.getSupportPtcls();
  dbVector& suppShapes = gPoint.getShapeFuncs();

  // Loop over its supporting particles.
  for(int i=0;i<suppPtcls.size();i++)

    // Loop over all rotational degrees of freedom.
    for(int j=usedDim,k=0;j<usedDOF;j++,k++)
      rotation[k] += suppShapes[i]*allDOF[suppPtcls[i]*usedDOF+j];

#ifdef _calculationDebugMode_
#endif

}
