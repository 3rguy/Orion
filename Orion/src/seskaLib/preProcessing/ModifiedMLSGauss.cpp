#include "ModifiedMLSGauss.h"

ModifiedMLSGauss::ModifiedMLSGauss(InputFileData* InputData,
				   std::ofstream& logFile) 
    : MLSGaussIntegral(InputData,logFile),
      MLSShapeFuncSet(InputData,logFile),
      ParticleDistribution(InputData,logFile),
      BackgroundMesh(InputData,logFile) {}



/**********************************************************************/
/**********************************************************************/
// Calculate for all Gauss points and all particles a vector containing 
// the calculated shape functions and their derivatives for all 
// supporting particles.
void ModifiedMLSGauss::setAllShapeFuncs(InputFileData* InputData,
					std::map<std::string,double>& calcData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile,
					PetscViewer& viewerMPI,
					PetscViewer& viewerSEQ) {

  using namespace std;

  if((int)InputData->getValue("shapefunctionType") != 2) {
    logFile <<"In ModifiedMLSGauss::setAllShapeFuncs only shapefunction "
	    <<"type 'EFG(=2)' is supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  // Create ghost particle on the deformation boundary (u=0/r=0 only!)
  createBoundGhostPtcls(InputData,modelData,logFile,viewerMPI,viewerSEQ);

  // Determine for all Gauss points and particles their supporting 
  // boundary ghost particles.
  setGhostPtcleConn(InputData,modelData,logFile,viewerMPI,viewerSEQ);

  /**********************************************************************/
  // Calculate for all particles a vector containing the calculated
  // shape functions and their derivations for all supporting particles.
  if((int)InputData->getValue("plotOnParticles") == 1)

    setShapeFuncsOnPtcls(InputData,modelData,logFile,
			 viewerMPI,viewerSEQ);

  // Calculate for all inner and boundary gauss points a vector 
  // containing the calculated shape functions and their derivations
  // for all supporting particles.
  setShapeFuncsOnGauss(InputData,modelData,logFile,viewerMPI,viewerSEQ);
  setShapeFuncsOnBGauss(InputData,modelData,logFile,viewerMPI,viewerSEQ);


  /**********************************************************************/
  // reset the original support lists

  particles.resize(particlesNum,Particle(usedDOF));

  // Loop over all particles
  for (int i=0;i<particles.size();i++) {

    intVector& suppGhostPtcls = particles[i].getSupportBoundGhostPtcls();
    suppGhostPtcls.resize(0);
  }

  // Loop over all local Gauss points.
  for (int i=0;i<localGaussPtsNum;i++) {

    GaussPoint& gPoint = gaussPoints[i];
    intVector& suppGhostPtcls = gPoint.getSupportBoundGhostPtcls();
    suppGhostPtcls.resize(0);
  }

  // Loop over all local boundary Gauss points.
  for (int i=0;i<localBGaussPtsNum;i++) {

    GaussPoint& gPoint = boundGaussPoints[i];
    intVector& suppGhostPtcls = gPoint.getSupportBoundGhostPtcls();
    suppGhostPtcls.resize(0);
  }

}

/**********************************************************************/
/**********************************************************************/
// Create ghost particle on the deformation boundary (u=0/r=0 only!)
void  ModifiedMLSGauss::createBoundGhostPtcls(InputFileData* InputData,
					      std::map<std::string,double>& modelData,
					      std::ofstream& logFile,
					      PetscViewer& viewerMPI,
					      PetscViewer& viewerSEQ) {

  using namespace std; 

  std::map<std::string,double>& backGroundMeshInfo =
    InputData->getBackGroundMeshInfo();

  int nodesPerSurfaceElem = (int)backGroundMeshInfo["nodesPerSurfaceElement"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  double multiplier = InputData->getValue("influenceRadiusMultiplier");

  dbMatrix& dispBoundConds = InputData->getSurfaceDispBoundConds();
  dbMatrix& rotBoundConds = InputData->getSurfaceRotBoundConds();

  dbMatrix defBoundConds = dispBoundConds;
  defBoundConds.insert(defBoundConds.begin(),rotBoundConds.begin(),
		       rotBoundConds.end());

  // Loop over all boundary elements applied a surface load.
  int elemID,idx,pos;
  double factor;
  intVector boundaryElemIdx;

  ElementTemplate* FEVolumeSet = new Cube8ElementTemplate();
  ElementTemplate* FESurfaceSet = new Rect4ElementTemplate();

  double surface;
  intVector nodes(nodesPerSurfaceElem);
  dbVector motherNormal(3);
  
  int currentIdx = particlesNum;
  int numOfGhostLayers = 3;

  int numOfNewPtcls = defBoundConds.size()*nodesPerSurfaceElem*numOfGhostLayers;

  particles.resize(particlesNum + numOfNewPtcls,Particle(usedDOF));

  // deformation boundary constraints
  for(int i=0;i<defBoundConds.size();i++) {

    // Check if a boundary condition is already applied to current element.
    elemID = (int)defBoundConds[i][0]-1;

    pos = findIntVecPos(elemID,0,boundaryElemIdx.size(),
			boundaryElemIdx);

    // no boundary condition is already applied
    if(pos == -1) {

      // loop over element's nodes
      for(int j=0;j<nodesPerSurfaceElem;j++)
	nodes[j] = newIdx[(int)defBoundConds[i][j+1]-1]+1;
      
      // loop over element's nodes
      for(int j=0;j<nodesPerSurfaceElem;j++) {

	idx = nodes[j]-1;
	Particle& motherPtcle = particles[idx];
	dbVector& motherCoords = motherPtcle.getCoords();
	
	// determine the surface normal vector
	FEVolumeSet->surfaceMetric(particles,nodes,
				   FESurfaceSet->nodalCoords[j],
				   motherNormal,surface,logFile);     

	dbVector& motherRadii = motherPtcle.getRadii();

	// loop over the desired ghost particle layers
	for(int k=1;k<=numOfGhostLayers;k++) {

	  factor = k/multiplier;

	  Particle& ptcle = particles[currentIdx];
	  int& motherIdx = ptcle.getMotherPtcle();
	  dbVector& coords = ptcle.getCoords();
	  dbVector& radii = ptcle.getRadii();

	  motherIdx = idx;
	  coords.resize(motherCoords.size());
	  radii = motherRadii;

	  for(int l=0;l<motherCoords.size();l++) 

	    coords[l] = motherCoords[l] + 
	      motherNormal[l]*motherRadii[l]*factor;


	  currentIdx++;
	}

	pushBackVector(boundaryElemIdx,elemID);
      }
      
    }

  }

  particles.resize(currentIdx,Particle(usedDOF));

  delete FEVolumeSet,FESurfaceSet;

#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"************** boundary ghost particles **************"<<endl;
  for(int i=particlesNum;i<particles.size();i++) {
    Particle& ptcle = particles[i];
    int& motherPtcle = ptcle.getMotherPtcle();
    dbVector& motherCoords = particles[motherPtcle].getCoords();
    dbVector& coords = ptcle.getCoords();
    dbVector& radii = ptcle.getRadii();
    logFile<<"BoundGhostPtcle "<<i<<" coords: ";
    for(int j=0;j<coords.size();j++)
      logFile<<coords[j]<<" ";
    logFile<<"radii: ";
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" ";
    logFile<<"motherPtcle "<<motherPtcle<<" coords: ";
    for(int j=0;j<motherCoords.size();j++)
      logFile<<motherCoords[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine for all Gauss points and particles their supporting 
// boundary ghost particles.
void ModifiedMLSGauss::setGhostPtcleConn(InputFileData* InputData,
					 std::map<std::string,double>& modelData,
					 std::ofstream& logFile,
					 PetscViewer& viewerMPI,
					 PetscViewer& viewerSEQ) {

  using namespace std;

  int gaussPtcleConnect = 
    (int)InputData->getValue("gaussParticleConnectivity");

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  dbVector radius(3);
  int supportCounts;


  int maxSupport = 0;
  int minSupport = 0;

  // Loop over all local gauss points.
  for (int i=0;i<localGaussPtsNum;i++) {

    GaussPoint& gPoint = gaussPoints[i];
    supportCounts = 0;
    intVector& suppGhostPtcls = gPoint.getSupportBoundGhostPtcls();

    // Loop over all ghost particles.
    for(int j=particlesNum;j<particles.size();j++) {
      Particle& ptcle = particles[j];
      
      radius[0] = fabs(ptcle.getCoord(0)-gPoint.getCoord(0));
      radius[1] = fabs(ptcle.getCoord(1)-gPoint.getCoord(1));
      radius[2] = fabs(ptcle.getCoord(2)-gPoint.getCoord(2));
      
      if(radius[0] <= ptcle.getRadius(0) && 
	 radius[1] <= ptcle.getRadius(1) &&
	 radius[2] <= ptcle.getRadius(2)) {

	if(suppGhostPtcls.size() > supportCounts)
	  suppGhostPtcls[supportCounts] = j;
	else
	  suppGhostPtcls.push_back(j);

	supportCounts++;
      }
    }

    if(supportCounts < suppGhostPtcls.size())
      suppGhostPtcls.resize(supportCounts);

  }


#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"********* supporting boundary ghost particles ********"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    GaussPoint& gPoint = gaussPoints[i];
    logFile<<i<<".) GAUSS POINT "<<gPoint.getGlobalID()<<": ";
    intVector& suppGhostPtcls = gPoint.getSupportBoundGhostPtcls();
    for(int j=0;j<suppGhostPtcls.size();j++) {
      logFile<<suppGhostPtcls[j]<<" ";
    }
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // Loop over all boundary gauss points.
  for (int i=0;i<localBGaussPtsNum;i++) {
    GaussPoint& gPoint = boundGaussPoints[i];

    supportCounts = 0;
    intVector& suppGhostPtcls = gPoint.getSupportBoundGhostPtcls();

    // Loop over all ghost particles.
    for(int j=particlesNum;j<particles.size();j++) {
      Particle& ptcle = particles[j];

      radius[0] = fabs(ptcle.getCoord(0)-gPoint.getCoord(0));
      radius[1] = fabs(ptcle.getCoord(1)-gPoint.getCoord(1));
      radius[2] = fabs(ptcle.getCoord(2)-gPoint.getCoord(2));

      if(radius[0] <= ptcle.getRadius(0) && 
	 radius[1] <= ptcle.getRadius(1) &&
	 radius[2] <= ptcle.getRadius(2)) {

	if(suppGhostPtcls.size() > supportCounts)
	  suppGhostPtcls[supportCounts] = j;
	else
	  suppGhostPtcls.push_back(j);

	supportCounts++;
      }
    }

    if(supportCounts < suppGhostPtcls.size())
      suppGhostPtcls.resize(supportCounts);

  }

#ifdef _geometryDebugMode_
  logFile<<"******************************************************"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    logFile<<"Bound GAUSSPOINT "<<i<<": ";
    GaussPoint& gPoint = boundGaussPoints[i];
    intVector& suppGhostPtcls = gPoint.getSupportBoundGhostPtcls();
    for(int j=0;j<suppGhostPtcls.size();j++) {
      logFile<<suppGhostPtcls[j]<<" ";
    }
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // Loop over all particles.
  for(int i=0;i<particlesNum;i++) {

    supportCounts = 0;
    intVector& suppPtcls = particles[i].getSupportPtcls();
    intVector& suppGhostPtcls = particles[i].getSupportBoundGhostPtcls();
    
    // Loop over all ghost particles.
    for(int j=particlesNum;j<particles.size();j++) {
      
      radius[0] = fabs(particles[j].getCoord(0)-particles[i].getCoord(0));
      radius[1] = fabs(particles[j].getCoord(1)-particles[i].getCoord(1));
      radius[2] = fabs(particles[j].getCoord(2)-particles[i].getCoord(2));
      
      if(radius[0] <= particles[j].getRadius(0) && 
	 radius[1] <= particles[j].getRadius(1) &&
	 radius[2] <= particles[j].getRadius(2)) {
	
	if(suppGhostPtcls.size() > supportCounts)
	  suppGhostPtcls[supportCounts] = j;
	else
	  suppGhostPtcls.push_back(j);
	
	supportCounts++;
      }
      
    }

    if(supportCounts < suppGhostPtcls.size())
      suppGhostPtcls.resize(supportCounts);

    if(supportCounts + suppPtcls.size() > maxSupport)
      maxSupport = supportCounts + suppPtcls.size();
  
    if(i == 0)
      minSupport = supportCounts + suppPtcls.size();

    else if(supportCounts + suppPtcls.size() < minSupport)
      minSupport = supportCounts + suppPtcls.size();

  }

  
  logFile<<"minPtcleSupport = "<<minSupport<<" (incl. ghost particles)"<<endl;
  logFile<<"maxPtcleSupport = "<<maxSupport<<" (incl. ghost particles)"<<endl;
  cout<<"ptcleSupport: "<<minSupport<<" - "<<maxSupport<<" (incl. ghost particles)"<<endl;
  
#ifdef _geometryDebugMode_
  logFile<<"********** neighbour particles of all particles *********"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppGhostPtcls = particles[i].getSupportBoundGhostPtcls();
    logFile<<"Particle "<<i<<": ";
   for(int j=0;j<suppGhostPtcls.size();j++) {
      logFile<<suppGhostPtcls[j]<<" ";
    }
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all gauss points a vector containing the calculated
// shape functions and their derivations for all supporting particles.
void ModifiedMLSGauss::setShapeFuncsOnGauss(InputFileData* InputData,
					    std::map<std::string,double>& modelData,
					    std::ofstream& logFile,
					    PetscViewer& viewerMPI,
					    PetscViewer& viewerSEQ) {

  using namespace std;

  int derivationOrder = 
    (int)modelData["shapesDerivationOrderOnIntPoints"];

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  dbMatrix secondDerivShapes(6);

  int allSupportSize;
  int supportSize;
  intVector allSupportPtcls;
  
  // Loop over all local gauss points to set their shape function 
  // vectors.
  int rank; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int basisTermNum,pos;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"******************* gausspoint shapes ********************"<<endl;
#endif

  for(int i=0;i<localGaussPtsNum;i++) {
    intVector& suppGhostPtcls = gaussPoints[i].getSupportBoundGhostPtcls();
    intVector& suppPtcls = gaussPoints[i].getSupportPtcls();

    allSupportPtcls = suppPtcls;
    allSupportPtcls.insert(allSupportPtcls.end(),suppGhostPtcls.begin(),
			   suppGhostPtcls.end());

    supportSize = suppPtcls.size();
    allSupportSize = allSupportPtcls.size();
    
    if(allSupportSize > shapeFuncs.size()) {
      shapeFuncs.resize(allSupportSize);

      if(derivationOrder > 0)

	for(int j=0;j<3;j++)
	  firstDerivShapes[j].resize(allSupportSize);

      if(derivationOrder > 1)

	for(int j=0;j<6;j++)
	  secondDerivShapes[j].resize(allSupportSize);

    }
    
#ifdef _geometryDebugMode_
    logFile<<"Gauss point "<<gaussPoints[i].getGlobalID()<<endl;
#endif

    double& x = gaussPoints[i].getCoord(0);
    double& y = gaussPoints[i].getCoord(1);
    double& z = gaussPoints[i].getCoord(2);

    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,basisTermNum,modelData,
		     logFile,viewerSEQ);
    
      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];
      }

      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);

    }
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,firstDerivShapes,modelData,
		     logFile,viewerSEQ);

      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];

	for(int k=0;k<firstDerivShapes.size();k++)
	  firstDerivShapes[k][pos] += firstDerivShapes[k][j];

      }

      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      gaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else if(derivationOrder == 2) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,firstDerivShapes,
		     secondDerivShapes,modelData,logFile,viewerSEQ);

      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];

	for(int k=0;k<firstDerivShapes.size();k++)
	  firstDerivShapes[k][pos] += firstDerivShapes[k][j];

	for(int k=0;k<secondDerivShapes.size();k++)
	  secondDerivShapes[k][pos] += secondDerivShapes[k][j];

      }

      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      gaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);
      gaussPoints[i].setSecondDerivShapes(supportSize,secondDerivShapes);

    }
    else {
      logFile <<"In ModifiedMLSGauss::setShapeFuncsOnGauss shape function "
	      <<"derivations higher than second order are not supported!"<<endl;
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
	  logFile<<gaussPoints[i].getSupportPtcle(j)<<" "<<firstDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
    if(derivationOrder > 1) {
      logFile<<"second order derivations"<<endl;
      dbMatrix& secondDerivs = gaussPoints[i].getSecondDerivShapes();
      for(int k=0;k<6;k++) {
	for(int j=0;j<secondDerivs[k].size();j++)
	  logFile<<gaussPoints[i].getSupportPtcle(j)<<" "<<secondDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
  }
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
}

/**********************************************************************/
/**********************************************************************/
// Calculate for all particles a vector containing the calculated
// shape functions and their derivations for all supported particles.
void ModifiedMLSGauss::setShapeFuncsOnPtcls(InputFileData* InputData,
					    std::map<std::string,double>& modelData,
					    std::ofstream& logFile,
					    PetscViewer& viewerMPI,
					    PetscViewer& viewerSEQ) {

  using namespace std;

  int rank,size; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"********************* particle shapes ********************"<<endl;
#endif

  int derivationOrder = 
    (int)modelData["shapesDerivationOrderOnPtcls"];

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  dbMatrix secondDerivShapes(6);
  int basisTermNum,pos,supportSize;

  int allSupportSize;
  intVector allSupportPtcls;

  int startIdx = exclusiveLocalPtcls.front();
  int endIdx = exclusiveLocalPtcls.back()+1;
  int ptclePortion = exclusiveLocalPtcls.size();
  
  // Loop over processor's particles portion.
  for(int i=startIdx;i<endIdx;i++) {

    intVector& suppGhostPtcls = particles[i].getSupportBoundGhostPtcls();
    intVector& suppPtcls = particles[i].getSupportPtcls();

    allSupportPtcls = suppPtcls;
    allSupportPtcls.insert(allSupportPtcls.end(),suppGhostPtcls.begin(),
			   suppGhostPtcls.end());

    supportSize = suppPtcls.size();
    allSupportSize = allSupportPtcls.size();
    
    if(allSupportSize > shapeFuncs.size()) {
      shapeFuncs.resize(allSupportSize);

      if(derivationOrder > 0)

	for(int j=0;j<3;j++)
	  firstDerivShapes[j].resize(allSupportSize);

      if(derivationOrder > 1)

	for(int j=0;j<6;j++)
	  secondDerivShapes[j].resize(allSupportSize);

    }
    
#ifdef _geometryDebugMode_
    logFile<<"PARTICLE "<<i<<endl;
#endif

    double& x = particles[i].getCoord(0);
    double& y = particles[i].getCoord(1);
    double& z = particles[i].getCoord(2);


    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,basisTermNum,modelData,
		     logFile,viewerSEQ);
    
      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];
      }

      // Store the calculated shape function set.
      particles[i].setShapeFuncs(supportSize,shapeFuncs);

    }
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,firstDerivShapes,modelData,
		     logFile,viewerSEQ);

      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];

	for(int k=0;k<firstDerivShapes.size();k++)
	  firstDerivShapes[k][pos] += firstDerivShapes[k][j];

      }

      // Store the calculated shape function set.
      particles[i].setShapeFuncs(supportSize,shapeFuncs);
      particles[i].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else if(derivationOrder == 2) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,firstDerivShapes,
		     secondDerivShapes,modelData,logFile,viewerSEQ);

      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];

	for(int k=0;k<firstDerivShapes.size();k++)
	  firstDerivShapes[k][pos] += firstDerivShapes[k][j];

	for(int k=0;k<secondDerivShapes.size();k++)
	  secondDerivShapes[k][pos] += secondDerivShapes[k][j];

      }

      // Store the calculated shape function set.
      particles[i].setShapeFuncs(supportSize,shapeFuncs);
      particles[i].setFirstDerivShapes(supportSize,firstDerivShapes);
      particles[i].setSecondDerivShapes(supportSize,secondDerivShapes);

    }
    else {
      logFile <<"In ModifiedMLSGauss::setShapeFuncsOnPtcls shape function\n "
	      <<"derivations higher than second order are not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }


  }

  /*********************************************************************/
  // Determine on all processor for all particles their supporting
  // particle numbers. Further determine for each particle its root
  // processor.
  intVector localSuppPtclsCounts(ptclePortion);
  int m=0;

  for(int i=startIdx;i<endIdx;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    localSuppPtclsCounts[m] = suppPtcls.size();
    m++;
  }
 
  int* recvCounts = new int[size];
  int* displs = new int[size];

  MPI_Allgather(&ptclePortion,1,MPI_INT,recvCounts,1,MPI_INT,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  intVector suppPtclsCounts(particlesNum);

  MPI_Allgatherv(&localSuppPtclsCounts[0],recvCounts[rank],MPI_INT,
		 &suppPtclsCounts[0],recvCounts,displs,MPI_INT,
		 MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"**************** supported particle counts ***************"<<endl;
  for(int i=0;i<particlesNum;i++) 
    logFile<<i<<" "<<suppPtclsCounts[i]<<endl;
  logFile<<"******************* root list ****************************"<<endl;
  for(int i=0;i<particlesNum;i++) 
    logFile<<i<<" "<<ptcleRootList[i]<<endl;
#endif


  // Assemble on all processors on all particles their supporting particles 
  // shape function sets and store them.
  for(int i=0;i<particlesNum;i++) {

    if(ptcleRootList[i] != rank)
      particles[i].setShapeFuncsSize(suppPtclsCounts[i],derivationOrder);

    dbVector& suppShapes = particles[i].getShapeFuncs();
    MPI_Bcast(&suppShapes[0],suppPtclsCounts[i],
	      MPI_DOUBLE,ptcleRootList[i],MPI_COMM_WORLD);

    if(derivationOrder > 0) {
      dbMatrix& firstDerivShapes = particles[i].getFirstDerivShapes();

      for(int j=0;j<firstDerivShapes.size();j++)
	
	MPI_Bcast(&firstDerivShapes[j][0],suppPtclsCounts[i],
		  MPI_DOUBLE,ptcleRootList[i],MPI_COMM_WORLD);

    }

    if(derivationOrder > 1) {
      dbMatrix& secondDerivShapes = particles[i].getSecondDerivShapes();

      for(int j=0;j<secondDerivShapes.size();j++)
	
	MPI_Bcast(&secondDerivShapes[j][0],suppPtclsCounts[i],
		  MPI_DOUBLE,ptcleRootList[i],MPI_COMM_WORLD);

    }

  }

  delete[] displs,recvCounts;

#ifdef _geometryDebugMode_
  double PUM;
  logFile<<"******************************************************"<<endl;
  logFile<<"*********** calculated shape sets on ptcls ***********"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    dbVector& suppShapes = particles[i].getShapeFuncs();
    dbMatrix& firstDerivShapes = particles[i].getFirstDerivShapes();
    dbMatrix& secondDerivShapes = particles[i].getSecondDerivShapes();
    PUM = 0;
    logFile<<"Ptcle "<<i<<":: ";
    for(int j=0;j<suppPtcls.size();j++) {
      PUM += suppShapes[j];
      logFile<<suppPtcls[j]<<": "
	     <<suppShapes[j]<<" | ";
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
	  logFile<<suppPtcls[j]<<": "<<firstDerivs[k][j]<<" | ";
	logFile<<endl;
      }
    }
    if(derivationOrder > 1) {
      logFile<<"second order derivations"<<endl;
      dbMatrix& secondDerivs = gaussPoints[i].getSecondDerivShapes();
      for(int k=0;k<6;k++) {
	for(int j=0;j<secondDerivs[k].size();j++)
	  logFile<<suppPtcls[j]<<": "<<secondDerivs[k][j]<<" | ";
	logFile<<endl;
      }
    }
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all boundary gauss points a vector containing the 
// calculated shape functions for all supported particles.
void ModifiedMLSGauss::setShapeFuncsOnBGauss(InputFileData* InputData,
					     std::map<std::string,double>& modelData,
					     std::ofstream& logFile,
					     PetscViewer& viewerMPI,
					     PetscViewer& viewerSEQ) {

  using namespace std;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"*************** boundary gausspoint shapes ***************"<<endl;
#endif

  unsigned int derivationOrder = (unsigned int)
    modelData["shapesDerivationOrderOnBoundIntPoints"];

  int basisTermNum,pos,supportSize;

  int allSupportSize;
  intVector allSupportPtcls;

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  dbMatrix secondDerivShapes(6);

  int m=0;

  // Loop over all local boundary gauss points to set their shape 
  // function vectors.
  for(int i=0;i<localBGaussPtsNum;i++) {

    intVector& suppGhostPtcls = boundGaussPoints[i].getSupportBoundGhostPtcls();
    intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls();

    allSupportPtcls = suppPtcls;
    allSupportPtcls.insert(allSupportPtcls.end(),suppGhostPtcls.begin(),
			   suppGhostPtcls.end());

    supportSize = suppPtcls.size();
    allSupportSize = allSupportPtcls.size();
    
    if(allSupportSize > shapeFuncs.size()) {
      shapeFuncs.resize(allSupportSize);

      if(derivationOrder > 0)

	for(int j=0;j<3;j++)
	  firstDerivShapes[j].resize(allSupportSize);

      if(derivationOrder > 1)

	for(int j=0;j<6;j++)
	  secondDerivShapes[j].resize(allSupportSize);

    }

    double& x = boundGaussPoints[i].getCoord(0);
    double& y = boundGaussPoints[i].getCoord(1);
    double& z = boundGaussPoints[i].getCoord(2);

#ifdef _geometryDebugMode_
    logFile<<"Boundary gauss point "<<i<<" coords: "<<x<<" "<<y<<" "<<z<<endl;
#endif

    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,basisTermNum,modelData,
		     logFile,viewerSEQ);

      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];
      }

      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);

    }
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,firstDerivShapes,modelData,
		     logFile,viewerSEQ);

      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);

#ifdef _geometryDebugMode_
	logFile<<"ghost ptcle "<<allSupportPtcls[j]<<" shape "<<shapeFuncs[j]
	       <<" mother ptcle "<<motherPtcle<<" shape "<<shapeFuncs[pos]
	       <<" coords: "<<particles[motherPtcle].getCoord(0)<<" "
	       <<particles[motherPtcle].getCoord(1)<<" "
	       <<particles[motherPtcle].getCoord(2)
	       <<" shape pos "<<pos<<" ";//endl;
#endif

	shapeFuncs[pos] += shapeFuncs[j];

#ifdef _geometryDebugMode_
	logFile<<"new mother shape "<<shapeFuncs[pos]<<endl;
#endif

	for(int k=0;k<firstDerivShapes.size();k++)
	  firstDerivShapes[k][pos] += firstDerivShapes[k][j];

      }

      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      boundGaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else if(derivationOrder == 2) {

      calcShapeFuncs(InputData,allSupportSize,allSupportPtcls,
		     particles,x,y,z,shapeFuncs,firstDerivShapes,
		     secondDerivShapes,modelData,logFile,viewerSEQ);

      // add shapefunction ordinate of ghost particle to its
      // mother particle 
      for(int j=supportSize;j<allSupportSize;j++) {
	int& motherPtcle = particles[allSupportPtcls[j]].getMotherPtcle();
	pos = findIntVecPos(motherPtcle,0,suppPtcls.size(),
			    suppPtcls);
	shapeFuncs[pos] += shapeFuncs[j];

	for(int k=0;k<firstDerivShapes.size();k++)
	  firstDerivShapes[k][pos] += firstDerivShapes[k][j];

	for(int k=0;k<secondDerivShapes.size();k++)
	  secondDerivShapes[k][pos] += secondDerivShapes[k][j];
      }

      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      boundGaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);
      boundGaussPoints[i].setSecondDerivShapes(supportSize,secondDerivShapes);

    }
    else {
      logFile <<"In ModifiedMLSGauss::setShapeFuncsOnBGauss shape function "
	      <<"derivations higher than second order are not supported!"<<endl;

      MPI_Abort(MPI_COMM_WORLD,1);
    }

#ifdef _geometryDebugMode_
      for(int j=0;j<allSupportSize;j++) {
	logFile<<allSupportPtcls[j]<<": "<<shapeFuncs[j]<<" | ";
      }
      logFile<<endl;
      if(derivationOrder > 0) {
	dbMatrix& firstDerivs = boundGaussPoints[i].getFirstDerivShapes();
	logFile<<"first order derivations"<<endl;
	for(int k=0;k<3;k++) {
	  for(int j=0;j<firstDerivs[k].size();j++)
	    logFile<<allSupportPtcls[j]<<" "<<firstDerivShapes[k][j]<<" |";
	  logFile<<endl;
	}
      }
      if(derivationOrder > 1) {
	dbMatrix& secondDerivs = boundGaussPoints[i].getSecondDerivShapes();
	logFile<<"second order derivations"<<endl;
	for(int k=0;k<6;k++) {
	  for(int j=0;j<secondDerivs[k].size();j++)
	    logFile<<allSupportPtcls[j]<<" "<<secondDerivShapes[k][j]<<" |";
	  logFile<<endl;
	}
      }
#endif

  }

#ifdef _geometryDebugMode_
  double PUM;
  logFile<<"######################################################"<<endl;
  logFile<<"** calculated shapes on boundary integration points **"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    PUM = 0;
    intVector& suppPtcls = boundGaussPoints[i].getSupportPtcls(); 
    logFile<<"Bound GAUSSPOINT "<<i<<": support = "
	   <<boundGaussPoints[i].getSupportCounts()<<" ";
    for(int j=0;j<boundGaussPoints[i].getSupportCounts();j++) {
      PUM += boundGaussPoints[i].getShapeFunc(j);
      logFile<<suppPtcls[j]<<": "<<boundGaussPoints[i].getShapeFunc(j)<<" | ";
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
	  logFile<<suppPtcls[j]<<" "<<firstDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
    if(derivationOrder > 1) {
      dbMatrix& secondDerivs = boundGaussPoints[i].getSecondDerivShapes();
      logFile<<"second order derivations"<<endl;
      for(int k=0;k<6;k++) {
	for(int j=0;j<secondDerivs[k].size();j++)
	  logFile<<suppPtcls[j]<<" "<<secondDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
  }
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
}
