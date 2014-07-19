#include "FEMGeometry.h"

FEMGeometry::FEMGeometry(InputFileData* InputData,
			 std::map<std::string,double>& modelData,
			 std::ifstream& meshFile,
			 std::ofstream& logFile,
			 PetscViewer& viewerMPI,
			 PetscViewer& viewerSEQ) :
  globalGaussPtsNum(0) {

  using namespace std;


  /*********************************************************************/
  // Read mesh file.
  readMeshFile(InputData,modelData,meshFile,logFile);

  // Create boundary elements.
  createBoundElems(InputData,modelData,logFile);

  // Ensure that all elements have the necessary approximation and
  // integration tools assigned.
  setApproxAndIntTools(InputData,modelData,logFile);


  /*********************************************************************/
  int integrationMethod = (int)modelData["integrationMethod"];

  // Gauss quadrature is chosen 
  if(integrationMethod == 1) {

    // Calculate all volume Gauss points.
    setVolumeGaussPoints(InputData,modelData,logFile,viewerMPI);

    // Calculate the boundary gauss points.
    setBoundGaussPoints(InputData,modelData,logFile,viewerMPI);

    // Calculate for all particles their particle weight (Nystroem 
    // Integration).
    setAllPtclsWeights(InputData,logFile);

    // Compute the shape function derivatives at all Gauss points with 
    // respect to global coordinates.
    if((int)InputData->getValue("shapefunctionType") == 7)

      setGlobalFEShapeFuncDerivs(InputData,modelData,logFile,viewerMPI);

  }

  //---------------------------------------------------------------------
  // Particle integration is chosen 
  else if(integrationMethod == 2) {

    // Calculate for all particles their particle weight 
    // (Nystroesm Integration) and surface and volume loads.
    setPtcleWeightsLoads(InputData,modelData,logFile,viewerMPI);

    // Determine weights and surface normals for these particles at 
    // which the boundary has to be integrated and store any applied 
    // loads.
    setBoundPtcleWeightsLoads(InputData,modelData,logFile,viewerMPI);

  }
}

FEMGeometry::~FEMGeometry() {

  for(int i=0;i<lineElemTemplates.size();i++)

    delete lineElemTemplates[i];

  for(int i=0;i<surfaceElemTemplates.size();i++)

    delete surfaceElemTemplates[i];

  for(int i=0;i<volumeElemTemplates.size();i++)

    delete volumeElemTemplates[i];

  for(int i=0;i<lineGaussPtTemplates.size();i++)

    delete lineGaussPtTemplates[i];

  for(int i=0;i<surfaceGaussPtTemplates.size();i++)

    delete surfaceGaussPtTemplates[i];

  for(int i=0;i<volumeGaussPtTemplates.size();i++)

    delete volumeGaussPtTemplates[i];
}

/***********************************************************************/
/***********************************************************************/
// Read a complete Finite Elements's mesh from the mesh file.
void FEMGeometry::readMeshFile(InputFileData* InputData,
			       std::map<std::string,double>& modelData,
			       std::ifstream& meshFile,
			       std::ofstream& logFile) {
  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  int n,nodeID,nodesNum;
  double coord1,coord2,coord3;
  double* dummyPointer;
  string name;


  /*********************************************************************/
  // Read the nodes' coordinates.
  meshFile >> name;
  meshFile >> nodesNum;

  if(name != "COORDINATES") {
    logFile<<"Mesh file 'mesh.dat' is invalid!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Contains the connectivity sorted node identifier to the read ones.
  oldIdx = intVector(nodesNum);

  // Contains the connectivity read node identifier to the sorted ones.
  newIdx = intVector(nodesNum);

  particles = vector<Particle>(nodesNum,Particle(usedDOF));
  dbMatrix ptcls(3,dbVector(nodesNum));

  for(int i=0;i<nodesNum;i++){    
    meshFile >> nodeID >> name >> ptcls[0][i] >> ptcls[1][i] >> ptcls[2][i];

    if(meshFile.eof()) {
      logFile<<"Mesh file 'mesh.dat' doesn't contain enough data!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

#ifdef _FEdebugMode_
  logFile<<"###################################################"<<endl;
  logFile<<"******************* unsorted nodes *****************"<<endl;
  for(int i=0;i<nodesNum;i++)
    logFile<<i<<" "<<ptcls[0][i]<<" "<<ptcls[1][i]<<" "<<ptcls[2][i]<<endl;
#endif


  /*********************************************************************/
  // Sort the particles in order to get a small stiffness matrix 
  // bandwidth.
  int sortCoord = (int)InputData->getValue("sortingCoordinate");

  // Run the sorting algorithm.
  if(sortCoord != 0) {

    for(int i=0;i<nodesNum;i++) { 
      oldIdx[i] = i;
    }   

    // Sort the particles in 3 dimensions.
    int numOfPlanes = 
      (int)InputData->getValue("numberOfSortingPlanes");

    sortParticles(ptcls,oldIdx,sortCoord,numOfPlanes,logFile);
    
    // Store the sorted particle coordinates in the particles vector and
    // set the index sets.
    for(int i=0;i<nodesNum;i++) {
      particles[i].setCoords(ptcls[0][i],ptcls[1][i],
			     ptcls[2][i]);

      newIdx[oldIdx[i]] = i;
      particles[i].setID(oldIdx[i]);
    }

  }
  // Don't run the sorting algorithm.
  else {

    // Store the node coordinates in the particles vector.
    for(int i=0;i<nodesNum;i++) {  
      particles[i].setCoords(ptcls[0][i],ptcls[1][i],
			     ptcls[2][i]);
      particles[i].setID(i);
    }

    // Store the connectivity unsorted node index to sorted node index.
    for(int i=0;i<nodesNum;i++) {
      newIdx[i] = i;
      oldIdx[i] = i;
    }

  }


#ifdef _FEdebugMode_
  logFile<<"************** MESH index -> SESKA index *************"<<endl;
  for(int i=0;i<nodesNum;i++)
    logFile<<i<<" -> "<<newIdx[i]<<" | "<<i+1<<" -> "<<newIdx[i]+1<<endl;
  logFile<<"************** SESKA index -> MESH index *************"<<endl;
  for(int i=0;i<nodesNum;i++)
    logFile<<i<<" -> "<<oldIdx[i]<<" | "<<i+1<<" -> "<<oldIdx[i]+1<<endl;
  logFile<<"******************* sorted nodes *********************"<<endl;
  for(int i=0;i<nodesNum;i++)
    logFile<<"SESKA-idx "<<i<<" MESH-idx "<<oldIdx[i]<<" coords: "
	   <<ptcls[0][i]<<" "
	   <<ptcls[1][i]<<" "
	   <<ptcls[2][i]<<endl;
  logFile<<"****************** sorted particles ******************"<<endl;
  for(int i=0;i<nodesNum;i++)
    logFile<<"SESKA-idx "<<i<<" MESH-Idx "<<particles[i].getID()
	   <<" coords: "
	   <<particles[i].getCoord(0)<<" "
	   <<particles[i].getCoord(1)<<" "
	   <<particles[i].getCoord(2)<<endl;
#endif


  /*********************************************************************/
  // Read and the element connectivity list.
  int globalElemNum,nodesPerElem;
  int referenceElemType;

  meshFile >> name;
  meshFile >> globalElemNum;

  if(globalElemNum == 0) {
    logFile<<"In FEMGeometry::readMeshFile no volume elements in 'mesh.dat'.\n" 
	   <<"Check elemType - must be hexahedral or tetrahedral!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  vector<FEMElement> allNodesElements(globalElemNum,FEMElement(usedDOF));

  int elemID,node;
  map<string,double> params;
  intVector data(4);

  // Loop over all global elements.
  for(int i=0;i<globalElemNum;i++) {

    int& elemType = allNodesElements[i].getElemType();
    int& elemOrder = allNodesElements[i].getElemOrder();
    int& matID = allNodesElements[i].getMaterialID();


    meshFile >> elemID >> name;

    // store element's material ID, type (tetrahedra,hexahedra), 
    // order (linear,quadratic)
    meshFile >> matID >> elemType >> elemOrder >> name;


    if(matID > (int)InputData->getNumberOfMats() || matID == 0) {
      logFile<<"In FEMGeometry::readMeshFile matID in mesh file 'mesh.dat'\n" 
	     <<"does not match materials specified in input file 'input.dat'!"
	     <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    matID -= 1; 

    // elemType is uniform thoughout the problem domain
    if(i == 0) 
      referenceElemType = elemType;

    else if(referenceElemType != elemType) {
      logFile<<"In FEMGeometry::readMeshFile elemType must be uniform\n" 
	     <<"thoughout the problem domain!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }


    // set nodes/particles
    data[0] = elemType;
    data[1] = elemOrder;

    getFEMMeshData(data,params);
    nodesPerElem = (int)params["nodesPerVolumeElement"];

    intVector& nodes = allNodesElements[i].getNodes();
    nodes.resize(nodesPerElem);

    // store the current element's nodes and their material ID
    for(int j=0;j<nodes.size();j++) {
      meshFile >> node;

      if(sortCoord != 0)
	nodes[j] = newIdx[node-1]+1;
      else
	nodes[j] = node;

      // store the current element's nodes material IDs
      int& ptcleMatID = particles[nodes[j]-1].getMaterialID();
      ptcleMatID = matID;

    }


    if(abs(elemType) > 2) {
      logFile<<"In FEMGeometry::readMeshFile elemType "<<elemType<<" "
	     <<"specified in file 'mesh.dat' is not supported!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    if(meshFile.eof()) {
      logFile<<"Mesh file 'mesh.dat' doesn't contain enough data!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

  /*********************************************************************/
  // Store information concerning the background mesh.
  std::map<std::string,double>& backGroundMeshInfo = 
    InputData->getBackGroundMeshInfo();
  
  backGroundMeshInfo["numberOfVolumeElements"] = (double)globalElemNum; 
  
#ifdef _FEdebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"*************** backgroundmesh info ******************"<<endl;
  logFile<<"nodesPerVolumeElement = "
	 <<backGroundMeshInfo["nodesPerVolumeElement"]<<endl;
  logFile<<"nodesPerSurfaceElement = "
	 <<backGroundMeshInfo["nodesPerSurfaceElement"]<<endl;
  logFile<<"nodesPerLineElement = "
	 <<backGroundMeshInfo["nodesPerLineElement"]<<endl;
  logFile<<"gaussPointsPerVolumeElement = "
	 <<backGroundMeshInfo["gaussPointsPerVolumeElement"]<<endl;
  logFile<<"gaussPointsPerSurfaceElement = "
	 <<backGroundMeshInfo["gaussPointsPerSurfaceElement"]<<endl;
  logFile<<"gaussPointsPerLineElement = "
	 <<backGroundMeshInfo["gaussPointsPerLineElement"]<<endl;
  logFile<<"numberOfVolumeElements = "
	 <<backGroundMeshInfo["numberOfVolumeElements"]<<endl;
#endif

  /*********************************************************************/
  // Write the FEM reference mesh.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  ofstream femMeshFile;

  if(rank == 0) {

    femMeshFile.open("fem.msh");
    femMeshFile.precision(12);
    femMeshFile.setf(ios_base::scientific,ios_base::floatfield);
    
    if((int)backGroundMeshInfo["elemType"] == 1)
      femMeshFile<<"MESH  dimension 3  ElemType Tetrahedra Nnode "
		 <<(int)backGroundMeshInfo["nodesPerVolumeElement"]<<endl;
    
    else if((int)backGroundMeshInfo["elemType"] == 2)
      femMeshFile<<"MESH  dimension 3  ElemType Hexahedra Nnode "
		 <<(int)backGroundMeshInfo["nodesPerVolumeElement"]<<endl;
    
    // Write nodes coordinates of reference mesh.
    femMeshFile<<"Coordinates"<<endl;
    
    for(int i=0;i<nodesNum;i++)
      femMeshFile<<i+1<<" "<<particles[newIdx[i]].getCoord(0)<<" "
		 <<particles[newIdx[i]].getCoord(1)<<" "
		 <<particles[newIdx[i]].getCoord(2)<<endl;
    
    femMeshFile<<"end coordinates\n"<<endl;
    
    femMeshFile<<"Elements"<<endl;
    
    // Loop over all global elements.
    for(int i=0;i<globalElemNum;i++) {
      int& elemType = allNodesElements[i].getElemType();
      
      int& matID = allNodesElements[i].getMaterialID();
      intVector& nodes = allNodesElements[i].getNodes();
      
      femMeshFile<<i+1<<" ";
      
      for(int j=0;j<nodes.size();j++)
	femMeshFile<<oldIdx[nodes[j]-1]+1<<" ";
      
      femMeshFile<<matID<<endl;
    }
    
    femMeshFile<<"end elements"<<endl;
  }

  
#ifdef _FEdebugMode_
  logFile<<"************** all unsorted elements *****************"<<endl;
  for(int i=0;i<globalElemNum;i++) {
    logFile<<i<<" ";
    intVector& nodes = allNodesElements[i].getNodes();
    int& matID = allNodesElements[i].getMaterialID();
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<matID<<endl;
  }
#endif

  /*********************************************************************/
  // Sort the nodes according to the chosen coordinate and then store
  // a portion of elements to each processor.
  storeElements(InputData,allNodesElements,modelData,logFile);


  logFile<<"**************** local elements **********************"<<endl;
  for(int i=0;i<nodesElements.size();i++) {
    logFile<<i<<" ";
    intVector& nodes = nodesElements[i].getNodes();
    int& elemType = nodesElements[i].getElemType();
    int& matID = nodesElements[i].getMaterialID();
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<" material = "<<matID<<" elemType = "<<elemType<<endl;
  }

  logFile<<"************************ nodes **********************"<<endl;
    for(int i=0;i<nodesNum;i++) {
      logFile<<"newIdx "<<i<<" oldIdx "<<particles[i].getID()<<" coords: "
  	   <<particles[i].getCoord(0)<<" "
  	   <<particles[i].getCoord(1)<<" "
  	   <<particles[i].getCoord(2)<<" material = "
  	   <<particles[i].getMaterialID()<<endl;
    }

#ifdef _FEdebugMode_
  logFile<<"**************** local elements **********************"<<endl;
  for(int i=0;i<nodesElements.size();i++) {
    logFile<<i<<" ";
    intVector& nodes = nodesElements[i].getNodes();
    int& elemType = nodesElements[i].getElemType();
    int& matID = nodesElements[i].getMaterialID();
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<" material = "<<matID<<" elemType = "<<elemType<<endl;
  }
  logFile<<"************************ nodes **********************"<<endl;
  for(int i=0;i<nodesNum;i++) {
    logFile<<"newIdx "<<i<<" oldIdx "<<particles[i].getID()<<" coords: "
	   <<particles[i].getCoord(0)<<" "
	   <<particles[i].getCoord(1)<<" "
	   <<particles[i].getCoord(2)<<" material = "
	   <<particles[i].getMaterialID()<<endl;
  }
  logFile<<"******* body force applied to local elements *********"<<endl;
  logFile<<"SIZE "<<bodyForceElemIdx.size()<<endl;
  for(int i=0;i<bodyForceElemIdx.size();i++) {
    int& currentIdx = bodyForceElemIdx[i][0];
    blVector& affectedDOF = nodesElements[currentIdx].getBodyForceDOF();
    dbVector& conditions = nodesElements[currentIdx].getBodyForce();
    logFile<<"VOLUME ELEMENT "<<currentIdx<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"******* body moment applied to local elements *********"<<endl;
  logFile<<"SIZE "<<bodyMomentElemIdx.size()<<endl;
  for(int i=0;i<bodyMomentElemIdx.size();i++) {
    int& currentIdx = bodyMomentElemIdx[i][0];
    blVector& affectedDOF = nodesElements[currentIdx].getBodyMomentDOF();
    dbVector& conditions = nodesElements[currentIdx].getBodyMoment();
    logFile<<"VOLUME ELEMENT "<<currentIdx<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"******* body electric charge applied to local elements *********"<<endl;
  logFile<<"SIZE "<<bodyElectricChargeElemIdx.size()<<endl;
  for(int i=0;i<bodyElectricChargeElemIdx.size();i++) {
    int& currentIdx = bodyElectricChargeElemIdx[i][0];
    blVector& affectedDOF = nodesElements[currentIdx].getBodyElectricChargeDOF();
    dbVector& conditions = nodesElements[currentIdx].getBodyElectricCharge();
    logFile<<"VOLUME ELEMENT "<<currentIdx<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Sort the nodes according to the chosen coordinate and then store
// a portion of elements to each processor.
void FEMGeometry::storeElements(InputFileData* InputData,
				std::vector<FEMElement>& allNodesElements,
				std::map<std::string,double>& modelData,
				std::ofstream& logFile) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int usedDims = (int)modelData["usedDimensions"];
  int sortCoord = (int)InputData->getValue("sortingCoordinate");

  int globalElemNum = allNodesElements.size();
  intVector oldGlobElemIdx(globalElemNum);
  intVector newGlobElemIdx(globalElemNum);
  dbVector coords(globalElemNum);

  // Elements are going to get sorted.
  if(sortCoord != 0) {

    // Store form each element the sorting coordinate of the first node.
    for(int i=0;i<globalElemNum;i++) {
      oldGlobElemIdx[i] = i;
      intVector& nodes = allNodesElements[i].getNodes();
      coords[i] = particles[oldIdx[nodes[0]-1]].getCoord(sortCoord-1);
    }

    // Sort the global elements - oldGlobElemIdx: new global index -> old 
    // global index.
    sortValuesIdx(coords,oldGlobElemIdx,0,globalElemNum-1);

    // Store the connectivity old global index -> new global index.
    for(int i=0;i<globalElemNum;i++) 
      newGlobElemIdx[oldGlobElemIdx[i]] = i;

  }
  else {

    matchElemNodeOrder(InputData,allNodesElements,oldGlobElemIdx,
		       newGlobElemIdx,modelData,logFile);

//    for(int i=0;i<globalElemNum;i++)
//      oldGlobElemIdx[i] = newGlobElemIdx[i] = i;

  }

  // arrange the elements in the new ordering and store only the
  // local elements

  rearrangeElements(InputData,allNodesElements,oldGlobElemIdx,
		    modelData,logFile);
  
  // connectivity new global to new local element indices

  globalLocalElemIdx.resize(newLocalElemIdx.size());

  for(int i=0;i<newLocalElemIdx.size();i++) {
    
    if(newLocalElemIdx[i] >= 0)
      globalLocalElemIdx[newGlobElemIdx[i]] = newLocalElemIdx[i];

    else
      
      globalLocalElemIdx[newGlobElemIdx[i]] = -1;

  } 

#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"*************** new element ordering ****************"<<endl;
  logFile<<"******* connectivity old -> new index ***************"<<endl;
  for(int i=0;i<globalElemNum;i++)
    logFile<<"old global "<<i<<" -> new global "<<newGlobElemIdx[i]
	   <<" / new local "<<newLocalElemIdx[i]<<endl;
  logFile<<"***** connectivity new global -> new local index ****"<<endl;
  for(int i=0;i<globalElemNum;i++)
    logFile<<"new global "<<i
	   <<" -> new local "<<globalLocalElemIdx[i]<<endl;
  logFile<<"******* connectivity new -> old index ***************"<<endl;
  for(int i=0;i<globalElemNum;i++)
    logFile<<i<<" -> old global "<<oldGlobElemIdx[i]<<endl;
  logFile<<"********* root list ************************"<<endl;
  for(int i=0;i<elementRootList.size();i++) {
    logFile<<i<<".) "<<elementRootList[i]<<endl;
  }
#endif

  /*********************************************************************/
  // Store the body force.

  // Determine the local portion.
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int number = (int)ceil((double)globalElemNum/size);

  int localElemNum = nodesElements.size();
  int currentIdx,dof,elem,position;
  intVector elemIdx;
  intVector indexSet(2);

  dbMatrix& inputForce = InputData->getBodyForceLoads();

  // loop over all body force conditions
  for(int i=0;i<inputForce.size();i++) {

    // Only loads at local elements are stored.
    if(newGlobElemIdx[(int)inputForce[i][0]-1] >= rank*number 
       && newGlobElemIdx[(int)inputForce[i][0]-1] < rank*number+localElemNum) {

      // Store these local element's identifiers which have a body force 
      // applied.
      currentIdx = newLocalElemIdx[(int)inputForce[i][0]-1];

      position = findIntVecPos(currentIdx,0,elemIdx.size(),elemIdx);

      // check current element has already got a body force
      if(position != -1) {
	logFile<<"In function FEMGeometry::storeElements repetetive "
	       <<"body force conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
      }
      else {
	
	// Store the the body force.
	blVector& affectedDOF = nodesElements[currentIdx].getBodyForceDOF();
	dbVector& bodyForce = nodesElements[currentIdx].getBodyForce();
	affectedDOF.resize(usedDims);
	bodyForce.resize(usedDims);
	
	intVector& nodes = nodesElements[currentIdx].getNodes();

	for(int j=nodes.size()+1;j<inputForce[i].size();j+=2) {
	  dof = (int)inputForce[i][j]-1;
	  affectedDOF[dof] = true;
	  bodyForce[dof] = inputForce[i][j+1];
	}
	
	indexSet[0] = currentIdx;
	indexSet[1] = 0; // load ID
	pushBackVector(bodyForceElemIdx,indexSet);
	elemIdx.push_back(currentIdx);
      }
      
    }
    
  }

  // --------------------------------------------------------------------
  // body moment conditions
  dbMatrix& inputMoment = InputData->getBodyMomentLoads();

  // loop over all body force conditions
  for(int i=0;i<inputMoment.size();i++) {

    // Only loads at local elements are stored.
    if(newGlobElemIdx[(int)inputMoment[i][0]-1] >= rank*number 
       && newGlobElemIdx[(int)inputMoment[i][0]-1] < rank*number+localElemNum) {

      // Store these local element's identifiers which have a body moment
      // applied.
      currentIdx = newLocalElemIdx[(int)inputMoment[i][0]-1];

      position = findIntVecPos(currentIdx,0,elemIdx.size(),elemIdx);

      // check current element has already got a body force
      if(position != -1) {
	logFile<<"In function FEMGeometry::storeElements repetetive "
	       <<"body moment conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
      }
      else {
	
	// Store the the body moment.
	blVector& affectedDOF = nodesElements[currentIdx].getBodyMomentDOF();
	dbVector& bodyMoment = nodesElements[currentIdx].getBodyMoment();
	affectedDOF.resize(usedDims);
	bodyMoment.resize(usedDims);
	
	intVector& nodes = nodesElements[currentIdx].getNodes();

	for(int j=nodes.size()+1;j<inputMoment[i].size();j+=2) {
	  dof = (int)inputMoment[i][j]-1;
	  affectedDOF[dof] = true;
	  bodyMoment[dof] = inputMoment[i][j+1];
	}
	
	indexSet[0] = currentIdx;
	indexSet[1] = 0; // load ID
	pushBackVector(bodyMomentElemIdx,indexSet);
	elemIdx.push_back(currentIdx);
      }
      
    }
    
  }

  // ---------------------------------------------------------------------
  // Store the body electric charge.

  dbMatrix& inputElecCharge = InputData->getBodyElectricChargeLoads();

  // loop over all body electric charge conditions
  for(int i=0;i<inputElecCharge.size();i++) {

    // Only loads at local elements are stored.
    if(newGlobElemIdx[(int)inputElecCharge[i][0]-1] >= rank*number 
       && newGlobElemIdx[(int)inputElecCharge[i][0]-1] < rank*number+localElemNum) {

      // Store these local element's identifiers which have a body electric charge
      // applied.
      currentIdx = newLocalElemIdx[(int)inputElecCharge[i][0]-1];

      position = findIntVecPos(currentIdx,0,elemIdx.size(),elemIdx);

      // check current element has already got a body charge
      if(position != -1) {
	logFile<<"In function FEMGeometry::storeElements repetetive "
	       <<"body electric charge conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
      }
      else {
	
	// Store the the body force.
	blVector& affectedDOF = nodesElements[currentIdx].getBodyElectricChargeDOF();
	dbVector& bodyCharge = nodesElements[currentIdx].getBodyElectricCharge();
	affectedDOF.resize(1);
	bodyCharge.resize(1);
	
	intVector& nodes = nodesElements[currentIdx].getNodes();

	for(int j=nodes.size()+1;j<inputElecCharge[i].size();j+=2) {
	  dof = (int)inputElecCharge[i][j]-1;
	  affectedDOF[dof] = true;
	  bodyCharge[dof] = inputElecCharge[i][j+1];
	}
	
	indexSet[0] = currentIdx;
	indexSet[1] = 0; // load ID
	pushBackVector(bodyElectricChargeElemIdx,indexSet);
	elemIdx.push_back(currentIdx);
      }
      
    }
    
  }

}

/**********************************************************************/
/**********************************************************************/
// Sort the elements to have a similar geometrical ordering as the 
// particles.
void FEMGeometry::matchElemNodeOrder(InputFileData* InputData,
				     std::vector<FEMElement>& allNodesElements,
				     std::vector<int>& oldGlobElemIdx,
				     std::vector<int>& newGlobElemIdx,
				     std::map<std::string,double>& modelData,
				     std::ofstream& logFile) {

  using namespace std;

  int rank,size; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  int particlesNum = particles.size();


  /**********************************************************************/
  // set up particle - element connectivity

  int ptcle;

  // Loop over all global elements.
  for(int i=0;i<allNodesElements.size();i++) {
    FEMElement& elem = allNodesElements[i];
    int& globalID = elem.getGlobalID();
    intVector& nodes = elem.getNodes();
    
    for(int j=0;j<nodes.size();j++) {
      ptcle = oldIdx[nodes[j]-1];
      particles[ptcle].setElems(i);
    }

  }

#ifdef _FEdebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"************* node-elem connectivity ****************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& connectedElems = particles[i].getElems();
    logFile<<"NODE "<<i<<" connected elems: ";
    for(int j=0;j<connectedElems.size();j++)
      logFile<<connectedElems[j]<<" ";
    logFile<<endl;
  }
#endif


  /**********************************************************************/
  // sort elements in accordance with ordering of the particles

  int pos;

  allocateArray(oldGlobElemIdx,allNodesElements.size());
  allocateArray(newGlobElemIdx,allNodesElements.size());
  int m=0;

  for(int i=0;i<particlesNum;i++) {

    intVector& elems = particles[i].getElems();

    for(int j=0;j<elems.size();j++) {

      pos = findIntVecPos(elems[j],0,m,oldGlobElemIdx);

      // check current element is already stored
      if(pos == -1) {
	oldGlobElemIdx[m] = elems[j];
	m++;
      }

    }

    resizeArray(elems,0);
  }

  // Store the connectivity old global index -> new global index.
  for(int i=0;i<allNodesElements.size();i++) 
    newGlobElemIdx[oldGlobElemIdx[i]] = i;


#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"*************** new element ordering ****************"<<endl;
  logFile<<"******* connectivity old -> new index ***************"<<endl;
  for(int i=0;i<newGlobElemIdx.size();i++)
    logFile<<i<<" -> new global "<<newGlobElemIdx[i]<<endl;
  logFile<<"******* connectivity new -> old index ***************"<<endl;
  for(int i=0;i<oldGlobElemIdx.size();i++)
    logFile<<i<<" -> old global "<<oldGlobElemIdx[i]<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// arrange the elements in the new ordering and store only the
// local elements
void FEMGeometry::rearrangeElements(InputFileData* InputData,
				    std::vector<FEMElement>& allNodesElements,
				    intVector& oldGlobElemIdx,
				    std::map<std::string,double>& modelData,
				    std::ofstream& logFile) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int globalElemNum = allNodesElements.size();

  // determine a local portion of elements
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  
  // Determine the local portion.
  int number = (int)ceil((double)globalElemNum/size);
  
  int localElemNum;

  if(number*(rank+1) <= globalElemNum) {
    localElemNum = number;
  }
  else if(number*rank >= globalElemNum) {
    localElemNum = 0;
    logFile<<"In function FEMGeometry::rearrangeElements "<<localElemNum
	   <<" local volume elements!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  else {
    localElemNum = globalElemNum - number*rank;
  }

  nodesElements = vector<FEMElement>(localElemNum,FEMElement(usedDOF));

  /*********************************************************************/
  // restructure the element ordering
  newLocalElemIdx.resize(globalElemNum);
  elementRootList.resize(globalElemNum);

  intVector localElementRootList(globalElemNum);

  int n = 0;

  // Store the nodes and the loads of an elements.
  for(int i=0;i<globalElemNum;i++) {
    
    // Store the nodes and the material ID.
    if(i >= rank*number && i<rank*number+localElemNum) {
      nodesElements[n].setGlobalID(i);
    
      localElementRootList[i] = rank;
  
      int& elemType = nodesElements[n].getElemType();
      int& elemOrder = nodesElements[n].getElemOrder();
      int& matID = nodesElements[n].getMaterialID();
      intVector& nodes = nodesElements[n].getNodes();

      nodes = allNodesElements[oldGlobElemIdx[i]].getNodes();
      matID = allNodesElements[oldGlobElemIdx[i]].getMaterialID();
      elemType = allNodesElements[oldGlobElemIdx[i]].getElemType();
      elemOrder = allNodesElements[oldGlobElemIdx[i]].getElemOrder();
  
      // clear storage of the global element 
      intVector& globalNodes = allNodesElements[oldGlobElemIdx[i]].getNodes();
      globalNodes.resize(0);

      // Store the connectivity old global index -> new local index.
      newLocalElemIdx[oldGlobElemIdx[i]] = n;
      
      n++;
    }
    else
      newLocalElemIdx[oldGlobElemIdx[i]] = -1;

  }

  MPI_Allreduce(&localElementRootList[0],&elementRootList[0],
		globalElemNum,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

}

/**********************************************************************/
/**********************************************************************/
// Transfer element info from volume to surface or line element.
bool FEMGeometry::swapVolToFaceElemInfo(InputFileData* InputData,
					FEMElement& sElem,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile) {
  
  using namespace std;

  string mode = "arbitrary-subvector";

  bool flag;
  intVector localvElemInfo(4);

  intVector& sNodes = sElem.getNodes();

  int localvElem = -1;

#ifdef _FEdebugMode_
  /*logFile<<"surface-Elem-nodes : ";
    for(int k=0;k<sNodes.size();k++)
    logFile<<sNodes[k]<<" ";
    logFile<<endl;*/
#endif

  // loop over all volume elements in order to find the volume element 
  // to which the line element belongs
  for(int j=0;j<nodesElements.size();j++) {
    intVector& vNodes = nodesElements[j].getNodes();

    flag = compareIntVecs(mode,sNodes,vNodes,logFile);

#ifdef _FEdebugMode_
    /* logFile<<"-> vElem "<<j<<": ";
       for(int k=0;k<vNodes.size();k++)
       logFile<<vNodes[k]<<" ";
       logFile<<endl;
       if(flag)
       logFile<<"check"<<endl; */
#endif

    if(flag) {

      localvElem = j;

      localvElemInfo[0] = nodesElements[j].getMaterialID();
      localvElemInfo[1] = nodesElements[j].getElemType();
      localvElemInfo[2] = nodesElements[j].getElemOrder();
      localvElemInfo[3] = nodesElements[j].getGlobalID(); 


      break;
    }

  }

  // determine the root process which is hosting the data
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  intVector globalvElem(size);

  MPI_Allgather(&localvElem,1,MPI_INT,&globalvElem[0],1,MPI_INT,
		MPI_COMM_WORLD);

  int root = -1;

  for(int i=0;i<size;i++)

    if(globalvElem[i] >= 0) {
      root = i;
      break;

    }

#ifdef _FEdebugMode_
  /*for(int i=0;i<globalvElem.size();i++)
    logFile<<"globalvElem["<<i<<"]="<<globalvElem[i]<<endl;
    logFile<<"root="<<root<<endl;*/
#endif


  // exchange data amongst threads
  if(root >= 0) {

    MPI_Bcast(&localvElemInfo[0],4,MPI_INT,root,MPI_COMM_WORLD);  

    int& elemMatID = sElem.getMaterialID();
    int& elemType = sElem.getElemType();
    int& elemOrder = sElem.getElemOrder();
    int& volElemID = sElem.getMotherElementID();
      
    elemMatID = localvElemInfo[0];
    elemType = localvElemInfo[1];
    elemOrder = localvElemInfo[2];
    volElemID = localvElemInfo[3];

    flag = true;
  }
  else 

    flag = false;


  return flag;
}

/**********************************************************************/
/**********************************************************************/
// Create boundary elements.
void FEMGeometry::createBoundElems(InputFileData* InputData,
				   std::map<std::string,double>& modelData,
				   std::ofstream& logFile) {

  using namespace std;

  // Create line elements and store their line load.
  createLineElems(InputData,modelData,logFile);

  // Create surface elements and store their surface load.
  createSurfaceElems(InputData,modelData,logFile);

}

/**********************************************************************/
/**********************************************************************/
// Create surface elements and store their surface load.
void FEMGeometry::createSurfaceElems(InputFileData* InputData,
				     std::map<std::string,double>& modelData,
				     std::ofstream& logFile) {

  using namespace std;

  std::map<std::string,double>& backGroundMeshInfo =
    InputData->getBackGroundMeshInfo();

  bool dirichletIntegration = 
    (bool)modelData["dirichletBoundaryIntegration"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];
  int dispDOF = (int)modelData["displacementDegreesOfFreedom"];
  int rotDOF = (int)modelData["rotationDegreesOfFreedom"];
  int electricDOF = (int)modelData["electricDegreesOfFreedom"];
  int depolarisationDOF = (int)modelData["depolarisationDegreesOfFreedom"];
  int microDOF = (int)modelData["microDegreesOfFreedom"];
  int usedDims = (int)modelData["usedDimensions"];
  int nodesPerSurfaceElem = (int)backGroundMeshInfo["nodesPerSurfaceElement"];

  dbMatrix& dispBoundConds = InputData->getSurfaceDispBoundConds();
  dbMatrix& rotBoundConds = InputData->getSurfaceRotBoundConds();
  dbMatrix& electricBoundConds = InputData->getSurfaceElectricBoundConds();
  dbMatrix& depolarisationBoundConds = InputData->getSurfaceDepolarisationBoundConds();
  dbMatrix& microBoundConds = InputData->getSurfaceMicroBoundConds();

  dbMatrix& tractionLoads = InputData->getTractionLoads();
  dbMatrix& pressureLoads = InputData->getSurfacePressureLoads();
  dbMatrix& momentLoads = InputData->getSurfaceMomentLoads();
  dbMatrix& electricChargeLoads = InputData->getSurfaceElectricChargeLoads();

  /**********************************************************************/
  // delete not used boundary conditions

  if(rotDOF == 0) {

    rotBoundConds.resize(0);
    momentLoads.resize(0);
  }

  if(electricDOF == 0) {

    electricBoundConds.resize(0);
    electricChargeLoads.resize(0);
  }

  if(depolarisationDOF == 0) {

    depolarisationBoundConds.resize(0);
  }

  if(microDOF == 0)
    microBoundConds.resize(0);


  dbMatrix defBoundConds = dispBoundConds;
  defBoundConds.insert(defBoundConds.begin(),rotBoundConds.begin(),
		       rotBoundConds.end());


  /**********************************************************************/
  // Create the vector containing all boundary elements.
  int surfaceElemNum = tractionLoads.size() + pressureLoads.size()+
    momentLoads.size() + electricChargeLoads.size() + defBoundConds.size();

  surfaceElemNum += depolarisationBoundConds.size() + electricBoundConds.size() + microBoundConds.size();

  surfaceNodesElems = 
    vector<FEMElement>(surfaceElemNum,FEMElement(usedDOF));

  // Loop over all boundary elements applied a surface load.
  int currentIdx,dof,elem,pos,vElem;
  intVector boundaryElemIdx;
  intVector indexSet(2);

  // Loop over all elements traction loads are applied.
  for(int i=0;i<tractionLoads.size();i++) {

    elem = (int)tractionLoads[i][0]-1;

    pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			boundaryElemIdx);

    if(pos != -1) {
      logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
	     <<"traction conditions!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {

      currentIdx = boundaryElemIdx.size();

      surfaceNodesElems[currentIdx].setGlobalID(elem);

      // Store the the traction load.
      blVector& boundDOF = surfaceNodesElems[currentIdx].getTractionDOF();
      dbVector& traction = surfaceNodesElems[currentIdx].getTraction();
      boundDOF.resize(usedDims);
      traction.resize(usedDims);
    
      for(int j=nodesPerSurfaceElem+1;j<tractionLoads[i].size();j+=2) {
	dof = (int)tractionLoads[i][j]-1;
	boundDOF[dof] = 1;
	traction[dof] = tractionLoads[i][j+1];
      }

      // store the element's nodes
      intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
      nodes.resize(nodesPerSurfaceElem);

      for(int j=0;j<nodesPerSurfaceElem;j++) 
	nodes[j] = newIdx[(int)tractionLoads[i][j+1]-1]+1;


      // set the element's material, type (tetrahedra,cube) 
      // and order (linear/quadratic)

      if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				modelData,logFile)) {
	logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
	       <<"volume element can be found for surface element '"<<elem<<"'!"
	       <<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }


      indexSet[0] = currentIdx;
      indexSet[1] = 0; // loading ID

      pushBackVector(tractionElemIdx,indexSet);
      pushBackVector(boundaryElemIdx,elem);
    }

  }

  /*********************************************************************/
  // Loop over all elements pressure loads are applied.
  for(int i=0;i<pressureLoads.size();i++) {

    // Check if a load is already applied to current elment.
    elem = (int)pressureLoads[i][0]-1;

    pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			boundaryElemIdx);

    // no surface loading is already applied
    if(pos == -1) {
      currentIdx = boundaryElemIdx.size();

      surfaceNodesElems[currentIdx].setGlobalID(elem);
 
      // Store the the pressure load.
      double& pressure = surfaceNodesElems[currentIdx].getSurfacePressure();
      pressure = pressureLoads[i][nodesPerSurfaceElem+1];

      // store the element's nodes
      intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
      nodes.resize(nodesPerSurfaceElem);

      for(int j=0;j<nodesPerSurfaceElem;j++) 
	nodes[j] = newIdx[(int)pressureLoads[i][j+1]-1]+1;
     

      // set the element's material, type (tetrahedra,cube) 
      // and order (linear/quadratic)

      if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				modelData,logFile)) {
	logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
	       <<"volume element can be found for surface element '"<<elem<<"'!"
	       <<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }


      indexSet[0] = currentIdx;
      indexSet[1] = 0; // loading ID

      pushBackVector(surfacePressureElemIdx,indexSet);
      pushBackVector(boundaryElemIdx,elem);
    }
    else {
      currentIdx = pos;

      // Check if a surface pressure load is already applied to current elment.
      pos = findIntMatPos(currentIdx,0,surfacePressureElemIdx.size(),0,
			  surfacePressureElemIdx);
      
      if(pos != -1) {
	logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
	       <<"surface pressure conditions!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      else {

	// Store the the pressure load.
	double& pressure = surfaceNodesElems[currentIdx].getSurfacePressure();
	pressure = pressureLoads[i][nodesPerSurfaceElem+1];

	indexSet[0] = currentIdx;
	indexSet[1] = 0; // loading ID

	pushBackVector(surfacePressureElemIdx,indexSet);
      }

    }

  }

  /*********************************************************************/
  // Loop over all elements surface moment loads are applied.
  for(int i=0;i<momentLoads.size();i++) {

    elem = (int)momentLoads[i][0]-1;

    pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			boundaryElemIdx);
    
    if(pos != -1) {
      logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
	     <<"surface moment conditions!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {

      currentIdx = boundaryElemIdx.size();

      surfaceNodesElems[currentIdx].setGlobalID(elem);

      // Store the the moment load.
      blVector& boundDOF = surfaceNodesElems[currentIdx].getSurfaceMomentDOF();
      dbVector& moment = surfaceNodesElems[currentIdx].getSurfaceMoment();
      boundDOF.resize(usedDims);
      moment.resize(usedDims);
    
      for(int j=nodesPerSurfaceElem+1;j<momentLoads[i].size();j+=2) {
	dof = (int)momentLoads[i][j]-1;
	boundDOF[dof] = 1;
	moment[dof] = momentLoads[i][j+1];
      }

      // store the element's nodes
      intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
      nodes.resize(nodesPerSurfaceElem);

      for(int j=0;j<nodesPerSurfaceElem;j++) 
	nodes[j] = newIdx[(int)momentLoads[i][j+1]-1]+1;


      // set the element's material, type (tetrahedra,cube) 
      // and order (linear/quadratic)

      if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				modelData,logFile)) {
	logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
	       <<"volume element can be found for surface element '"<<elem<<"'!"
	       <<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }


      indexSet[0] = currentIdx;
      indexSet[1] = 0; // loading ID

      pushBackVector(surfaceMomentElemIdx,indexSet);
      pushBackVector(boundaryElemIdx,elem);
    }

  }

  /*********************************************************************/
  // Loop over all elements surface electric charge loads are applied.
  for(int i=0;i<electricChargeLoads.size();i++) {

    elem = (int)electricChargeLoads[i][0]-1;

    pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			boundaryElemIdx);
    
    if(pos != -1) {
      logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
	     <<"surface electric charge conditions!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {

      currentIdx = boundaryElemIdx.size();

      surfaceNodesElems[currentIdx].setGlobalID(elem);

      // Store the the electricCharge load.
      blVector& boundDOF = surfaceNodesElems[currentIdx].getSurfaceElectricChargeDOF();
      dbVector& loads = surfaceNodesElems[currentIdx].getSurfaceElectricCharge();
      boundDOF.resize(1);
      loads.resize(1);
    
      for(int j=nodesPerSurfaceElem+1;j<electricChargeLoads[i].size();j+=2) {
	dof = (int)electricChargeLoads[i][j]-1;
	boundDOF[dof] = 1;
	loads[dof] = electricChargeLoads[i][j+1];
      }

      // store the element's nodes
      intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
      nodes.resize(nodesPerSurfaceElem);

      for(int j=0;j<nodesPerSurfaceElem;j++) 
	nodes[j] = newIdx[(int)electricChargeLoads[i][j+1]-1]+1;


      // set the element's material, type (tetrahedra,cube) 
      // and order (linear/quadratic)

      if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				modelData,logFile)) {
	logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
	       <<"volume element can be found for surface element '"<<elem<<"'!"
	       <<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }


      indexSet[0] = currentIdx;
      indexSet[1] = 0; // loading ID

      pushBackVector(surfaceElectricChargeElemIdx,indexSet);
      pushBackVector(boundaryElemIdx,elem);
    }

  }

  /*********************************************************************/
  // Loop over all elements deformation boundary conditions are applied,
  // if integration over deformation boundary area is required. 
  if(dirichletIntegration) {

    // displacement boundary constraints
    for(int i=0;i<dispBoundConds.size();i++) {

      if(dispDOF <= 0) {
	logFile<<"In function FEMGeometry::createSurfaceElems "
	       <<"displacement boundary conditions are not supported!"<<endl;
	break;
	//MPI_Abort(MPI_COMM_WORLD,1);
      }


      // Check if a load is already applied to current element.
      elem = (int)dispBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading is already applied
      if(pos == -1) {

	currentIdx = boundaryElemIdx.size();
	surfaceNodesElems[currentIdx].setGlobalID(elem);
	
	// Store the displacement boundary conditions.
	blVector& boundDOF = 
	  surfaceNodesElems[currentIdx].getSurfaceDeformationBoundDOF();
	dbVector& boundConds = 
	  surfaceNodesElems[currentIdx].getSurfaceDeformationBoundConds();
	boundDOF.resize(defDOF);
	boundConds.resize(defDOF);
      
	// Loop over all affected degrees of current element.
	for(int j=nodesPerSurfaceElem+1;j<dispBoundConds[i].size();j+=2) {
	  dof = (int)dispBoundConds[i][j]-1;
	  boundDOF[dof] = 1;
	  boundConds[dof] = dispBoundConds[i][j+1];
	}

	// store the element's nodes
	intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
	nodes.resize(nodesPerSurfaceElem);

	for(int j=0;j<nodesPerSurfaceElem;j++) 
	  nodes[j] = newIdx[(int)dispBoundConds[i][j+1]-1]+1;
       

	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
	
	if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID
	
	pushBackVector(surfaceDispBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }
      // a loading is already applied
      else {
	currentIdx = pos; 

	// Check if a surface displacement boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(currentIdx,0,surfaceDispBoundElemIdx.size(),0,
			    surfaceDispBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
		 <<"displacement boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the displacement boundary conditions.
	  blVector& boundDOF = 
	    surfaceNodesElems[currentIdx].getSurfaceDeformationBoundDOF();
	  dbVector& boundConds = 
	    surfaceNodesElems[currentIdx].getSurfaceDeformationBoundConds();
	  boundDOF.resize(defDOF);
	  boundConds.resize(defDOF);
	  
	  // Loop over all affected degrees of current element.
	  for(int j=nodesPerSurfaceElem+1;j<dispBoundConds[i].size();j+=2) {
	    dof = (int)dispBoundConds[i][j]-1;
	    boundDOF[dof] = 1;
	    boundConds[dof] = dispBoundConds[i][j+1];
	  }
	  
	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(surfaceDispBoundElemIdx,indexSet);
	}
      
      }
      
    }

    // ------------------------------------------------------------------
    // rotation surface constraints
    for(int i=0;i<rotBoundConds.size();i++) {

      if(rotDOF <= 0) {
	logFile<<"In function FEMGeometry::createSurfaceElems "
	       <<"rotation boundary conditions are not supported!"<<endl;
	break;
	//MPI_Abort(MPI_COMM_WORLD,1);
      }

      // Check if any loading or displacement constraint is already 
      // applied to current element.
      elem = (int)rotBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no condition is already applied
      if(pos == -1) {

	currentIdx = boundaryElemIdx.size();
	surfaceNodesElems[currentIdx].setGlobalID(elem);
	
	// Store the rotation boundary conditions.
	blVector& boundDOF = 
	  surfaceNodesElems[currentIdx].getSurfaceDeformationBoundDOF();
	dbVector& boundConds = 
	  surfaceNodesElems[currentIdx].getSurfaceDeformationBoundConds();

	if(boundDOF.size() < defDOF) {
	  boundDOF.resize(defDOF);
	  boundConds.resize(defDOF);
	}
      
	// Loop over all affected degrees of current element.
	for(int j=nodesPerSurfaceElem+1;j<rotBoundConds[i].size();j+=2) {
	  dof = (int)rotBoundConds[i][j]-1;
	  boundDOF[dispDOF+dof] = 1;
	  boundConds[dispDOF+dof] = rotBoundConds[i][j+1];
	}

	// store the element's nodes
	intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
	nodes.resize(nodesPerSurfaceElem);

	for(int j=0;j<nodesPerSurfaceElem;j++) 
	  nodes[j] = newIdx[(int)rotBoundConds[i][j+1]-1]+1;
       

	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
	
	if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID
	
	pushBackVector(surfaceRotBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }
      // a loading or displacement constraint is already applied
      else {
	currentIdx = pos; 

	// Check if a surface rotation boundary conditions is already 
	// applied to current element.
	pos = findIntMatPos(currentIdx,0,surfaceRotBoundElemIdx.size(),0,
			    surfaceRotBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
		 <<"rotation boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the rotation boundary conditions.
	  blVector& boundDOF = 
	    surfaceNodesElems[currentIdx].getSurfaceDeformationBoundDOF();
	  dbVector& boundConds = 
	    surfaceNodesElems[currentIdx].getSurfaceDeformationBoundConds();

	  if(boundDOF.size() < defDOF) {
	    boundDOF.resize(defDOF);
	    boundConds.resize(defDOF);
	  }
	  
	  // Loop over all affected degrees of current element.
	  for(int j=nodesPerSurfaceElem+1;j<rotBoundConds[i].size();j+=2) {
	    dof = (int)rotBoundConds[i][j]-1;
	    boundDOF[dispDOF+dof] = 1;
	    boundConds[dispDOF+dof] = rotBoundConds[i][j+1];
	  }
	  
	  indexSet[0] = currentIdx; // element's position in surfaceNodesElems
	  indexSet[1] = 0; // condition ID

	  pushBackVector(surfaceRotBoundElemIdx,indexSet);
	}
      
      }
      
    }

    // ------------------------------------------------------------------

    // electric boundary constraints
    for(int i=0;i<electricBoundConds.size();i++) {


      // Check if any loading or deformation constraint applied to current 
      // element.
      elem = (int)electricBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading is already applied
      if(pos == -1) {

	currentIdx = boundaryElemIdx.size();
	surfaceNodesElems[currentIdx].setGlobalID(elem);
	
	// Store the electric boundary conditions.
	blVector& boundDOF = 
	  surfaceNodesElems[currentIdx].getSurfaceElectricBoundDOF();
	dbVector& boundConds = 
	  surfaceNodesElems[currentIdx].getSurfaceElectricBoundConds();
	boundDOF.resize(electricDOF);
	boundConds.resize(electricDOF);
      
	// Loop over all affected degrees of current element.
	for(int j=nodesPerSurfaceElem+1;j<electricBoundConds[i].size();j+=2) {
	  dof = (int)electricBoundConds[i][j]-1;

	  if(dof < boundDOF.size()) {
	    boundDOF[dof] = 1;
	    boundConds[dof] = electricBoundConds[i][j+1];
	  }
	  else {
	    boundDOF.resize(electricDOF+usedDims);
	    boundConds.resize(electricDOF+usedDims);

	    boundDOF[dof] = 1;
	    boundConds[dof] = electricBoundConds[i][j+1];
	  }

	}

	// store the element's nodes
	intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
	nodes.resize(nodesPerSurfaceElem);

	for(int j=0;j<nodesPerSurfaceElem;j++) 
	  nodes[j] = newIdx[(int)electricBoundConds[i][j+1]-1]+1;
       

	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
	
	if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID
	
	pushBackVector(surfaceElectricBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }
      // a loading is already applied
      else {
	currentIdx = pos; 

	// Check if a surface electric boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(currentIdx,0,surfaceElectricBoundElemIdx.size(),0,
			    surfaceElectricBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
		 <<"electric boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the electric boundary conditions.
	  blVector& boundDOF = 
	    surfaceNodesElems[currentIdx].getSurfaceElectricBoundDOF();
	  dbVector& boundConds = 
	    surfaceNodesElems[currentIdx].getSurfaceElectricBoundConds();
	  boundDOF.resize(electricDOF);
	  boundConds.resize(electricDOF);
	  
	  // Loop over all affected degrees of current element.
	  for(int j=nodesPerSurfaceElem+1;j<electricBoundConds[i].size();j+=2) {
	    dof = (int)electricBoundConds[i][j]-1;

	    if(dof < boundDOF.size()) {
	      boundDOF[dof] = 1;
	      boundConds[dof] = electricBoundConds[i][j+1];
	    }
	    else {
	      boundDOF.resize(electricDOF+usedDims);
	      boundConds.resize(electricDOF+usedDims);
	      
	      boundDOF[dof] = 1;
	      boundConds[dof] = electricBoundConds[i][j+1];
	    }

	  }
	  
	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(surfaceElectricBoundElemIdx,indexSet);
	}
      
      }
      
    }

    // ------------------------------------------------------------------

    // depolarisation boundary constraints
    for(int i=0;i<depolarisationBoundConds.size();i++) {


      // Check if any loading or deformation constraint applied to current 
      // element.
      elem = (int)depolarisationBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading is already applied
      if(pos == -1) {

	currentIdx = boundaryElemIdx.size();
	surfaceNodesElems[currentIdx].setGlobalID(elem);
	
	// Store the depolarisation boundary conditions.
	blVector& boundDOF = 
	  surfaceNodesElems[currentIdx].getSurfaceDepolarisationBoundDOF();
	dbVector& boundConds = 
	  surfaceNodesElems[currentIdx].getSurfaceDepolarisationBoundConds();
	boundDOF.resize(depolarisationDOF);
	boundConds.resize(depolarisationDOF);
      
	// Loop over all affected degrees of current element.
	for(int j=nodesPerSurfaceElem+1;j<depolarisationBoundConds[i].size();j+=2) {
	  dof = (int)depolarisationBoundConds[i][j]-1;

	  if(dof < boundDOF.size()) {
	    boundDOF[dof] = 1;
	    boundConds[dof] = depolarisationBoundConds[i][j+1];
	  }
	  else {
	    boundDOF.resize(depolarisationDOF+usedDims);
	    boundConds.resize(depolarisationDOF+usedDims);

	    boundDOF[dof] = 1;
	    boundConds[dof] = depolarisationBoundConds[i][j+1];
	  }

	}

	// store the element's nodes
	intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
	nodes.resize(nodesPerSurfaceElem);

	for(int j=0;j<nodesPerSurfaceElem;j++) 
	  nodes[j] = newIdx[(int)depolarisationBoundConds[i][j+1]-1]+1;
       

	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
	
	if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID
	
	pushBackVector(surfaceDepolarisationBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }
      // a loading is already applied
      else {
	currentIdx = pos; 

	// Check if a surface depolarisation boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(currentIdx,0,surfaceDepolarisationBoundElemIdx.size(),0,
			    surfaceDepolarisationBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
		 <<"depolarisation boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the depolarisation boundary conditions.
	  blVector& boundDOF = 
	    surfaceNodesElems[currentIdx].getSurfaceDepolarisationBoundDOF();
	  dbVector& boundConds = 
	    surfaceNodesElems[currentIdx].getSurfaceDepolarisationBoundConds();
	  boundDOF.resize(depolarisationDOF);
	  boundConds.resize(depolarisationDOF);
	  
	  // Loop over all affected degrees of current element.
	  for(int j=nodesPerSurfaceElem+1;j<depolarisationBoundConds[i].size();j+=2) {
	    dof = (int)depolarisationBoundConds[i][j]-1;

	    if(dof < boundDOF.size()) {
	      boundDOF[dof] = 1;
	      boundConds[dof] = depolarisationBoundConds[i][j+1];
	    }
	    else {
	      boundDOF.resize(depolarisationDOF+usedDims);
	      boundConds.resize(depolarisationDOF+usedDims);
	      
	      boundDOF[dof] = 1;
	      boundConds[dof] = depolarisationBoundConds[i][j+1];
	    }

	  }
	  
	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(surfaceDepolarisationBoundElemIdx,indexSet);
	}
      
      }
      
    }

    // ------------------------------------------------------------------

    // micro boundary constraints
    for(int i=0;i<microBoundConds.size();i++) {


      // Check if any loading or deformation constraint applied to current 
      // element.
      elem = (int)microBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading is already applied
      if(pos == -1) {

	currentIdx = boundaryElemIdx.size();
	surfaceNodesElems[currentIdx].setGlobalID(elem);
	
	// Store the micro boundary conditions.
	blVector& boundDOF = 
	  surfaceNodesElems[currentIdx].getSurfaceMicroBoundDOF();
	dbVector& boundConds = 
	  surfaceNodesElems[currentIdx].getSurfaceMicroBoundConds();
	boundDOF.resize(microDOF);
	boundConds.resize(microDOF);
      
	// Loop over all affected degrees of current element.
	for(int j=nodesPerSurfaceElem+1;j<microBoundConds[i].size();j+=2) {
	  dof = (int)microBoundConds[i][j]-1;
	  boundDOF[dof] = 1;
	  boundConds[dof] = microBoundConds[i][j+1];
	}

	// store the element's nodes
	intVector& nodes = surfaceNodesElems[currentIdx].getNodes();
	nodes.resize(nodesPerSurfaceElem);

	for(int j=0;j<nodesPerSurfaceElem;j++) 
	  nodes[j] = newIdx[(int)microBoundConds[i][j+1]-1]+1;
       

	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
	
	if(!swapVolToFaceElemInfo(InputData,surfaceNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createSurfaceElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID
	
	pushBackVector(surfaceMicroBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }
      // a loading is already applied
      else {
	currentIdx = pos; 

	// Check if a surface micro boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(currentIdx,0,surfaceMicroBoundElemIdx.size(),0,
			    surfaceMicroBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createSurfaceElems repetetive "
		 <<"micro boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the micro boundary conditions.
	  blVector& boundDOF = 
	    surfaceNodesElems[currentIdx].getSurfaceMicroBoundDOF();
	  dbVector& boundConds = 
	    surfaceNodesElems[currentIdx].getSurfaceMicroBoundConds();
	  boundDOF.resize(microDOF);
	  boundConds.resize(microDOF);
	  
	  // Loop over all affected degrees of current element.
	  for(int j=nodesPerSurfaceElem+1;j<microBoundConds[i].size();j+=2) {
	    dof = (int)microBoundConds[i][j]-1;
	    boundDOF[dof] = 1;
	    boundConds[dof] = microBoundConds[i][j+1];
	  }
	  
	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(surfaceMicroBoundElemIdx,indexSet);
	}
      
      }
      
    }

  }


  surfaceNodesElems.resize(boundaryElemIdx.size(),FEMElement(usedDOF));

#ifdef _FEdebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"****************** surface elements ******************"<<endl;
  logFile<<"******************************************************"<<endl;
  logFile<<"*********************** nodes ************************"<<endl;
  for(int i=0;i<surfaceNodesElems.size();i++) {
    logFile<<i<<".) SURFACE ELEMENT "
	   <<surfaceNodesElems[i].getGlobalID()<<": ";
    intVector& nodes = surfaceNodesElems[i].getNodes();
    for(int j=0;j<nodes.size();j++) {
      logFile<<nodes[j]<<" ";
    }
    int& elemMatID = surfaceNodesElems[i].getMaterialID();
    int& elemType = surfaceNodesElems[i].getElemType();
    int& elemOrder = surfaceNodesElems[i].getElemOrder();
    int& volID = surfaceNodesElems[i].getMotherElementID();
    logFile<<"elemType "<<elemType<<" elemOrder "<<elemOrder
	   <<" volID "<<volID<<" matID "<<elemMatID<<endl;
  }
  logFile<<"***************** traction loads *********************"<<endl;
  for(int i=0;i<tractionElemIdx.size();i++) {
    int& idx = tractionElemIdx[i][0];
    int& volID = surfaceNodesElems[idx].getMotherElementID();
    blVector& affectedDOF = surfaceNodesElems[idx].getTractionDOF();
    dbVector& conditions = surfaceNodesElems[idx].getTraction();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<" volID "<<volID<<endl;
  }
  logFile<<"************** surface pressure loads ****************"<<endl;
  for(int i=0;i<surfacePressureElemIdx.size();i++) {
    int& idx = surfacePressureElemIdx[i][0];
    double& condition = surfaceNodesElems[idx].getSurfacePressure();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()
	   <<": surface pressure = "<<condition<<endl;
  }
  logFile<<"************** surface moment loads ******************"<<endl;
  for(int i=0;i<surfaceMomentElemIdx.size();i++) {
    int& idx = surfaceMomentElemIdx[i][0];
    blVector& affectedDOF = surfaceNodesElems[idx].getSurfaceMomentDOF();
    dbVector& conditions = surfaceNodesElems[idx].getSurfaceMoment();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"*********** surface electric charge loads *************"<<endl;
  for(int i=0;i<surfaceElectricChargeElemIdx.size();i++) {
    int& idx = surfaceElectricChargeElemIdx[i][0];
    blVector& affectedDOF = surfaceNodesElems[idx].getSurfaceElectricChargeDOF();
    dbVector& conditions = surfaceNodesElems[idx].getSurfaceElectricCharge();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"********* displacement boundary conditions ***********"<<endl;
  for(int i=0;i<surfaceDispBoundElemIdx.size();i++) {
    int& idx = surfaceDispBoundElemIdx[i][0];
    int& volID = surfaceNodesElems[idx].getMotherElementID();
    blVector& affectedDOF = 
      surfaceNodesElems[idx].getSurfaceDeformationBoundDOF();
    dbVector& conditions = 
      surfaceNodesElems[idx].getSurfaceDeformationBoundConds();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<" volID "<<volID<<endl;
  }
  logFile<<"************* rotation boundary conditions ***********"<<endl;
  for(int i=0;i<surfaceRotBoundElemIdx.size();i++) {
    int& idx = surfaceRotBoundElemIdx[i][0];
    blVector& affectedDOF = 
      surfaceNodesElems[idx].getSurfaceDeformationBoundDOF();
    dbVector& conditions = 
      surfaceNodesElems[idx].getSurfaceDeformationBoundConds();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"************* electric boundary conditions ***********"<<endl;
  for(int i=0;i<surfaceElectricBoundElemIdx.size();i++) {
    int& idx = surfaceElectricBoundElemIdx[i][0];
    blVector& affectedDOF = 
      surfaceNodesElems[idx].getSurfaceElectricBoundDOF();
    dbVector& conditions = 
      surfaceNodesElems[idx].getSurfaceElectricBoundConds();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"************* depolarisation boundary conditions ***********"<<endl;
  for(int i=0;i<surfaceDepolarisationBoundElemIdx.size();i++) {
    int& idx = surfaceDepolarisationBoundElemIdx[i][0];
    blVector& affectedDOF = 
      surfaceNodesElems[idx].getSurfaceDepolarisationBoundDOF();
    dbVector& conditions = 
      surfaceNodesElems[idx].getSurfaceDepolarisationBoundConds();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
  logFile<<"************* micro boundary conditions ***********"<<endl;
  for(int i=0;i<surfaceMicroBoundElemIdx.size();i++) {
    int& idx = surfaceMicroBoundElemIdx[i][0];
    blVector& affectedDOF = 
      surfaceNodesElems[idx].getSurfaceMicroBoundDOF();
    dbVector& conditions = 
      surfaceNodesElems[idx].getSurfaceMicroBoundConds();
    logFile<<"SURFACE ELEMENT "<<surfaceNodesElems[idx].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<endl;
  }
#endif
  
}

/**********************************************************************/
/**********************************************************************/
// Create line elements and store their line load.
void FEMGeometry::createLineElems(InputFileData* InputData,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile) {

  using namespace std;

  std::map<std::string,double>& backGroundMeshInfo =
    InputData->getBackGroundMeshInfo();

  bool dirichletIntegration = 
    (bool)modelData["dirichletBoundaryIntegration"];
  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];
  int dispDOF = (int)modelData["displacementDegreesOfFreedom"];
  int rotDOF = (int)modelData["rotationDegreesOfFreedom"];
  int electricDOF = (int)modelData["electricDegreesOfFreedom"];
  int depolarisationDOF = (int)modelData["depolarisationDegreesOfFreedom"];
  int nodesPerLineElem = (int)backGroundMeshInfo["nodesPerLineElement"];

  dbMatrix& lineForceLoads = InputData->getLineForceLoads();
  dbMatrix& lineMomentLoads = InputData->getLineMomentLoads();

  dbMatrix& dispBoundConds = InputData->getLineDispBoundConds();
  dbMatrix& rotBoundConds = InputData->getLineRotBoundConds();
  dbMatrix& electricBoundConds = InputData->getLineElectricBoundConds();
  dbMatrix& depolarisationBoundConds = InputData->getLineDepolarisationBoundConds();


  /**********************************************************************/
  // delete not used boundary conditions

  if(rotDOF == 0) {

    rotBoundConds.resize(0);
    lineMomentLoads.resize(0);
  }

  if(electricDOF == 0)

    electricBoundConds.resize(0);


  if(depolarisationDOF == 0)

    depolarisationBoundConds.resize(0);


  /**********************************************************************/

  int lineElemNum = lineForceLoads.size() + dispBoundConds.size()+
    rotBoundConds.size() + electricBoundConds.size() + depolarisationBoundConds.size();

  if(lineMomentLoads.size() > 0) {
    logFile<<"In function FEMGeometry::createLineElems "
	   <<"line moments are not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Create the vector containing all line boundary elements.
  lineNodesElems = 
    vector<FEMElement>(lineElemNum,FEMElement(usedDOF));

  // Loop over all boundary elements applied a line load.
  int currentIdx,dof,elem,pos,vElem;
  intVector boundaryElemIdx;

  int idxNormal = nodesPerLineElem+1;
  int startIdxConds = nodesPerLineElem+4;

  intVector indexSet(2);

  // Loop over all elements line loads are applied 
  for(int i=0;i<lineForceLoads.size();i++) {

    elem = (int)lineForceLoads[i][0]-1;

    pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			boundaryElemIdx);

    if(pos != -1) {
      logFile<<"In function FEMGeometry::createSurfaceElems repetetive\n" 
	     <<"line force conditions!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {

      currentIdx = boundaryElemIdx.size();

      lineNodesElems[currentIdx].setGlobalID(elem);

      // Store the the line load.
      blVector& boundDOF = lineNodesElems[currentIdx].getLineForceDOF();
      dbVector& lineForce = lineNodesElems[currentIdx].getLineForce();
      boundDOF.resize(usedDims);
      lineForce.resize(usedDims);
    
      for(int j=startIdxConds;j<lineForceLoads[i].size();j+=2) {
	dof = (int)lineForceLoads[i][j]-1;
	boundDOF[dof] = true;
	lineForce[dof] = lineForceLoads[i][j+1];
      }
      
      // Store the normal vector.
      dbVector& normal = lineNodesElems[currentIdx].getSurfaceNormal();
      normal.resize(usedDims);

      normal[0] = lineForceLoads[i][idxNormal];
      normal[1] = lineForceLoads[i][idxNormal+1];
      normal[2] = lineForceLoads[i][idxNormal+2];

      // Store the node numbers.
      intVector& nodes = lineNodesElems[currentIdx].getNodes();
      nodes.resize(nodesPerLineElem);

      for(int j=0;j<nodesPerLineElem;j++) 
	nodes[j] = newIdx[(int)lineForceLoads[i][j+1]-1]+1;


      // set the element's material, type (tetrahedra,cube) 
      // and order (linear/quadratic)
      
      if(!swapVolToFaceElemInfo(InputData,lineNodesElems[currentIdx],
				modelData,logFile)) {
	logFile<<"In function FEMGeometry::createLineElems no corresponding\n" 
	       <<"volume element can be found for surface element '"<<elem<<"'!"
	       <<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }


      indexSet[0] = currentIdx;
      indexSet[1] = 0; // loading ID

      pushBackVector(lineForceElemIdx,indexSet);
      pushBackVector(boundaryElemIdx,elem);
    }
 
  }

  /*********************************************************************/
  // Loop over all line elements deformation boundary conditions are
  // applied, if and integration of Dirichlet boundary conditions is 
  // desired.
  if(dirichletIntegration) {
    
    // loop over all displacement boundary conditions
    for(int i=0;i<dispBoundConds.size();i++) {
      
      if(dispDOF <= 0) {
	logFile<<"In function FEMGeometry::createLineElems "
	       <<"displacement boundary conditions are not supported!"<<endl;
	break;
	//MPI_Abort(MPI_COMM_WORLD,1);
      }

      // Check if a load is already applied to current elment.
      elem = (int)dispBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading is already applied
      if(pos == -1) {
	currentIdx = boundaryElemIdx.size();
	lineNodesElems[currentIdx].setGlobalID(elem);

        // Store the displacement boundary conditions.
        blVector& boundDOF = 
	  lineNodesElems[currentIdx].getLineDeformationBoundDOF();
        dbVector& boundConds = 
	  lineNodesElems[currentIdx].getLineDeformationBoundConds();
	boundDOF.resize(defDOF);
	boundConds.resize(defDOF);
     
        // Loop over all affected degrees of current element.
        for(int j=startIdxConds;j<dispBoundConds[i].size();j+=2) {
	  dof = (int)dispBoundConds[i][j]-1;
	  boundDOF[dof] = true;
	  boundConds[dof] = dispBoundConds[i][j+1];
        }

	// Store the normal vector.
	dbVector& normal = lineNodesElems[currentIdx].getSurfaceNormal();
	normal.resize(usedDims);
	
	normal[0] = dispBoundConds[i][idxNormal];
	normal[1] = dispBoundConds[i][idxNormal+1];
	normal[2] = dispBoundConds[i][idxNormal+2];
	
        // Store the node numbers.
        intVector& nodes = lineNodesElems[currentIdx].getNodes();
        nodes.resize(nodesPerLineElem);
      
        for(int j=0;j<nodes.size();j++) 
	  nodes[j] = newIdx[(int)dispBoundConds[i][j+1]-1]+1;


	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
      
	if(!swapVolToFaceElemInfo(InputData,lineNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createLineElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID

	pushBackVector(lineDispBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }

      // a loading is already applied
      else { 
	currentIdx = pos; 
      
	// Check if a line displacement boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(elem,0,lineDispBoundElemIdx.size(),0,
			    lineDispBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createLineElems repetetive "
		 <<"line displacement boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the displacement boundary conditions.
	  blVector& boundDOF = 
	    lineNodesElems[currentIdx].getLineDeformationBoundDOF();
	  dbVector& boundConds = 
	    lineNodesElems[currentIdx].getLineDeformationBoundConds();

	  if(boundDOF.size() < defDOF) {
	    boundDOF.resize(defDOF);
	    boundConds.resize(defDOF);
	  }
	  
	  // Loop over all affected degrees of current element.
	  for(int j=startIdxConds;j<dispBoundConds[i].size();j+=2) {
	    dof = (int)dispBoundConds[i][j]-1;
	    boundDOF[dof] = 1;
	    boundConds[dof] = dispBoundConds[i][j+1];
	  }

	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(lineDispBoundElemIdx,indexSet);
	}
      
      }
      
    }

    // ------------------------------------------------------------------
    // loop over all rotation boundary conditions
    for(int i=0;i<rotBoundConds.size();i++) {
      
      if(rotDOF <= 0) {
	logFile<<"In function FEMGeometry::createLineElems "
	       <<"rotation boundary conditions are not supported!"<<endl;
	break;
	//MPI_Abort(MPI_COMM_WORLD,1);
      }

      // Check if a load is already applied to current elment.
      elem = (int)rotBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading or displacement boundary is already applied
      if(pos == -1) {
	currentIdx = boundaryElemIdx.size();
	lineNodesElems[currentIdx].setGlobalID(elem);

        // Store the rotation boundary conditions.
        blVector& boundDOF = 
	  lineNodesElems[currentIdx].getLineDeformationBoundDOF();
        dbVector& boundConds = 
	  lineNodesElems[currentIdx].getLineDeformationBoundConds();
	boundDOF.resize(defDOF);
	boundConds.resize(defDOF);
     
        // Loop over all affected degrees of current element.
        for(int j=startIdxConds;j<rotBoundConds[i].size();j+=2) {
	  dof = (int)rotBoundConds[i][j]-1;
	  boundDOF[dispDOF+dof] = 1;
	  boundConds[dispDOF+dof] = rotBoundConds[i][j+1];
        }

	// Store the normal vector.
	dbVector& normal = lineNodesElems[currentIdx].getSurfaceNormal();
	normal.resize(usedDims);
	
	normal[0] = rotBoundConds[i][idxNormal];
	normal[1] = rotBoundConds[i][idxNormal+1];
	normal[2] = rotBoundConds[i][idxNormal+2];
	
        // Store the node numbers.
        intVector& nodes = lineNodesElems[currentIdx].getNodes();
        nodes.resize(nodesPerLineElem);
      
        for(int j=0;j<nodes.size();j++) 
	  nodes[j] = newIdx[(int)rotBoundConds[i][j+1]-1]+1;


	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
      
	if(!swapVolToFaceElemInfo(InputData,lineNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createLineElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID

	pushBackVector(lineRotBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }

      // a loading or displacement boundary condition is already applied
      else { 
	currentIdx = pos; 
      
	// Check if a line rotation boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(elem,0,lineRotBoundElemIdx.size(),0,
			    lineRotBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createLineElems repetetive "
		 <<"line rotation boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the rotation boundary conditions.
	  blVector& boundDOF = 
	    lineNodesElems[currentIdx].getLineDeformationBoundDOF();
	  dbVector& boundConds = 
	    lineNodesElems[currentIdx].getLineDeformationBoundConds();
	  
	  if(boundDOF.size() < defDOF) {
	    boundDOF.resize(defDOF);
	    boundConds.resize(defDOF);
	  }
	  
	  // Loop over all affected degrees of current element.
	  for(int j=startIdxConds;j<rotBoundConds[i].size();j+=2) {
	    dof = (int)rotBoundConds[i][j]-1;
	    boundDOF[dispDOF+dof] = 1;
	    boundConds[dispDOF+dof] = rotBoundConds[i][j+1];
	  }

	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(lineRotBoundElemIdx,indexSet);
	}
      
      }
      
    }

    // ------------------------------------------------------------------

    // electric boundary constraints
    for(int i=0;i<electricBoundConds.size();i++) {


      // Check if any loading or deformation constraint applied to current 
      // element.
      elem = (int)electricBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading is already applied
      if(pos == -1) {

	currentIdx = boundaryElemIdx.size();
	lineNodesElems[currentIdx].setGlobalID(elem);
	
	// Store the electric boundary conditions.
	blVector& boundDOF = 
	  lineNodesElems[currentIdx].getLineElectricBoundDOF();
	dbVector& boundConds = 
	  lineNodesElems[currentIdx].getLineElectricBoundConds();
	boundDOF.resize(electricDOF);
	boundConds.resize(electricDOF);
      
	// Loop over all affected degrees of current element.
	for(int j=nodesPerLineElem+1;j<electricBoundConds[i].size();j+=2) {
	  dof = (int)electricBoundConds[i][j]-1;
	  boundDOF[dof] = 1;
	  boundConds[dof] = electricBoundConds[i][j+1];
	}

	// Store the normal vector.
	dbVector& normal = lineNodesElems[currentIdx].getSurfaceNormal();
	normal.resize(usedDims);
	
	normal[0] = electricBoundConds[i][idxNormal];
	normal[1] = electricBoundConds[i][idxNormal+1];
	normal[2] = electricBoundConds[i][idxNormal+2];

	// store the element's nodes
	intVector& nodes = lineNodesElems[currentIdx].getNodes();
	nodes.resize(nodesPerLineElem);

	for(int j=0;j<nodesPerLineElem;j++) 
	  nodes[j] = newIdx[(int)electricBoundConds[i][j+1]-1]+1;
       

	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
      
	if(!swapVolToFaceElemInfo(InputData,lineNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createLineElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID
	
	pushBackVector(lineElectricBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }
      // a loading is already applied
      else {
	currentIdx = pos; 

	// Check if a line electric boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(currentIdx,0,lineElectricBoundElemIdx.size(),0,
			    lineElectricBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createLineElems repetetive "
		 <<"electric boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the electric boundary conditions.
	  blVector& boundDOF = 
	    lineNodesElems[currentIdx].getLineElectricBoundDOF();
	  dbVector& boundConds = 
	    lineNodesElems[currentIdx].getLineElectricBoundConds();
	  boundDOF.resize(electricDOF);
	  boundConds.resize(electricDOF);
	  
	  // Loop over all affected degrees of current element.
	  for(int j=nodesPerLineElem+1;j<electricBoundConds[i].size();j+=2) {
	    dof = (int)electricBoundConds[i][j]-1;
	    boundDOF[dof] = 1;
	    boundConds[dof] = electricBoundConds[i][j+1];
	  }
	  
	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(lineElectricBoundElemIdx,indexSet);
	}
      
      }
      
    }
    
    // ------------------------------------------------------------------

    // depolarisation boundary constraints
    for(int i=0;i<depolarisationBoundConds.size();i++) {


      // Check if any loading or deformation constraint applied to current 
      // element.
      elem = (int)depolarisationBoundConds[i][0]-1;

      pos = findIntVecPos(elem,0,boundaryElemIdx.size(),
			  boundaryElemIdx);

      // no loading is already applied
      if(pos == -1) {

	currentIdx = boundaryElemIdx.size();
	lineNodesElems[currentIdx].setGlobalID(elem);
	
	// Store the depolarisation boundary conditions.
	blVector& boundDOF = 
	  lineNodesElems[currentIdx].getLineDepolarisationBoundDOF();
	dbVector& boundConds = 
	  lineNodesElems[currentIdx].getLineDepolarisationBoundConds();
	boundDOF.resize(depolarisationDOF);
	boundConds.resize(depolarisationDOF);
      
	// Loop over all affected degrees of current element.
	for(int j=nodesPerLineElem+1;j<depolarisationBoundConds[i].size();j+=2) {
	  dof = (int)depolarisationBoundConds[i][j]-1;
	  boundDOF[dof] = 1;
	  boundConds[dof] = depolarisationBoundConds[i][j+1];
	}

	// Store the normal vector.
	dbVector& normal = lineNodesElems[currentIdx].getSurfaceNormal();
	normal.resize(usedDims);
	
	normal[0] = depolarisationBoundConds[i][idxNormal];
	normal[1] = depolarisationBoundConds[i][idxNormal+1];
	normal[2] = depolarisationBoundConds[i][idxNormal+2];

	// store the element's nodes
	intVector& nodes = lineNodesElems[currentIdx].getNodes();
	nodes.resize(nodesPerLineElem);

	for(int j=0;j<nodesPerLineElem;j++) 
	  nodes[j] = newIdx[(int)depolarisationBoundConds[i][j+1]-1]+1;
       

	// set the element's material, type (tetrahedra,cube) 
	// and order (linear/quadratic)
      
	if(!swapVolToFaceElemInfo(InputData,lineNodesElems[currentIdx],
				  modelData,logFile)) {
	  logFile<<"In function FEMGeometry::createLineElems no corresponding\n" 
		 <<"volume element can be found for surface element '"<<elem<<"'!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}


	indexSet[0] = currentIdx;
	indexSet[1] = 0; // condition ID
	
	pushBackVector(lineDepolarisationBoundElemIdx,indexSet);
	pushBackVector(boundaryElemIdx,elem);
      }
      // a loading is already applied
      else {
	currentIdx = pos; 

	// Check if a line depolarisation boundary conditions is already 
	// applied to current elment.
	pos = findIntMatPos(currentIdx,0,lineDepolarisationBoundElemIdx.size(),0,
			    lineDepolarisationBoundElemIdx);
      
	if(pos != -1) {
	  logFile<<"In function FEMGeometry::createLineElems repetetive "
		 <<"depolarisation boundary conditions!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

	else {

	  // Store the depolarisation boundary conditions.
	  blVector& boundDOF = 
	    lineNodesElems[currentIdx].getLineDepolarisationBoundDOF();
	  dbVector& boundConds = 
	    lineNodesElems[currentIdx].getLineDepolarisationBoundConds();
	  boundDOF.resize(depolarisationDOF);
	  boundConds.resize(depolarisationDOF);
	  
	  // Loop over all affected degrees of current element.
	  for(int j=nodesPerLineElem+1;j<depolarisationBoundConds[i].size();j+=2) {
	    dof = (int)depolarisationBoundConds[i][j]-1;
	    boundDOF[dof] = 1;
	    boundConds[dof] = depolarisationBoundConds[i][j+1];
	  }
	  
	  indexSet[0] = currentIdx;
	  indexSet[1] = 0; // condition ID

	  pushBackVector(lineDepolarisationBoundElemIdx,indexSet);
	}
      
      }
      
    }

  }

  lineNodesElems.resize(boundaryElemIdx.size(),FEMElement(usedDOF));


#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"******************** line elements ******************"<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"dispDOF = "<<dispDOF<<endl;
  logFile<<"*********************** nodes ***********************"<<endl;
  for(int i=0;i<lineNodesElems.size();i++) {
    logFile<<i<<".) LINE ELEMENT "
	   <<lineNodesElems[i].getGlobalID()<<": ";
    intVector& nodes = lineNodesElems[i].getNodes();
    for(int j=0;j<nodes.size();j++) {
      logFile<<nodes[j]<<" ";
    }
    int& elemMatID = lineNodesElems[i].getMaterialID();
    int& elemType = lineNodesElems[i].getElemType();
    int& elemOrder = lineNodesElems[i].getElemOrder();
    logFile<<"elemType "<<elemType<<" elemOrder "<<elemOrder
	   <<" matID "<<elemMatID<<endl;
  }
  logFile<<"*************** line force loads *********************"<<endl;
  for(int i=0;i<lineForceElemIdx.size();i++) {
    int& idx = lineForceElemIdx[i][0];
    blVector& affectedDOF = lineNodesElems[idx].getLineForceDOF();
    dbVector& conditions = lineNodesElems[idx].getLineForce();
    dbVector& normal = lineNodesElems[idx].getSurfaceNormal();
    logFile<<"LINE ELEMENT "<<lineNodesElems[i].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<endl;
  }

  logFile<<"********** displacement boundary conditions **********"<<endl;
  for(int i=0;i<lineDispBoundElemIdx.size();i++) {
    int& idx = lineDispBoundElemIdx[i][0];
    blVector& affectedDOF = 
      lineNodesElems[idx].getLineDeformationBoundDOF();
    dbVector& conditions = 
      lineNodesElems[idx].getLineDeformationBoundConds();
    dbVector& normal = lineNodesElems[idx].getSurfaceNormal();
    logFile<<"LINE ELEMENT "<<lineNodesElems[i].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<endl;
  }
  logFile<<"************ rotation boundary conditions ************"<<endl;
  for(int i=0;i<lineRotBoundElemIdx.size();i++) {
    int& idx = lineRotBoundElemIdx[i][0];
    blVector& affectedDOF = 
      lineNodesElems[idx].getLineDeformationBoundDOF();
    dbVector& conditions = 
      lineNodesElems[idx].getLineDeformationBoundConds();
    dbVector& normal = lineNodesElems[idx].getSurfaceNormal();
    logFile<<"LINE ELEMENT "<<lineNodesElems[i].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<endl;
  }
  logFile<<"************ electric boundary conditions ************"<<endl;
  for(int i=0;i<lineElectricBoundElemIdx.size();i++) {
    int& idx = lineElectricBoundElemIdx[i][0];
    blVector& affectedDOF = 
      lineNodesElems[idx].getLineElectricBoundDOF();
    dbVector& conditions = 
      lineNodesElems[idx].getLineElectricBoundConds();
    dbVector& normal = lineNodesElems[idx].getSurfaceNormal();
    logFile<<"LINE ELEMENT "<<lineNodesElems[i].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<endl;
  }

  logFile<<"************ depolarisation boundary conditions ************"<<endl;
  for(int i=0;i<lineDepolarisationBoundElemIdx.size();i++) {
    int& idx = lineDepolarisationBoundElemIdx[i][0];
    blVector& affectedDOF = 
      lineNodesElems[idx].getLineDepolarisationBoundDOF();
    dbVector& conditions = 
      lineNodesElems[idx].getLineDepolarisationBoundConds();
    dbVector& normal = lineNodesElems[idx].getSurfaceNormal();
    logFile<<"LINE ELEMENT "<<lineNodesElems[i].getGlobalID()<<": ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Rearrange the element-node configuration and set corresponding FEM-
// approximation set and Gauss quadrature set..
void FEMGeometry::rearrangeMeshConfiguration(InputFileData* InputData,
					     std::map<std::string,double>& modelData,
					     std::ofstream& logFile,
					     PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  bool flag;
  string mode = "arbitrary-subvector";
  intVector vElemIdx;
  intMatrix lElemTOvElemIdx,sElemTOvElemIdx;
  intVector IDs(2);

#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"********** rearrange mesh configuration *************"<<endl;
  logFile<<"*****************************************************"<<endl;
#endif


  // loop over all surface elements
  for(int i=0;i<surfaceNodesElems.size();i++) {
    intVector& sNodes = surfaceNodesElems[i].getNodes();

#ifdef _FEdebugMode_
    logFile<<"sElem "<<i<<": ";
    for(int k=0;k<sNodes.size();k++)
      logFile<<sNodes[k]<<" ";
    logFile<<endl;
#endif

    // loop over all volume elements in order to find the volume element 
    // to which the surface element belongs
    for(int j=0;j<nodesElements.size();j++) {
      intVector& vNodes = nodesElements[j].getNodes();

      flag = compareIntVecs(mode,sNodes,vNodes,logFile);

#ifdef _FEdebugMode_
	logFile<<"-> vElem "<<j<<": ";
	for(int k=0;k<vNodes.size();k++)
	  logFile<<vNodes[k]<<" ";
	logFile<<endl;
	if(flag)
	  logFile<<"check"<<endl;
#endif

      if(flag) {
	IDs[0] = i;
	IDs[1] = j;
	pushBackVector(vElemIdx,j);
	pushBackVector(sElemTOvElemIdx,IDs);
	break;
      }

    }

  }

  // loop over all line elements
  for(int i=0;i<lineNodesElems.size();i++) {
    intVector& lNodes = lineNodesElems[i].getNodes();

#ifdef _FEdebugMode_
	logFile<<"lElem "<<i<<": ";
	for(int k=0;k<lNodes.size();k++)
	  logFile<<lNodes[k]<<" ";
	logFile<<endl;
#endif

    // loop over all volume elements in order to find the volume element 
    // to which the line element belongs
    for(int j=0;j<nodesElements.size();j++) {
      intVector& vNodes = nodesElements[j].getNodes();

      flag = compareIntVecs(mode,lNodes,vNodes,logFile);

#ifdef _FEdebugMode_
	logFile<<"-> vElem "<<j<<": ";
	for(int k=0;k<vNodes.size();k++)
	  logFile<<vNodes[k]<<" ";
	logFile<<endl;
	if(flag)
	  logFile<<"check"<<endl;
#endif

      if(flag) {
	IDs[0] = i;
	IDs[1] = j;
	pushBackVector(vElemIdx,j);
	pushBackVector(lElemTOvElemIdx,IDs);
	break;
      }

    }

  }
  

  // remove duplex entries
  sortIntVector(vElemIdx,0,vElemIdx.size()-1);

  for(int i=1;i<vElemIdx.size();i++)
    
    if(vElemIdx[i-1] == vElemIdx[i]) {
      vElemIdx.erase(vElemIdx.begin()+i);
      --i;
      }


#ifdef _FEdebugMode_
  int position;
  logFile<<"************* affected volume elements **************"<<endl;
  for(int i=0;i<vElemIdx.size();i++)
    logFile<<i<<".) VOLUME ELEMENT "<<vElemIdx[i]<<endl;
  logFile<<"********** surface to volume element conn ***********"<<endl;
  for(int i=0;i<sElemTOvElemIdx.size();i++) {
    int& sID = sElemTOvElemIdx[i][0];
    int& vID = sElemTOvElemIdx[i][1];
    logFile<<i<<".) SURFACE ELEMENT "<<surfaceNodesElems[sID].getGlobalID()
	   <<" -> volume element "<<vID<<endl;
    position = findIntVecPos(vID,0,vElemIdx.size(),
			     vElemIdx);
    if(position == -1) {
      logFile<<"In FEMGeometry::rearrangeMeshConfiguration "
	     <<"volume-surface-element\n connectivity list "
	     <<"is incorrect!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  logFile<<"************ line to volume element conn ************"<<endl;
  for(int i=0;i<lElemTOvElemIdx.size();i++) {
    int& lID = lElemTOvElemIdx[i][0];
    int& vID = lElemTOvElemIdx[i][1];
    logFile<<i<<".) LINE ELEMENT "<<lineNodesElems[lID].getGlobalID()
	   <<" -> volume element "<<vID<<endl;
    position = findIntVecPos(vID,0,vElemIdx.size(),
			     vElemIdx);
    if(position == -1) {
      logFile<<"In FEMGeometry::rearrangeMeshConfiguration "
	     <<"volume-line-element\n connectivity list "
	     <<"is incorrect!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
#endif

  /**********************************************************************/
  // set the backgroundmesh FEM parameter and tools for all elements
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int pos;
  intVector localRootList;

  lineNodesElems.size() > surfaceNodesElems.size() ? 
    localRootList.resize(lineNodesElems.size()) : 
    localRootList.resize(surfaceNodesElems.size());

  localRootList.assign(localRootList.size(),-1);

  intVector surfaceRootList(surfaceNodesElems.size());

  // loop over all (local) surface elements to sycronize their element
  // parameters with these of the volume element they belong to
  for(int i=0;i<sElemTOvElemIdx.size();i++) {
    
    int& sID = sElemTOvElemIdx[i][0];
    int& vID = sElemTOvElemIdx[i][1];

    int& elemMatID = surfaceNodesElems[sID].getMaterialID();
    int& elemType = surfaceNodesElems[sID].getElemType();
    int& elemOrder = surfaceNodesElems[sID].getElemOrder();
    
    elemMatID = nodesElements[vID].getMaterialID();
    elemType = nodesElements[vID].getElemType();
    elemOrder = nodesElements[vID].getElemOrder();

    localRootList[sID] = rank;
    
  }

  // set up a processor root list for the surface elements
  intVector recvdIDs(size);

  for(int i=0;i<surfaceRootList.size();i++) {

    MPI_Allgather(&localRootList[i],1,MPI_INT,&recvdIDs[0],1,MPI_INT,
		  MPI_COMM_WORLD);

    for(int j=0;j<recvdIDs.size();j++)

      if(recvdIDs[j] > -1) {
	surfaceRootList[i] = recvdIDs[j];
	break;
      }

  }

  // data exchange between processors
  for(int i=0;i<surfaceRootList.size();i++) {
    int& elemMatID = surfaceNodesElems[i].getMaterialID();
    int& elemType = surfaceNodesElems[i].getElemType();
    int& elemOrder = surfaceNodesElems[i].getElemOrder();
    MPI_Bcast(&elemMatID,1,MPI_INT,surfaceRootList[i],MPI_COMM_WORLD);  
    MPI_Bcast(&elemType,1,MPI_INT,surfaceRootList[i],MPI_COMM_WORLD);  
    MPI_Bcast(&elemOrder,1,MPI_INT,surfaceRootList[i],MPI_COMM_WORLD); 
  }  

#ifdef _FEdebugMode_
    logFile<<"**************** surface elements ******************"<<endl;
    for(int i=0;i<surfaceNodesElems.size();i++) {
      int& elemMatID = surfaceNodesElems[i].getMaterialID();
      int& elemType = surfaceNodesElems[i].getElemType();
      int& elemOrder = surfaceNodesElems[i].getElemOrder();
      logFile<<"SURFACE ELEMENT "<<i<<" elemType "<<elemType
	     <<" elemOrder "<<elemOrder<<" matID "<<elemMatID
	     <<" root "<<surfaceRootList[i]<<endl;
    }
#endif

  // --------------------------------------------------------------------
  // loop over all (local) line elements to sycronize their element
  // parameters with these of the volume element they belong to

  localRootList.assign(localRootList.size(),-1);
  intVector lineRootList(lineNodesElems.size());

  for(int i=0;i<lElemTOvElemIdx.size();i++) {
   
    int& lID = lElemTOvElemIdx[i][0];
    int& vID = lElemTOvElemIdx[i][1];
    int& elemMatID = lineNodesElems[lID].getMaterialID();
    int& elemType = lineNodesElems[lID].getElemType();
    int& elemOrder = lineNodesElems[lID].getElemOrder();
    
    elemMatID = nodesElements[vID].getMaterialID();
    elemType = nodesElements[vID].getElemType();
    elemOrder = nodesElements[vID].getElemOrder();

    localRootList[lID] = rank;
 
  }

  // set up a processor root list for the line elements
  for(int i=0;i<lineRootList.size();i++) {

    MPI_Allgather(&localRootList[i],1,MPI_INT,&recvdIDs[0],1,MPI_INT,
		  MPI_COMM_WORLD);

    for(int j=0;j<recvdIDs.size();j++)

      if(recvdIDs[j] > -1) {
	lineRootList[i] = recvdIDs[j];
	break;
      }

  }

  // data exchange between processors
  for(int i=0;i<lineRootList.size();i++) {   
    int& elemMatID = lineNodesElems[i].getMaterialID();
    int& elemType = lineNodesElems[i].getElemType();
    int& elemOrder = lineNodesElems[i].getElemOrder();
    MPI_Bcast(&elemMatID,1,MPI_INT,lineRootList[i],MPI_COMM_WORLD);  
    MPI_Bcast(&elemType,1,MPI_INT,lineRootList[i],MPI_COMM_WORLD);  
    MPI_Bcast(&elemOrder,1,MPI_INT,lineRootList[i],MPI_COMM_WORLD); 
  }  

#ifdef _FEdebugMode_
    logFile<<"**************** line elements ******************"<<endl;
    for(int i=0;i<lineNodesElems.size();i++) {
      int& elemMatID = lineNodesElems[i].getMaterialID();
      int& elemType = lineNodesElems[i].getElemType();
      int& elemOrder = lineNodesElems[i].getElemOrder();
      logFile<<"LINE ELEMENT "<<i<<" elemType "<<elemType
	     <<" elemOrder "<<elemOrder<<" matID "<<elemMatID
	     <<" root "<<lineRootList[i]<<endl;
    }
#endif

  // Ensure that all elements have the necessary approximation and
  // integration tools assigned.
  setApproxAndIntTools(InputData,modelData,logFile);

  /*********************************************************************/
  // increase the element order on the boundary provided by the 
  // background mesh by one degree
  std::map<std::string,double>& backGroundMeshInfo =
    InputData->getBackGroundMeshInfo();

  int increaseElemOrder = 
    (int)backGroundMeshInfo["increaseElemOrderOnBoundary"];

  if(increaseElemOrder == 1) {

    // loop over all affected volume element on the boundary
    //    for(int i=0;i<vElemIdx.size();i++) {
      // exchange data between processors

    //      FEMElement& elem = nodesElements[i];

    //      int& elemOrder = elem.getElemOrder();
    //      ++elemOrder;

    //    }

    // loop over all surface elements
    //for(int i=0;i<surfaceNodesElems.size();i++) {

    //  FEMElement& elem = surfaceNodesElems[i];

    //  int& elemOrder = elem.getElemOrder();
    //  elemOrder = nodesElements[sElemTOvElemIdx[i]].getElemOrder();
     
    //  int& elemType = elem.getElemType();
    //  elemType = nodesElements[sElemTOvElemIdx[i]].getElemType();

    //}

    // loop over all line elements
    //for(int i=0;i<lineNodesElems.size();i++) {

    //  FEMElement& elem = lineNodesElems[i];

    // int& elemOrder = elem.getElemOrder();
    // elemOrder = nodesElements[lElemTOvElemIdx[i]].getElemOrder();
     
    //int& elemType = elem.getElemType();
    //elemType = nodesElements[lElemTOvElemIdx[i]].getElemType();

    //}

 
    // ------------------------------------------------------------------
    // ------------------------------------------------------------------
    // Compute the new nodes

    int numOfNewNodes,node,oldNodesPerVolElem;
    intVector data(2);
    map<string,double> oldParams;
    vector<map<int,int> > shapeSetOnCoordsID;

    dbMatrix3 shapeSetOnCoords;
   
    int shapeOrdsPos;
    FEMElement newElem(usedDOF);

    int currentNode = particles.size();
    int numOfOldPtcls = particles.size();

    vector<Particle> newParticles;

    // loop over all affected volume elements on the boundary
    for(int i=0;i<vElemIdx.size();i++) {
      int& elemID = vElemIdx[i];

      int& elemType = nodesElements[elemID].getElemType();
      int& oldElemOrder = nodesElements[elemID].getElemOrder();

      intVector& nodes = nodesElements[elemID].getNodes();
      oldNodesPerVolElem = nodes.size();

      int& newElemType = newElem.getElemType();
      int& newElemOrder = newElem.getElemOrder();
      newElemType = elemType;
      newElemOrder = oldElemOrder+1;

      ElementTemplate* newFEVolumeSet = 
	getVolumeElementTemplate(InputData,newElem,logFile);
      dbMatrix& newNodalCoords = newFEVolumeSet->nodalCoords;
	
      numOfNewNodes = newNodalCoords.size() - nodes.size();
      newParticles.resize(currentNode + numOfNewNodes,Particle(usedDOF));


      // check if the shape function ordinates are already calculated
     if(shapeSetOnCoordsID.size() <= elemType)
       shapeSetOnCoordsID.resize(elemType+1);

      map<int,int>::iterator p = 
	shapeSetOnCoordsID[elemType].find(oldElemOrder);

      // calculate the volume shape function ordinates at the new nodes
      if(p == shapeSetOnCoordsID[elemType].end()) {

	shapeOrdsPos = shapeSetOnCoords.size();
	shapeSetOnCoordsID[elemType][oldElemOrder] = shapeOrdsPos;
	shapeSetOnCoords.resize(shapeOrdsPos+1);
	shapeSetOnCoords[shapeOrdsPos] = 
	  dbMatrix(numOfNewNodes,dbVector(oldNodesPerVolElem));

	ElementTemplate* oldFEVolumeSet = 
	  nodesElements[elemID].getVolumeElementTemplate();


	// loop over the newly to be created nodes
	for(int i=0,m=nodes.size();i<numOfNewNodes;i++,m++)
	  
	  for(int j=0;j<oldNodesPerVolElem;j++)
	    
	    shapeSetOnCoords[shapeOrdsPos][i][j] = 
	      oldFEVolumeSet->N(j,newNodalCoords[m]);
	
      }
      else 

	shapeOrdsPos = p->second;
      

      // loop over the newly to be created nodes and calculate
      // their coordinates
      for(int i=0;i<numOfNewNodes;i++) {
	dbVector& newCoords = newParticles[currentNode+i].getCoords();

	for(int j=0;j<usedDims;j++) {

	  for(int k=0;k<oldNodesPerVolElem;k++) {
	    node = nodes[k]-1;
	
	    newCoords[j] += 
	      shapeSetOnCoords[shapeOrdsPos][i][k]*
	      particles[node].getCoord(j);	

	  }

	}

      }

      // extend the element node connectivity list
      for(int i=nodes.size();i<newNodalCoords.size();i++) {

	++currentNode;
	pushBackVector(nodes,currentNode);

      }

      oldElemOrder = newElemOrder;

    }

    // ------------------------------------------------------------------
    // exchange newly created particle coordinates between the processors
    int localNumOfNewPtcls = newParticles.size();
    int numOfNewPtcls;

    MPI_Allreduce(&localNumOfNewPtcls,&numOfNewPtcls,1,
		  MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    // determine a processor root list for all new particles
    int ptcleStartIdx,ptcleEndIdx;
    int defaultPortion = (int)ceil((double)numOfNewPtcls/size);

    if(defaultPortion*rank < numOfNewPtcls && 
       defaultPortion*(rank+1) <= numOfNewPtcls) {
      ptcleStartIdx = defaultPortion*rank;
      ptcleEndIdx = defaultPortion*(rank+1);
    }
    else if(defaultPortion*rank <= numOfNewPtcls && 
	    defaultPortion*(rank+1) >= numOfNewPtcls) {
      ptcleStartIdx = defaultPortion*rank;
      ptcleEndIdx = numOfNewPtcls;
    }
    else {
      ptcleStartIdx = ptcleEndIdx = numOfNewPtcls;
    }
  
    int ptclePortion = ptcleEndIdx-ptcleStartIdx;  

    int root = 0;
    intVector ptcleRootList(numOfNewPtcls);

    for(int i=0;i<numOfNewPtcls;i++) {

      // next processors turn
      if(i >= (root+1)*defaultPortion) {
	++root;
      }

      // set the processor ID
      ptcleRootList.push_back(root);

    }

    particles.resize(particles.size()+numOfNewPtcls,Particle(usedDOF));
    
    // loop over all new particles and exchange the particle coordinates
    for(int i=0;i<numOfNewPtcls;i++) {

      dbVector& coords = particles[numOfOldPtcls+i].getCoords();
    
      if(ptcleRootList[i] == rank) {
	dbVector& newCoords = newParticles[i].getCoords();
	coords = newCoords;
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&coords[0],3,MPI_DOUBLE,ptcleRootList[i],MPI_COMM_WORLD);
    }

 
#ifdef _FEdebugMode_
    logFile<<"**************** particle coords *******************"<<endl;
    for(int i=0;i<particles.size();i++) {
      dbVector& coords = particles[i].getCoords();
      logFile<<"PARTICLE "<<i<<": ";
      for(int j=0;j<coords.size();j++)
	logFile<<coords[j]<<" ";
      logFile<<endl;
    }
    logFile<<"****************** volume elements *****************"<<endl;
    logFile<<"----------------------------------------------------"<<endl;
    logFile<<"***************** volume shapes ********************"<<endl;
    for(int m=0;m<shapeSetOnCoords.size();m++) {
      logFile<<"volume shape set "<<m<<":"<<endl;
      for(int i=0;i<shapeSetOnCoords[m].size();i++) {
	logFile<<"node "<<i<<": ";
	for(int j=0;j<shapeSetOnCoords[m][i].size();j++)
	  logFile<<shapeSetOnCoords[m][i][j]<<" ";
	logFile<<endl;
      }
      logFile<<"--------------------------------------------------"<<endl;
    }
    for(int i=0;i<nodesElements.size();i++) {
      intVector& nodes = nodesElements[i].getNodes();
      logFile<<"Volume-element "<<nodesElements[i].getGlobalID()<<": "
	     <<"type = "<<nodesElements[i].getElemType()<<" "
	     <<"order = "<<nodesElements[i].getElemOrder()<<" "
	     <<"nodes: "<<endl;
      for(int j=0;j<nodes.size();j++) {
	logFile<<j<<".) "<<nodes[j]<<" coords: ";
	for(int k=0;k<usedDims;k++)
	  logFile<<particles[nodes[j]-1].getCoord(k)<<" ";
	logFile<<endl;
      }
    }
#endif

    // ------------------------------------------------------------------
    // update surface node element connectivity list
    int numOfNewVolElemNodes,pos;
    dbMatrix volCoords;
    dbVector newCoords(usedDims);

    shapeSetOnCoordsID.resize(0);
    shapeSetOnCoords.resize(0);

#ifdef _FEdebugMode_
    logFile<<"***************************************************"<<endl;
    logFile<<"************** surface elements *******************"<<endl;
#endif

    // loop over all affected surface elements
    for(int i=0;i<sElemTOvElemIdx.size();i++) {
      int& sElemID = sElemTOvElemIdx[i][0];
      int& vElemID = sElemTOvElemIdx[i][1];

      int& elemType = surfaceNodesElems[sElemID].getElemType();
      int& oldElemOrder = surfaceNodesElems[sElemID].getElemOrder();

      int& newElemType = newElem.getElemType();
      int& newElemOrder = newElem.getElemOrder();
      newElemType = elemType;
      newElemOrder = oldElemOrder+1;

      // determine parameter 'nodesPerVolumeElement' for the old
      // element order
      data[0] = elemType;
      data[1] = oldElemOrder;
      getFEMMeshData(data,oldParams);
     
      intVector& vNodes = nodesElements[vElemID].getNodes();
      oldNodesPerVolElem = (int)oldParams["nodesPerVolumeElement"];
      numOfNewVolElemNodes = vNodes.size() - oldNodesPerVolElem;

      if(numOfNewVolElemNodes > volCoords.size())
	volCoords.resize(numOfNewVolElemNodes);

      ElementTemplate* newFESurfaceSet = 
	getSurfaceElementTemplate(InputData,newElem,logFile);
      dbMatrix& newNodalCoords = newFESurfaceSet->nodalCoords;

      intVector& nodes = surfaceNodesElems[sElemID].getNodes();
      numOfNewNodes = newNodalCoords.size() - nodes.size();

#ifdef _FEdebugMode_
      logFile<<"SURFACE ELEMENT "<<i<<": nodes: ";
      for(int j=0;j<nodes.size();j++)
	logFile<<nodes[j]<<" ";
      logFile<<endl;
      logFile<<"elemType = "<<elemType<<" oldElemOrder = "<<oldElemOrder
	     <<" newElemOrder = "<<newElemOrder<<" numOfNewNodes = "
	     <<numOfNewNodes<<endl;
      logFile<<"-> vElemID = "<<vElemID<<" ";
      logFile<<" oldNodesPerVolElem = "<<oldNodesPerVolElem
	     <<" numOfNewVolElemNodes = "<<numOfNewVolElemNodes<<" "<<endl;
      for(int j=0;j<vNodes.size();j++) {
	logFile<<"node "<<vNodes[j]<<" coords: ";
	for(int k=0;k<usedDims;k++)
	  logFile<<particles[vNodes[j]-1].getCoord(k)<<" ";
	logFile<<endl;
      }
#endif 

      // check if the shape function ordinates are already calculated
      if(shapeSetOnCoordsID.size() <= elemType)
	shapeSetOnCoordsID.resize(elemType+1);

      map<int,int>::iterator p = 
	shapeSetOnCoordsID[elemType].find(oldElemOrder);

      // calculate the volume shape function ordinates at the new nodes
      if(p == shapeSetOnCoordsID[elemType].end()) {

	shapeOrdsPos = shapeSetOnCoords.size();
	shapeSetOnCoordsID[elemType][oldElemOrder] = shapeOrdsPos;
	shapeSetOnCoords.resize(shapeOrdsPos+1);
	shapeSetOnCoords[shapeOrdsPos] = 
	  dbMatrix(numOfNewNodes,dbVector(oldNodesPerVolElem));

	ElementTemplate* oldFEVolumeSet = 
	  surfaceNodesElems[sElemID].getVolumeElementTemplate();
	
	// loop over the newly to be created surface nodes
	for(int i=0,m=nodes.size();i<numOfNewNodes;i++,m++)
	  
	  // loop over all volume nodes
	  for(int j=0;j<oldNodesPerVolElem;j++)
	    
	    shapeSetOnCoords[shapeOrdsPos][i][j] = 
	      oldFEVolumeSet->N(j,newNodalCoords[m]);
	
      }
      else 

	shapeOrdsPos = p->second;
      

      // loop over the newly to be created surface nodes and calculate
      // their coordinates
      for(int i=0;i<numOfNewNodes;i++) {
	clearArray(newCoords);

	for(int j=0;j<usedDims;j++) {

	  for(int k=0;k<oldNodesPerVolElem;k++) {
	    node = vNodes[k]-1;
	
	    newCoords[j] += 
	      shapeSetOnCoords[shapeOrdsPos][i][k]*
	      particles[node].getCoord(j);	

	  }

	}

#ifdef _FEdebugMode_
	logFile<<"***********************************************"<<endl;
	logFile<<"** new node "<<i<<": coords: ";
	for(int j=0;j<newCoords.size();j++)
	  logFile<<newCoords[j]<<" ";
	logFile<<endl;
#endif

	// loop over the new volume element nodes and assemble their
	// coordinates
	for(int j=oldNodesPerVolElem,k=0;j<vNodes.size();j++,k++) {
	  node = vNodes[j]-1;
	  volCoords[k] = particles[node].getCoords();
	}

	// find the volume element node which corresponds to the surface 
	// element node
	mode = "strict";
	pos = findDoubleMatPos(newCoords,0,numOfNewVolElemNodes,mode,
			       1.0e-04,volCoords,logFile);

	if(pos != -1)
	  pushBackVector(nodes,vNodes[pos]);

	else {
	  logFile<<"In FEMGeometry::rearrangeMeshConfiguration "
		 <<"volume element node was not found!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

      }

      oldElemOrder = newElemOrder; 

    }

    // exchange data between processors
    //intVector data(2);
    map<string,double> params;

    for(int i=0;i<surfaceRootList.size();i++) { 
      int& elemType = surfaceNodesElems[i].getElemType();
      int& elemOrder = surfaceNodesElems[i].getElemOrder();

      data[0] = elemType;
      data[1] = elemOrder;
      getFEMMeshData(data,params);

      MPI_Bcast(&elemOrder,1,MPI_INT,surfaceRootList[i],MPI_COMM_WORLD);

      intVector& nodes = surfaceNodesElems[i].getNodes();
      MPI_Bcast(&nodes[0],(int)params["nodesPerSurfaceElement"],MPI_INT,
		surfaceRootList[i],MPI_COMM_WORLD);
    } 

#ifdef _FEdebugMode_
    logFile<<"**************** surface elements ******************"<<endl;
    logFile<<"----------------------------------------------------"<<endl;
    logFile<<"***************** volume shapes ********************"<<endl;
    for(int m=0;m<shapeSetOnCoords.size();m++) {
      logFile<<"element type "<<m<<":"<<endl;
      for(int i=0;i<shapeSetOnCoords[m].size();i++) {
	logFile<<"node "<<i<<": ";
	for(int j=0;j<shapeSetOnCoords[m][i].size();j++)
	  logFile<<shapeSetOnCoords[m][i][j]<<" ";
	logFile<<endl;
      }
      logFile<<"--------------------------------------------------"<<endl;
    }
    for(int i=0;i<surfaceNodesElems.size();i++) {
      intVector& nodes = surfaceNodesElems[i].getNodes();
      logFile<<"Surface-element "<<i<<"("<<surfaceNodesElems[i].getGlobalID()<<"): "
	     <<"type = "<<surfaceNodesElems[i].getElemType()<<" "
	     <<"order = "<<surfaceNodesElems[i].getElemOrder()<<" "
	     <<"nodes: "<<endl;
      for(int j=0;j<nodes.size();j++)
	logFile<<nodes[j]<<" ";
      logFile<<endl; 
    }
#endif

    // ------------------------------------------------------------------
    // ------------------------------------------------------------------
    // Ensure that all elements have the necessary approximation and
    // integration tools assigned.
    setApproxAndIntTools(InputData,modelData,logFile);

  }

  else if(increaseElemOrder > 1) {
    logFile<<"In FEMGeometry::rearrangeMeshConfiguration "
	   <<"boundary element order can not be increased\n"
	   <<"by "<<increaseElemOrder<<" degrees!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

/*
  * loop over all line and surface elements and determine those volume elements they are part of
    -> intVector volumeElementIdx
  * determine for each volume element its 6 faces
    -> calculate for each face ksi = -1/1, eta = -1/1 and zeta = -1/1 the four set of nodal coordinates
       -> dbMatrix elementCoords(6,dbVector(4))
    -> compare for each face the just calculated nodal coordinates with the ones 
       in particles[elemNodeIdx].getCoords() in order to determine the global nodal indices
       -> matrix which contains for each face the corresponding global nodal indices 
          intMatrix elementFaces(6,intVector(4))
  * determine for each volume element its neighbours and by which face they are connected to
    -> loop over all elements faces which are still not set (=> blMatrix setFaces(surfVolElems,blVector(6)))
       and then loop over all surface-volume-elements and for each surface-volume element loop over 
       all still not set face (=> setFaces) and compare the nodal coordinates 
       -> if there is a match store the current element its neighbour and the face ID and do the same
          for the neighbour in intMatrix3 volumeElementNeighbours(surfVolElems,intMatrix(*,intVector(2))),
          which contains the element indices and the face IDs
  * determine array setNodes(numOfNodes) 
    -> assign(numOfNodes,true)
    -> loop over all all surface-volume elements and set setNodes[surfVolElemIdx] = false
  * loop over all surface-volume elements and add extra nodes
    -> add internal nodes
    -> check blVector setNodes and add face node if necessary and store its coordinates in 
       allNewNodeCoords(addedFaceNodes,dbVector(4)) which contains the coords and the global nodal indices
  * add nodes of all surface boundary elements
    -> loop over all surface boundaryelements and determine the additional nodal coordinates
       dbVector newNodalCoords(5) and find the corresponding nodal index in allNewNodeCoords
       and store these indices
  * add nodes of all line boundary elements 
    -> in the same as done for the surface boundary elements
*/
     



}

/************************************************************************/
/************************************************************************/
// Ensure that all elements have the necessary approximation and
// integration tools assigned.
void FEMGeometry::setApproxAndIntTools(InputFileData* InputData,
				       std::map<std::string,double>& modelData,
				       std::ofstream& logFile) {
  
  using namespace std;

  std::map<std::string,double>& backGroundMeshInfo =
    InputData->getBackGroundMeshInfo();

#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"************* set approximation tools ***************"<<endl;
  logFile<<"*****************************************************"<<endl;
#endif

  // Assign to all elements the necessary FEM approximation tools.

  // loop over all volume elements
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];
    elem.setVolumeElementTemplate(getVolumeElementTemplate(InputData,elem,logFile));
  }

  // loop over all surface elements
  for(int i=0;i<surfaceNodesElems.size();i++) {
    FEMElement& elem = surfaceNodesElems[i];
    elem.setVolumeElementTemplate(getVolumeElementTemplate(InputData,elem,logFile));
    elem.setSurfaceElementTemplate(getSurfaceElementTemplate(InputData,elem,logFile));
  }

  // loop over all line elements
  for(int i=0;i<lineNodesElems.size();i++) {
    FEMElement& elem = lineNodesElems[i];
    elem.setVolumeElementTemplate(getVolumeElementTemplate(InputData,elem,logFile));
    elem.setLineElementTemplate(getLineElementTemplate(InputData,elem,logFile));
  }

  // Assign to all elements the FEM the necessary Gauss quadrature tools.
  int integrationMethod = (int)modelData["integrationMethod"];
  
  if(integrationMethod == 1) {
    
    for(int i=0;i<nodesElements.size();i++) {
      FEMElement& elem = nodesElements[i];
      intVector& intPts = elem.getVolumeIntegrationPts();
      
      intPts.resize((int)backGroundMeshInfo["gaussPointsPerVolumeElement"]);
      elem.setVolumeGaussSet(getVolumeGaussSet(InputData,elem,logFile));
    }
    
    for(int i=0;i<surfaceNodesElems.size();i++) {
      FEMElement& elem = surfaceNodesElems[i];
      intVector& intPts = elem.getSurfaceIntegrationPts();
      
      intPts.resize((int)backGroundMeshInfo["gaussPointsPerSurfaceElement"]);
      elem.setSurfaceGaussSet(getSurfaceGaussSet(InputData,elem,logFile));
    }
    
    for(int i=0;i<lineNodesElems.size();i++) {
      FEMElement& elem = lineNodesElems[i];
      intVector& intPts = elem.getLineIntegrationPts();
      
      intPts.resize((int)backGroundMeshInfo["gaussPointsPerLineElement"]);
      elem.setLineGaussSet(getLineGaussSet(InputData,elem,logFile));
    }

    /********************************************************************/
    // change the order of the Gauss quadrature for selected elements
    //intMatrix& gaussPtsPerVolumeElem = InputData->getAllGaussPtsPerVolumeElem();
    //intMatrix& gaussPtsPerSurfaceElem = InputData->getAllGaussPtsPerSurfaceElem();
    //intMatrix& gaussPtsPerLineElem = InputData->getAllGaussPtsPerLineElem();

    //for(int i=0;i<gaussPtsPerVolumeElem.size();i++) {
    //  int& elemID = gaussPtsPerVolumeElem[i][0];
    //  intVector& intPts = nodesElements[elemID].getVolumeIntegrationPts();

    //      intPts.resize(gaussPtsPerVolumeElem[i][1]);
    //}

    //for(int i=0;i<gaussPtsPerSurfaceElem.size();i++) {
    //  int& elemID = gaussPtsPerSurfaceElem[i][0];
    //  intVector& intPts = surfaceNodesElems[elemID].getSurfaceIntegrationPts();

    //      intPts.resize(gaussPtsPerSurfaceElem[i][1]);
    //}

    //for(int i=0;i<gaussPtsPerLineElem.size();i++) {
    //  int& elemID = gaussPtsPerLineElem[i][0];
    //  intVector& intPts = lineNodesElems[elemID].getLineIntegrationPts();

    //  intPts.resize(gaussPtsPerLineElem[i][1]);
    //}

  }

  /*********************************************************************/
  // calculate for each element the shape function ordinates at all
  // element's Gauss integration points

  // loop over all volume elements
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];
    elem.setVolumeShapeFuncOrds(getVolumeFEShapeFunctions(InputData,elem,logFile));
  }

  // loop over all surface elements
  for(int i=0;i<surfaceNodesElems.size();i++) {
    FEMElement& elem = surfaceNodesElems[i];
    elem.setSurfaceShapeFuncOrds(getSurfaceFEShapeFunctions(InputData,elem,logFile));
  }

  // loop over all line elements
  for(int i=0;i<lineNodesElems.size();i++) {
    FEMElement& elem = lineNodesElems[i];
    elem.setLineShapeFuncOrds(getLineFEShapeFunctions(InputData,elem,logFile));
  }

  // calculate for each element the shape function derivative ordinates 
  // at all element's Gauss integration points

  if((int)InputData->getValue("shapefunctionType") == 7) {
    
    // loop over all volume elements
    for(int i=0;i<nodesElements.size();i++) {
      FEMElement& elem = nodesElements[i];
      elem.setVolumeShapeFuncDerivOrds(getVolumeFEShapeFuncDerivs(InputData,elem,logFile));
    }

  }
  
}

/**********************************************************************/
/**********************************************************************/
// Determine for all local volume Gauss points their coordinates.
void FEMGeometry::setVolumeGaussPoints(InputFileData* InputData,
				       std::map<std::string,double>& modelData,
				       std::ofstream& logFile,
				       PetscViewer& viewerMPI) {
  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  // Determine the number of global and local Gauss points and their 
  // global IDs.

  // determine the local number of Gauss points and assign the local Gauss 
  // point IDs to the corresponding elements.
  int localGaussPtsNum = 0;
  elemGaussIdx.resize(nodesElements.size());

  for(int i=0;i<nodesElements.size();i++) {
    intVector& intPts = nodesElements[i].getVolumeIntegrationPts();
    elemGaussIdx[i].resize(intPts.size());

    for(int j=0;j<intPts.size();j++) {
      intPts[j] = localGaussPtsNum;
      elemGaussIdx[i][j] = localGaussPtsNum;
      localGaussPtsNum++;
    }
 
  }

  gaussPoints.resize(localGaussPtsNum);

  // determine the number of global Gauss points
  int globalStartIdx,rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  intVector recvCounts(size);

  MPI_Allgather(&localGaussPtsNum,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		MPI_COMM_WORLD);

  globalGaussPtsNum = 0;
  
  for(int i=0;i<recvCounts.size();i++) {

    if(i == rank)
      globalStartIdx = globalGaussPtsNum;
    
    globalGaussPtsNum += recvCounts[i];
  }

  int globalEndIdx = globalStartIdx + recvCounts[rank];

  // set global Gauss points indices
  for(int i=globalStartIdx,j=0;i<globalEndIdx;i++,j++)
    gaussPoints[j].setGlobalID(i);

#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"***************** VOLUME GAUSS POINTS ***************"<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"globalGaussPtsNum = "<<globalGaussPtsNum<<endl;
  logFile<<"localGaussPtsNum = "<<localGaussPtsNum<<endl;
  logFile<<"globalStartIdx = "<<globalStartIdx<<" globalEndIdx = "
	 <<globalEndIdx<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"******************* elemGaussIdx ********************"<<endl;
  for(int i=0;i<elemGaussIdx.size();i++) {
    logFile<<"ELEM "<<i<<": ";
    for(int j=0;j<elemGaussIdx[i].size();j++)
      logFile<<elemGaussIdx[i][j]<<" ";
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // Determine the coordinates of all local Gauss points.  
  int node;
  
  int idx = 0;

  // Loop over all local elements.
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];

    int& elemMatID = nodesElements[i].getMaterialID();
    intVector& intPts = nodesElements[i].getVolumeIntegrationPts();
    intVector& nodes = elem.getNodes();
    int nodesPerElem = nodes.size();

    // Determine the set of shape functions.
    dbMatrix& sFuncs = elem.getVolumeShapeFuncOrds();

    // Loop over all gauss coordinates of a single element to determine
    // its gauss points.
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      dbVector& gCoords = gaussPoints[idx].getCoords();      
      gCoords.resize(usedDims);
      
      for(int k=0;k<nodes.size();k++) {
	node = nodes[k]-1;
	
	gCoords[0] += 
	  sFuncs[j][k]*particles[node].getCoord(0);	
	gCoords[1] += 
	  sFuncs[j][k]*particles[node].getCoord(1);
	gCoords[2] += 
	  sFuncs[j][k]*particles[node].getCoord(2);
      }

      // set material ID
      int& matID = gaussPoints[idx].getMaterialID();
      matID = elemMatID;

      // set location of the current Gauss point
      intVector& elemInfo = gaussPoints[idx].getElementInfo();
      elemInfo.resize(2);
      elemInfo[0] = i; // element ID
      elemInfo[1] = j; // Gauss point ID within the element

      //idx++;
    } 
    
  }

  
  
#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"***************** VOLUME GAUSS POINTS ***************"<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"******************* template element ****************"<<endl;
  std::map<std::string,double>& backGroundMeshInfo = 
    InputData->getBackGroundMeshInfo();
  dbVector gCoords(usedDOF);
  int ndsPerElem = (int)backGroundMeshInfo["nodesPerVolumeElement"];
  int gPtsPerElem = (int)backGroundMeshInfo["gaussPointsPerVolumeElement"];
  FEMElement elem(usedDOF);
  intVector& intPts = elem.getVolumeIntegrationPts();
  intVector& nds = elem.getNodes();
  int& elemType = elem.getElemType();
  int& elemOrder = elem.getElemOrder();
  intPts.resize(gPtsPerElem);
  nds.resize(ndsPerElem);
  elemType = (int)backGroundMeshInfo["elemType"];
  elemOrder = (int)backGroundMeshInfo["elemOrder"];
  ElementTemplate* FESet = 
    getVolumeElementTemplate(InputData,elem,logFile);
  GaussPointSet* GaussSet = 
    getVolumeGaussSet(InputData,elem,logFile);
  dbMatrix& sFuncs = getVolumeFEShapeFunctions(InputData,elem,logFile);
  for(int i=0;i<gPtsPerElem;i++) {
    clearArray(gCoords);
    for(int j=0;j<ndsPerElem;j++)
      for(int k=0;k<usedDOF;k++)
	gCoords[k] += sFuncs[i][j]*FESet->nodalCoords[j][k];
    logFile<<"point "<<i<<": "<<GaussSet->coord[i][0]<<" "
	   <<GaussSet->coord[i][1]<<" "
	   <<GaussSet->coord[i][2]<<" ?= "<<endl;
    logFile<<gCoords[0]<<" "<<gCoords[1]<<" "<<gCoords[2]<<endl;
    logFile<<endl;
  }
  logFile<<"******************* gauss coords ********************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    int& matID = gaussPoints[i].getMaterialID();
    intVector& elemInfo = gaussPoints[i].getElementInfo();
    logFile<<i<<" "<<gaussPoints[i].getGlobalID()<<" "
	   <<gaussPoints[i].getCoord(0)<<" "
	   <<gaussPoints[i].getCoord(1)<<" "
	   <<gaussPoints[i].getCoord(2)<<" matID "<<matID
	   <<" element "<<elemInfo[0]<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"******************* elemGaussCoords *****************"<<endl;
  for(int i=0;i<elemGaussIdx.size();i++) {
    logFile<<"ELEM "<<i<<": "<<endl;
    for(int j=0;j<elemGaussIdx[i].size();j++) {
      GaussPoint& gPoint = gaussPoints[elemGaussIdx[i][j]];
      logFile<<"GPt "<<elemGaussIdx[i][j]<<" (globalID "
	     <<gPoint.getGlobalID()<<") "
	     <<gPoint.getCoord(0)<<" "
	     <<gPoint.getCoord(1)<<" "
	     <<gPoint.getCoord(2)<<endl;
    }
  }
#endif

    
  /*********************************************************************/
  // Store the element body force to the corresponding gauss points.
  int point;
  intVector indexSet(3);
  //bodyForceGaussPtsIdx;
  
  for(int i=0;i<bodyForceElemIdx.size();i++) {
    int& currentIdx = bodyForceElemIdx[i][0];
    FEMElement& elem = nodesElements[currentIdx];

    intVector& intPts = elem.getVolumeIntegrationPts();
    blVector& elemDOF = elem.getBodyForceDOF();
    dbVector& elemConds = elem.getBodyForce();
    
    for(int j=0;j<intPts.size();j++) {
      point = intPts[j];

      // Store the the body force.
      blVector& affectedDOF = gaussPoints[point].getBodyForceDOF();
      dbVector& conditions = gaussPoints[point].getBodyForce();
 
      affectedDOF = elemDOF;
      conditions = elemConds;
      
      indexSet[0] = point; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // body force ID

      pushBackVector(bodyForceGaussPtsIdx,indexSet); 
    }
    
  }
   
  /*********************************************************************/
  // Store the element body moment to the corresponding gauss points.
  for(int i=0;i<bodyMomentElemIdx.size();i++) {
    int& currentIdx = bodyMomentElemIdx[i][0];
    FEMElement& elem = nodesElements[currentIdx];

    intVector& intPts = elem.getVolumeIntegrationPts();
    blVector& elemDOF = elem.getBodyMomentDOF();
    dbVector& elemConds = elem.getBodyMoment();
    
    for(int j=0;j<intPts.size();j++) {
      point = intPts[j];
      
      // Store the the body moment.
      blVector& affectedDOF = gaussPoints[point].getBodyMomentDOF();
      dbVector& conditions = gaussPoints[point].getBodyMoment();
 
      affectedDOF = elemDOF;
      conditions = elemConds;
      
      indexSet[0] = point; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // body moment ID

      pushBackVector(bodyMomentGaussPtsIdx,indexSet); 
    }
    
  }

  for(int i=0;i<bodyElectricChargeElemIdx.size();i++) {
    int& currentIdx = bodyElectricChargeElemIdx[i][0];
    FEMElement& elem = nodesElements[currentIdx];

    intVector& intPts = elem.getVolumeIntegrationPts();
    blVector& elemDOF = elem.getBodyElectricChargeDOF();
    dbVector& elemConds = elem.getBodyElectricCharge();
    
    for(int j=0;j<intPts.size();j++) {
      point = intPts[j];

      // Store the the body force.
      blVector& affectedDOF = gaussPoints[point].getBodyElectricChargeDOF();
      dbVector& conditions = gaussPoints[point].getBodyElectricCharge();
 
      affectedDOF = elemDOF;
      conditions = elemConds;
      
      indexSet[0] = point; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // body force ID

      pushBackVector(bodyElectricChargeGaussPtsIdx,indexSet); 
    }
    
  }

  // Calculate for all local volume Gauss points their weights.
  setVolumeGaussWeights(InputData,modelData,logFile);
  
#ifdef _FEdebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"************** gauss weights ************************"<<endl;
  for(int i=0;i<localGaussPtsNum;i++)
    logFile<<i<<" "<<"GAUSS POINT "<<gaussPoints[i].getGlobalID()<<": "
	   <<gaussPoints[i].getWeight()<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"**************** body forces on gauss ***************"<<endl;
  for(int i=0;i<bodyForceGaussPtsIdx.size();i++) {
    point = bodyForceGaussPtsIdx[i][0];
    int& weightID = bodyForceGaussPtsIdx[i][1];
    int& loadID = bodyForceGaussPtsIdx[i][2];
    dbVector& loads = gaussPoints[point].getBodyForce(loadID);
    blVector& affectedDOF = gaussPoints[point].getBodyForceDOF(loadID);
    double& weight = gaussPoints[point].getIntWeight(weightID);
    logFile<<"GAUSS POINT "<<gaussPoints[i].getGlobalID()<<": ";
    logFile<<gaussPoints[i].getCoord(0)<<" "
	   <<gaussPoints[i].getCoord(1)<<" "
	   <<gaussPoints[i].getCoord(2)<<" ";
    logFile<<"body force: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"**************** body moments on gauss ***************"<<endl;
  for(int i=0;i<bodyMomentGaussPtsIdx.size();i++) {
    point = bodyMomentGaussPtsIdx[i][0];
    int& weightID = bodyMomentGaussPtsIdx[i][1];
    int& loadID = bodyMomentGaussPtsIdx[i][2];
    dbVector& loads = gaussPoints[point].getBodyMoment(loadID);
    blVector& affectedDOF = gaussPoints[point].getBodyMomentDOF(loadID);
    double& weight = gaussPoints[point].getIntWeight(weightID);
    logFile<<"GAUSS POINT "<<gaussPoints[i].getGlobalID()<<": ";
    logFile<<gaussPoints[i].getCoord(0)<<" "
	   <<gaussPoints[i].getCoord(1)<<" "
	   <<gaussPoints[i].getCoord(2)<<" ";
    logFile<<"body moment: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<endl;
  } 
  logFile<<"*****************************************************"<<endl;
  logFile<<"*************** body electric charge on gauss *******"<<endl;
  for(int i=0;i<bodyElectricChargeGaussPtsIdx.size();i++) {
    point = bodyElectricChargeGaussPtsIdx[i][0];
    int& weightID = bodyElectricChargeGaussPtsIdx[i][1];
    int& loadID = bodyElectricChargeGaussPtsIdx[i][2];
    dbVector& loads = gaussPoints[point].getBodyElectricCharge(loadID);
    blVector& affectedDOF = gaussPoints[point].getBodyElectricChargeDOF(loadID);
    double& weight = gaussPoints[point].getIntWeight(weightID);
    logFile<<"GAUSS POINT "<<gaussPoints[i].getGlobalID()<<": ";
    logFile<<gaussPoints[i].getCoord(0)<<" "
	   <<gaussPoints[i].getCoord(1)<<" "
	   <<gaussPoints[i].getCoord(2)<<" ";
    logFile<<"body force: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate boundary Gauss points.
void FEMGeometry::setBoundGaussPoints(InputFileData* InputData,
				      std::map<std::string,double>& modelData,
				      std::ofstream& logFile,
				      PetscViewer& viewerMPI) {
  using namespace std;

  int boundaryEnforcement = (int)modelData["boundaryEnforcementMethod"];
  bool lineElemConversion = 
    (bool)InputData->getValue("lineToSurfaceElementConversion");

  // Compute all surface Gauss points and store their surface load.
  setSurfaceGaussPoints(InputData,modelData,logFile,viewerMPI);

  // Compute all line Gauss points and store their line load.
  if(boundaryEnforcement == 2 && lineElemConversion)

    setLineSurfaceGPoints(InputData,modelData,logFile,viewerMPI);

  else

    setLineGaussPoints(InputData,modelData,logFile,viewerMPI);

  // Store the point integration points coordinate's, any deformation 
  // boundary condition and any point loads if applied.
  setBoundPointIntegration(InputData,modelData,logFile,viewerMPI);

}

/**********************************************************************/
/**********************************************************************/
// Calculate the surface gauss points coordinate's and store the 
// surface load if applied.
void FEMGeometry::setSurfaceGaussPoints(InputFileData* InputData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile,
					PetscViewer& viewerMPI) {
  
  using namespace std;
  
  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];
  vector<FEMElement>& bElems = surfaceNodesElems;
  
  /*********************************************************************/
  // Determine a local portion of boundary Gauss points and the 
  // corresponding rootList.
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int elemStartIdx,elemEndIdx;
  int surfaceElemNum = surfaceNodesElems.size();
  int defaultPortion = (int)ceil((double)surfaceElemNum/size);

  if(defaultPortion*rank < surfaceElemNum && 
     defaultPortion*(rank+1) <= surfaceElemNum) {
    elemStartIdx = defaultPortion*rank;
    elemEndIdx = defaultPortion*(rank+1);
  }
  else if(defaultPortion*rank <= surfaceElemNum && 
	  defaultPortion*(rank+1) >= surfaceElemNum) {
    elemStartIdx = defaultPortion*rank;
    elemEndIdx = surfaceElemNum;
  }
  else {
    elemStartIdx = elemEndIdx = surfaceElemNum;
  }
  
  int elemPortion = elemEndIdx-elemStartIdx;  
  
  int root = 0;
  int gaussStartIdx = 0;
  int localGaussPtPortion = 0;

  // determine a Gauss point rootList
  intVector rootList;
  
  for(int i=0;i<surfaceElemNum;i++) {
    intVector& intPts = 
      bElems[i].getSurfaceIntegrationPts();

    for(int j=0;j<intPts.size();j++) {
    
      // next processors turn
      if(i >= (root+1)*defaultPortion) {

	// next processor Gauss point start index
	if(root+1 == rank)
	  gaussStartIdx = rootList.size();

	++root;
      }
    
      if(root == rank)
	localGaussPtPortion++;

      intPts[j] = rootList.size();

      // set the processor ID
      rootList.push_back(root);

    }

  }

  int surfaceGaussPtsNum = rootList.size();
  int gaussEndIdx = gaussStartIdx + localGaussPtPortion;

  globalBGaussPtsNum = surfaceGaussPtsNum;
  boundGaussPoints.resize(globalBGaussPtsNum);
  
#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"************** SURFACE GAUSS POINTS *****************"<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"elemPortion = "<<elemPortion
	 <<" elemStartIdx = "<<elemStartIdx
	 <<" elemEndIdx = "<<elemEndIdx
	 <<" surfaceGaussPtsNum "<<surfaceGaussPtsNum
	 <<" localGaussPtPortion = "<<localGaussPtPortion<<endl;
  logFile<<"******** element Gauss point connectivity ************"<<endl;
  for(int i=0;i<surfaceElemNum;i++) {
    intVector& intPts =
      bElems[i].getSurfaceIntegrationPts();
    logFile<<"Element "<<bElems[i].getGlobalID()<<": ";
    for(int j=0;j<intPts.size();j++) 
      logFile<<intPts[j]<<" ";
    logFile<<endl;
  }
  logFile<<"********************* rootList **********************"<<endl;
  for(int i=0;i<rootList.size();i++)
    logFile<<i<<".) "<<rootList[i]<<endl;
  logFile<<"gaussStartIdx = "<<gaussStartIdx<<endl;
  logFile<<"gaussEndIdx = "<<gaussEndIdx<<endl;
#endif
  
  /*********************************************************************/
  // Determine the coordinates of a local portion of boundary Gauss
  // points.
  int node;
  int idx;
  
  // Loop over the local portion of boundary elements.
  for(int i=elemStartIdx;i<elemEndIdx;i++) {
    intVector& intPts = 
      bElems[i].getSurfaceIntegrationPts();
    intVector& nodes = bElems[i].getNodes();
    int nodesPerElem = nodes.size();
    
    // Determine the set of FEM shape functions.
    dbMatrix& sFuncs = bElems[i].getSurfaceShapeFuncOrds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      dbVector& gCoords = boundGaussPoints[idx].getCoords();      
      gCoords.resize(usedDims);
      
      for(int k=0;k<nodes.size();k++) {
	node = nodes[k]-1;
	
	gCoords[0] += 
	  sFuncs[j][k]*particles[node].getCoord(0);	
	gCoords[1] += 
	  sFuncs[j][k]*particles[node].getCoord(1);
	gCoords[2] += 
	  sFuncs[j][k]*particles[node].getCoord(2);
      }

      // set location of the current Gauss point
      intVector& elemInfo = boundGaussPoints[idx].getElementInfo();
      elemInfo.resize(3);
      elemInfo[0] = i; // element ID
      elemInfo[1] = j; // Gauss point ID within the element
      elemInfo[2] = bElems[i].getMotherElementID(); // corresponding global volume element ID

      //idx++;
    } 

  }

  // --------------------------------------------------------------------
  // Store material ID.
  for(int i=0;i<surfaceNodesElems.size();i++) {

    int& elemMatID = bElems[i].getMaterialID();
    intVector& intPts = 
      bElems[i].getSurfaceIntegrationPts();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      int& matID =  boundGaussPoints[idx].getMaterialID();
      matID = elemMatID;
    }

  }

#ifdef _FEdebugMode_
  logFile<<"******************* template element ****************"<<endl;
  std::map<std::string,double>& backGroundMeshInfo = 
    InputData->getBackGroundMeshInfo();
  dbVector gCoords(2);
  int ndsPerElem = (int)backGroundMeshInfo["nodesPerSurfaceElement"];
  int gPtsPerElem = (int)backGroundMeshInfo["gaussPointsPerSurfaceElement"];
  FEMElement elem(usedDOF);
  intVector& intPts = elem.getSurfaceIntegrationPts();
  intVector& nds = elem.getNodes();
  int& elemType = elem.getElemType();
  int& elemOrder = elem.getElemOrder();
  intPts.resize(gPtsPerElem);
  nds.resize(ndsPerElem);
  elemType = (int)backGroundMeshInfo["elemType"];
  elemOrder = (int)backGroundMeshInfo["elemOrder"];
  elem.setSurfaceElementTemplate(getSurfaceElementTemplate(InputData,elem,logFile));
  elem.setSurfaceGaussSet(getSurfaceGaussSet(InputData,elem,logFile));
  ElementTemplate* FESet = elem.getSurfaceElementTemplate();
  GaussPointSet* GaussSet = elem.getSurfaceGaussSet();
  dbMatrix& sFuncs = getSurfaceFEShapeFunctions(InputData,elem,logFile);
  for(int i=0;i<gPtsPerElem;i++) {
    clearArray(gCoords);
    for(int j=0;j<ndsPerElem;j++)
      for(int k=0;k<2;k++)
	gCoords[k] += sFuncs[i][j]*FESet->nodalCoords[j][k];
    logFile<<"point "<<i<<": "<<GaussSet->coord[i][0]<<" "
	   <<GaussSet->coord[i][1]<<" ?= "<<gCoords[0]<<" "<<gCoords[1]<<endl;
  }
  logFile<<"**************** local gauss coords ******************"<<endl;
  logFile<<"Bound Gauss point portion: "
	 <<localGaussPtPortion<<" = "<<gaussEndIdx-gaussStartIdx<<endl;
  logFile<<"gaussStartIdx = "<<gaussStartIdx<<" gaussEndIdx = "
	 <<gaussEndIdx<<"; surfaceGaussPtsNum = "<<surfaceGaussPtsNum<<endl;
  for(int i=gaussStartIdx;i<gaussEndIdx;i++) {
    int& matID =  boundGaussPoints[i].getMaterialID();
    dbVector& gCoords = boundGaussPoints[i].getCoords();      
    intVector& elemInfo = boundGaussPoints[i].getElementInfo();
    logFile<<i<<".) Bound GAUSSPOINT: ";
    for(int j=0;j<gCoords.size();j++)
      logFile<<gCoords[j]<<" ";
    logFile<<"; s-elem "<<elemInfo[0]<<", local gPoint ID "<<elemInfo[1]
	   <<" volID "<<elemInfo[2]<<" matID "<<matID<<endl;
  }
#endif
  
  // Assemble all boundary gauss points coordinates from all processors,
  // store these.

  // Loop over the local calculated portion of Gauss points.
  for(int i=0;i<surfaceGaussPtsNum;i++) {
    boundGaussPoints[i].setGlobalID(i);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&boundGaussPoints[i].getCoord(0),usedDims,MPI_DOUBLE,
	      rootList[i],MPI_COMM_WORLD);

    intVector& elemInfo = boundGaussPoints[i].getElementInfo();

    if(elemInfo.size() < 3)
      elemInfo.resize(3);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&elemInfo[0],3,MPI_INT,rootList[i],MPI_COMM_WORLD);
  }
  
#ifdef _FEdebugMode_
  logFile<<"************* global bound gauss points **************"<<endl;
  for(int i=0;i<boundGaussPoints.size();i++) {
    dbVector& gCoords = boundGaussPoints[i].getCoords();      
    intVector& elemInfo = boundGaussPoints[i].getElementInfo();
    logFile<<i<<".) Bound GAUSSPOINT "
	   <<boundGaussPoints[i].getGlobalID()<<": ";
    for(int j=0;j<gCoords.size();j++)
      logFile<<gCoords[j]<<" ";
    logFile<<"; s-elem "<<elemInfo[0]<<", local gPoint ID "<<elemInfo[1]
	   <<" volID "<<elemInfo[2]<<endl;
  }
#endif


  /**********************************************************************/
  // Store the traction loads on all boundary gauss points and the
  // (local) indices of these boundary Gauss points.
  intVector indexSet(3);

  for(int i=0;i<tractionElemIdx.size();i++) {
    int& currentElem = tractionElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getTractionDOF();
    dbVector& elemConds = bElems[currentElem].getTraction();
    
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      dbVector& conditions = boundGaussPoints[idx].getTraction();
      blVector& affectedDOF = boundGaussPoints[idx].getTractionDOF();

      affectedDOF = elemDOF;
      conditions = elemConds;
	
      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // load ID
      pushBackVector(tractionBoundGaussPtsIdx,indexSet);
    }
    
  }

  /*********************************************************************/
  // Store the surface pressure loads on all boundary gauss points and 
  // the (local) indices of these boundary Gauss points.
  for(int i=0;i<surfacePressureElemIdx.size();i++) {
    int& currentElem = surfacePressureElemIdx[i][0];

    double& elemCond = bElems[currentElem].getSurfacePressure();
    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      
      double& condition = boundGaussPoints[idx].getSurfacePressure();

      condition = elemCond;
      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // load ID
      pushBackVector(surfacePressureBoundGaussPtsIdx,indexSet);
    }
    
  }
  
  /*********************************************************************/
  // Store the surface moment loads on all boundary gauss points and the
  // (local) indices of these boundary Gauss points.

  for(int i=0;i<surfaceMomentElemIdx.size();i++) {
    int& currentElem = surfaceMomentElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getSurfaceMomentDOF();
    dbVector& elemConds = bElems[currentElem].getSurfaceMoment();
    
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      dbVector& conditions = boundGaussPoints[idx].getSurfaceMoment();
      blVector& affectedDOF = boundGaussPoints[idx].getSurfaceMomentDOF();

      affectedDOF = elemDOF;
      conditions = elemConds;
	
      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // load ID
      pushBackVector(surfaceMomentBoundGaussPtsIdx,indexSet);
    }
    
  }

  /*********************************************************************/
  // Store the surface electric charge loads on all boundary gauss 
  // points and the (local) indices of these boundary Gauss points.

  for(int i=0;i<surfaceElectricChargeElemIdx.size();i++) {
    int& currentElem = surfaceElectricChargeElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getSurfaceElectricChargeDOF();
    dbVector& elemConds = bElems[currentElem].getSurfaceElectricCharge();
    
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      dbVector& conditions = boundGaussPoints[idx].getSurfaceElectricCharge();
      blVector& affectedDOF = boundGaussPoints[idx].getSurfaceElectricChargeDOF();

      affectedDOF = elemDOF;
      conditions = elemConds;
	
      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // load ID
      pushBackVector(surfaceElectricChargeBoundGaussPtsIdx,indexSet);
    }
    
  }

  /*********************************************************************/
  // Store all (local) deformation boundary gauss point indices and 
  // their boundary conditions.

  // displacement boundary constraints on surfaces
  for(int i=0;i<surfaceDispBoundElemIdx.size();i++) {
    int& currentElem = surfaceDispBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getSurfaceDeformationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getSurfaceDeformationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      
      blVector& affectedDOF = boundGaussPoints[idx].getDeformationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDeformationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;
      

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID
      pushBackVector(surfaceDispBoundGaussPtsIdx,indexSet);      
    }
    
  }

  // --------------------------------------------------------------------
  // rotation boundary constraints on surfaces
  for(int i=0;i<surfaceRotBoundElemIdx.size();i++) {
    int& currentElem = surfaceRotBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getSurfaceDeformationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getSurfaceDeformationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      
      blVector& affectedDOF = boundGaussPoints[idx].getDeformationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDeformationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;
      

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID
      pushBackVector(surfaceRotBoundGaussPtsIdx,indexSet);      
    }
    
  }

  /*********************************************************************/
  // electric boundary constraints on surfaces
  for(int i=0;i<surfaceElectricBoundElemIdx.size();i++) {
    int& currentElem = surfaceElectricBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getSurfaceElectricBoundDOF();
    dbVector& elemConds = bElems[currentElem].getSurfaceElectricBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      
      blVector& affectedDOF = boundGaussPoints[idx].getElectricBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getElectricBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;
      

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID
      pushBackVector(surfaceElectricBoundGaussPtsIdx,indexSet);      
    }
    
  }

  /*********************************************************************/
  // depolarisation boundary constraints on surfaces
  for(int i=0;i<surfaceDepolarisationBoundElemIdx.size();i++) {
    int& currentElem = surfaceDepolarisationBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getSurfaceDepolarisationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getSurfaceDepolarisationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      
      blVector& affectedDOF = boundGaussPoints[idx].getDepolarisationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDepolarisationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;
      

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID
      pushBackVector(surfaceDepolarisationBoundGaussPtsIdx,indexSet);      
    }
    
  }

  /*********************************************************************/
  // micro boundary constraints on surfaces
  for(int i=0;i<surfaceMicroBoundElemIdx.size();i++) {
    int& currentElem = surfaceMicroBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getSurfaceMicroBoundDOF();
    dbVector& elemConds = bElems[currentElem].getSurfaceMicroBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      
      blVector& affectedDOF = boundGaussPoints[idx].getMicroBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getMicroBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;
      

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID
      pushBackVector(surfaceMicroBoundGaussPtsIdx,indexSet);      
    }
    
  }


  // Calculate for all boundary gauss points their weights and their
  // surface normal vector.
  setSGaussWeightsNormals(InputData,modelData,logFile);

#ifdef _FEdebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"************ traction boundary gauss points *********"<<endl;
  for(int i=0;i<tractionBoundGaussPtsIdx.size();i++) {
    idx = tractionBoundGaussPtsIdx[i][0];
    int& weightID = tractionBoundGaussPtsIdx[i][1];
    int& loadID = tractionBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = boundGaussPoints[idx].getTractionDOF(loadID);
    dbVector& loads = boundGaussPoints[idx].getTraction(loadID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"traction: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"****** surface pressure boundary gauss points *******"<<endl;
  for(int i=0;i<surfacePressureBoundGaussPtsIdx.size();i++) {
    idx = surfacePressureBoundGaussPtsIdx[i][0];
    int& weightID = surfacePressureBoundGaussPtsIdx[i][1];
    int& loadID = surfacePressureBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    double& pressure = boundGaussPoints[idx].getSurfacePressure(loadID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"surface pressure = "<<pressure<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"******* surface moment boundary gauss points ********"<<endl;
  for(int i=0;i<surfaceMomentBoundGaussPtsIdx.size();i++) {
    idx = surfaceMomentBoundGaussPtsIdx[i][0];
    int& weightID = surfaceMomentBoundGaussPtsIdx[i][1];
    int& loadID = surfaceMomentBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = boundGaussPoints[idx].getSurfaceMomentDOF(loadID);
    dbVector& loads = boundGaussPoints[idx].getSurfaceMoment(loadID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"surface moment: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"*** surface electric charge boundary gauss points ***"<<endl;
  for(int i=0;i<surfaceElectricChargeBoundGaussPtsIdx.size();i++) {
    idx = surfaceElectricChargeBoundGaussPtsIdx[i][0];
    int& weightID = surfaceElectricChargeBoundGaussPtsIdx[i][1];
    int& loadID = surfaceElectricChargeBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = boundGaussPoints[idx].getSurfaceElectricChargeDOF(loadID);
    dbVector& loads = boundGaussPoints[idx].getSurfaceElectricCharge(loadID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"surface electric charge: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"******** displacement boundary gauss points *********"<<endl;
  for(int i=0;i<surfaceDispBoundGaussPtsIdx.size();i++) {
    idx = surfaceDispBoundGaussPtsIdx[i][0];
    int& weightID = surfaceDispBoundGaussPtsIdx[i][1];
    int& conditionID = surfaceDispBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDeformationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDeformationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* rotation boundary gauss points ************"<<endl;
  for(int i=0;i<surfaceRotBoundGaussPtsIdx.size();i++) {
    idx = surfaceRotBoundGaussPtsIdx[i][0];
    int& weightID = surfaceRotBoundGaussPtsIdx[i][1];
    int& conditionID = surfaceRotBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDeformationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDeformationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* electric boundary gauss points ************"<<endl;
  for(int i=0;i<surfaceElectricBoundGaussPtsIdx.size();i++) {
    idx = surfaceElectricBoundGaussPtsIdx[i][0];
    int& weightID = surfaceElectricBoundGaussPtsIdx[i][1];
    int& conditionID = surfaceElectricBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getElectricBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getElectricBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* depolarisation boundary gauss points ************"<<endl;
  for(int i=0;i<surfaceDepolarisationBoundGaussPtsIdx.size();i++) {
    idx = surfaceDepolarisationBoundGaussPtsIdx[i][0];
    int& weightID = surfaceDepolarisationBoundGaussPtsIdx[i][1];
    int& conditionID = surfaceDepolarisationBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDepolarisationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDepolarisationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* micro boundary gauss points ************"<<endl;
  for(int i=0;i<surfaceMicroBoundGaussPtsIdx.size();i++) {
    idx = surfaceMicroBoundGaussPtsIdx[i][0];
    int& weightID = surfaceMicroBoundGaussPtsIdx[i][1];
    int& conditionID = surfaceMicroBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getMicroBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getMicroBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
#endif
  
}
  
/**********************************************************************/
/**********************************************************************/
// Calculate the surface gauss points coordinate's and store the 
// surface load if applied.
void FEMGeometry::setLineGaussPoints(InputFileData* InputData,
				     std::map<std::string,double>& modelData,
				     std::ofstream& logFile,
				     PetscViewer& viewerMPI) {
  
  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];
  vector<FEMElement>& bElems = lineNodesElems;
  
  /*********************************************************************/
  // Determine a local portion of line Gauss points and the 
  // corresponding rootList.
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int elemStartIdx,elemEndIdx;
  int lineElemNum = lineNodesElems.size();
  int defaultPortion = (int)ceil((double)lineElemNum/size);

  if(defaultPortion*rank < lineElemNum && 
     defaultPortion*(rank+1) <= lineElemNum) {
    elemStartIdx = defaultPortion*rank;
    elemEndIdx = defaultPortion*(rank+1);
  }
  else if(defaultPortion*rank <= lineElemNum && 
	  defaultPortion*(rank+1) >= lineElemNum) {
    elemStartIdx = defaultPortion*rank;
    elemEndIdx = lineElemNum;
  }
  else {
    elemStartIdx = elemEndIdx = lineElemNum;
  }
  
  int elemPortion = elemEndIdx-elemStartIdx;  
  
  int root = 0;
  int lineGaussStartIdx = boundGaussPoints.size();
  int gaussStartIdx = lineGaussStartIdx;
  int localGaussPtPortion = 0;

  // determine a Gauss point rootList
  intVector rootList;
  
  for(int i=0;i<lineElemNum;i++) {
   
    intVector& intPts =
      bElems[i].getLineIntegrationPts();

    for(int j=0;j<intPts.size();j++) {
    
      // next processors turn
      if(i >= (root+1)*defaultPortion) {

	// next processor Gauss point start index
	if(root+1 == rank)
	  gaussStartIdx += rootList.size();

	root++;
      }
    
      if(root == rank)
	localGaussPtPortion++;
      
      intPts[j] = lineGaussStartIdx + rootList.size();

      // set the processor ID
      rootList.push_back(root);

    }

  }

  int lineGaussPtsNum = rootList.size();
  int gaussEndIdx = gaussStartIdx + localGaussPtPortion;

  globalBGaussPtsNum += lineGaussPtsNum;
  boundGaussPoints.resize(globalBGaussPtsNum);

#ifdef _FEdebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"************** LINE GAUSS POINTS *********************"<<endl;
  logFile<<"******************************************************"<<endl;
  logFile<<"elemPortion = "<<elemPortion
	 <<" elemStartIdx = "<<elemStartIdx
	 <<" elemEndIdx = "<<elemEndIdx
	 <<" lineGaussPtsNum "<<lineGaussPtsNum
	 <<" localGaussPtPortion = "<<localGaussPtPortion<<endl;
  logFile<<"******** element Gauss point connectivity ************"<<endl;
  for(int i=0;i<lineElemNum;i++) {
    intVector& intPts =
      bElems[i].getLineIntegrationPts();
    logFile<<"Element "<<bElems[i].getGlobalID()<<": ";
    for(int j=0;j<intPts.size();j++) 
      logFile<<intPts[j]<<" ";
    logFile<<endl;
  }
  logFile<<"********************* rootList ***********************"<<endl;
  for(int i=0;i<rootList.size();i++)
    logFile<<i<<".) "<<rootList[i]<<endl;
  logFile<<"gaussStartIdx = "<<gaussStartIdx<<endl;
  logFile<<"gaussEndIdx = "<<gaussEndIdx<<endl;
#endif

  /*********************************************************************/
  // Determine the coordinates of a local portion of boundary Gauss
  // points.
  int node;
  int idx;

  // Loop over the local portion of boundary elements.
  for(int i=elemStartIdx;i<elemEndIdx;i++) {
    intVector& intPts = 
      bElems[i].getLineIntegrationPts();
    intVector& nodes = bElems[i].getNodes();
    int nodesPerElem = nodes.size();
    
    // Get the set of FEM shape functions.
    dbMatrix& sFuncs = bElems[i].getLineShapeFuncOrds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      dbVector& gCoords = boundGaussPoints[idx].getCoords();      
      gCoords.resize(usedDims);
      
      for(int k=0;k<nodes.size();k++) {
	node = nodes[k]-1;
	
	gCoords[0] += 
	  sFuncs[j][k]*particles[node].getCoord(0);	
	gCoords[1] += 
	  sFuncs[j][k]*particles[node].getCoord(1);
	gCoords[2] += 
	  sFuncs[j][k]*particles[node].getCoord(2);
      }

      // set location of the current Gauss point
      intVector& elemInfo = boundGaussPoints[idx].getElementInfo();
      elemInfo.resize(3);
      elemInfo[0] = i; // element ID
      elemInfo[1] = j; // Gauss point ID within the element
      elemInfo[2] = bElems[i].getMotherElementID(); // corresponding volume element

      //idx++;
    } 

  }

  // --------------------------------------------------------------------
  // Loop over the local portion of boundary elements.
  for(int i=0;i<lineNodesElems.size();i++) {
    int& elemMatID = bElems[i].getMaterialID();
    intVector& intPts = 
      bElems[i].getLineIntegrationPts();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      int& matID = boundGaussPoints[idx].getMaterialID();
      matID = elemMatID;
    }

  }

#ifdef _FEdebugMode_
  logFile<<"******************* template element ****************"<<endl;
  std::map<std::string,double>& backGroundMeshInfo = 
    InputData->getBackGroundMeshInfo();
  double gCoord;
  int ndsPerElem = (int)backGroundMeshInfo["nodesPerLineElement"];
  int gPtsPerElem = (int)backGroundMeshInfo["gaussPointsPerLineElement"];
  FEMElement elem(usedDOF);
  intVector& intPts = elem.getLineIntegrationPts();
  intVector& nds = elem.getNodes();
  int& elemType = elem.getElemType();
  int& elemOrder = elem.getElemOrder();
  intPts.resize(gPtsPerElem);
  nds.resize(ndsPerElem);
  elemType = (int)backGroundMeshInfo["elemType"];
  elemOrder = (int)backGroundMeshInfo["elemOrder"];
  elem.setLineElementTemplate(getLineElementTemplate(InputData,elem,logFile));
  elem.setLineGaussSet(getLineGaussSet(InputData,elem,logFile));
  ElementTemplate* FESet = elem.getLineElementTemplate();
  GaussPointSet* GaussSet = elem.getLineGaussSet();
  dbMatrix& sFuncs = getLineFEShapeFunctions(InputData,elem,logFile);
  for(int i=0;i<gPtsPerElem;i++) {
    gCoord = 0;
    for(int j=0;j<ndsPerElem;j++)
	gCoord += sFuncs[i][j]*FESet->nodalCoords[j][0];
    logFile<<"point "<<i<<": "<<GaussSet->coord[i][0]<<" ?= "
	   <<gCoord<<endl;
  }
  logFile<<"**************** local gauss coords ******************"<<endl;
  logFile<<"Bound Gauss point portion: "
	 <<localGaussPtPortion<<" = "<<gaussEndIdx-gaussStartIdx<<endl;
  logFile<<"gaussStartIdx = "<<gaussStartIdx<<" gaussEndIdx = "
	 <<gaussEndIdx<<"; lineGaussPtsNum = "<<lineGaussPtsNum<<endl;
  for(int i=gaussStartIdx;i<gaussEndIdx;i++) {
    int& matID = boundGaussPoints[i].getMaterialID();
    intVector& elemInfo = boundGaussPoints[i].getElementInfo();
    dbVector& gCoords = boundGaussPoints[i].getCoords();      
    logFile<<i<<".) BOUND GAUSSPOINT: ";
    for(int j=0;j<gCoords.size();j++)
      logFile<<gCoords[j]<<" ";
    logFile<<"; l-elem "<<elemInfo[0]<<", local gPoint ID "<<elemInfo[1]
	   <<" volID "<<elemInfo[2]<<" matID "<<matID<<endl;
  }
#endif

  // Assemble all boundary gauss points coordinates from all processors,
  // store these.
  
  // Loop over the local calculated portion of Gauss points.
  for(int i=lineGaussStartIdx,j=0;i<globalBGaussPtsNum;i++,j++) {
    boundGaussPoints[i].setGlobalID(i);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&boundGaussPoints[i].getCoord(0),usedDims,MPI_DOUBLE,
	      rootList[j],MPI_COMM_WORLD);

    intVector& elemInfo = boundGaussPoints[i].getElementInfo();

    if(rootList[j] != rank)
      elemInfo.resize(3);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&elemInfo[0],3,MPI_INT,rootList[j],MPI_COMM_WORLD);
  }

  /*********************************************************************/
  // Store the line loads and the normal vector on all boundary gauss 
  // points, and the (local) indices of these boundary Gauss points.
  intVector indexSet(3);
 
  for(int i=0;i<lineForceElemIdx.size();i++) {
    int& currentElem = lineForceElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getLineIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineForceDOF();
    dbVector& elemConds = bElems[currentElem].getLineForce();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      // set the loading
      blVector& affectedDOF = boundGaussPoints[idx].getLineForceDOF();
      dbVector& conditions = boundGaussPoints[idx].getLineForce();

      affectedDOF = elemDOF;
      conditions = elemConds;

      // set the surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());
	
      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // load ID

      pushBackVector(lineForceBoundGaussPtsIdx,indexSet);
    }

  }

  //---------------------------------------------------------------------
  // Store all (local) displacement boundary gauss point indices and 
  // boundary conditions.
  for(int i=0;i<lineDispBoundElemIdx.size();i++) {
    int& currentElem = lineDispBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getLineIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineDeformationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getLineDeformationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      // set conditions
      blVector& affectedDOF = boundGaussPoints[idx].getDeformationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDeformationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;

      // set surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID

      pushBackVector(lineDispBoundGaussPtsIdx,indexSet);
    }

  }

  // --------------------------------------------------------------------
  // Store all (local) rotation boundary gauss point indices and 
  // boundary conditions.

  for(int i=0;i<lineRotBoundElemIdx.size();i++) {
    int& currentElem = lineRotBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getLineIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineDeformationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getLineDeformationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      // set conditions
      blVector& affectedDOF = boundGaussPoints[idx].getDeformationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDeformationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;

      // set surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID

      pushBackVector(lineRotBoundGaussPtsIdx,indexSet);
    }

  }

  /*********************************************************************/
  // Store all (local) electric boundary gauss point indices and 
  // boundary conditions.

  for(int i=0;i<lineElectricBoundElemIdx.size();i++) {
    int& currentElem = lineElectricBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getLineIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineElectricBoundDOF();
    dbVector& elemConds = bElems[currentElem].getLineElectricBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      // set conditions
      blVector& affectedDOF = boundGaussPoints[idx].getElectricBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getElectricBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;

      // set surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID

      pushBackVector(lineElectricBoundGaussPtsIdx,indexSet);
    }

  }


  /*********************************************************************/
  // Store all (local) depolarisation boundary gauss point indices and 
  // boundary conditions.

  for(int i=0;i<lineDepolarisationBoundElemIdx.size();i++) {
    int& currentElem = lineDepolarisationBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getLineIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineDepolarisationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getLineDepolarisationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

      // set conditions
      blVector& affectedDOF = boundGaussPoints[idx].getDepolarisationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDepolarisationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;

      // set surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID

      pushBackVector(lineDepolarisationBoundGaussPtsIdx,indexSet);
    }

  }

  /*********************************************************************/
  // Calculate for all boundary gauss points their weights.

  setLineGaussPtsWeights(InputData,modelData,logFile);

#ifdef _FEdebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"globalBGaussPtsNum = "<<globalBGaussPtsNum
	 <<" = surfaceGaussPtsNum + lineGaussPtsNum = "
	 <<lineGaussStartIdx+lineGaussPtsNum<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"********** line force boundary gauss points *********"<<endl;
  for(int i=0;i<lineForceBoundGaussPtsIdx.size();i++) {
    idx = lineForceBoundGaussPtsIdx[i][0];
    int& weightID = lineForceBoundGaussPtsIdx[i][1];
    int& loadID = lineForceBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = boundGaussPoints[idx].getLineForceDOF(loadID);
    dbVector& loads = boundGaussPoints[idx].getLineForce(loadID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"lineForce: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* displacement boundary gauss points *********"<<endl;
  for(int i=0;i<lineDispBoundGaussPtsIdx.size();i++) {
    idx = lineDispBoundGaussPtsIdx[i][0];
    int& weightID = lineDispBoundGaussPtsIdx[i][1];
    int& conditionID = lineDispBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDeformationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDeformationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* rotation boundary gauss points *********"<<endl;
  for(int i=0;i<lineRotBoundGaussPtsIdx.size();i++) {
    idx = lineRotBoundGaussPtsIdx[i][0];
    int& weightID = lineRotBoundGaussPtsIdx[i][1];
    int& conditionID = lineRotBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDeformationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDeformationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"*********** electric boundary gauss points **********"<<endl;
  for(int i=0;i<lineElectricBoundGaussPtsIdx.size();i++) {
    idx = lineElectricBoundGaussPtsIdx[i][0];
    int& weightID = lineElectricBoundGaussPtsIdx[i][1];
    int& conditionID = lineElectricBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getElectricBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getElectricBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"*********** depolarisation boundary gauss points **********"<<endl;
  for(int i=0;i<lineDepolarisationBoundGaussPtsIdx.size();i++) {
    idx = lineDepolarisationBoundGaussPtsIdx[i][0];
    int& weightID = lineDepolarisationBoundGaussPtsIdx[i][1];
    int& conditionID = lineDepolarisationBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDepolarisationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDepolarisationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate from line elements surface gauss points coordinate's and 
// integration weight. Store the line force as surface load if applied.
void FEMGeometry::setLineSurfaceGPoints(InputFileData* InputData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile,
					PetscViewer& viewerMPI) {
  
  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];

  std::map<std::string,double>& backGroundMeshInfo = 
    InputData->getBackGroundMeshInfo();
  int nodesPerSurfaceElem = 
    (int)backGroundMeshInfo["nodesPerSurfaceElement"];
  int gaussPtsPerSurfaceElem = 
    (int)backGroundMeshInfo["gaussPointsPerSurfaceElement"];
  vector<FEMElement>& bElems = lineNodesElems;
  
  /*********************************************************************/
  // Determine a local portion of line Gauss points and the 
  // corresponding rootList.
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int elemStartIdx,elemEndIdx;
  int lineElemNum = lineNodesElems.size();
  int defaultPortion = (int)ceil((double)lineElemNum/size);

  if(defaultPortion*rank < lineElemNum && 
     defaultPortion*(rank+1) <= lineElemNum) {
    elemStartIdx = defaultPortion*rank;
    elemEndIdx = defaultPortion*(rank+1);
  }
  else if(defaultPortion*rank <= lineElemNum && 
	  defaultPortion*(rank+1) >= lineElemNum) {
    elemStartIdx = defaultPortion*rank;
    elemEndIdx = lineElemNum;
  }
  else {
    elemStartIdx = elemEndIdx = lineElemNum;
  }
  
  int elemPortion = elemEndIdx-elemStartIdx;  
  
  int root = 0;
  int lineGaussStartIdx = boundGaussPoints.size();
  int gaussStartIdx = lineGaussStartIdx;
  int localGaussPtPortion = 0;

  // determine a Gauss point rootList
  intVector rootList;

  for(int i=0;i<lineElemNum;i++) {
   
    intVector& intPts =
      bElems[i].getSurfaceIntegrationPts();
    intPts.resize(gaussPtsPerSurfaceElem);

    for(int j=0;j<intPts.size();j++) {
    
      // next processors turn
      if(i >= (root+1)*defaultPortion) {

	// next processor Gauss point start index
	if(root+1 == rank)
	  gaussStartIdx += rootList.size();

	root++;
      }
    
      if(root == rank)
	localGaussPtPortion++;
      
      intPts[j] = lineGaussStartIdx + rootList.size();

      // set the processor ID
      rootList.push_back(root);

    }

  }

  int gaussPtsNum = rootList.size();
  int gaussEndIdx = gaussStartIdx + localGaussPtPortion;

  globalBGaussPtsNum += gaussPtsNum;
  boundGaussPoints.resize(globalBGaussPtsNum);

#ifdef _FEdebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"************** LINE-SURFACE GAUSS POINTS *************"<<endl;
  logFile<<"******************************************************"<<endl;
  logFile<<"elemPortion = "<<elemPortion
	 <<" elemStartIdx = "<<elemStartIdx
	 <<" elemEndIdx = "<<elemEndIdx<<endl;
  logFile<<"************* elem gauss connectivity ****************"<<endl;
  for(int i=elemStartIdx;i<elemEndIdx;i++) {
    logFile<<"Line element "<<i<<" : ";
    intVector& intPts =
      bElems[i].getSurfaceIntegrationPts();
    for(int j=0;j<intPts.size();j++)
      logFile<<intPts[j]<<" ";
    logFile<<endl;
  }
  logFile<<" localGaussPtPortion = "<<localGaussPtPortion
	 <<" globalBGaussPtsNum = "<<globalBGaussPtsNum
	 <<" lineGaussStartIdx "<<lineGaussStartIdx
	 <<" gaussPtsNum "<<gaussPtsNum<<endl;
  logFile<<"************* temporary parallel distribution ********"<<endl;
  logFile<<"gaussStartIdx = "<<gaussStartIdx<<endl;
  logFile<<"gaussEndIdx = "<<gaussEndIdx<<endl;
  logFile<<"********************* rootList ***********************"<<endl;
  for(int i=0;i<rootList.size();i++)
    logFile<<i<<".) "<<rootList[i]<<endl;
  logFile<<"******************************************************"<<endl;
#endif

  double lineWidth = InputData->getValue("lineWidth");
  intMatrix e = getPermutations(3);

  int ptcle1,ptcle2;
  double absValue;
  dbVector arcVec(usedDims);
  dbVector tangVec(usedDims);

  vector<Particle> 
    newParticles(bElems.size()*nodesPerSurfaceElem,Particle(usedDOF));

  int currentPtcle = 0;

  // Loop over the local portion of boundary elements to convert
  // the line elements to surface elements and create the needed new
  // particles.
  for(int i=elemStartIdx;i<elemEndIdx;i++) {

    // calculate a tangential vector to surface element to be created
    dbVector& normalVec = bElems[i].getSurfaceNormal();    
    intVector lineNodes = bElems[i].getNodes();

    ptcle1 = lineNodes[0] - 1;
    ptcle2 = lineNodes[lineNodes.size()-1] - 1;

#ifdef _FEdebugMode_
    logFile<<"Line element "<<i<<" : "<<endl;
    logFile<<"node1 "<<ptcle1+1<<" coords ";
    for(int k=0;k<normalVec.size();k++)
      logFile<<particles[ptcle1].getCoord(k)<<" ";
    logFile<<endl;
    logFile<<"node2 "<<ptcle2+1<<" coords ";
    for(int k=0;k<normalVec.size();k++)
      logFile<<particles[ptcle2].getCoord(k)<<" ";
    logFile<<endl;
#endif

    // convert the line element to a vector
    for(int k=0;k<normalVec.size();k++)

      arcVec[k] = particles[ptcle1].getCoord(k) - 
	particles[ptcle2].getCoord(k);

    // tangential vector is the crossproduct of normal and line element 
    // vector
    clearArray(tangVec);

    for(int k=0;k<e.size();k++)

      tangVec[e[k][2]] += e[k][3]*normalVec[e[k][0]]*arcVec[e[k][1]];

    absValue = 0;

    for(int k=0;k<tangVec.size();k++)
      absValue += tangVec[k]*tangVec[k];

    absValue = sqrt(absValue);

    // give the tangential vector as length 'lineWidth/2'
    for(int k=0;k<tangVec.size();k++)

      tangVec[k] = tangVec[k]/absValue*lineWidth/2.0;

    intVector& surfaceNodes = bElems[i].getNodes();
    surfaceNodes.resize(nodesPerSurfaceElem);

    // compute coordinates of the surface element nodes
    for(int k=0;k<surfaceNodes.size();k++) {

      surfaceNodes[k] = currentPtcle + 1;

      dbVector& newCoords = newParticles[currentPtcle].getCoords();
      newCoords.resize(usedDims);

      dbVector& oldCoords1 = particles[ptcle1].getCoords();
      dbVector& oldCoords2 = particles[ptcle2].getCoords();

      switch(k) {

      case 0:
      
	for(int m=0;m<usedDims;m++)
	  newCoords[m] = oldCoords1[m] + tangVec[m];

	break;

      case 1:
      
	for(int m=0;m<usedDims;m++)
	  newCoords[m] = oldCoords1[m] - tangVec[m];

	break;

      case 2:
      
	for(int m=0;m<usedDims;m++)
	  newCoords[m] = oldCoords2[m] - tangVec[m];

	break;

      case 3:

	for(int m=0;m<usedDims;m++)
	  newCoords[m] = oldCoords2[m] + tangVec[m];

	break;

      default:

	logFile<<"In FEMGeometry::setLineSurfaceGPoints node ID "<<k
	       <<" is not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
	
      }

      currentPtcle++;

    }

#ifdef _FEdebugMode_
    logFile<<"normal: ";
    for(int k=0;k<normalVec.size();k++)
      logFile<<normalVec[k]<<" ";
    logFile<<endl;
    logFile<<"arc: ";
    for(int k=0;k<arcVec.size();k++)
      logFile<<arcVec[k]<<" ";
    logFile<<endl;
    logFile<<"tangent: ";
    for(int k=0;k<tangVec.size();k++)
      logFile<<tangVec[k]<<" ";
    logFile<<endl;
    logFile<<"-----------------------------------------------------"<<endl;
    for(int k=0;k<surfaceNodes.size();k++) {
      dbVector& newCoords = newParticles[surfaceNodes[k]-1].getCoords();
      logFile<<"node "<<surfaceNodes[k]<<": ";
      for(int m=0;m<newCoords.size();m++)
	logFile<<newCoords[m]<<" ";
      logFile<<endl;
    }
    logFile<<"-----------------------------------------------------"<<endl;
#endif

  }

  /*********************************************************************/
  // Loop over the local portion of boundary elements and compute the
  // surface Gauss point coodinates and their integration weight
  int integrationMethod = (int)modelData["integrationMethod"];

  if(integrationMethod != 1) {

    logFile<<"In FEMGeometry::setLineSurfaceGPoints integration method\n"
	   <<"is not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);

  }

  int node;

  for(int i=elemStartIdx;i<elemEndIdx;i++) {

    FEMElement& elem = bElems[i];

    intVector& intPts =
      elem.getSurfaceIntegrationPts();

    // set surface element approximation tools
    elem.setVolumeElementTemplate(getVolumeElementTemplate(InputData,elem,logFile));
    elem.setSurfaceElementTemplate(getSurfaceElementTemplate(InputData,elem,logFile));
    elem.setSurfaceGaussSet(getSurfaceGaussSet(InputData,elem,logFile));
    elem.setSurfaceShapeFuncOrds(getSurfaceFEShapeFunctions(InputData,elem,logFile));

    intVector& nodes = elem.getNodes();
    
    // Get the set of FEM shape functions.
    dbMatrix& sFuncs = elem.getSurfaceShapeFuncOrds();

    for(int j=0;j<intPts.size();j++) {
      int& idx = intPts[j];
      dbVector& gCoords = boundGaussPoints[idx].getCoords();      
      gCoords.resize(usedDims);
      
      for(int k=0;k<nodes.size();k++) {
	node = nodes[k]-1;
	
	gCoords[0] += 
	  sFuncs[j][k]*newParticles[node].getCoord(0);	
	gCoords[1] += 
	  sFuncs[j][k]*newParticles[node].getCoord(1);
	gCoords[2] += 
	  sFuncs[j][k]*newParticles[node].getCoord(2);
      }

      // set location of the current Gauss point
      intVector& elemInfo = boundGaussPoints[idx].getElementInfo();
      elemInfo.resize(3);
      elemInfo[0] = i; // element ID
      elemInfo[1] = j; // Gauss point ID within the element
      elemInfo[2] = elem.getMotherElementID(); // corresponding volume element

    } 

 
  }

#ifdef _FEdebugMode_
  logFile<<"**************** local gauss coords ******************"<<endl;
  logFile<<"Bound Gauss point portion: "
	 <<localGaussPtPortion<<" = "<<gaussEndIdx-gaussStartIdx<<endl;
  logFile<<"gaussStartIdx = "<<gaussStartIdx<<" gaussEndIdx = "
	 <<gaussEndIdx<<"; gaussPtsNum = "<<gaussPtsNum<<endl;
  for(int i=gaussStartIdx;i<gaussEndIdx;i++) {
    intVector& elemInfo = boundGaussPoints[i].getElementInfo();
    dbVector& gCoords = boundGaussPoints[i].getCoords();      
    logFile<<i<<".) BOUND GAUSSPOINT: ";
    for(int j=0;j<gCoords.size();j++)
      logFile<<gCoords[j]<<" ";
    logFile<<"; l-elem "<<elemInfo[0]<<", local gPoint ID "<<elemInfo[1]
	   <<" volID "<<elemInfo[2]<<endl;
  }
#endif

  // Assemble all boundary gauss points coordinates from all processors,
  // store these.
  
  // Loop over the local calculated portion of Gauss points.
  for(int i=lineGaussStartIdx,j=0;i<globalBGaussPtsNum;i++,j++) {
    boundGaussPoints[i].setGlobalID(i);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&boundGaussPoints[i].getCoord(0),usedDims,MPI_DOUBLE,
	      rootList[j],MPI_COMM_WORLD);

    intVector& elemInfo = boundGaussPoints[i].getElementInfo();

    if(rootList[j] != rank)
      elemInfo.resize(3);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&elemInfo[0],3,MPI_INT,rootList[j],MPI_COMM_WORLD);
  }

  /*********************************************************************/
  // Store the line loads and the normal vector on all boundary gauss 
  // points, and the (local) indices of these boundary Gauss points.
  intVector indexSet(3);
 
  for(int i=0;i<lineForceElemIdx.size();i++) {
    int& currentElem = lineForceElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineForceDOF();
    dbVector& elemConds = bElems[currentElem].getLineForce();

    for(int j=0;j<intPts.size();j++) {
      int& idx = intPts[j];

      // set the loading
      blVector& affectedDOF = boundGaussPoints[idx].getTractionDOF();
      dbVector& conditions = boundGaussPoints[idx].getTraction();

      affectedDOF = elemDOF;
      conditions = elemConds;

      for(int k=0;k<conditions.size();k++)

	conditions[k] /= lineWidth;

      // set the surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());
	
      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // load ID

      pushBackVector(lineForceBoundGaussPtsIdx,indexSet);
    }

  }

  //---------------------------------------------------------------------
  // Store all (local) displacement boundary gauss point indices and 
  // boundary conditions.
  for(int i=0;i<lineDispBoundElemIdx.size();i++) {
    int& currentElem = lineDispBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineDeformationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getLineDeformationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      int& idx = intPts[j];

      // set conditions
      blVector& affectedDOF = boundGaussPoints[idx].getDeformationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDeformationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;

      // set surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID

      pushBackVector(lineDispBoundGaussPtsIdx,indexSet);
    }

  }

  // --------------------------------------------------------------------
  // Store all (local) rotation boundary gauss point indices and 
  // boundary conditions.
  //if(lineRotBoundElemIdx.size() > 0) {
  //logFile<<"In function FEMGeometry::setLineGaussPoints "
  //   <<"line rotation boundary constraints are not supported!"<<endl;
  //MPI_Abort(MPI_COMM_WORLD,1);
  //}

  for(int i=0;i<lineRotBoundElemIdx.size();i++) {
    int& currentElem = lineRotBoundElemIdx[i][0];

    intVector& intPts = 
      bElems[currentElem].getSurfaceIntegrationPts();
    blVector& elemDOF = bElems[currentElem].getLineDeformationBoundDOF();
    dbVector& elemConds = bElems[currentElem].getLineDeformationBoundConds();

    for(int j=0;j<intPts.size();j++) {
      int& idx = intPts[j];

      // set conditions
      blVector& affectedDOF = boundGaussPoints[idx].getDeformationBoundDOF();
      dbVector& conditions = boundGaussPoints[idx].getDeformationBoundConds();

      affectedDOF = elemDOF;
      conditions = elemConds;

      // set surface normal vector
      boundGaussPoints[idx].
	setSurfaceNormal(bElems[currentElem].getSurfaceNormal());

      indexSet[0] = idx; // Gauss point ID
      indexSet[1] = 0; // weight ID
      indexSet[2] = 0; // condition ID

      pushBackVector(lineRotBoundGaussPtsIdx,indexSet);
    }

  }

  /*********************************************************************/
  // determine the integration weight of a local portion of boundary
  // Gauss points
  double surface;
  dbVector normalVec(3);

  int m=0;
      
  // loop over a portion of surface boundary Gauss points
  for(int i=gaussStartIdx;i<gaussEndIdx;i++) {
    int& currentPoint = i;

    // Calculate surface normal vector and its length - the surface 
    // of the volume element.
    double& weight = boundGaussPoints[currentPoint].getWeight();
    intVector& elemInfo = boundGaussPoints[currentPoint].getElementInfo();

    int& elemID = elemInfo[0];
    int& gaussID = elemInfo[1];

    FEMElement& elem = bElems[elemID];
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* SurfaceGaussSet = elem.getSurfaceGaussSet();
    intVector& nodes = elem.getNodes();

    FEVolumeSet->surfaceMetric(newParticles,nodes,
			       SurfaceGaussSet->coord[gaussID],
			       normalVec,surface,logFile);

    weight = fabs(surface)*SurfaceGaussSet->weight[gaussID];
      
  }

#ifdef _FEdebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"************* surface element weights ***************"<<endl;
  double area,weight,sumWeights,sumAreas;
  double sumElemAreas = 0;
  double sumElemWeights = 0;
  for(int i=gaussStartIdx;i<gaussEndIdx;i++) {
    sumAreas = sumWeights = 0;
    int& currentPoint = i;
    intVector& elemInfo = boundGaussPoints[currentPoint].getElementInfo();
    int& elemID = elemInfo[0];
    int& gaussID = elemInfo[1];
    FEMElement& elem = bElems[elemID];
    intVector& nodes = elem.getNodes();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* SurfaceGaussSet = elem.getSurfaceGaussSet();
    logFile<<"Element "<<elemID<<" nodes: ";
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<endl;
    logFile<<"delta areas/weights:"<<endl;
    FEVolumeSet->surfaceMetric(newParticles,nodes,
			       SurfaceGaussSet->coord[gaussID],
			       normalVec,surface,logFile);
    area = surface;
    weight = area*SurfaceGaussSet->weight[gaussID];
    sumElemAreas += area;
    sumElemWeights += weight;
    logFile<<area<<" / "<<weight<<endl;
  }
  logFile<<"System area/weight = "<<sumElemAreas<<" / "
	 <<sumElemWeights<<endl;
#endif
 
  /*********************************************************************/ 
  // Assemble all boundary gauss points weights from all processors 
  // and store them.
  
  // Loop over all boundary integration points.
  for(int i=lineGaussStartIdx,j=0;i<globalBGaussPtsNum;i++,j++) {

    MPI_Bcast(&boundGaussPoints[i].getWeight(),1,MPI_DOUBLE,
	      rootList[j],MPI_COMM_WORLD);
  }

#ifdef _FEdebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"globalBGaussPtsNum = "<<globalBGaussPtsNum<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"********** line force boundary gauss points *********"<<endl;
  for(int i=0;i<lineForceBoundGaussPtsIdx.size();i++) {
    int& idx = lineForceBoundGaussPtsIdx[i][0];
    int& weightID = lineForceBoundGaussPtsIdx[i][1];
    int& loadID = lineForceBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = boundGaussPoints[idx].getTractionDOF(loadID);
    dbVector& loads = boundGaussPoints[idx].getTraction(loadID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"Traction: ";
    for(int j=0;j<loads.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<loads[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* displacement boundary gauss points *********"<<endl;
  for(int i=0;i<lineDispBoundGaussPtsIdx.size();i++) {
    int& idx = lineDispBoundGaussPtsIdx[i][0];
    int& weightID = lineDispBoundGaussPtsIdx[i][1];
    int& conditionID = lineDispBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDeformationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDeformationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"********* rotation boundary gauss points *********"<<endl;
  for(int i=0;i<lineRotBoundGaussPtsIdx.size();i++) {
    int& idx = lineRotBoundGaussPtsIdx[i][0];
    int& weightID = lineRotBoundGaussPtsIdx[i][1];
    int& conditionID = lineRotBoundGaussPtsIdx[i][2];
    double& weight = boundGaussPoints[idx].getIntWeight(weightID);
    blVector& affectedDOF = 
      boundGaussPoints[idx].getDeformationBoundDOF(conditionID);
    dbVector& conditions = 
      boundGaussPoints[idx].getDeformationBoundConds(conditionID);
    dbVector& normal = boundGaussPoints[idx].getSurfaceNormal();
    logFile<<i<<".) BOUND GAUSSPOINT "<<boundGaussPoints[idx].getGlobalID()<<": ";
    logFile<<boundGaussPoints[idx].getCoord(0)<<" "
	   <<boundGaussPoints[idx].getCoord(1)<<" "
	   <<boundGaussPoints[idx].getCoord(2)<<" ";
    logFile<<"boundary conditions: ";
    for(int j=0;j<conditions.size();j++)
      if(affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<conditions[j]<<" ";
    logFile<<"normal: ";
    for(int j=0;j<normal.size();j++)
      logFile<<normal[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calucate for all local volume Gauss points their weight.
void FEMGeometry::setVolumeGaussWeights(InputFileData* InputData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile) {


  using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  double minWeight,localMinWeight;
  
  // Loop over all local elements.
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];

    intVector& intPts = elem.getVolumeIntegrationPts();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* VolumeGaussSet = elem.getVolumeGaussSet();
    intVector& nodes = elem.getNodes();

    for(int j=0;j<intPts.size();j++) {
      int& idx = intPts[j]; 
      double& weight = gaussPoints[idx].getWeight();

      weight = 
	FEVolumeSet->metricFactorVolume(particles,nodes,
					VolumeGaussSet->coord[j],
					logFile) 
	* VolumeGaussSet->weight[j];

      if(i == 0 && j ==0)
	localMinWeight = weight;

      else if(weight < minWeight)
	localMinWeight = weight;

    }

  }

  MPI_Allreduce(&localMinWeight,&minWeight,1,
		MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  
  logFile<<"minimum volume integration weight = "<<minWeight<<endl;


#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"***** volume gauss weights - element volumes ********"<<endl;
  double sumVols,sumWeights,volume,weight;
  double sumElemVols = 0;
  double sumElemWeights = 0;
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];
    intVector& nodes = elem.getNodes();
    intVector& intPts = elem.getVolumeIntegrationPts();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* VolumeGaussSet = elem.getVolumeGaussSet();
    sumVols = sumWeights = 0;
    logFile<<"Element "<<i<<" volumes/weights: \n";
    for(int j=0;j<intPts.size();j++) {
      int& idx = intPts[j];
      volume = 
	FEVolumeSet->metricFactorVolume(particles,nodes,
					VolumeGaussSet->coord[j],
					logFile);
      weight = volume*VolumeGaussSet->weight[j];
      sumVols += volume;
      sumWeights += weight;
      logFile<<volume<<" / "<<weight<<endl;
    }
    logFile<<"sum volumes/weights: "<<sumVols<<" / "<<sumWeights<<endl;
    sumElemVols += sumVols;
    sumElemWeights += sumWeights;
  }
  logFile<<"system volume/weight sum = "<<sumElemVols<<" / "
	 <<sumElemWeights<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all surface Gauss points their weight and their 
// surface normal vector.
void FEMGeometry::setSGaussWeightsNormals(InputFileData* InputData,
					  std::map<std::string,double>& modelData,
					  std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];

  // Determine the a local portion of boundary elements for calculating.
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  int startIdx,endIdx;

  allSurfaceBoundGaussPtsIdx =
    getAllSurfaceGaussPtsIdx(logFile);

  intMatrix& surfaceGaussPtsIdx = allSurfaceBoundGaussPtsIdx;

  int surfaceGaussPtsNum = surfaceGaussPtsIdx.size();
  int gaussPortion = (int)ceil((double)surfaceGaussPtsNum/size);

  if(gaussPortion*rank < surfaceGaussPtsNum && 
     gaussPortion*(rank+1) <= surfaceGaussPtsNum) {
    startIdx = gaussPortion*rank;
    endIdx = gaussPortion*(rank+1);
  }
  else if(gaussPortion*rank <= surfaceGaussPtsNum && 
	  gaussPortion*(rank+1) >= surfaceGaussPtsNum) {
    startIdx = gaussPortion*rank;
    endIdx = surfaceGaussPtsNum;
  }
  else {
    startIdx = endIdx = surfaceGaussPtsNum;
  }
  
  gaussPortion = endIdx-startIdx;
    
  double surface;
  intVector localRootList(surfaceGaussPtsNum);
  dbVector normalVec(usedDims);

  int m=0;

  // loop over a portion of surface boundary Gauss points
  for(int i=startIdx;i<endIdx;i++) {
    int& currentPoint = surfaceGaussPtsIdx[i][0];
    localRootList[i] = rank;

    // Calculate surface normal vector and its length - the surface 
    // of the volume element.
    double& weight = boundGaussPoints[currentPoint].getWeight();
    intVector& elemInfo = boundGaussPoints[currentPoint].getElementInfo();

    int& elemID = elemInfo[0];
    int& gaussID = elemInfo[1];

    FEMElement& elem = surfaceNodesElems[elemID];
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* SurfaceGaussSet = elem.getSurfaceGaussSet();
    intVector& nodes = elem.getNodes();

    FEVolumeSet->surfaceMetric(particles,nodes,
			       SurfaceGaussSet->coord[gaussID],
			       normalVec,surface,logFile);

    weight = fabs(surface)*SurfaceGaussSet->weight[gaussID];
    boundGaussPoints[currentPoint].setSurfaceNormal(normalVec);
      
  }
      
#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"************* local surface normals *****************"<<endl;
  for(int i=startIdx;i<endIdx;i++) {
    int& currentPoint = surfaceGaussPtsIdx[i][0];
    logFile<<"BOUND GAUSS POINT "<<currentPoint<<" ";
    dbVector& surfaceNormal = 
      boundGaussPoints[currentPoint].getSurfaceNormal();
    logFile<<"surface normal: "<<surfaceNormal[0]<<" "
	   <<surfaceNormal[1]<<" "<<surfaceNormal[2]<<endl;
  }
  logFile<<"*****************************************************"<<endl;
  logFile<<"******* surface element weights and normals *********"<<endl;
  double area,weight,sumWeights,sumAreas;
  double sumElemAreas = 0;
  double sumElemWeights = 0;
  for(int i=startIdx;i<endIdx;i++) {
    sumAreas = sumWeights = 0;
    int& currentPoint = surfaceGaussPtsIdx[i][0];
    intVector& elemInfo = boundGaussPoints[currentPoint].getElementInfo();
    int& elemID = elemInfo[0];
    int& gaussID = elemInfo[1];
    FEMElement& elem = surfaceNodesElems[elemID];
    intVector& nodes = elem.getNodes();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* SurfaceGaussSet = elem.getSurfaceGaussSet();
    logFile<<"Element "<<elemID<<" nodes: ";
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<endl;
    FEVolumeSet->surfaceMetric(particles,nodes,
			       SurfaceGaussSet->coord[gaussID],
			       normalVec,surface,logFile);
    area = surface;
    weight = area*SurfaceGaussSet->weight[gaussID];
    sumElemAreas += area;
    sumElemWeights += weight;
    logFile<<"normal:"<<endl;
    for(int k=0;k<normalVec.size();k++)
      logFile<<normalVec[k]<<" ";
    logFile<<" delta areas/weights:"<<endl;
    logFile<<area<<" / "<<weight<<endl;
  }
  logFile<<"System area/weight = "<<sumElemAreas<<" / "
	 <<sumElemWeights<<endl;
#endif
 
  /*********************************************************************/ 
  // Assemble all boundary gauss points weights and normal vectors 
  // from all processors and store them.

  // Determine for all boundary integration points which processor has 
  // calculated them.
  intVector rootList(surfaceGaussPtsNum);

  if(surfaceGaussPtsNum > 0)
    MPI_Allreduce(&localRootList[0],&rootList[0],surfaceGaussPtsNum,
		  MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  // Loop over all boundary integration points.
  for(int i=0;i<surfaceGaussPtsNum;i++) {
    int& currentPoint = surfaceGaussPtsIdx[i][0];

    // weight
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&boundGaussPoints[currentPoint].getWeight(),1,MPI_DOUBLE,
	      rootList[i],MPI_COMM_WORLD);

    // surface normal
    dbVector& surfaceNormal = 
      boundGaussPoints[currentPoint].getSurfaceNormal();

    if(surfaceNormal.size() < usedDims)
      surfaceNormal.resize(usedDims);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&surfaceNormal[0],3,MPI_DOUBLE,rootList[i],MPI_COMM_WORLD);

  }
  
  /**********************************************************************/
  // plot computed surface normals

  if((bool)InputData->getValue("plotSurfaceNormalDistribution")) {

    ofstream normalMeshFile;
    normalMeshFile.open("normal.msh");
    normalMeshFile<<"MESH  dimension 3  ElemType Point Nnode 1"<<endl;
    normalMeshFile<<"Coordinates"<<endl;
    for(int i=0;i<surfaceGaussPtsNum;i++) {
      int& currentPoint = surfaceGaussPtsIdx[i][0];
      dbVector& coords = boundGaussPoints[currentPoint].getCoords();
      normalMeshFile<<i+1<<" ";
      for(int j=0;j<coords.size();j++)
	normalMeshFile<<coords[j]<<" ";
      normalMeshFile<<endl;
    }
    normalMeshFile<<"end coordinates"<<endl;
    normalMeshFile<<endl;
    normalMeshFile<<"Elements"<<endl;
    for(int i=0;i<surfaceGaussPtsNum;i++)
      normalMeshFile<<i+1<<" "<<i+1<<"  1"<<endl;
    normalMeshFile<<"end elements"<<endl;
    std::ofstream normalResFile;
    normalResFile.open("normal.res");
    normalResFile.precision(12);
    normalResFile.setf(ios_base::scientific,ios_base::floatfield);
    normalResFile<<"GiD Post Results File 1.0"<<endl;
    normalResFile<<"Result \"normal_vector\" \" \" 1.0"
		 <<"  Vector  OnNodes"<<endl;
    normalResFile<<"ComponentNames ";
    for(int i=0;i<3;i++)
      normalResFile<<"\"DOF_"<<i+1<<"\", ";
    normalResFile<<endl;
    normalResFile<<"Values"<<endl;

    for(int i=0;i<surfaceGaussPtsNum;i++) {
      int& currentPoint = surfaceGaussPtsIdx[i][0];
      dbVector& surfaceNormal = 
	boundGaussPoints[currentPoint].getSurfaceNormal();
      normalResFile<<i+1<<" ";
      for(int j=0;j<surfaceNormal.size();j++)
	normalResFile<<surfaceNormal[j]<<" ";
      normalResFile<<endl;
      normalResFile.flush();
    }
    normalResFile<<"End Values"<<endl;
  }

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all line gauss points their weight.
void FEMGeometry::setLineGaussPtsWeights(InputFileData* InputData,
					 std::map<std::string,double>& modelData,
					 std::ofstream& logFile) {

  using namespace std;

  // Determine the a local portion of boundary elements for calculating.
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  allLineBoundGaussPtsIdx =
    getAllLineGaussPtsIdx(logFile);

  intMatrix& lineGaussPtsIdx = allLineBoundGaussPtsIdx;

  int lineGaussPtsNum = lineGaussPtsIdx.size();

  int startIdx,endIdx;
  int gaussPortion = (int)ceil((double)lineGaussPtsNum/size);

  if(gaussPortion*rank < lineGaussPtsNum && 
     gaussPortion*(rank+1) <= lineGaussPtsNum) {
    startIdx = gaussPortion*rank;
    endIdx = gaussPortion*(rank+1);
  }
  else if(gaussPortion*rank <= lineGaussPtsNum && 
	  gaussPortion*(rank+1) >= lineGaussPtsNum) {
    startIdx = gaussPortion*rank;
    endIdx = lineGaussPtsNum;
  }
  else {
    startIdx = endIdx = lineGaussPtsNum;
  }
  
  gaussPortion = endIdx-startIdx;

  double length;
  intVector localRootList(lineGaussPtsNum);

  // loop over a portion of line boundary Gauss points
  for(int i=startIdx;i<endIdx;i++) {
    int& currentPoint = lineGaussPtsIdx[i][0];
    localRootList[i] = rank;

    // Calculate surface normal vector and its length - the surface 
    // of the volume element.
    double& weight = boundGaussPoints[currentPoint].getWeight();
    intVector& elemInfo = boundGaussPoints[currentPoint].getElementInfo();

    int& elemID = elemInfo[0]; // line element ID
    int& gaussID = elemInfo[1]; // Gauss point ID within the element

    FEMElement& elem = lineNodesElems[elemID];

    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* LineGaussSet = elem.getLineGaussSet();
    intVector& nodes = elem.getNodes();
 
    FEVolumeSet->lineMetric(particles,nodes,
			    LineGaussSet->coord[gaussID],
			    length,logFile);

    weight = fabs(length)*LineGaussSet->weight[gaussID];

  }

#ifdef _FEdebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"************** line element weights  *****************"<<endl;
  double weight;
  double sumElemLengths = 0;
  double sumElemWeights = 0;
  for(int i=startIdx;i<endIdx;i++) {
    int& currentPoint = lineGaussPtsIdx[i][0];
    intVector& elemInfo = boundGaussPoints[currentPoint].getElementInfo();
    int& elemID = elemInfo[0];
    int& gaussID = elemInfo[1];
    FEMElement& elem = lineNodesElems[elemID];
    intVector& nodes = elem.getNodes();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* LineGaussSet = elem.getLineGaussSet();
    logFile<<"Element "<<elemID<<" nodes: ";
    for(int j=0;j<nodes.size();j++)
      logFile<<nodes[j]<<" ";
    logFile<<endl;
    logFile<<"delta lengths/weights:"<<endl;
    FEVolumeSet->lineMetric(particles,nodes,
			    LineGaussSet->coord[gaussID],
			    length,logFile);
    weight = fabs(length)*LineGaussSet->weight[gaussID];
    logFile<<length<<" / "<<weight<<endl;
    sumElemLengths += length;
    sumElemWeights += weight;
  }
  logFile<<"System length/weight = "<<sumElemLengths<<" / "
	 <<sumElemWeights<<endl;
#endif

  /*********************************************************************/ 
  // Assemble all boundary gauss points weights and normal vectors 
  // from all processors and store them.

  // Determine for all boundary integration points which processor has 
  // calculated them.
  intVector rootList(lineGaussPtsNum);

  if(lineGaussPtsNum > 0)
    MPI_Allreduce(&localRootList[0],&rootList[0],lineGaussPtsNum,
		  MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  // Loop over all line deformation boundary integration points.
  double lineWidth = InputData->getValue("lineWidth");

  for(int i=0;i<lineGaussPtsIdx.size();i++) {
    int& currentPoint = lineGaussPtsIdx[i][0];
    double& weight = boundGaussPoints[currentPoint].getWeight();

    if(rootList[i] == rank) 
      weight *= lineWidth;
  }

  // Loop over all boundary integration points.
  for(int i=0;i<lineGaussPtsNum;i++) {
    int& currentPoint = lineGaussPtsIdx[i][0];

    // weight
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&boundGaussPoints[currentPoint].getWeight(),1,MPI_DOUBLE,
	      rootList[i],MPI_COMM_WORLD);

  }

#ifdef _FEdebugMode_
  logFile<<"************** local root list  *****************"<<endl;
  for(int i=0;i<localRootList.size();i++)
    logFile<<"LINE GAUSSPOINT "<<lineGaussPtsIdx[i][0]<<" : "
	   <<localRootList[i]<<endl;
  logFile<<"************** global root list  *****************"<<endl;
  for(int i=0;i<rootList.size();i++)
    logFile<<"LINE GAUSSPOINT "<<lineGaussPtsIdx[i][0]<<" : "
	   <<rootList[i]<<endl;
  logFile<<"************** final weights  *****************"<<endl;
  for(int i=startIdx;i<endIdx;i++) {
    int& currentPoint = lineGaussPtsIdx[i][0];
    double& weight = boundGaussPoints[currentPoint].getWeight();
    logFile<<"LINE GAUSSPOINT "<<lineGaussPtsIdx[i][0]<<" : "
	   <<weight<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Calculate for all particles their particle weight
// (needed for Nystroem Integration and window function norming).
void FEMGeometry::setAllPtclsWeights(InputFileData* InputData,
				     std::ofstream& logFile) {
  
  using namespace std;

#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
#endif

  int ptcle;
  int nodesNum = particles.size();

  double weight;
  dbVector localPtcleWeights(nodesNum);
  
  // Loop over all local elements.
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];

#ifdef _FEdebugMode_
    logFile<<"----------------------------------------------------"<<endl; 
    logFile<<"ELEMENT "<<i<<endl;
#endif

    int& globalID = elem.getGlobalID();
    int& elemType = elem.getElemType();
    intVector& nodes = elem.getNodes();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    
    for(int j=0;j<nodes.size();j++) {
      ptcle = nodes[j]-1;

      // tetrahedral
      if(elemType == 1)

	weight = 
	  FEVolumeSet->metricFactorVolume(particles,nodes,
					  FEVolumeSet->nodalCoords[j],
					  logFile);

	

      // hexahedral
      else if(elemType == 2) 


	weight = 
	  FEVolumeSet->nodalWeights[j]*
	  FEVolumeSet->metricFactorVolume(particles,nodes,
					  FEVolumeSet->nodalCoords[j],
					  logFile);


      
      else {
	logFile<<"In FEMGeometry::setAllPtclsWeights element type "
	       <<elemType<<"is not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      localPtcleWeights[ptcle] += weight;
      particles[ptcle].setElems(globalID);

#ifdef _FEdebugMode_
      logFile<<"PARTICLE "<<ptcle<<" weight portion "<<weight<<endl;
#endif

    }

  }

#ifdef _FEdebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"************* local particle weights ****************"<<endl;
  for(int i=0;i<nodesNum;i++) {
    intVector& connectedElems = particles[i].getElems();
    logFile<<i<<" "<<localPtcleWeights[i]<<"; connected elems: ";
    for(int j=0;j<connectedElems.size();j++)
      logFile<<connectedElems[j]<<" ";
    logFile<<endl;
  }
#endif

  /*********************************************************************/
  // Assemble the particle weights of all processors.
  dbVector globalWeights(nodesNum);
   
  MPI_Allreduce(&localPtcleWeights[0],&globalWeights[0],nodesNum,
		MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int sendBuf; 
  intVector recvCounts(size);
  intVector displs(size);

  int globalElemSize,m;
  intVector globalElems;

  for(int i=0;i<nodesNum;i++) {

    // particle weights
    dbVector& ptcleWeights = particles[i].getAllIntWeights();
    ptcleWeights.resize(1);
    ptcleWeights[0] = globalWeights[i];

    // weight function integral (not implemented yet -> so use particle weight instead
    //double& windowFuncIntegral = particles[i].getWindowFuncIntegral();
    //windowFuncIntegral = ptcleWeights[0];

    // particle element connectivity lists
    intVector& localElems = particles[i].getElems();
    sendBuf = localElems.size();
    MPI_Allgather(&sendBuf,1,MPI_INT,&recvCounts[0],1,
		  MPI_INT,MPI_COMM_WORLD);

    displs[0] = 0;

    for(int j=1;j<size;j++)
      displs[j] = displs[j-1]+recvCounts[j-1];

    globalElemSize = 0;

    for(int j=0;j<size;j++)
      globalElemSize += recvCounts[j];


    if(globalElemSize == 0) {
      logFile<<"In FEMGeometry::setAllPtclsWeights particle "<<i
	     <<" is not member of any element. Check your mesh!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    if(globalElems.size() < globalElemSize)
      globalElems.resize(globalElemSize);

    MPI_Allgatherv(&localElems[0],recvCounts[rank],MPI_INT,
		   &globalElems[0],&recvCounts[0],&displs[0],MPI_INT,
		   MPI_COMM_WORLD);

    // sort and store the complete list for the current particle
    sortIntVector(globalElems,0,globalElemSize-1);

    localElems.resize(globalElemSize);

    m=1;
    localElems[0] = globalElems[0];

    for(int j=1;j<globalElemSize;j++)

      // skip already stored elements
      if(globalElems[j-1] != globalElems[j]) {

	localElems[m] = globalElems[j];
	m++;

      }

  }

#ifdef _FEdebugMode_
  double systemVolume = 0;
  logFile<<"*****************************************************"<<endl;
  logFile<<"****************** particle weights *****************"<<endl;
  for(int i=0;i<nodesNum;i++) {
    intVector& connectedElems = particles[i].getElems();
    dbVector& ptcleWeights = particles[i].getAllIntWeights();
    systemVolume += ptcleWeights[0];
    logFile<<i<<" "<<ptcleWeights[0]<<"; connected elems: ";
    for(int j=0;j<connectedElems.size();j++)
      logFile<<connectedElems[j]<<" ";
    logFile<<endl;
  }
  logFile<<"System Volum = "<<systemVolume<<endl;
#endif

}

/********************************************************************/
/********************************************************************/
// Calculate for all particles their particle weight 
// (Nystroem Integration) and surface and volume loads.
void FEMGeometry::setPtcleWeightsLoads(InputFileData* InputData,
				       std::map<std::string,double>& modelData,
				       std::ofstream& logFile,
				       PetscViewer& viewerMPI) {
  
  using namespace std;
  
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int usedDims = (int)modelData["usedDimensions"];
  
#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"******* calculation of particle weights & loads *****"<<endl;
#endif
  
  int integrationMethod = (int)modelData["integrationMethod"];
  
  // number of which conditions get their own integration weight 
  // (entire body integration,body force)
  int numOfEntries = 2;
  
  int currentNode,position;
  double weight;
  
  intVector loadPtcleIdx;
  
  // Loop over all local elements.
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];
  
    intVector& nodes = elem.getNodes();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    
    // Loop over the element's nodes.
    for(int j=0;j<nodes.size();j++) {
      currentNode = nodes[j]-1;
      
      weight = 
	FEVolumeSet->metricFactorVolume(particles,nodes,
					FEVolumeSet->nodalCoords[j],
					logFile);
      
      dbVector& ptcleWeights =  particles[currentNode].getAllIntWeights();
      ptcleWeights.resize(numOfEntries);
      
      // weight for the entire body integration
      ptcleWeights[0] += weight;
    
      // ----------------------------------------------------------------
      // weight for the body force only
      position = findIntMatPos(i,0,bodyForceElemIdx.size(),0,
			       bodyForceElemIdx);
      
      // On current element is a body force applied
      if(position != -1) {
      
	blVector& elemDOF = elem.getBodyForceDOF();
	dbVector& elemConds = elem.getBodyForce();
	
	position = findIntVecPos(currentNode,0,loadPtcleIdx.size(),
				 loadPtcleIdx);
	
	blMatrix& affectedDOF = particles[currentNode].getAllBodyForceDOF();
	dbMatrix& conditions = particles[currentNode].getAllBodyForceLoads();
	
	// for current node is still not a body force applied
	if(position == -1) {
	
	  // Store the the body force.  
	  allocateArray(affectedDOF,1);
	  allocateArray(conditions,1);
	  
	  affectedDOF[0] = elemDOF;
	  conditions[0] = elemConds;
	  
	  // Store the corresponding integration weight.
	  ptcleWeights[1] = weight;
	  
	  loadPtcleIdx.push_back(currentNode);
	}
      
	// for current node is already a body force applied.
	else {
	      
	  for(int l=0;l<elemConds.size();l++) {
	    
	    if((elemDOF[l] != affectedDOF[0][l]) || 
	       (conditions[0][l] != elemConds[l])) {
	        
	      logFile<<"In FEMGeometry::setPtcleWeightsLoads non-"
		     <<"constant body force loads are not supported!"
		     <<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	     
	    }
	    
	  }
	  
	  ptcleWeights[1] += weight;
	  
	}
	
      }
      // ----------------------------------------------------------------
      // body moments
      if(bodyMomentElemIdx.size() > 0) {
	logFile<<"In function FEMGeometry::setPtcleWeightsLoads "
	       <<"body moments are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
    }
    
  }
  
#ifdef _FEdebugMode_
  logFile<<"************ local particle weights ******************"<<endl;
  for(int i=0;i<particles.size();i++) {
    dbVector& ptcleWeights = particles[i].getAllIntWeights();
    logFile<<"PARTICLE "<<i<<": "<<ptcleWeights[0]<<endl;
  }
  logFile<<"********* local body loading particle indices ********"<<endl;
  for(int i=0;i<loadPtcleIdx.size();i++)
    logFile<<i<<" -> particle "<<loadPtcleIdx[i]<<endl;
  logFile<<"*********** local body forces on particles ***********"<<endl;
  for(int i=0;i<loadPtcleIdx.size();i++) {
    logFile<<"PARTICLE "<<loadPtcleIdx[i]<<": ";
    blMatrix& affectedDOF = 
      particles[loadPtcleIdx[i]].getAllBodyForceDOF();
    dbMatrix& loads = particles[loadPtcleIdx[i]].getAllBodyForceLoads();
    for(int j=0;j<loads.size();j++) {
      logFile<<"body force ID "<<j<<" *********"<<endl;
      for(int k=0;k<loads[j].size();k++)
	logFile<<"affected: "<<affectedDOF[j][k]<<"load: "
	       <<loads[j][k]<<" ";
      logFile<<endl;
    }
  }
#endif
  
  /*********************************************************************/
  // Assemble the entire body integration particle weights of all 
  // processors.
  int nodesNum = particles.size();
  dbVector localPtcleWeights(nodesNum);
  dbVector globalPtcleWeights(nodesNum);
  
  for(int i=0;i<nodesNum;i++) {
    dbVector& ptcleWeights = particles[i].getAllIntWeights();
    localPtcleWeights[i] = ptcleWeights[0];
  }
  
  MPI_Allreduce(&localPtcleWeights[0],&globalPtcleWeights[0],nodesNum,
		MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  for(int i=0;i<nodesNum;i++) {
    dbVector& ptcleWeights = particles[i].getAllIntWeights();
    ptcleWeights[0] = globalPtcleWeights[i];
  }

#ifdef _FEdebugMode_
  double systemVolume = 0;
  logFile<<"****** entire body integration particle weights ******"<<endl;
  for(int i=0;i<nodesNum;i++) {
    systemVolume += particles[i].getIntWeight(0);
    logFile<<"PARTICLE "<<i<<": "<<particles[i].getIntWeight(0)<<endl;
  }
  logFile<<"System Volum = "<<systemVolume<<endl;
#endif

  //---------------------------------------------------------------------
  // Assemble body force integration particle weights of all 
  // processors.
  for(int i=0;i<nodesNum;i++) {
    dbVector& ptcleWeights = particles[i].getAllIntWeights();
    localPtcleWeights[i] = ptcleWeights[1];
  }
  
  MPI_Allreduce(&localPtcleWeights[0],&globalPtcleWeights[0],nodesNum,
		MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  intVector inputDOF(usedDOF);
  intVector outputDOF(usedDOF);
  dbVector localLoads;
  bodyForcePtcleIdx.resize(nodesNum,intVector(numOfEntries));
  
  int m = 0;
  
  for(int i=0;i<nodesNum;i++) {
    
    if(globalPtcleWeights[i] > 0) {
      
      // body force weights
      dbVector& ptcleWeights = particles[i].getAllIntWeights();
      ptcleWeights[1] = globalPtcleWeights[i];
      
      // list of particles carrying a body force loading
      bodyForcePtcleIdx[m][0] = i; // particle ID
      bodyForcePtcleIdx[m][1] = 1; // weight ID
      
      // particle's DOF with a body force loading
      blMatrix& globalDOF = particles[i].getAllBodyForceDOF();

      if(globalDOF.size() < 1) {
	globalDOF.resize(1);

	for(int j=0;j<globalDOF.size();j++)
	  globalDOF[j].resize(usedDOF);
	
      }
      
      for(int j=0;j<globalDOF[0].size();j++)
	inputDOF[j] = (int)globalDOF[0][j];

      MPI_Allreduce(&inputDOF[0],&outputDOF[0],usedDOF,MPI_INT,
      		    MPI_MAX,MPI_COMM_WORLD);
      
      for(int j=0;j<globalDOF[0].size();j++)
	globalDOF[0][j] = (bool)outputDOF[j];

      // particle's body force loading
      dbMatrix& globalLoads = particles[i].getAllBodyForceLoads();
      
      if(globalLoads.size() < 1) {
	globalLoads.resize(1);
	
	for(int j=0;j<globalLoads.size();j++)
	  globalLoads[j].resize(usedDOF);
	
      }

      localLoads = globalLoads[0];

      MPI_Allreduce(&localLoads[0],&globalLoads[0][0],usedDOF,MPI_DOUBLE,
		    MPI_SUM,MPI_COMM_WORLD);
      m++;
    }
    
  }
  
#ifdef _FEdebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"********** global body forces on particles ***********"<<endl;
  for(int i=0;i<bodyForcePtcleIdx.size();i++) {
    int& ptcle = bodyForcePtcleIdx[i][0];
    dbVector& ptcleWeights = particles[ptcle].getAllIntWeights();
    blMatrix& globalDOF = particles[ptcle].getAllBodyForceDOF();
    dbMatrix& globalLoads = particles[ptcle].getAllBodyForceLoads();
    logFile<<i<<".) PARTICLE "<<bodyForcePtcleIdx[i][0]<<" weight ID "
	   <<bodyForcePtcleIdx[i][1]<<" | ";
    for(int j=0;j<globalLoads.size();j++) {
      for(int k=0;k<globalLoads[j].size();k++) {
	if(globalDOF[j][k])
	  logFile<<globalLoads[j][k]<<" ";
      }
    }
    logFile<<"; weight = "<<ptcleWeights[1]<<endl;
  }
#endif
  
}

/**********************************************************************/
/**********************************************************************/
// Calculate the weights for boundary particles and store loads if
// applied.
  void FEMGeometry::setBoundPtcleWeightsLoads(InputFileData* InputData,
					    std::map<std::string,double>& modelData,
					    std::ofstream& logFile,
					    PetscViewer& viewerMPI) {

  using namespace std;

  // Determine the integration weights for these particles at which 
  // a line boundary has to be integrated and store the line loads if any 
  // are applied.
  setLinePtcleWeightsLoads(InputData,modelData,logFile,viewerMPI);

  // Determine weights and surface normals for these particles at which 
  // the boundary has to be integrated and store the surface loads if any 
  // are applied.
  setSurfacePtcleWeightsLoads(InputData,modelData,logFile);

  // Store the point integration points coordinate's, any deformation 
  // boundary condition and any point loads if applied.
  setBoundPointIntegration(InputData,modelData,logFile,viewerMPI);
}

/**********************************************************************/
/**********************************************************************/
// Store the boundary point integration points essential boundary
// conditions and any point forces if applied.
void FEMGeometry::setBoundPointIntegration(InputFileData* InputData,
					   std::map<std::string,double>& modelData,
					   std::ofstream& logFile,
					   PetscViewer& viewerMPI) {
  
  using namespace std;
  
  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];
  int dispDOF = (int)modelData["displacementDegreesOfFreedom"];
  dbMatrix& pointForces = InputData->getPointForceLoads();
  dbMatrix& pointDispBounds = InputData->getPointDispBoundConds();
  
  double pointWidth = InputData->getValue("pointWidth");

  bool conditionSet,normalSet;
  int idx,position,ptcle;
  double ptcleWeight;
  intVector indexSet(4);
  dbVector normal(3);

  intVector boundPtcleIdx;
  intMatrix storagePositions(pointForces.size()+pointDispBounds.size());

  // Loop over all point loads.
  for(int i=0;i<pointForces.size();i++) {
    
    ptcle = newIdx[(int)pointForces[i][0]-1];
  
    dbVector& weights = particles[ptcle].getAllIntWeights();
    dbMatrix& conditions = particles[ptcle].getAllPointForceLoads();
    blMatrix& affectedDOF = particles[ptcle].getAllPointForceDOF();
    dbMatrix& surfaceNormals = particles[ptcle].getAllSurfaceNormals();
    
    position = findIntVecPos(ptcle,0,boundPtcleIdx.size(),boundPtcleIdx);

    // no bound point force applied yet
    if(position == -1) {

      // set identifier
      indexSet[0] = ptcle;

      // set weight
      ptcleWeight = pointWidth;
      pushBackVector(weights,ptcleWeight);
      indexSet[1] = 0;

      // set loads    
      conditions.resize(1);
      conditions[0].resize(usedDOF);
      affectedDOF.resize(1);
      affectedDOF[0].resize(usedDOF);
      
      for(int k=4;k<pointForces[i].size();k+=2) {
	affectedDOF[0][(int)pointForces[i][k]-1] = true;
	conditions[0][(int)pointForces[i][k]-1] = pointForces[i][k+1];
      }

      indexSet[2] = 0;

      // set surface normal
      normal[0] = pointForces[i][1];
      normal[1] = pointForces[i][2];
      normal[2] = pointForces[i][3];
      indexSet[3] = 0;

      pushBackVector(surfaceNormals,normal);

      pushBackVector(pointForceBoundPtcleIdx,indexSet);
      boundPtcleIdx.push_back(ptcle);
    }

    // current particle has got point forces already
    else {

      conditionSet = false;

      // check the already applied conditions
      for(int j=0;j<storagePositions[position].size();j++) {
	idx = storagePositions[position][j];

	int& weightID = pointForceBoundPtcleIdx[idx][1];
	int& conditionID = pointForceBoundPtcleIdx[idx][2];
	int& surfaceID = pointForceBoundPtcleIdx[idx][3];
      
	// particle has got point force data for current surface already
	if(surfaceNormals[surfaceID][0] == pointForces[i][0] && 
	   surfaceNormals[surfaceID][1] == pointForces[i][1] && 
	   surfaceNormals[surfaceID][2] == pointForces[i][2]) {
	  
	  // add the new conditions
	  for(int k=4;k<pointForces[i].size();k+=2) {
	    affectedDOF[conditionID][(int)pointForces[i][k]-1] = true;
	    conditions[conditionID][(int)pointForces[i][k]-1] += pointForces[i][k+1];
	  }

	  conditionSet = true;
	}

      }
	
      // current particle has not got point force data for current surface 
      // already
      if(!conditionSet) {

	// set identifier
	indexSet[0] = ptcle;

	// set loads    
	idx = conditions.size();
	conditions.resize(idx+1);
	conditions[idx].resize(usedDOF);
	affectedDOF.resize(idx+1);
	affectedDOF[idx].resize(usedDOF);
      
	for(int k=4;k<pointForces[i].size();k+=2) {
	  affectedDOF[idx][(int)pointForces[i][k]-1] = true;
	  conditions[idx][(int)pointForces[i][k]-1] = pointForces[i][k+1];
	}

	indexSet[2] = idx;
      
	// set surface normal and the integration weight
	
	// check first if the normal vector is already set
	normalSet = false;

	for(int j=0;j<surfaceNormals.size();j++)

	  if(surfaceNormals[j][0] == pointForces[i][0] && 
	     surfaceNormals[j][1] == pointForces[i][1] && 
	     surfaceNormals[j][2] == pointForces[i][2]) {

	    indexSet[1] = j;
	    indexSet[3] = j;
	    normalSet = true;
	  }
        
	if(!normalSet) {
	  normal[0] = pointForces[i][1];
	  normal[1] = pointForces[i][2];
	  normal[2] = pointForces[i][3];

	  indexSet[1] = surfaceNormals.size();
	  indexSet[3] = surfaceNormals.size();

	  pushBackVector(surfaceNormals,normal);

	  ptcleWeight = 1.0;
	  pushBackVector(weights,ptcleWeight);
	}

	pushBackVector(pointForceBoundPtcleIdx,indexSet);
	storagePositions[position].
	  push_back(pointForceBoundPtcleIdx.size());
      }

    }

  }

  /*********************************************************************/
  // Store all point moment loads.
  dbMatrix& pointMoments = InputData->getPointMomentLoads();
 
  if(pointMoments.size() > 0) {
    logFile<<"In function FEMGeometry::setBoundPointIntegration "
	   <<"point moments are not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /*********************************************************************/
  // Store all point deformation boundary conditions.

  for(int i=0;i<storagePositions.size();i++)
    storagePositions[i].resize(0);

  // Loop over all displacement point boundary conditions.
  for(int i=0;i<pointDispBounds.size();i++) {

    if(dispDOF <= 0) {
      logFile<<"In function FEMGeometry::setBoundPointIntegration "
	     <<"displacement boundary conditions are not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    ptcle = newIdx[(int)pointDispBounds[i][0]-1];

    dbVector& weights = particles[ptcle].getAllIntWeights();
    dbMatrix& conditions = particles[ptcle].getAllDeformationBoundConds();
    blMatrix& affectedDOF = particles[ptcle].getAllDeformationBoundDOF();
    dbMatrix& surfaceNormals = particles[ptcle].getAllSurfaceNormals();

    position = findIntVecPos(ptcle,0,boundPtcleIdx.size(),boundPtcleIdx);

    // no boundary particle integration data applied yet
    if(position == -1) {

      // set identifier
      indexSet[0] = ptcle;

      // set weight
      ptcleWeight = pointWidth;
      pushBackVector(weights,ptcleWeight);
      indexSet[1] = 0;

      // set essential boundary conditions  
      conditions.resize(1);
      conditions[0].resize(defDOF);
      affectedDOF.resize(1);
      affectedDOF[0].resize(defDOF);
      
      for(int k=4;k<pointDispBounds[i].size();k+=2) {
	affectedDOF[0][(int)pointDispBounds[i][k]-1] = true;
	conditions[0][(int)pointDispBounds[i][k]-1] = pointDispBounds[i][k+1];
      }

      indexSet[2] = 0;

      // set surface normal
      normal[0] = pointDispBounds[i][1];
      normal[1] = pointDispBounds[i][2];
      normal[2] = pointDispBounds[i][3];
      indexSet[3] = 0;

      pushBackVector(surfaceNormals,normal);

      pushBackVector(pointDispBoundPtcleIdx,indexSet);
      boundPtcleIdx.push_back(ptcle);
    }


    // current particle has got particle integration data already
    else {

      conditionSet = false;

      // check the already applied conditions
      for(int j=0;j<storagePositions[position].size();j++) {
	idx = storagePositions[position][j];

	int& weightID = pointDispBoundPtcleIdx[idx][1];
	int& conditionID = pointDispBoundPtcleIdx[idx][2];
	int& surfaceID = pointDispBoundPtcleIdx[idx][3];
      
	// particle has got particle integration data for current surface 
	// already
	if(surfaceNormals[surfaceID][0] == pointDispBounds[i][0] && 
	   surfaceNormals[surfaceID][1] == pointDispBounds[i][1] && 
	   surfaceNormals[surfaceID][2] == pointDispBounds[i][2]) {
	  
	  // add the new conditions
	  for(int k=4;k<pointDispBounds[i].size();k+=2) {
	    idx = (int)pointDispBounds[i][k]-1;

	    // particle carries for a DOF no conditions yet
	    if(!affectedDOF[conditionID][idx]) {
	      affectedDOF[conditionID][idx] = true;
	      conditions[conditionID][idx] = pointDispBounds[i][k+1];
	    }
	    // particle carries for a DOF different conditions which is not admissible
	    else if(affectedDOF[conditionID][idx] &&
		    conditions[conditionID][idx] != pointDispBounds[i][k+1]) {
	      logFile<<"In FEMGeometry::setBoundPointIntegration non-constant\n"
		     <<"conditions are not supported!"
		     << endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	    }
	    
	  }

	  conditionSet = true;
	}

      }
	
      // current particle has not got any point integration data for 
      // current surface yet
      if(!conditionSet) {

	// set identifier
	indexSet[0] = ptcle;

	// set weight
	indexSet[1] = weights.size(); // weight ID
	ptcleWeight = pointWidth;
	pushBackVector(weights,ptcleWeight);

	// set conditions  
	idx = conditions.size();
	conditions.resize(idx+1);
	conditions[idx].resize(defDOF);
	affectedDOF.resize(idx+1);
	affectedDOF[idx].resize(defDOF);
      
	for(int k=4;k<pointDispBounds[i].size();k+=2) {
	  affectedDOF[idx][(int)pointDispBounds[i][k]-1] = true;
	  conditions[idx][(int)pointDispBounds[i][k]-1] = pointDispBounds[i][k+1];
	}

	indexSet[2] = idx; // condition ID

	// set surface normal and the integration weight
	
	// check first if the normal vector is already set
	normalSet = false;

	for(int j=0;j<surfaceNormals.size();j++)

	  if(surfaceNormals[j][0] == pointDispBounds[i][0] && 
	     surfaceNormals[j][1] == pointDispBounds[i][1] && 
	     surfaceNormals[j][2] == pointDispBounds[i][2]) {

	    indexSet[3] = j; // surface normal ID;
	    normalSet = true;
	  }
        
	if(!normalSet) {
	  normal[0] = pointDispBounds[i][1];
	  normal[1] = pointDispBounds[i][2];
	  normal[2] = pointDispBounds[i][3];

	  indexSet[3] = surfaceNormals.size(); // surface normal ID

	  pushBackVector(surfaceNormals,normal);
	}

	pushBackVector(pointDispBoundPtcleIdx,indexSet);
	storagePositions[position].
	  push_back(pointDispBoundPtcleIdx.size());
      }

    }

  }

  /*********************************************************************/
  // Store all point rotation boundary constraints.
  dbMatrix& pointRotBounds = InputData->getPointRotBoundConds();
 
  if(pointRotBounds.size() > 0) {
    logFile<<"In function FEMGeometry::setBoundPointIntegration "
	   <<"point rotation constraints are not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"******************** point forces *******************"<<endl;
  for(int i=0;i<pointForceBoundPtcleIdx.size();i++) {
    int& ptcle = pointForceBoundPtcleIdx[i][0];
    int& weightID = pointForceBoundPtcleIdx[i][1];
    int& loadID = pointForceBoundPtcleIdx[i][2];
    int& surfaceID = pointForceBoundPtcleIdx[i][3];
    double& weight = particles[ptcle].getIntWeight(weightID);
    dbVector& pointForce = particles[ptcle].getPointForce(loadID);
    blVector& affectedDOF = particles[ptcle].getPointForceDOF(loadID);
    dbVector& surfaceNormal = particles[ptcle].getSurfaceNormal(surfaceID);
    logFile<<"PARTICLE "<<ptcle<<": coords: "
	   <<particles[ptcle].getCoord(0)<<" "
	   <<particles[ptcle].getCoord(1)<<" "
	   <<particles[ptcle].getCoord(2)<<" ";
    logFile<<"normal: ";
    for(int j=0;j<surfaceNormal.size();j++)
      logFile<<surfaceNormal[j]<<" ";
    logFile<<"load: ";
    for(int j=0;j<pointForce.size();j++)
      if((bool)affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<pointForce[j]<<" ";
    logFile<<"weight "<<weight<<endl;
  }
  logFile<<"******************************************************"<<endl;
  logFile<<"******* point displacement boundary conditions *******"<<endl;
  for(int i=0;i<pointDispBoundPtcleIdx.size();i++) {
    int& ptcle = pointDispBoundPtcleIdx[i][0]; 
    int& weightID = pointDispBoundPtcleIdx[i][1]; 
    int& condID = pointDispBoundPtcleIdx[i][2];  
    int& surfaceID = pointDispBoundPtcleIdx[i][3];
    double& weight = particles[ptcle].getIntWeight(weightID);
    dbVector& boundConds = particles[ptcle].getDeformationBoundConds(condID);
    blVector& affectedDOF = particles[ptcle].getDeformationBoundDOF(condID);
    dbVector& surfaceNormal = particles[ptcle].getSurfaceNormal(surfaceID);
    logFile<<"PARTICLE "<<ptcle<<": coords: "
	   <<particles[ptcle].getCoord(0)<<" "
	   <<particles[ptcle].getCoord(1)<<" "
	   <<particles[ptcle].getCoord(2)<<" ";
    logFile<<"normal: ";
    for(int j=0;j<surfaceNormal.size();j++)
      logFile<<surfaceNormal[j]<<" ";
    logFile<<"bound cond: ";
    for(int j=0;j<boundConds.size();j++)
      if((bool)affectedDOF[j])
	logFile<<"DOF "<<j<<" condition "<<boundConds[j]<<" ";
    logFile<<"weight "<<weight<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine the integration weights for these particles at which 
// a line boundary has to be integrated and store the line loads if any 
// are applied.
void FEMGeometry::setLinePtcleWeightsLoads(InputFileData* InputData,
					   std::map<std::string,double>& modelData,
					   std::ofstream& logFile,
					   PetscViewer& viewerMPI) {

  using namespace std;

  // Compute all line Gauss points and store their line load.
  if(lineNodesElems.size() > 0) {
    logFile<<"In FEMGeometry::setLinePtcleWeightsLoads line boundary \n" 
	   <<"conditions are currently not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/**********************************************************************/
/**********************************************************************/
// Determine weights and surface normals for these particles at which 
// the boundary has to be integrated and store the surface loads if any 
// are applied.
void FEMGeometry::setSurfacePtcleWeightsLoads(InputFileData* InputData,
					      std::map<std::string,double>& modelData,
					      std::ofstream& logFile) {
  
  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int usedDims = (int)modelData["usedDimensions"];
  int defDOF = (int)modelData["deformationDegreesOfFreedom"];
  vector<FEMElement>& bElems = surfaceNodesElems;
  
  int currentNode,position;
  double surface;
  dbVector normal(3);

  intMatrix storagePositions;

  // get the point loading and point deformation boundary particles 
  intVector boundPtcleIdx(pointForceBoundPtcleIdx.size() +
			  pointDispBoundPtcleIdx.size());

  int m=0;
  
  for(int i=0;i<pointForceBoundPtcleIdx.size();i++) {
    boundPtcleIdx[m] = pointForceBoundPtcleIdx[i][0];
    m++;
  }

  for(int i=0;i<pointDispBoundPtcleIdx.size();i++) {
    boundPtcleIdx[m] = pointDispBoundPtcleIdx[i][0];
    m++;
  }

#ifdef _FEdebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"*********** FORCE BOUNDARY INTEGRAL PARTICLES ********"<<endl;
#endif

  // Store the surface loads and for the corresponding particles the
  // particle integration weight.

  /*********************************************************************/
  // traction loads

  // Determine the a local portion of traction boundary elements for 
  // calculating.
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  int startIdx,endIdx;
  int numOfCondElems = tractionElemIdx.size();
  int elemPortion = (int)ceil((double)numOfCondElems/size);
  
  if(elemPortion*rank < numOfCondElems && 
     elemPortion*(rank+1) <= numOfCondElems) {
    startIdx = elemPortion*rank;
    endIdx = elemPortion*(rank+1);
  }
  else if(elemPortion*rank <= numOfCondElems && 
	  elemPortion*(rank+1) >= numOfCondElems) {
    startIdx = elemPortion*rank;
    endIdx = numOfCondElems;
  }
  else {
    startIdx = endIdx = numOfCondElems;
  }
  
  elemPortion = endIdx-startIdx;  
  

  bool conditionSet,normalSet;
  int idx;
  intVector indexSet(4);

  // --------------------------------------------------------------------
  // Loop over the local portion of traction boundary elements.
  for(int i=startIdx;i<endIdx;i++) {
    FEMElement& elem = surfaceNodesElems[tractionElemIdx[i][0]];

    int& elemType = elem.getElemType();
    intVector& nodes = elem.getNodes();
    blVector& elemDOF = elem.getTractionDOF();
    dbVector& elemConds = elem.getTraction();

    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    ElementTemplate* FESurfaceSet = elem.getSurfaceElementTemplate();

    // Loop over the element's nodes.
    for(int j=0;j<nodes.size();j++) {
      
      currentNode = nodes[j]-1;
      
      dbMatrix& surfaceNormals = particles[currentNode].getAllSurfaceNormals();
      dbVector& weights = particles[currentNode].getAllIntWeights();
      dbMatrix& conditions = particles[currentNode].getAllTractionLoads();
      blMatrix& affectedDOF = particles[currentNode].getAllTractionDOF();
	
      // Calculate surface normal vector and its length, which is the surface 
      // of the volume element.
      
      FEVolumeSet->surfaceMetric(particles,nodes,
				 FESurfaceSet->nodalCoords[j],
				 normal,surface,logFile);      

      position = findIntVecPos(currentNode,0,boundPtcleIdx.size(),
			       boundPtcleIdx);
      
      // current particle has not got any loading for any surface yet
      if(position == -1) {

	// set identifier
	indexSet[0] = currentNode;

	// set weight
	pushBackVector(weights,surface);
	indexSet[1] = 0;
	
	// set the traction conditions
	allocateArray(affectedDOF,1);
	allocateArray(conditions,1);
	  
	affectedDOF[0] = elemDOF;
	conditions[0] = elemConds;
	
	indexSet[2] = 0;

	// set surface normal
	indexSet[3] = 0;

	pushBackVector(surfaceNormals,normal);


	pushBackVector(tractionBoundPtcleIdx,indexSet);
	boundPtcleIdx.push_back(currentNode);
	storagePositions.resize(boundPtcleIdx.size());

      }
	  
      // current particle has got particle integration data already
      else {
	
	conditionSet = false;

	// check the already applied conditions
	for(int j=0;j<storagePositions[position].size();j++) {
	  idx = storagePositions[position][j];

	  int& weightID = tractionBoundPtcleIdx[idx][1];
	  int& conditionID = tractionBoundPtcleIdx[idx][2];
	  int& surfaceID = tractionBoundPtcleIdx[idx][3];
      
	  // particle has got traction condition data for current surface 
	  // already
	  if(surfaceNormals[surfaceID][0] == normal[0] && 
	     surfaceNormals[surfaceID][1] == normal[1] && 
	     surfaceNormals[surfaceID][2] == normal[2]) {
	  
	    // check for non-constant conditions
	    for(int l=0;l<elemConds.size();l++) {

	      // particle carries for a DOF different conditions which is not admissible
	      if((affectedDOF[conditionID][l] != elemDOF[l]) || 
		 (conditions[conditionID][l] != elemConds[l])) {
		logFile<<"In FEMGeometry::setSurfacePtcleWeightsLoads non-constant\n"
		       <<"conditions are not supported!"
		       << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	      }

	    }

	    // add weight
	    weights[weightID] += surface;

	    conditionSet = true;
	  }

	}
	
	// current particle has not got any traction condition
	// data for current surface yet
	if(!conditionSet) {

	  // set identifier
	  indexSet[0] = currentNode;

	  // set weight
	  indexSet[1] = weights.size(); // weight ID
	  pushBackVector(weights,surface);

	  // set conditions  
	  idx = conditions.size();

	  affectedDOF.resize(idx+1);
	  conditions.resize(idx+1);

	  affectedDOF[idx] = elemDOF;
	  conditions[idx] = elemConds;

	  indexSet[2] = idx; // condition ID

	  // set surface normal

	  // check first if the normal vector is already set
	  normalSet = false;

	  for(int j=0;j<surfaceNormals.size();j++)

	    // other particle integration data than traction
	    // are already applied
	    if(surfaceNormals[j][0] == normal[0] && 
	       surfaceNormals[j][1] == normal[1] && 
	       surfaceNormals[j][2] == normal[2]) {

	      indexSet[3] = j; // surface normal ID;
	      normalSet = true;
	    }
        
	  if(!normalSet) {
	    indexSet[3] = surfaceNormals.size(); // surface normal ID

	    pushBackVector(surfaceNormals,normal);
	  }

	  storagePositions[position].
	    push_back(tractionBoundPtcleIdx.size());
	  pushBackVector(tractionBoundPtcleIdx,indexSet);
	}

      }

    }
    
  }

  /*********************************************************************/
  // surface pressure loads
  for(int i=0;i<storagePositions.size();i++)
    storagePositions[i].resize(0);

  // Determine the a local portion of surface pressure boundary elements 
  // for calculating.
  
  numOfCondElems = surfacePressureElemIdx.size();
  elemPortion = (int)ceil((double)numOfCondElems/size);
  
  if(elemPortion*rank < numOfCondElems && 
     elemPortion*(rank+1) <= numOfCondElems) {
    startIdx = elemPortion*rank;
    endIdx = elemPortion*(rank+1);
  }
  else if(elemPortion*rank <= numOfCondElems && 
	  elemPortion*(rank+1) >= numOfCondElems) {
    startIdx = elemPortion*rank;
    endIdx = numOfCondElems;
  }
  else {
    startIdx = endIdx = numOfCondElems;
  }
  
  elemPortion = endIdx-startIdx;  
  
  // Loop over the local portion of boundary elements.
  for(int i=startIdx;i<endIdx;i++) {
    FEMElement& elem = surfaceNodesElems[surfacePressureElemIdx[i][0]];

    int& elemType = elem.getElemType();
    intVector& nodes = elem.getNodes();
    double& pressure = elem.getSurfacePressure();

    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    ElementTemplate* FESurfaceSet = elem.getSurfaceElementTemplate();
    
    // Loop over the element's nodes.
    for(int j=0;j<nodes.size();j++) {
      currentNode = nodes[j]-1;
      
      dbMatrix& surfaceNormals = particles[currentNode].getAllSurfaceNormals();
      dbVector& weights = particles[currentNode].getAllIntWeights();
      dbVector& conditions = particles[currentNode].getAllSurfacePressureLoads();
	
      // Calculate surface normal vector and its length, which is the surface 
      // of the volume element.

      FEVolumeSet->surfaceMetric(particles,nodes,
				 FESurfaceSet->nodalCoords[j],
				 normal,surface,logFile);      

      position = findIntVecPos(currentNode,0,boundPtcleIdx.size(),
			       boundPtcleIdx);
      
      // current particle has not got any loading for any surface yet
      if(position == -1) {

	// set identifier
	indexSet[0] = currentNode;

	// set weight
	pushBackVector(weights,surface);
	indexSet[1] = 0;

	// set the surface pressure conditions
	pushBackVector(conditions,pressure);	
	
	indexSet[2] = 0;

	// set surface normal
	indexSet[3] = 0;

	pushBackVector(surfaceNormals,normal);
	
	pushBackVector(surfacePressureBoundPtcleIdx,indexSet);
	boundPtcleIdx.push_back(currentNode);
	storagePositions.resize(boundPtcleIdx.size());

      }
	  
      // current particle has got particle integration data already
      else {
	
	conditionSet = false;

	// check the already applied conditions
	for(int j=0;j<storagePositions[position].size();j++) {
	  idx = storagePositions[position][j];

	  int& weightID = surfacePressureBoundPtcleIdx[idx][1];
	  int& conditionID = surfacePressureBoundPtcleIdx[idx][2];
	  int& surfaceID = surfacePressureBoundPtcleIdx[idx][3];
      
	  // particle has got surface pressure data for current surface 
	  // already
	  if(surfaceNormals[surfaceID][0] == normal[0] && 
	     surfaceNormals[surfaceID][1] == normal[1] && 
	     surfaceNormals[surfaceID][2] == normal[2]) {

	    // particle carries for a DOF different conditions which is not admissible
	    if(conditions[conditionID] != pressure) {
	      logFile<<"In FEMGeometry::setSurfacePtcleWeightsLoads non-constant\n"
		     <<"conditions are not supported!"
		     << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	    }

	    // add weight
	    weights[weightID] += fabs(surface);

	    conditionSet = true;
	  }

	}
	
	// current particle has not got any surface pressure condition 
	// data for current surface yet
	if(!conditionSet) {

	  // set identifier
	  indexSet[0] = currentNode;

	  // set weight
	  indexSet[1] = weights.size(); // weight ID
	  pushBackVector(weights,surface);

	  // set the surface pressure conditions
	  indexSet[3] = conditions.size(); // condition ID

	  pushBackVector(conditions,pressure);	

	  // set surface normal and the integration weight
	  
	  // check first if the normal vector is already set
	  normalSet = false;

	  for(int j=0;j<surfaceNormals.size();j++)

	    // other particle integration data than traction
	    // are already applied
	    if(surfaceNormals[j][0] == normal[0] && 
	       surfaceNormals[j][1] == normal[1] && 
	       surfaceNormals[j][2] == normal[2]) {

	      indexSet[2] = j; // surface normal ID;
	      normalSet = true;
	    }
        
	  if(!normalSet) {
	    indexSet[2] = surfaceNormals.size(); // surface normal ID

	    pushBackVector(surfaceNormals,normal);
	  }

	  storagePositions[position].
	    push_back(surfacePressureBoundPtcleIdx.size());
	  pushBackVector(surfacePressureBoundPtcleIdx,indexSet);
	}

      }

    }
    
  }

  /*********************************************************************/
  // surface moment loads
  if(surfaceMomentElemIdx.size() > 0) {
    logFile<<"In function FEMGeometry::surfacePtcleWeightsLoads "
	   <<"surface moment loads are not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }


#ifdef _FEdebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"***** DISPLACEMENT BOUNDARY INTEGRAL PARTICLES *******"<<endl;
#endif

  for(int i=0;i<storagePositions.size();i++)
    storagePositions[i].resize(0);

  // Determine the a local portion of deformation boundary elements 
  // for calculating.  

  numOfCondElems = surfaceDispBoundElemIdx.size();
  elemPortion = (int)ceil((double)numOfCondElems/size);
  
  if(elemPortion*rank < numOfCondElems && 
     elemPortion*(rank+1) <= numOfCondElems) {
    startIdx = elemPortion*rank;
    endIdx = elemPortion*(rank+1);
  }
  else if(elemPortion*rank <= numOfCondElems && 
	  elemPortion*(rank+1) >= numOfCondElems) {
    startIdx = elemPortion*rank;
    endIdx = numOfCondElems;
  }
  else {
    startIdx = endIdx = numOfCondElems;
  }
  
  elemPortion = endIdx-startIdx;  

  // Loop over the local portion of boundary elements.
  for(int i=startIdx;i<endIdx;i++) {
    FEMElement& elem = surfaceNodesElems[surfaceDispBoundElemIdx[i][0]];

    int& elemType = elem.getElemType();
    intVector& nodes = elem.getNodes();
    blVector& elemDOF = elem.getSurfaceDeformationBoundDOF();
    dbVector& elemConds = elem.getSurfaceDeformationBoundConds();

    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    ElementTemplate* FESurfaceSet = elem.getSurfaceElementTemplate();
    
    // Loop over the element's nodes.
    for(int j=0;j<nodes.size();j++) {
      
      currentNode = nodes[j]-1;
      
      dbMatrix& surfaceNormals = particles[currentNode].getAllSurfaceNormals();
      dbVector& weights = particles[currentNode].getAllIntWeights();
      dbMatrix& conditions = particles[currentNode].getAllDeformationBoundConds();
      blMatrix& affectedDOF = particles[currentNode].getAllDeformationBoundDOF();
	
      // Calculate surface normal vector and its length, which is the surface 
      // of the volume element.

      FEVolumeSet->surfaceMetric(particles,nodes,
				 FESurfaceSet->nodalCoords[j],
				 normal,surface,logFile);      

      position = findIntVecPos(currentNode,0,boundPtcleIdx.size(),
			       boundPtcleIdx);
      
      // current particle has not got any loading for any surface yet
      if(position == -1) {

	// set identifier
	indexSet[0] = currentNode;

	// set weight
	pushBackVector(weights,surface);
	indexSet[1] = 0;

	// set the surface deformation boundary conditions
	allocateArray(affectedDOF,1);
	allocateArray(conditions,1);
	  
	affectedDOF[0] = elemDOF;
	conditions[0] = elemConds;
	
	indexSet[2] = 0;
	
	// set surface normal
	indexSet[3] = 0;

	pushBackVector(surfaceNormals,normal);

	pushBackVector(surfaceDispBoundPtcleIdx,indexSet);
	boundPtcleIdx.push_back(currentNode);
	storagePositions.resize(boundPtcleIdx.size());

      }
	  
      // current particle has got particle integration data already
      else {
	
	conditionSet = false;

	// check the already applied conditions
	for(int j=0;j<storagePositions[position].size();j++) {
	  idx = storagePositions[position][j];

	  int& weightID = surfaceDispBoundPtcleIdx[idx][1];
	  int& conditionID = surfaceDispBoundPtcleIdx[idx][2];
	  int& surfaceID = surfaceDispBoundPtcleIdx[idx][3];
      
	  // particle has got displacement boundary conditions for current surface 
	  // already
	  if(surfaceNormals[surfaceID][0] == normal[0] && 
	     surfaceNormals[surfaceID][1] == normal[1] && 
	     surfaceNormals[surfaceID][2] == normal[2]) {
	  
	    // Loop over all deformation DOF.
	    for(int l=0;l<elemDOF.size();l++) {

	      if((affectedDOF[conditionID][l] != elemDOF[l]) ||
		 (conditions[conditionID][l] != elemConds[l])) {
		logFile<<"In FEMGeometry::setSurfacePtcleWeightsLoads non-constant\n"
		       <<"conditions are not supported!"
		       << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	      }

	    }

	    // add weight
	    weights[weightID] += surface;

	    conditionSet = true;
	  }

	}
	
	// current particle has NOT got any deformation boundary condition 
	// for current surface yet
	if(!conditionSet) {

	  // set identifier
	  indexSet[0] = currentNode;

	  // set weight
	  indexSet[1] = weights.size(); // weight ID
	  pushBackVector(weights,surface);

	  // set the deformation boundary conditions  
	  idx = conditions.size();

	  affectedDOF.resize(idx+1);
	  conditions.resize(idx+1);

	  affectedDOF[idx] = elemDOF;
	  conditions[idx] = elemConds;

	  indexSet[2] = idx; // condition ID

	  // set surface normal

	  // check first if the normal vector is already set
	  normalSet = false;

	  for(int j=0;j<surfaceNormals.size();j++)

	    // other particle integration data than traction
	    // are already applied
	    if(surfaceNormals[j][0] == normal[0] && 
	       surfaceNormals[j][1] == normal[1] && 
	       surfaceNormals[j][2] == normal[2]) {

	      indexSet[3] = j; // surface normal ID;
	      normalSet = true;
	    }
        
	  if(!normalSet) {
	    indexSet[3] = surfaceNormals.size(); // surface normal ID

	    pushBackVector(surfaceNormals,normal);
	  }

	  storagePositions[position].
	    push_back(surfaceDispBoundPtcleIdx.size());
	  pushBackVector(surfaceDispBoundPtcleIdx,indexSet);
	}

      }

    }
    
  }

  /*********************************************************************/
  // surface rotation boundary constraints
  if(surfaceRotBoundElemIdx.size() > 0) {
    logFile<<"In function FEMGeometry::surfacePtcleWeightsLoads "
	   <<"surface rotation constraints are not supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }


#ifdef _FEdebugMode_
  dbVector conditionSurface(3);
  logFile<<"************ traction boundary particles ************"<<endl;
  for(int i=0;i<tractionBoundPtcleIdx.size();i++) {
    currentNode = tractionBoundPtcleIdx[i][0];
    int& weightID = tractionBoundPtcleIdx[i][1];
    int& conditionID = tractionBoundPtcleIdx[i][2];
    int& surfaceID = tractionBoundPtcleIdx[i][3];
    logFile<<i<<") BOUND INT PARTICLE "<<currentNode<<": "
	   <<particles[currentNode].getCoord(0)<<" "
	   <<particles[currentNode].getCoord(1)<<" "
	   <<particles[currentNode].getCoord(2)<<" ";
    dbVector& vec = particles[currentNode].getSurfaceNormal(surfaceID);
    double& weight = particles[currentNode].getIntWeight(weightID);
    dbVector& traction = particles[currentNode].getTraction(conditionID);
    blVector& affectedDOF = particles[currentNode].getTractionDOF(conditionID);
    logFile<<"normal "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<" ";
    logFile<<"traction: ";
    for(int j=0;j<affectedDOF.size();j++)
      if(affectedDOF[j])
    logFile<<traction[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
    conditionSurface[0] += weight;
  }
  logFile<<"************ pressure boundary particles *************"<<endl;
  for(int i=0;i<surfacePressureBoundPtcleIdx.size();i++) {
    currentNode = surfacePressureBoundPtcleIdx[i][0];
    int& weightID = surfacePressureBoundPtcleIdx[i][1];
    int& conditionID = surfacePressureBoundPtcleIdx[i][2];
    int& surfaceID = surfacePressureBoundPtcleIdx[i][3];
    logFile<<i<<") BOUND INT PARTICLE "<<currentNode<<": "
	   <<particles[currentNode].getCoord(0)<<" "
	   <<particles[currentNode].getCoord(1)<<" "  
	   <<particles[currentNode].getCoord(2)<<" ";  
    dbVector& vec = particles[currentNode].getSurfaceNormal(surfaceID);
    double& weight = particles[currentNode].getIntWeight(weightID);
    double& pressure = particles[currentNode].getSurfacePressure(conditionID);
    logFile<<"normal "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<" ";
    logFile<<"pressure = "<<pressure<<" weight = "<<weight<<endl;
    conditionSurface[1] += weight;
  }
  logFile<<"********** displacement boundary particles ************"<<endl;
  for(int i=0;i<surfaceDispBoundPtcleIdx.size();i++) {
    currentNode = surfaceDispBoundPtcleIdx[i][0];
    int& weightID = surfaceDispBoundPtcleIdx[i][1];
    int& conditionID = surfaceDispBoundPtcleIdx[i][2];
    int& surfaceID = surfaceDispBoundPtcleIdx[i][3];
    logFile<<i<<") BOUND INT PARTICLE "<<currentNode<<": "
	   <<particles[currentNode].getCoord(0)<<" "
	   <<particles[currentNode].getCoord(1)<<" "
	   <<particles[currentNode].getCoord(2)<<" ";
    dbVector& vec = particles[currentNode].getSurfaceNormal(surfaceID);
    double& weight = particles[currentNode].getIntWeight(weightID);
    dbVector& conditions = particles[currentNode].getDeformationBoundConds(conditionID);
    blVector& affectedDOF = particles[currentNode].getDeformationBoundDOF(conditionID);
    logFile<<"normal "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<" ";
    logFile<<"deformation bounds: ";
    for(int j=0;j<affectedDOF.size();j++)
      if(affectedDOF[j])
    logFile<<conditions[j]<<" ";
    logFile<<"weight = "<<weight<<endl;
    conditionSurface[2] += weight;
  }
  logFile<<"TRACTION SURFACE = "<<conditionSurface[0]<<endl;
  logFile<<"PRESSURE SURFACE = "<<conditionSurface[1]<<endl;
  logFile<<"DEFORMATION SURFACE = "<<conditionSurface[2]<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Return all surface integration Gauss point indices.
intMatrix& FEMGeometry::getAllSurfaceGaussPtsIdx(std::ofstream& logFile) {

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
intMatrix& FEMGeometry::getAllLineGaussPtsIdx(std::ofstream& logFile) {

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
// Return a pointer to tools and properties of a volume standard element
ElementTemplate* FEMGeometry::getVolumeElementTemplate(InputFileData* InputData,
						       FEMElement& elem,
						       std::ofstream& logFile) { 

  using namespace std;

  int& elemType = elem.getElemType();

  intVector data(2);
  data[0] = elemType;
  data[1] = elem.getElemOrder();

  map<string,double> params;
  getFEMMeshData(data,params);

  int nodesPerElem = (int)params["nodesPerVolumeElement"];

  if(elemType >= volumeElemTemplateID.size())
    volumeElemTemplateID.resize(elemType+1);

  int elemID = nodesPerElem;

  map<int,int>::iterator p = 
    volumeElemTemplateID[elemType].find(elemID);

  // Key was found.
  if(p != volumeElemTemplateID[elemType].end())
  
    return volumeElemTemplates[p->second];

  // Key was not found. 
  else {

    int position = volumeElemTemplates.size();
    volumeElemTemplateID[elemType][elemID] = position;
    volumeElemTemplates.resize(position+1);

    // Set a pointer to all necessary FEM shape generation data.
    switch (elemType) {

      // tetrahedral
    case 1:
      
      switch(nodesPerElem) {

      case 4:
	volumeElemTemplates[position] = new Tetra4ElementTemplate();
	break;

      default:
	logFile<<"In function FEMGeometry::getVolumeElementTemplate element "
	       <<nodesPerElem<<" nodes per element are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      break;
      
      // hexahedral
    case 2:

      switch(nodesPerElem) {

      case 8:
	volumeElemTemplates[position] = new Cube8ElementTemplate();
	break;

     case 27:
	volumeElemTemplates[position] = new Cube27ElementTemplate();
	break;

      default:
	logFile<<"In function FEMGeometry::getVolumeElementTemplate element "
	       <<nodesPerElem<<" nodes per element are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
      break;
      
    default:
      logFile<<"In function FEMGeometry::getVolumeElementTemplate element "
	     <<"type "<<elemType<<"\n is not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    return volumeElemTemplates[position];
  }

}

/************************************************************************/
// Return a pointer to tools and properties of a surface standard element
ElementTemplate* FEMGeometry::getSurfaceElementTemplate(InputFileData* InputData,
							FEMElement& elem,
							std::ofstream& logFile) { 
  
  using namespace std;

  int& elemType = elem.getElemType();

  intVector data(2);
  data[0] = elemType;
  data[1] = elem.getElemOrder();

  map<string,double> params;
  getFEMMeshData(data,params);

  int nodesPerElem = (int)params["nodesPerSurfaceElement"];


  if(elemType >= surfaceElemTemplateID.size())
    surfaceElemTemplateID.resize(elemType+1);

  int elemID = nodesPerElem;

  map<int,int>::iterator p = 
    surfaceElemTemplateID[elemType].find(elemID);

  // Key was found.
  if(p != surfaceElemTemplateID[elemType].end())

    return surfaceElemTemplates[p->second];

  // Key was not found. 
  else {

    int position = surfaceElemTemplates.size();
    surfaceElemTemplateID[elemType][elemID] = position;
    surfaceElemTemplates.resize(position+1);

    // Set a pointer to all necessary FEM shape generation data.
    switch (elemType) {

      // tetrahedral
    case 1:
      
      switch(nodesPerElem) {

      case 3:
	surfaceElemTemplates[position] = new Tria3ElementTemplate();

	break;

      default:
	logFile<<"In function FEMGeometry::getSurfaceElementTemplate element "
	       <<nodesPerElem<<" nodes per element are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      break;
      
      // hexahedral
    case 2:

      switch(nodesPerElem) {

      case 4:
	surfaceElemTemplates[position] = new Rect4ElementTemplate();
	break;

      case 9:
	surfaceElemTemplates[position] = new Rect9ElementTemplate();
	break;

      default:
	logFile<<"In function FEMGeometry::getSurfaceElementTemplate element "
	       <<nodesPerElem<<" nodes per element are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
      break;
      
    default:
      logFile<<"In function FEMGeometry::getSurfaceElementTemplate element "
	     <<"type "<<elemType<<"\n is not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    return surfaceElemTemplates[position];
  }

}

/************************************************************************/
// Return a pointer to tools and properties of a line standard element
ElementTemplate* FEMGeometry::getLineElementTemplate(InputFileData* InputData,
						     FEMElement& elem,
						     std::ofstream& logFile) { 
  
  using namespace std;

  int& elemType = elem.getElemType();

  intVector data(2);
  data[0] = elemType;
  data[1] = elem.getElemOrder();

  map<string,double> params;
  getFEMMeshData(data,params);

  int nodesPerElem = (int)params["nodesPerLineElement"];


  if(elemType >= lineElemTemplateID.size())
    lineElemTemplateID.resize(elemType+1);

  int elemID = nodesPerElem;

  map<int,int>::iterator p = 
    lineElemTemplateID[elemType].find(elemID);

  // Key was found.
  if(p != lineElemTemplateID[elemType].end())

    return lineElemTemplates[p->second];

  // Key was not found. 
  else {

    int position = lineElemTemplates.size();
    lineElemTemplateID[elemType][elemID] = position;
    lineElemTemplates.resize(position+1);

    // Set a pointer to all necessary FEM shape generation data.
    switch(nodesPerElem) {
      
    case 2:
      lineElemTemplates[position] = new Line2ElementTemplate();
      break;
      
    case 3:
      lineElemTemplates[position] = new Line3ElementTemplate();
      break;
      
    default:
      logFile<<"In function FEMGeometry::getLineElementTemplate element "
	     <<nodesPerElem<<" nodes per element are not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    return lineElemTemplates[position];
  }

}

/************************************************************************/
/************************************************************************/
dbMatrix& FEMGeometry::getVolumeFEShapeFunctions(InputFileData* InputData,
						 FEMElement& elem,
						 std::ofstream& logFile) { 

  using namespace std;

  int& elemType = elem.getElemType();
  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getVolumeIntegrationPts();

  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();

  if(elemType >= volumeFEShapeFuncSetID.size())
    volumeFEShapeFuncSetID.resize(elemType+1);

  int shapeFuncSetID = (nodesPerElem+1000)*gaussPtsPerElem;

  map<int,int>::iterator p = 
    volumeFEShapeFuncSetID[elemType].find(shapeFuncSetID);

  // Key was found.
  if(p != volumeFEShapeFuncSetID[elemType].end())

    return volumeFEShapeFuncSets[p->second];

 
  // Key was not found. 
  else {

    int position = volumeFEShapeFuncSets.size();

    volumeFEShapeFuncSetID[elemType][shapeFuncSetID] = position;
    allocateArray(volumeFEShapeFuncSets,position+1);

    setVolumeFEShapeFunctions(InputData,elem,
			      volumeFEShapeFuncSets[position],
			      logFile);
 
    return volumeFEShapeFuncSets[position];
  }

}

/************************************************************************/
/************************************************************************/
dbMatrix3& FEMGeometry::getVolumeFEShapeFuncDerivs(InputFileData* InputData,
						   FEMElement& elem,
						   std::ofstream& logFile) { 
  
  using namespace std;

  int& elemType = elem.getElemType();
  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getVolumeIntegrationPts();

  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();

  if(elemType >= volumeFEShapeFuncDerivSetID.size())
    volumeFEShapeFuncDerivSetID.resize(elemType+1);

  int shapeFuncSetID = (nodesPerElem+1000)*gaussPtsPerElem;

  map<int,int>::iterator p = 
    volumeFEShapeFuncDerivSetID[elemType].find(shapeFuncSetID);

  // Key was found.
  if(p != volumeFEShapeFuncDerivSetID[elemType].end())

    return volumeFEShapeFuncDerivSets[p->second];

 
  // Key was not found. 
  else {

    int position = volumeFEShapeFuncDerivSets.size();

    volumeFEShapeFuncDerivSetID[elemType][shapeFuncSetID] = position;
    allocateArray(volumeFEShapeFuncDerivSets,position+1);

    setVolumeFEShapeFuncDerivs(InputData,elem,
			       volumeFEShapeFuncDerivSets[position],
			       logFile);
 
    return volumeFEShapeFuncDerivSets[position];
  }

}

/************************************************************************/
/************************************************************************/
dbMatrix& FEMGeometry::getSurfaceFEShapeFunctions(InputFileData* InputData,
						  FEMElement& elem,
						  std::ofstream& logFile) { 
  using namespace std;

  int& elemType = elem.getElemType();
  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getSurfaceIntegrationPts();
 
  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();

  if(elemType >= surfaceFEShapeFuncSetID.size())
    surfaceFEShapeFuncSetID.resize(elemType+1);

  int shapeFuncSetID = (nodesPerElem+1000)*gaussPtsPerElem;

  map<int,int>::iterator p = 
    surfaceFEShapeFuncSetID[elemType].find(shapeFuncSetID);

  // Key was found.
  if(p != surfaceFEShapeFuncSetID[elemType].end())

    return surfaceFEShapeFuncSets[p->second];

  // Key was not found. 
  else {

    int position = surfaceFEShapeFuncSets.size();

    surfaceFEShapeFuncSetID[elemType][shapeFuncSetID] = position;
    allocateArray(surfaceFEShapeFuncSets,position+1);

    setSurfaceFEShapeFunctions(InputData,elem,
			       surfaceFEShapeFuncSets[position],
			       logFile);

    return surfaceFEShapeFuncSets[position];
  }

}

/************************************************************************/
/************************************************************************/
dbMatrix& FEMGeometry::getLineFEShapeFunctions(InputFileData* InputData,
					       FEMElement& elem,
					       std::ofstream& logFile) { 
  using namespace std;

  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getLineIntegrationPts();
  int& elemType = elem.getElemType();

  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();

  if(elemType >= lineFEShapeFuncSetID.size())
    lineFEShapeFuncSetID.resize(elemType+1);

  int shapeFuncSetID = (nodesPerElem+1000)*gaussPtsPerElem;

  map<int,int>::iterator p = 
    lineFEShapeFuncSetID[elemType].find(shapeFuncSetID);

  // Key was found.
  if(p != lineFEShapeFuncSetID[elemType].end())

    return lineFEShapeFuncSets[p->second];

  // Key was not found. 
  else {

    int position = lineFEShapeFuncSets.size();

    lineFEShapeFuncSetID[elemType][shapeFuncSetID] = position;
    allocateArray(lineFEShapeFuncSets,position+1);

    setLineFEShapeFunctions(InputData,elem,
			    lineFEShapeFuncSets[position],
			    logFile);

    return lineFEShapeFuncSets[position];
  }

}

/************************************************************************/
/************************************************************************/
// Calculate a set of volume shape function used to determine the Gauss 
// points of a volume element.
void FEMGeometry::setVolumeFEShapeFunctions(InputFileData* InputData,
					    FEMElement& elem,
					    dbMatrix& sFuncs,
					    std::ofstream& logFile) { 

  using namespace std;

  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getVolumeIntegrationPts();
  int& elemType = elem.getElemType();
  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();
  
  ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
  GaussPointSet* VolumeGaussSet = elem.getVolumeGaussSet();
 
  allocateArray(sFuncs,gaussPtsPerElem,nodesPerElem);

  for(int i=0;i<gaussPtsPerElem;i++)

    for(int j=0;j<nodesPerElem;j++)

      sFuncs[i][j] = 
	FEVolumeSet->N(j,VolumeGaussSet->coord[i]);


#ifdef _FEdebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"********* volume element shape functions *************"<<endl;
  for(int i=0;i<sFuncs.size();i++)
    for(int j=0;j<sFuncs[i].size();j++)
      logFile<<"N"<<j<<"("<<i<<")="<<sFuncs[i][j]<<endl;
  double ord;
  for(int i=0;i<gaussPtsPerElem;i++) {
    logFile<<"GAUSS POINT "<<i<<": ";
    ord = 0;
    for(int j=0;j<nodesPerElem;j++) {
      logFile<<sFuncs[i][j]<<" ";
      ord += sFuncs[i][j];
    }
    if(ord != 1.0)
      logFile<<"<==> sum ord = "<<ord<<endl;
    else
      logFile<<endl;
  }
  logFile<<"---------------------------------------------------"<<endl;
  for(int i=0;i<nodesPerElem;i++) {
    for(int j=0;j<nodesPerElem;j++) {
      if(i == j)
	logFile<<"N"<<j<<"(x_"<<i<<") = "
	       <<FEVolumeSet->N(j,FEVolumeSet->nodalCoords[i])<<endl;
      else if(fabs(FEVolumeSet->N(j,FEVolumeSet->nodalCoords[i])) > 0)
	logFile<<"<==> N_"<<j<<"(x_"<<i<<") = "
	       <<FEVolumeSet->N(j,FEVolumeSet->nodalCoords[i])<<endl;
    }
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Calculate a set of volume shape function derivatives of a volume 
// element with respect to the local element coordinates.
void FEMGeometry::setVolumeFEShapeFuncDerivs(InputFileData* InputData,
					     FEMElement& elem,
					     dbMatrix3& sFuncs,
					     std::ofstream& logFile) { 

  using namespace std;

  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getVolumeIntegrationPts();
  int& elemType = elem.getElemType();
  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();
  
  ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
  GaussPointSet* VolumeGaussSet = elem.getVolumeGaussSet();
 
  allocateArray(sFuncs,3,gaussPtsPerElem,nodesPerElem);

  // dN_j/dx_k = dN_j/dL_r dL_r/dx_k

  for(int k=0;k<3;k++) 

    for(int i=0;i<gaussPtsPerElem;i++)
      
      for(int j=0;j<nodesPerElem;j++)
	
	sFuncs[k][i][j] = 
	  FEVolumeSet->dN(j,k,VolumeGaussSet->coord[i]);
  
  
#ifdef _FEdebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"********* local FE-shape function derivs *************"<<endl;
  for(int k=0;k<sFuncs.size();k++) {
    logFile<<"****************************************************"<<endl;
    logFile<<"dN/dx"<<k<<" ---------------------------------------"<<endl;
    for(int i=0;i<sFuncs[k].size();i++) {
      logFile<<"GPt "<<i<<": "<<endl;
      for(int j=0;j<sFuncs[k][i].size();j++)
	logFile<<"dN"<<j<<"("<<VolumeGaussSet->coord[i][0]
	       <<","<<VolumeGaussSet->coord[i][1]<<","
	       <<VolumeGaussSet->coord[i][2]<<")/dL"<<k<<"="
	       <<sFuncs[k][i][j]<<endl;
    }
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate a set of surface shape function used to determine the Gauss 
// points for a surface element.
void FEMGeometry::setSurfaceFEShapeFunctions(InputFileData* InputData,
					     FEMElement& elem,
					     dbMatrix& sFuncs,
					     std::ofstream& logFile) { 
  using namespace std;

  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getSurfaceIntegrationPts();
  int& elemType = elem.getElemType();
  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();

  ElementTemplate* FESurfaceSet = elem.getSurfaceElementTemplate();
  GaussPointSet* SurfaceGaussSet = elem.getSurfaceGaussSet();

  allocateArray(sFuncs,gaussPtsPerElem,nodesPerElem);

  for(int i=0;i<gaussPtsPerElem;i++)

    for(int j=0;j<nodesPerElem;j++)

      sFuncs[i][j] = 
	FESurfaceSet->N(j,SurfaceGaussSet->coord[i]);

#ifdef _FEdebugMode_
  logFile<<"********* surface element shape functions ************"<<endl;
  double ord;
  for(int i=0;i<gaussPtsPerElem;i++) {
    ord = 0;
    logFile<<"GAUSS POINT "<<i<<": ";
    for(int j=0;j<nodesPerElem;j++) {
      logFile<<sFuncs[i][j]<<" ";
      ord += sFuncs[i][j];
    }
    if(ord != 1.0)
      logFile<<"<==> sum ord = "<<ord<<endl;
    else
      logFile<<endl;
  }
  logFile<<"---------------------------------------------------"<<endl;
  for(int i=0;i<nodesPerElem;i++) {
    for(int j=0;j<nodesPerElem;j++) {
      if(i == j)
	logFile<<"N_"<<j<<"(x_"<<i<<") = "
	       <<FESurfaceSet->N(j,FESurfaceSet->nodalCoords[i])<<endl;
      else if(fabs(FESurfaceSet->N(j,FESurfaceSet->nodalCoords[i])) > 0)
	logFile<<"<==> N_"<<j<<"(x_"<<i<<") = "
	       <<FESurfaceSet->N(j,FESurfaceSet->nodalCoords[i])<<endl;
    }
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Calculate a set of line shape function used to determine the Gauss 
// points for a line element.
void FEMGeometry::setLineFEShapeFunctions(InputFileData* InputData,
					  FEMElement& elem,
					  dbMatrix& sFuncs,
					  std::ofstream& logFile) { 
  using namespace std;

  intVector& nodes = elem.getNodes();
  intVector& intPts = elem.getLineIntegrationPts();
  int& elemType = elem.getElemType();
  int nodesPerElem = nodes.size();
  int gaussPtsPerElem = intPts.size();

  ElementTemplate* FELineSet = elem.getLineElementTemplate();
  GaussPointSet* LineGaussSet = elem.getLineGaussSet();

  allocateArray(sFuncs,gaussPtsPerElem,nodesPerElem);

  for(int i=0;i<gaussPtsPerElem;i++)

    for(int j=0;j<nodesPerElem;j++)
      
      sFuncs[i][j] = 
	FELineSet->N(j,LineGaussSet->coord[i]);



#ifdef _FEdebugMode_
  double ord;
  logFile<<"*********** line element shape functions *************"<<endl;
  for(int i=0;i<gaussPtsPerElem;i++) {
    ord = 0;
    logFile<<"GAUSS POINT "<<i<<": ";
    for(int j=0;j<nodesPerElem;j++) {
      logFile<<sFuncs[i][j]<<" ";
      ord += sFuncs[i][j];
    }
    if(ord != 1.0)
      logFile<<"<==> sum ord = "<<ord<<endl;
    else
      logFile<<endl;
  }
  logFile<<"---------------------------------------------------"<<endl;
  for(int i=0;i<nodesPerElem;i++) {
    for(int j=0;j<nodesPerElem;j++) {
      if(i == j)
	logFile<<"N_"<<j<<"(x_"<<i<<") = "
	       <<FELineSet->N(j,FELineSet->nodalCoords[i])<<endl;
      else if(fabs(FELineSet->N(j,FELineSet->nodalCoords[i])) > 0)
	logFile<<"<==> N_"<<j<<"(x_"<<i<<") = "
	       <<FELineSet->N(j,FELineSet->nodalCoords[i])<<endl;
    }
  }
#endif

}

/************************************************************************/
/************************************************************************/
// get a pointer to the tools and properties of a volume Gauss quadrature
GaussPointSet* FEMGeometry::getVolumeGaussSet(InputFileData* InputData,
					      FEMElement& elem,
					      std::ofstream& logFile) { 

  using namespace std;

  intVector& intPts = elem.getVolumeIntegrationPts();
  int gaussPtsPerElem = intPts.size();
  int& elemType = elem.getElemType();

  if(elemType >= volumeGaussPtTemplateID.size())
    volumeGaussPtTemplateID.resize(elemType+1);

  int gaussID = gaussPtsPerElem;

  map<int,int>::iterator p = 
    volumeGaussPtTemplateID[elemType].find(gaussID);

  // Key was found.
  if(p != volumeGaussPtTemplateID[elemType].end())
   
    return volumeGaussPtTemplates[p->second];
  
  // Key was not found. 
  else {

    int position = volumeGaussPtTemplates.size();
    volumeGaussPtTemplateID[elemType][gaussID] = position;
    volumeGaussPtTemplates.resize(position+1);

    // Set a pointer to all necessary FEM shape generation data.
    switch (elemType) {

      // tetrahedral
    case 1:
      
      switch(gaussPtsPerElem) {

      case 1:
	volumeGaussPtTemplates[position] = new GaussSetTetra1();
	break;

      case 4:
	volumeGaussPtTemplates[position] = new GaussSetTetra4();
	break;

      case 5:
	volumeGaussPtTemplates[position] = new GaussSetTetra5();
	break;

      default:
	logFile<<"In function FEMGeometry::getVolumeGaussSet element\n "
	       <<gaussPtsPerElem<<" Gauss points per tetrahedral element\n " 
	       <<" are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      break;
      
      // hexahedral
    case 2:

      switch(gaussPtsPerElem) {

      case 1:
	volumeGaussPtTemplates[position] = new GaussSetCube1();
	break;

      case 8:
	volumeGaussPtTemplates[position] = new GaussSetCube8();
	break;

      case 27:
	volumeGaussPtTemplates[position] = new GaussSetCube27();
	break;

      case 64:
	volumeGaussPtTemplates[position] = new GaussSetCube64();
	break;

      case 125:
	volumeGaussPtTemplates[position] = new GaussSetCube125();
	break;

      default:
	logFile<<"In function FEMGeometry::getVolumeGaussSet element\n "
	       <<gaussPtsPerElem<<" Gauss points per hexahedral element\n " 
	       <<" are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
      break;
      
    default:
      logFile<<"In function FEMGeometry::getVolumeGaussSet element "
	     <<"type "<<elemType<<"\n is not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    return volumeGaussPtTemplates[position];
  }

}

/************************************************************************/
/************************************************************************/
// get a pointer to the tools and properties of a surface Gauss quadrature
GaussPointSet* FEMGeometry::getSurfaceGaussSet(InputFileData* InputData,
					       FEMElement& elem,
					       std::ofstream& logFile) { 

  using namespace std;

  intVector& intPts = elem.getSurfaceIntegrationPts();
  int gaussPtsPerElem = intPts.size();
  int& elemType = elem.getElemType();

  if(elemType >= surfaceGaussPtTemplateID.size())
    surfaceGaussPtTemplateID.resize(elemType+1);

  int gaussID = gaussPtsPerElem;

  map<int,int>::iterator p = 
    surfaceGaussPtTemplateID[elemType].find(gaussID);

  // Key was found.
  if(p != surfaceGaussPtTemplateID[elemType].end())

    return surfaceGaussPtTemplates[p->second];

  // Key was not found. 
  else {

    int position = surfaceGaussPtTemplates.size();
    surfaceGaussPtTemplateID[elemType][gaussID] = position;
    surfaceGaussPtTemplates.resize(position+1);

    // Set a pointer to all necessary FEM shape generation data.
    switch (elemType) {

      // triangle
    case 1:
      
      switch(gaussPtsPerElem) {

      case 1:
	surfaceGaussPtTemplates[position] = new GaussSetTria1();
	break;

      case 3:
	surfaceGaussPtTemplates[position] = new GaussSetTria3();
	break;

      case 4:
	surfaceGaussPtTemplates[position] = new GaussSetTria4();
	break;

      default:
	logFile<<"In function FEMGeometry::getSurfaceGaussSet element\n "
	       <<gaussPtsPerElem<<" Gauss points per triangle element\n " 
	       <<" are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      break;
      
      // rectangle
    case 2:

      switch(gaussPtsPerElem) {

      case 1:
	surfaceGaussPtTemplates[position] = new GaussSetRect1();
	break;

      case 4:
	surfaceGaussPtTemplates[position] = new GaussSetRect4();
	break;

      case 9:
	surfaceGaussPtTemplates[position] = new GaussSetRect9();
	break;

      case 16:
	surfaceGaussPtTemplates[position] = new GaussSetRect16();
	break;

      case 25:
	surfaceGaussPtTemplates[position] = new GaussSetRect25();
	break;

      default:
	logFile<<"In function FEMGeometry::getSurfaceGaussSet element\n "
	       <<gaussPtsPerElem<<" Gauss points per rectangle element\n " 
	       <<" are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
      break;
      
    default:
      logFile<<"In function FEMGeometry::getSurfaceGaussSet element "
	     <<"type "<<elemType<<"\n is not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    return surfaceGaussPtTemplates[position];
  }

}

/************************************************************************/
/************************************************************************/
// get a pointer to the tools and properties of a line Gauss quadrature
GaussPointSet* FEMGeometry::getLineGaussSet(InputFileData* InputData,
					    FEMElement& elem,
					    std::ofstream& logFile) { 

  using namespace std;

  intVector& intPts = elem.getLineIntegrationPts();
  int gaussPtsPerElem = intPts.size();

  int elemType = 1;

  if(elemType >= lineGaussPtTemplateID.size())
    lineGaussPtTemplateID.resize(elemType+1);

  int gaussID = gaussPtsPerElem;

  map<int,int>::iterator p = 
    lineGaussPtTemplateID[elemType].find(gaussID);

  // Key was found.
  if(p != lineGaussPtTemplateID[elemType].end())

    return lineGaussPtTemplates[p->second];

  // Key was not found. 
  else {

    int position = lineGaussPtTemplates.size();
    lineGaussPtTemplateID[elemType][gaussID] = position;
    lineGaussPtTemplates.resize(position+1);

    // Set a pointer to all necessary FEM shape generation data.
    switch (elemType) {

    case 1:
      
      switch(gaussPtsPerElem) {

      case 1:
	lineGaussPtTemplates[position] = new GaussSetLine1();
	break;

      case 2:
	lineGaussPtTemplates[position] = new GaussSetLine2();
	break;

      case 3:
	lineGaussPtTemplates[position] = new GaussSetLine3();
	break;

      case 4:
	lineGaussPtTemplates[position] = new GaussSetLine4();
	break;

      case 5:
	lineGaussPtTemplates[position] = new GaussSetLine5();
	break;

      default:
	logFile<<"In function FEMGeometry::getLineGaussSet element\n "
	       <<gaussPtsPerElem<<" Gauss points per line element\n " 
	       <<" are not supported!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      break;
      
    default:
      logFile<<"In function FEMGeometry::getLineGaussSet element "
	     <<"type "<<elemType<<"\n is not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    return lineGaussPtTemplates[position];
  }

}


/************************************************************************/
/************************************************************************/
// Compute shape function derivatives at all Gauss points with respect
// to global coordinates.
void FEMGeometry::setGlobalFEShapeFuncDerivs(InputFileData* InputData,
					     std::map<std::string,double>& modelData,
					     std::ofstream& logFile,
					     PetscViewer& viewerMPI) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];

  int idx;
  dbMatrix jacobian,jacobianInv;

#ifdef _FEdebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"********* compute FE-shape function derivs **********"<<endl;
#endif


  // dN^I/dx_r = dN^I/dL_s dL_s/dx_r

  // Loop over all local elements.
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];

#ifdef _FEdebugMode_
    logFile<<"****************************************************"<<endl;
    logFile<<"Element "<<i<<" ("<<elem.getGlobalID()<<")"<<endl;
#endif

    intVector& intPts = nodesElements[i].getVolumeIntegrationPts();
    ElementTemplate* FEVolumeSet = elem.getVolumeElementTemplate();
    GaussPointSet* VolumeGaussSet = elem.getVolumeGaussSet();
    intVector& nodes = elem.getNodes();
    int nodesPerElem = nodes.size();

    // Determine the set of shape functions.
    dbMatrix3& sFuncs = elem.getVolumeShapeFuncDerivOrds();

    // Loop over all gauss coordinates of a single element to determine
    // its gauss points.

    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];

#ifdef _FEdebugMode_
      logFile<<"--------------------------------------------------"<<endl;
      logFile<<"GPoint "<<j<<" ("<<gaussPoints[idx].getGlobalID()<<")"
	     <<endl;
#endif

      dbMatrix& feShapeDerivs = gaussPoints[idx].getFEShapeFuncDerivs();      
      allocateArray(feShapeDerivs,usedDims,nodes.size());

      // compute dL_s/dx_r
      FEVolumeSet->jacobian(particles,nodes,
			    VolumeGaussSet->coord[j],
			    jacobian,logFile);

      calcInvDoubleDense(jacobian,jacobianInv,logFile);

#ifdef _FEdebugMode_
      logFile<<"---------------- jacobi matrix inv ---------------"<<endl;
      for(int v=0;v<jacobianInv.size();v++)
	for(int w=0;w<jacobianInv[v].size();w++)
	  logFile<<"J["<<v<<"]["<<w<<"]="<<jacobianInv[v][w]<<endl;
#endif
      
      // loop over element nodes and compute dN^k/dx_r = dN^k/dL_s dL_s/dx_r
      for(int k=0;k<nodes.size();k++)

	// loop the global derivatives
	for(int r=0;r<usedDims;r++)

	  // loop over the local derivatives
	  for(int s=0;s<usedDims;s++)

	    feShapeDerivs[r][k] += sFuncs[s][j][k]*jacobianInv[s][r];
	

    }
    
  }

#ifdef _FEdebugMode_
  logFile<<"******************************************************"<<endl;
  logFile<<"************** FE-shape function derivs **************"<<endl;
  for(int i=0;i<nodesElements.size();i++) {
    FEMElement& elem = nodesElements[i];
    intVector& intPts = nodesElements[i].getVolumeIntegrationPts();
    for(int j=0;j<intPts.size();j++) {
      idx = intPts[j];
      logFile<<"--------------------------------------------------"<<endl;
      logFile<<"GPoint "<<j<<" ("<<gaussPoints[idx].getGlobalID()<<")"
	     <<endl;
      dbMatrix& feShapeDerivs = gaussPoints[idx].getFEShapeFuncDerivs(); 
      for(int v=0;v<feShapeDerivs.size();v++)
	for(int w=0;w<feShapeDerivs[v].size();w++)
	  logFile<<"dN["<<w<<"]/dx["<<v<<"]="<<feShapeDerivs[v][w]<<endl;
    }
  }
#endif

}
