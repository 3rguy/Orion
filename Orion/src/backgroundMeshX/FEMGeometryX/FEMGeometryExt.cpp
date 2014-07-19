#include "FEMGeometryExt.h"

/*!****************************************************************************/
/*!****************************************************************************/
FEMGeometryExt::FEMGeometryExt(InputFileData* InputData,
		std::map<std::string, double>& modelData, std::string& meshFileName,
		std::ofstream& logFile) :
		FEMGeoData(NULL),FEMInputData(NULL) {

	using namespace std;

	readMeshFile(InputData, modelData, meshFileName, logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
FEMGeometryExt::FEMGeometryExt(InputFileData* InputData,
		std::map<std::string, double>& modelData, std::string& meshFileName,
		std::string& inputFileName, std::ofstream& logFile) :
		FEMGeoData(NULL),FEMInputData(NULL) {

	using namespace std;

	int choice = (int) InputData->getValue("FEMGeometrySetupType");

	switch (choice) {
	case 1:
		readMeshFile(InputData, modelData, meshFileName, logFile);
		break;

	case 2:
		FEMGeoDataSetup(InputData, modelData, meshFileName, inputFileName,
				logFile);
		break;

	default:
		logFile
				<< "In FEMGeometryExt::FEMGeometryExt, FEMGeometrySetupType incorrect \n "
				<< "Valid options are: 1, 2" << endl;
		cout
				<< "In FEMGeometryExt::FEMGeometryExt, FEMGeometrySetupType incorrect \n "
				<< "Valid options are: 1, 2" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

FEMGeometryExt::~FEMGeometryExt() {

	for (int i = 0; i < lineElemTemplates.size(); i++)

		delete lineElemTemplates[i];

	for (int i = 0; i < surfaceElemTemplates.size(); i++)

		delete surfaceElemTemplates[i];

	for (int i = 0; i < volumeElemTemplates.size(); i++)

		delete volumeElemTemplates[i];

	for (int i = 0; i < lineGaussPtTemplates.size(); i++)

		delete lineGaussPtTemplates[i];

	for (int i = 0; i < surfaceGaussPtTemplates.size(); i++)

		delete surfaceGaussPtTemplates[i];

	for (int i = 0; i < volumeGaussPtTemplates.size(); i++)

		delete volumeGaussPtTemplates[i];
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read a complete Finite Elements's mesh from the mesh file.
void FEMGeometryExt::readMeshFile(InputFileData* InputData,
		std::map<std::string, double>& modelData, std::string& meshFileName,
		std::ofstream& logFile) {

	using namespace std;

	/*********************************************************************/
	// Read mesh file.
	ifstream meshFile(meshFileName.c_str());

	if (!meshFile) {
		logFile << "Can't open FEM mesh file: "<< meshFileName << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int usedDOF = InputData->getValue("usedDegreesOfFreedom");

	// -----------------------------------------------
	logFile <<"In FEMGeometryExt::readMeshFile "<< endl;
	// -----------------------------------------------

	int n, nodeID, nodesNum;
	double coord1, coord2, coord3;
	double* dummyPointer;
	string name;

	/*!********************************************************************/
	//! Read the nodes' coordinates.
	meshFile >> name;
	meshFile >> nodesNum;

	if (name != "COORDINATES") {
		logFile << "Mesh file 'mesh.dat' is invalid!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//! Define the particles vector using the nodal coordinates
	particles.resize(nodesNum, ParticleExt(usedDOF));

	dbMatrix ptcls(3, dbVector(nodesNum));

	//! Loop over the nodes number in the file to read their coordinates
	for (int i = 0; i < nodesNum; i++) {

		if (meshFile.eof()) {
			logFile << "Mesh file 'mesh.dat' doesn't contain enough data!"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		//! Read the nodeID, name and the coordinates in each axis
		meshFile >> nodeID >> name >> ptcls[0][i] >> ptcls[1][i] >> ptcls[2][i];

		//! Stores the coordinates in the particles vector
		particles[i].setCoords(ptcls[0][i], ptcls[1][i], ptcls[2][i]);

		//! Set the particles ID
		particles[i].setID(nodeID);

	}

#ifdef _FEdebugMode_
	logFile << "###################################################" << endl;
	logFile << "******************** Mesh Coords ******************" << endl;
	for (int i = 0; i < nodesNum; i++)
		logFile << i << " " << ptcls[0][i] << " " << ptcls[1][i] << " "
				<< ptcls[2][i] << endl;
	logFile << "****************** Particles ******************" << endl;
	for (int i = 0; i < nodesNum; i++)
		logFile << "MESH-Id " << particles[i].getID() << " coords: "
				<< particles[i].getCoord(0) << " " << particles[i].getCoord(1)
				<< " " << particles[i].getCoord(2) << endl;
#endif

	/*!************************************************************************/
	//! Read and the element connectivity list.
	int globalElemNum, nodesPerElem;
	int referenceElemType;

	//! Read the title name of the element connectivity list and
	//! the number of elements
	meshFile >> name;
	meshFile >> globalElemNum;

	if (globalElemNum == 0) {
		logFile
				<< "In FEMGeometryExt::readMeshFile no volume elements in 'mesh.dat'.\n"
				<< "Check elemType - must be hexahedral or tetrahedral!"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//! Defines a vector of FEMElement containing all the elements
	//! allNodesElements: Not sure of its purpose !
	//! nodesElements: Contains all the volume elements
	vector<FEMElementExt> allNodesElements(globalElemNum, FEMElementExt(usedDOF));
	nodesElements.resize(globalElemNum, FEMElementExt(usedDOF));

	int elemID, node;
	map<string, double> params;
	intVector data(4);

	//! Loop over all global elements to extract their material and
	//! element information
	for (int i = 0; i < globalElemNum; i++) {

		//! Set the element's type, order and material ID
		int& elemType = allNodesElements[i].getElemType();
		int& elemOrder = allNodesElements[i].getElemOrder();
		int& matID = allNodesElements[i].getMaterialID();

		meshFile >> elemID >> name;

		nodesElements[i].setGlobalID(elemID);

		//! store element's material ID, type (tetrahedra,hexahedra),
		//! order (linear,quadratic)
		meshFile >> matID >> elemType >> elemOrder >> name;

		nodesElements[i].setElemType(elemType);
		nodesElements[i].setElemOrder(elemOrder);
		nodesElements[i].setMaterialID(matID);

		if (matID > (int) InputData->getNumberOfMats() || matID == 0) {
			logFile
					<< "In FEMGeometryExt::readMeshFile matID in mesh file 'mesh.dat'\n"
					<< "does not match materials specified in input file 'input.dat'!"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		//! Auto-correct the material ID number
		matID -= 1;

		//! elemType is uniform through-out the problem domain
		if (i == 0)
			referenceElemType = elemType;

		else if (referenceElemType != elemType) {
			logFile << "In FEMGeometryExt::readMeshFile elemType must be uniform\n"
					<< "thoughout the problem domain!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		//! set element's read info so as to extract their corresponding details
		//! from getFEMMeshData using the map object "params"
		data[0] = elemType;
		data[1] = elemOrder;
		getFEMMeshData(data, params);
		nodesPerElem = (int) params["nodesPerVolumeElement"];

		intVector nodeVec(nodesPerElem);

		//! Read the nodes' ID
		for (int j = 0; j < nodesPerElem; j++) {
			meshFile >> nodeVec[j];

			//! For each nodes read, store its material ID
			particles[nodeVec[j] - 1].setMaterialID(matID);

			//! For each nodes read, store its element ID
			particles[nodeVec[j] - 1].getElems().resize(
					particles[nodeVec[j] - 1].getElems().size() + 1); //! resizing
			particles[nodeVec[j] - 1].getElems()[particles[nodeVec[j] - 1].getElems().size()
					- 1] = i; //! storing
		}

		//! Store the nodes of the element to the nodesElement vector
		nodesElements[i].setNodes(nodeVec);


		/*!************************ TESTING(start) ************************/

		intMatrix tempMat = decompVolumeToSurfaceElems(nodeVec,elemType,
				InputData,modelData,logFile);

		nodesElements[i].setSurfaceElems(
				surfaceIDGenerator(tempMat,InputData,modelData,logFile)) ;


		/*!************************* TESTING(end) *************************/


		//! Element greater than 2 are not yet supported
		if (abs(elemType) > 2) {
			logFile << "In FEMGeometryExt::readMeshFile elemType " << elemType
					<< " " << "specified in file 'mesh.dat' is not supported!"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		//! If file does not contain element's info, then program cannot proceed
		if (meshFile.eof()) {
			logFile << "Mesh file 'mesh.dat' doesn't contain enough data!"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

	}


#ifdef _FEdebugMode_
	logFile << "######################################################" << endl;
	logFile << "************** Surface Elements Info *****************" << endl;
	logFile << "Number of Surface Elements: " << surfaceNodesElems.size()
			<< endl;
	logFile << "Surface -> Node connectivity:" << endl;

	for (int l = 0; l < surfaceNodesElems.size(); l++) {
		logFile << surfaceNodesElems[l].getGlobalID() << ".) "
				<< surfaceNodesElems[l].getMaterialID() << " "
				<< surfaceNodesElems[l].getElemOrder() << " "
				<< surfaceNodesElems[l].getElemType() << ":";
		for (int m = 0; m < surfaceNodesElems[l].getNodes().size(); m++) {
			logFile << " " << surfaceNodesElems[l].getNodes()[m];
		}
		logFile << endl;
	}
#endif

#ifdef _FEdebugMode_

	logFile << "###################################################" << endl;
	logFile << "**************** Volume-Surface info **************" << endl;

	for (int i = 0; i < nodesElements.size(); i++) {
		logFile << nodesElements[i].getGlobalID() << ":";
		for (int j = 0; j < nodesElements[i].getSurfaceElems().size(); j++) {
			logFile << " " << nodesElements[i].getSurfaceElems()[j];
		}
		logFile << endl;
	}

#endif

	//! setting up the surface details in the volume element list
	//combiningSurfaceVolumeElems(InputData, modelData, logFile);

	/*********************************************************************/
	// Store information concerning the background mesh.
	std::map<std::string, double>& backGroundMeshInfo =
			InputData->getBackGroundMeshInfo();

	backGroundMeshInfo["numberOfVolumeElements"] = (double) globalElemNum;
	backGroundMeshInfo["numberOfSurfaceElements"] = (double) globalElemNum;

#ifdef _FEdebugMode_
	logFile << "######################################################" << endl;
	logFile << "*************** backgroundmesh info ******************" << endl;
	logFile << "nodesPerVolumeElement = "
			<< backGroundMeshInfo["nodesPerVolumeElement"] << endl;
	logFile << "nodesPerSurfaceElement = "
			<< backGroundMeshInfo["nodesPerSurfaceElement"] << endl;
	logFile << "nodesPerLineElement = "
			<< backGroundMeshInfo["nodesPerLineElement"] << endl;
	logFile << "gaussPointsPerVolumeElement = "
			<< backGroundMeshInfo["gaussPointsPerVolumeElement"] << endl;
	logFile << "gaussPointsPerSurfaceElement = "
			<< backGroundMeshInfo["gaussPointsPerSurfaceElement"] << endl;
	logFile << "gaussPointsPerLineElement = "
			<< backGroundMeshInfo["gaussPointsPerLineElement"] << endl;
	logFile << "numberOfVolumeElements = "
			<< backGroundMeshInfo["numberOfVolumeElements"] << endl;
#endif

	/*********************************************************************/
	// Write the FEM reference mesh.
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	ofstream femMeshFile;

	if (rank == 0) {

		femMeshFile.open("fem.msh");
		femMeshFile.precision(12);
		femMeshFile.setf(ios_base::scientific, ios_base::floatfield);

		if ((int) backGroundMeshInfo["elemType"] == 1)
			femMeshFile << "MESH  dimension 3  ElemType Tetrahedra Nnode "
					<< (int) backGroundMeshInfo["nodesPerVolumeElement"]
					<< endl;

		else if ((int) backGroundMeshInfo["elemType"] == 2)
			femMeshFile << "MESH  dimension 3  ElemType Hexahedra Nnode "
					<< (int) backGroundMeshInfo["nodesPerVolumeElement"]
					<< endl;

		// Write nodes coordinates of reference mesh.
		femMeshFile << "Coordinates" << endl;

		for (int i = 0; i < nodesNum; i++)
			femMeshFile << i + 1 << " " << particles[i].getCoord(0) << " "
					<< particles[i].getCoord(1) << " "
					<< particles[i].getCoord(2) << endl;

		femMeshFile << "end coordinates\n" << endl;

		femMeshFile << "Elements" << endl;

		// Loop over all global elements.
		for (int i = 0; i < globalElemNum; i++) {
			int& elemType = nodesElements[i].getElemType();

			int& matID = allNodesElements[i].getMaterialID();
			intVector& nodes = nodesElements[i].getNodes();

			femMeshFile << i + 1 << " ";

			for (int j = 0; j < nodes.size(); j++)
				femMeshFile << nodes[j] << " ";

			femMeshFile << matID << endl;
		}

		femMeshFile << "end elements" << endl;
	}

#ifdef _FEdebugMode_
	logFile << "************** all unsorted elements *****************" << endl;
	for (int i = 0; i < globalElemNum; i++) {
		logFile << i << " ";
		intVector& nodes = allNodesElements[i].getNodes();
		int& matID = allNodesElements[i].getMaterialID();
		for (int j = 0; j < nodes.size(); j++)
			logFile << nodes[j] << " ";
		logFile << matID << endl;
	}
#endif

#ifdef _FEdebugMode_
	logFile << "**************** local elements **********************" << endl;
	for (int i = 0; i < nodesElements.size(); i++) {
		logFile << i << " ";
		intVector& nodes = nodesElements[i].getNodes();
		int& elemType = nodesElements[i].getElemType();
		int& matID = nodesElements[i].getMaterialID();
		for (int j = 0; j < nodes.size(); j++)
			logFile << nodes[j] << " ";
		logFile << " material = " << matID << " elemType = " << elemType
				<< endl;
	}
	logFile << "************************ nodes **********************" << endl;
	for (int i = 0; i < nodesNum; i++) {
		logFile << "newIdx " << i << " oldIdx " << particles[i].getID()
				<< " coords: " << particles[i].getCoord(0) << " "
				<< particles[i].getCoord(1) << " " << particles[i].getCoord(2)
				<< " material = " << particles[i].getMaterialID() << endl;
	}

#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Define the surface elements in the volume element list.
void FEMGeometryExt::combiningSurfaceVolumeElems(InputFileData* InputData,
		std::map<std::string, double>& modelData, std::ofstream& logFile) {

	using namespace std;

	int position, counter;

	//! Loop over each surface element to extract nodes
	for (int i = 0; i < surfaceNodesElems.size(); i++) {
		intVector& surfNodesVec = surfaceNodesElems[i].getNodes();

		//! Loop over each volume element to extract nodes
		for (int j = 0; j < nodesElements.size(); j++) {
			intVector& volNodesVec = nodesElements[j].getNodes();

			nodesElements[j].getSurfaceElems().push_back(j);

			//! Reset counter
			counter = 0;

			//! Find if the surface nodes belongs to the volume nodes set
			for (int k = 0; k < surfNodesVec.size(); k++) {

				//! Return the position of the surface node in the volume nodes set.
				//! If it is not found, the findIntVecPos return -1.
				position = findIntVecPos(surfNodesVec[k], 0, volNodesVec.size(),
						volNodesVec);

				//! Increase the counter if the surface node has been found.
				if (position > -1)
					counter++;
			}

			//! If all surface nodes have been found in the volume nodes set,
			//! then surface element belongs to volume.
			if (counter == surfNodesVec.size())
				nodesElements[j].getSurfaceElems().push_back(
						surfaceNodesElems[i].getGlobalID());
		}
	}

#ifdef _FEdebugMode_

	logFile << "###################################################" << endl;
	logFile << "**************** Volume-Surface info **************" << endl;

	for (int i = 0; i < nodesElements.size(); i++) {
		logFile << nodesElements[i].getGlobalID() << ":";
		for (int j = 0; j < nodesElements[i].getSurfaceElems().size(); j++) {
			logFile << " " << nodesElements[i].getSurfaceElems()[j];
		}
		logFile << endl;
	}

#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Define the surface elements in the volume element list.
intMatrix FEMGeometryExt::decompVolumeToSurfaceElems(intVector& volNodes,
		int& elemType, InputFileData* InputData,
		std::map<std::string, double>& modelData, std::ofstream& logFile) {

	using namespace std;

	intMatrix surfaceElems;

	switch (elemType) {
	case 1:
		surfaceElems =
				decompTetraToTrianElems(volNodes, InputData, modelData, logFile);
		break;

	case 2:
		surfaceElems =
				decompHexaToQuadElems(volNodes, InputData, modelData, logFile);
		break;

	default:
		cout
				<< "In FEMGeometryExt::decompVolumeToSurfaceElems: Element"
						" type is not yet supported"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return surfaceElems;
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Forms the triangular surface elements of a tetrahedra element.
intMatrix FEMGeometryExt::decompTetraToTrianElems(intVector& volNodes,
		InputFileData* InputData, std::map<std::string, double>& modelData,
		std::ofstream& logFile) {

	using namespace std;

	logFile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" <<endl;
	printVector(volNodes,"Using volume nodes:",logFile);

	//! Number of surfaces defining tetrahedra = 4
	//! Number nodes of quadrilateral = 3
	intMatrix surfaceElems(4, intVector(3));
	int counter = 0;

	for(int i=0; i<volNodes.size() ;i++){
		int Node1 = volNodes[i];
		for(int j=i+1; j<volNodes.size() ;j++){
			int Node2 = volNodes[j];
			for(int k=j+1;k<volNodes.size();k++){
				if(counter <5){
				surfaceElems[counter][0]=Node1;
				surfaceElems[counter][1]=Node2;
				surfaceElems[counter][2]=volNodes[k];
				logFile<<"Plane formed: "<< Node1 <<" "<< Node2 <<" "<< volNodes[k] <<endl;
				counter++;
				}
				else
					break;
			}
			if(counter>4)break;
		}
		if(counter>4)break;
	}


	return surfaceElems;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Forms the quadrilateral surface elements of a hexahedra element.
intMatrix FEMGeometryExt::decompHexaToQuadElems(intVector& volNodes,
		InputFileData* InputData, std::map<std::string, double>& modelData,
		std::ofstream& logFile) {

	using namespace std;

	//! Number of surfaces defining hexahedra = 6
	//! Number nodes of quadrilateral = 4
	intMatrix surfaceElems(6,intVector(4));

	//! Copy the reference faces in the surface elements list
	int count=0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 4; j++) {
			surfaceElems[i][j] = volNodes[count];
			count++;
		}
	}

	//! Check if volNodes data is appropriate
	if (volNodes.size() != 8) {
		logFile
				<< "In FEMGeometryExt::decompHexaToQuadElems: Too few or too many"
						" nodes defining the volume element"<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	intMatrix faceMatrix(2, intVector(2));
	int Numface = 2;

	for (int i = 0; i < surfaceElems[0].size(); i++) {

			//!
			if (i < surfaceElems[0].size()-1 ) {

				//! Record the face nodes in the surfaceElems list
				surfaceElems[Numface][0] = surfaceElems[0][i];
				surfaceElems[Numface][1] = surfaceElems[0][i+1];
				surfaceElems[Numface][2] = surfaceElems[1][i+1];
				surfaceElems[Numface][3] = surfaceElems[1][i];

				Numface++;

			}else if(i == surfaceElems[0].size()-1){
				surfaceElems[Numface][0] = surfaceElems[0][i];
				surfaceElems[Numface][1] = surfaceElems[0][0];
				surfaceElems[Numface][2] = surfaceElems[1][0];
				surfaceElems[Numface][3] = surfaceElems[1][i];
		}

	}

#ifdef _FEdebugMode_

	logFile << "###################################################" << endl;
	logFile << "**************** Hexahedra-Quad info **************" << endl;

	logFile << "Reference Nodes:" << endl;
	for (int i = 0; i < volNodes.size(); i++) {
			logFile << volNodes[i] << " " << endl;
	}

	logFile << "Reconstructed Faces:" << endl;
	for (int k = 0; k < surfaceElems.size(); k++) {
		for (int l = 0; l < surfaceElems[k].size(); l++) {
			logFile << surfaceElems[k][l] << " ";
		}
		logFile << endl;
	}

#endif

	return surfaceElems;

}

/*!****************************************************************************/
/*!****************************************************************************/
intVector FEMGeometryExt::surfaceIDGenerator(intMatrix surfaceNodesList,
       		InputFileData* InputData,std::map<std::string, double>& modelData,
       		std::ofstream& logFile){

	using namespace std;

	intVector surfaceIDList(surfaceNodesList.size());

	int position=0;

	for(int i=0;i<surfaceNodesList.size();i++){
		intVector& surfElem = surfaceNodesList[i];
		//printVector(surfElem,"surfElem",logFile);

		position = findSurfElemInSurfList(surfElem,InputData,modelData,logFile);

		if(position == -1)
			surfaceIDList[i] = storeSurfaceElem(surfaceNodesList[i],InputData,modelData,logFile);
		else
			surfaceIDList[i] = surfaceNodesElems[position].getGlobalID();
	}

	return surfaceIDList;

}

/*!****************************************************************************/
/*!****************************************************************************/
int FEMGeometryExt::findSurfElemInSurfList(intVector& surfaceNodes,
		InputFileData* InputData, std::map<std::string, double>& modelData,
		std::ofstream& logFile) {

	using namespace std;

	int position = -1;
	int counter;

	for (int j = 0; j < surfaceNodesElems.size(); j++) {
		intVector& surfListElem = surfaceNodesElems[j].getNodes();
		//printVector(surfListElem,"surfListElem",logFile);

		counter = 0;
		for (int k = 0; k < surfaceNodes.size(); k++) {

			position = findIntVecPos(surfaceNodes[k], 0, surfListElem.size(),
					surfListElem);

			if (position == -1)
				break;
			else
				counter++;

		}

		if (counter == surfListElem.size()) {
			position = j;
			break;
		}
	}

	return position;
}

/*!****************************************************************************/
/*!****************************************************************************/
int FEMGeometryExt::storeSurfaceElem(intVector& surfaceNodes,
       		InputFileData* InputData,std::map<std::string, double>& modelData,
       		std::ofstream& logFile){

	using namespace std;

	int usedDOF = (int) modelData["usedDegreesOfFreedom"];
	FEMElementExt newSurfaceElem(usedDOF);

	newSurfaceElem.setNodes(surfaceNodes);

//	newSurfaceElem.setNodes(surfaceNodes);

	if(surfaceNodesElems.size() == 0)
		newSurfaceElem.setGlobalID(1);
	else
		newSurfaceElem.setGlobalID(surfaceNodesElems[surfaceNodesElems.size()-1].getGlobalID()+1);

	surfaceNodesElems.push_back(newSurfaceElem);

	return surfaceNodesElems[surfaceNodesElems.size()-1].getGlobalID();
}

/*!****************************************************************************/
/*!****************************************************************************/
void FEMGeometryExt::FEMGeoDataSetup(InputFileData* InputData,
		std::map<std::string, double>& modelData, std::string& meshFileName,
		std::string& inputFileName, std::ofstream& logFile) {

	using namespace std;

	// Read mesh file.
	ifstream meshFile(meshFileName.c_str());

	if (!meshFile) {
		logFile << "Can't open FEM mesh file: " << meshFileName << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	FEMInputData = new InputFileData(inputFileName, logFile);

	modelData["usedDegreesOfFreedom"] = 3;
	modelData["usedDimensions"] = 3;
	modelData["integrationMethod"] = 1; // Gauss integration scheme.

	PetscViewer viewerMPI, viewerSEQ;
	FEMGeoData = new FEMGeometry(FEMInputData, modelData, meshFile, logFile,
			viewerMPI, viewerSEQ);

//	delete FEMInputData;

	vector<Particle>& ptcls = FEMGeoData->getNodesVec();
	vector<FEMElement>& volElems = FEMGeoData->getNodesElemsVec();
	vector<FEMElement>& surfElems = FEMGeoData->getSurfaceNodesElemsVec();

	particles.resize(ptcls.size(), ParticleExt(modelData["usedDOF"]));
	for (int i = 0; i < ptcls.size(); i++) {
		particles[i] = ParticleExt(ptcls[i]);
		particles[i].setID(i + 1);
	}
//	ptcls.clear();

	nodesElements.resize(volElems.size(), FEMElementExt(modelData["usedDOF"]));
	for (int i = 0; i < volElems.size(); i++) {
		nodesElements[i] = FEMElementExt(volElems[i]);
	}
//	volElems.clear();

//	surfaceNodesElems.resize(volElems.size(),FEMElementExt(modelData["usedDOF"]));
	for (int i = 0; i < surfElems.size(); i++) {
//		FEMElementExt myElem(modelData["usedDOF"]) =
//		surfaceNodesElems[i] = FEMElementExt(surfElems[i]);

		surfaceNodesElems.resize(volElems.size(), FEMElementExt(surfElems[i]));

	}
//	surfElems.clear();

	for (int i = 0; i < nodesElements.size(); i++) {

		intVector& nodeVec = nodesElements[i].getNodes();
		int& elemType = nodesElements[i].getElemType();

		intMatrix tempMat = decompVolumeToSurfaceElems(nodeVec, elemType,
				InputData, modelData, logFile);

		nodesElements[i].setSurfaceElems(
				surfaceIDGenerator(tempMat, InputData, modelData, logFile));
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Find Point in geometry
int FEMGeometryExt::findPointInGeometry(dbVector& pointCoord,
		InputFileData* InputData, std::map<std::string, double>& modelData,
		std::ofstream& logFile) {

	using namespace std;

#ifdef _FEdebugMode_PointInPoly_
	logFile <<"################################################"<<endl;
	printVector(pointCoord,"Main Coordinate",logFile);
	logFile <<"################################################"<<endl;
#endif

	int volumeElem = -1;
	bool isInside = false;

	//! Loop over all volumes element
	for (int i = 0; i < nodesElements.size(); i++) {

#ifdef _FEdebugMode_PointInPoly_

		int location = 0;

		logFile <<"----------------------------------" <<endl;
		logFile <<"Searching in volume element: " << i <<endl;
		logFile <<"----------------------------------" <<endl;

		logFile << "with the following nodes: " <<endl;;

		for(int j=0; j< nodesElements[i].getNodes().size(); j++){
			logFile << nodesElements[i].getNodes()[j] << ": ";

			location = findPtcleIDInPtcleList(nodesElements[i].getNodes()[j], logFile);
			for(int k=0;k<particles[location].getCoords().size();k++){
				logFile << particles[location].getCoord(k) << " ";
			}
			logFile << endl;
		}
		logFile << endl;
		// --------------------------------------------------------------------
#endif

		// find if pointCoord is found in volume element:nodesElements[i]
		isInside = findPointInVolumeElem(pointCoord,
				nodesElements[i].getSurfaceElems(), nodesElements[i].getNodes(),
				InputData, modelData, logFile);

		if(isInside == true){
			volumeElem = i;
			break;
		}
	}

	return volumeElem;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Find if a point belongs to a volume.
bool FEMGeometryExt::findPointInVolumeElem(dbVector& pointCoord,
		intVector& surfaceElems, intVector& volumeNodes,
		InputFileData* InputData, std::map<std::string, double>& modelData,
		std::ofstream& logFile) {

	using namespace std;

	int surfaceLocation;
	bool isInside = false;
	intVector surfaceNodes;

	dbMatrix surfacePointsCoords(surfaceNodes.size());
	dbMatrix projSurfacePointsCoords;
	int location;

	dbVector newOriginCoord;
	int refNodeID, refNodeLocation;
	dbVector projRefNodeCoord;

	dbVector projPointCoord;

	int counter = 0;

	// int axisCheck = 0	// For x-axis
	// int axisCheck = 1	// For y-axis
	int axisCheck = 2;		// For z-axis

	oPType tolCheck = InputData->getValue("pointInPolyTolerance");

#ifdef _FEdebugMode_PointInPoly_
	printVector(surfaceElems, "Considering elements", logFile);
#endif

	// loop over each surface element
	for (int i = 0; i < surfaceElems.size(); i++) {

		//! Extract the surface particles
		surfaceLocation = findSurfaceElemIDInList(surfaceElems[i], logFile);
		if (surfaceLocation == -1) {
			logFile << "In FEMGeometryExt::findPointInVolume: Surface element "
					"cannot be found in list " << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		surfaceNodes = surfaceNodesElems[surfaceLocation].getNodes();


		//! Extract the coordinates of the surface points
		surfacePointsCoords.clear();
		for (int j = 0; j < surfaceNodes.size(); j++) {

			location = findPtcleIDInPtcleList(surfaceNodes[j], logFile);
			if (location == -1) {
				logFile << "In FEMGeometryExt::findPointInVolumeElem: Particle not "
						"found in particles list" << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			surfacePointsCoords.push_back(particles[location].getCoords());
		}


		//! Define a new origin for surface: Usually, the first coordinate is chosen
		newOriginCoord.clear();
		for (int n = 0; n < surfacePointsCoords[0].size(); n++) {
			newOriginCoord.push_back(surfacePointsCoords[0][n]);
		}

		//! new surfacePoints Coordinates with a new origin
		dbMatrix surfacePointsCoords_(surfacePointsCoords.size(),
				dbVector(surfacePointsCoords[0].size()));
		for (int l = 0; l < surfacePointsCoords.size(); l++) {
			for (int m = 0; m < surfacePointsCoords[l].size(); m++) {
				// if condition used to prevent computer accurracy error
				// such as: 2-2 = 1e-16
				if(abs(surfacePointsCoords[l][m]
						- newOriginCoord[m]) < 1e-14){
					surfacePointsCoords_[l][m] = 0;
				}
				else{
					surfacePointsCoords_[l][m] = surfacePointsCoords[l][m]
						- newOriginCoord[m];
				}
			}
		}

#ifdef _FEdebugMode_PointInPoly_
		logFile << "For surface element: "
				<< surfaceNodesElems[surfaceLocation].getGlobalID() << endl;
		printVector(surfaceNodes, "With nodes:", logFile);
		printMatrix(surfacePointsCoords, "of Coordinates:", logFile);

		logFile << "For surface element (New Origin): "
				<< surfaceNodesElems[surfaceLocation].getGlobalID() << endl;
		printVector(surfaceNodes, "With nodes:", logFile);
		printMatrix(surfacePointsCoords_, "of Coordinates:", logFile);
#endif

		//! Find the local basis which is on the surface of the element
		dbMatrix localBasis = setupLocalBasis(surfacePointsCoords, InputData,
				modelData, logFile);


		//! Project surface coordinates on local basis
		//! Use for checking -> can be deleted later-on
		innerTensorProduct(localBasis, surfacePointsCoords_,
				projSurfacePointsCoords, false, true, logFile);


#ifdef _FEdebugMode_PointInPoly_
		printMatrix(localBasis, "Basis of surface", logFile);
		printMatrix(projSurfacePointsCoords,
				"Projected Surface Coordinates", logFile);
#endif

		//! Select a volume element node not on the surface
		for (int k = 0; k < volumeNodes.size(); k++) {

			findIntVecPos(volumeNodes[k], 0, surfaceNodes.size(), surfaceNodes);

			if (findIntVecPos(volumeNodes[k], 0, surfaceNodes.size(),
					surfaceNodes) == -1) {
				refNodeID = volumeNodes[k];
				break;
			}
		}


#ifdef _FEdebugMode_PointInPoly_
		logFile << "Selected reference node: " << refNodeID << endl;
#endif
		//! Extract reference node and project on local basis
		refNodeLocation = findPtcleIDInPtcleList(refNodeID, logFile);
		dbVector refNodeCoord_;
		if (refNodeLocation != -1) {
			dbVector refNodeCoord = particles[refNodeLocation].getCoords();
			for (int m = 0; m < refNodeCoord.size(); m++) {
				if(abs(refNodeCoord[m] - newOriginCoord[m]) < 1e-14)
					refNodeCoord_.push_back(0);
				else
					refNodeCoord_.push_back(refNodeCoord[m] - newOriginCoord[m]);
			}
			innerTensorProduct(localBasis, refNodeCoord_, projRefNodeCoord,
					false, logFile);
#ifdef _FEdebugMode_PointInPoly_
			printVector(particles[refNodeLocation].getCoords(),
					"REF-Point Coordinates:", logFile);
			printVector(refNodeCoord_, "REF-Point Coordinates(newOrg):",
					logFile);
			printVector(projRefNodeCoord, "Projected REF-Point Coordinates:",
					logFile);
#endif
		} else {
			logFile
					<< "In FEMGeometryExt::findPointInVolume: Reference point not found in list"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		//! Project unknown point on local basis
		dbVector pointCoord_(pointCoord.size());
		for (int m = 0; m < pointCoord.size(); m++) {
			if(abs(pointCoord[m] - newOriginCoord[m]) < 1e-14)
				pointCoord_[m] = 0;
			else
				pointCoord_[m] = pointCoord[m] - newOriginCoord[m];
		}
		innerTensorProduct(localBasis, pointCoord_, projPointCoord, false,
				logFile);

#ifdef _FEdebugMode_PointInPoly_
		printVector(pointCoord, "Point Coordinates:", logFile);
		printVector(pointCoord_, "Point Coordinates(newOrg):", logFile);
		printVector(projPointCoord, "Projected Point Coordinates:",	logFile);
#endif


		//! Check if the sign of the 3rd axis(Normal to surface) is the same
		//! for both refNode and pointCoord
		int temp = projRefNodeCoord[axisCheck];
//		if ((projRefNodeCoord[axisCheck] > 0 && projPointCoord[axisCheck] > 0) ||
//				(projRefNodeCoord[axisCheck] < 0 && projPointCoord[axisCheck] < 0) ||
//				projPointCoord[axisCheck] == 0){

		if ((projRefNodeCoord[axisCheck] > 0 && projPointCoord[axisCheck] > -tolCheck)
				|| (projRefNodeCoord[axisCheck] < 0
						&& projPointCoord[axisCheck] < tolCheck)
				|| projPointCoord[axisCheck] == 0) {

			counter++;

#ifdef _FEdebugMode_PointInPoly_
			logFile << " --> Point maybe inside <--" << endl;
			logFile << "Counter Value: " << counter << endl;
#endif
		}else
			break;
		// MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (counter == surfaceElems.size()) {
		isInside = true;

#ifdef _FEdebugMode_PointInPoly_
		logFile << "--->>> Point is inside !! <<<---" << endl;
#endif

	}

	return isInside;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Find Particle ID in particles list.
int FEMGeometryExt::findPtcleIDInPtcleList(int& ptcleID, std::ofstream& logFile) {

	using namespace std;

	int position = -1;

	for (int i = 0; i < particles.size(); i++) {
		if (ptcleID == particles[i].getID()) {
			position = i;
			break;
		}
	}

	return position;
}

/*!****************************************************************************/
/*!****************************************************************************/
//!  Find surface ID in surfaceNodesElems list.
int FEMGeometryExt::findSurfaceElemIDInList(int& surfaceID,
		std::ofstream& logFile) {

	using namespace std;

	int position = -1;

	for (int i = 0; i < surfaceNodesElems.size(); i++) {
		if (surfaceID == surfaceNodesElems[i].getGlobalID()) {
			position = i;
			break;
		}
	}

	return position;
}

/*!****************************************************************************/
/*!****************************************************************************/
//!  Find volume ID in nodesElements list.
int FEMGeometryExt::findVolumeElemIDInVolList(int& volumeID,
		std::ofstream& logFile) {

	using namespace std;

	int position = -1;

	for (int i = 0; i < nodesElements.size(); i++) {
		if (volumeID == nodesElements[i].getGlobalID()) {
			position = i;
			break;
		}
	}

	return position;
}

/*!****************************************************************************/
/*!****************************************************************************/
//!  Find the volume elements neighbouring a particular volume element.
intVector FEMGeometryExt::findAdjacentVolElemsOfVolumeElement(int& volID,
		std::ofstream& logFile) {

	// Find volume position in vector of volume element list
	int sVolume = this->findVolumeElemIDInVolList(volID, logFile);

	// Extract the nodes of the volume element
	intVector& sVolumeNodes = nodesElements[sVolume].getNodes();

	// For each node, get the list of volume elements to which they belong to
	// and set up list of all the volume elements.
	intVector adjacentElements;
	bool isRepeated;
	for (int i = 0; i < sVolumeNodes.size(); i++) {

		intVector elements = particles[sVolumeNodes[i]].getElems();

		for (int j = 0; j < elements.size(); j++) {
			isRepeated = false;

			for (int k = 0; k < adjacentElements.size(); k++) {
				if (elements[j] == adjacentElements[k]) {
					isRepeated = true;
					break;
				}
			}

			if (isRepeated == false && elements[j] != volID) {
				resizeArray(adjacentElements, adjacentElements.size() + 1);
				adjacentElements[adjacentElements.size()] = elements[j];
			}
		}
	}

	return adjacentElements;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Setting up the local basis of the surface plane
dbMatrix FEMGeometryExt::setupLocalBasis(dbMatrix& surfacePointsCoords,
		InputFileData* InputData, std::map<std::string, double>& modelData,
		std::ofstream& logFile){

	using namespace std;

	dbVector normalSurfaceVec, parallelSurfaceVec, orthoSurfaceVec;

	//!---------------------------------------------------------------------
	//! Calculating the normal surface vector
	normalSurfaceVec = caclSurfaceNormal(surfacePointsCoords,InputData,logFile);
	dbVector normalSurfaceVec_ = normaliseVec(normalSurfaceVec,logFile);


	//!---------------------------------------------------------------------
	//! Finding a vector parallel to surface plane
	parallelSurfaceVec.resize(normalSurfaceVec.size());
	for(int j=0;j<surfacePointsCoords[0].size();j++){
		if(abs(surfacePointsCoords[1][j] - surfacePointsCoords[0][j])<1e14)
			parallelSurfaceVec[j] = 0;
		else
			parallelSurfaceVec[j] = surfacePointsCoords[1][j] - surfacePointsCoords[0][j];
	}
	dbVector parallelSurfaceVec_ = normaliseVec(parallelSurfaceVec,logFile);

	//!---------------------------------------------------------------------
	//! Finding a vector parallel to surface plane and orthogonal
	//! to parallelSurfaceVec
	crossProduct(normalSurfaceVec,parallelSurfaceVec,orthoSurfaceVec);
	dbVector orthoSurfaceVec_ = normaliseVec(orthoSurfaceVec,logFile);


	//!---------------------------------------------------------------------
	//! Setting up the basis matrix
	dbMatrix basisMatrix(normalSurfaceVec.size(),dbVector(normalSurfaceVec.size()));
	for(int k=0; k<normalSurfaceVec.size();k++){
		basisMatrix[0][k] = orthoSurfaceVec_[k];
		basisMatrix[1][k] = parallelSurfaceVec_[k];
		basisMatrix[2][k] = normalSurfaceVec_[k];
	}

	return basisMatrix;
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the surface normal from a set of nodal coordinates
dbVector FEMGeometryExt::caclSurfaceNormal(dbMatrix& surfacePointsCoords,
		InputFileData* InputData,
		std::ofstream& logFile){

	using namespace std;

	if(surfacePointsCoords.size() < 3){
		logFile <<"In FEMGeometryExt::caclSurfaceNormal: Too few coordinates"
				" to compute surface normal"<<endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	dbVector vecOne(surfacePointsCoords[0].size());
	dbVector vecTwo(surfacePointsCoords[0].size());

	for(int i=0; i<surfacePointsCoords[0].size();i++){

		if(abs(surfacePointsCoords[1][i] - surfacePointsCoords[0][i])<1e-14)
			vecOne[i] = 0;
		else
			vecOne[i] = surfacePointsCoords[1][i] - surfacePointsCoords[0][i];

		if(abs(surfacePointsCoords[2][i] - surfacePointsCoords[0][i])<1e-14)
			vecTwo[i] = 0;
		else
			vecTwo[i] = surfacePointsCoords[2][i] - surfacePointsCoords[0][i];
	}

	dbVector surfaceNormalVec;
	crossProduct(vecOne,vecTwo,surfaceNormalVec);

	return surfaceNormalVec;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Normalise any given vector
dbVector FEMGeometryExt::normaliseVec(dbVector& vec, std::ofstream& logFile){

	using namespace std;

	double normVal = computeNorm(vec,2,logFile);

	dbVector aVec(vec.size());
	for(int i=0; i<vec.size();i++){
		aVec[i] = vec[i]/normVal;
	}

	return aVec;
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Create a GiD mesh file
void FEMGeometryExt::writeMeshFile(InputFileData* InputData,
		std::ofstream& logFile){

	using namespace std;

	int globalElemNum = nodesElements.size();
	int nodesNum = particles.size();

	std::map<std::string, double>& backGroundMeshInfo =
						InputData->getBackGroundMeshInfo();

	/*********************************************************************/
		// Write the FEM reference mesh.
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		ofstream femMeshFile;

		if (rank == 0) {

			femMeshFile.open("fem2.msh");
			femMeshFile.precision(12);
			femMeshFile.setf(ios_base::scientific, ios_base::floatfield);

			if ((int) backGroundMeshInfo["elemType"] == 1)
				femMeshFile << "MESH  dimension 3  ElemType Tetrahedra Nnode "
						<< (int) backGroundMeshInfo["nodesPerVolumeElement"]
						<< endl;

			else if ((int) backGroundMeshInfo["elemType"] == 2)
				femMeshFile << "MESH  dimension 3  ElemType Hexahedra Nnode "
						<< (int) backGroundMeshInfo["nodesPerVolumeElement"]
						<< endl;

			// Write nodes coordinates of reference mesh.
			femMeshFile << "Coordinates" << endl;

			for (int i = 0; i < nodesNum; i++)
				femMeshFile << i + 1 << " " << particles[i].getCoord(0) << " "
						<< particles[i].getCoord(1) << " "
						<< particles[i].getCoord(2) << endl;

			femMeshFile << "end coordinates\n" << endl;

			femMeshFile << "Elements" << endl;

			// Loop over all global elements.
			for (int i = 0; i < globalElemNum; i++) {
				int& elemType = nodesElements[i].getElemType();

				int& matID = nodesElements[i].getMaterialID();
				intVector& nodes = nodesElements[i].getNodes();

				femMeshFile << i + 1 << " ";

				for (int j = 0; j < nodes.size(); j++)
					femMeshFile << nodes[j] << " ";

				femMeshFile << matID << endl;
			}

			femMeshFile << "end elements" << endl;
		}
}

void FEMGeometryExt::printVolumePtclsDetails(InputFileData* InputData,
			std::ofstream& logFile){

	using namespace std;

	logFile << "********** In FEMGeometryExt::printVolumePtclsDetails **********" << endl;

	for(int i = 0; i < nodesElements.size(); i++){
		logFile << "Volume[" << i << "] -> ";
		intVector& volPtcls = nodesElements[i].getNodes();
		for(int j = 0; j < volPtcls.size(); j++){
			logFile << volPtcls[j] << "(";
			dbVector& volPtclsCoords = particles[volPtcls[j]-1].getCoords();
			for(int k = 0; k < volPtclsCoords.size(); k++){
				logFile << volPtclsCoords[k] << ",";
			}
			logFile <<") | ";
		}
		logFile << endl;
	}

}
