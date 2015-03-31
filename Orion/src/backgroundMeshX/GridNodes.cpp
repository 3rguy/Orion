#include "GridNodes.h"

/*!****************************************************************************/
/*!****************************************************************************/
// Read a complete Finite Elements's mesh from the mesh file.
//GridNodes::GridNodes(std::vector<Data>& iDataList, InputFileData* InputData,
//		ofstream& logFile):gridGeometry(NULL),MeshlessData(NULL) {
//
//	std::map<std::string, oPType> modelData;
//	// modelData["usedDegreesOfFreedom"] = 3;
//	InputData->setValue("usedDegreesOfFreedom", 3);
//
//	nDims = InputData->getValue("dimension"); //temp default value
//
//	// Define all sample nodes
//	logFile << "***** Setting up Sample Nodes *****" << endl << endl;
//	setGridNodes(InputData, logFile);
//	logFile << "#######################################################" << endl;
//
//	for (int i = 0; i < iDataList.size(); i++) {
//		iDataList[i].readMeshDataFile(InputData, logFile);
//		calcInterpolants(iDataList[i], InputData, logFile);
//		interpolateDisp(iDataList[i], InputData, logFile);
//
//		printMatrix(iDataList[i].getDisplacement(),
//				"TRANSFORMED DISPLACEMENT MATRIX", logFile);
//
//		logFile << "#######################################################"
//				<< endl;
//
//		if (i == 1) {
//			cout << "The End" << endl;
//			logFile << "The End" << endl;
//			MPI_Abort(MPI_COMM_WORLD, 1);
//		}
//	}
//}

/*!****************************************************************************/
/*!****************************************************************************/
GridNodes::GridNodes(Data& iData, InputFileData* InputData, ofstream& logFile)
:gridGeometry(NULL),MeshlessData(NULL){

	logFile << "************************" << endl;
	logFile << "GridNodes::GridNodes(Single Data)" << endl;
	logFile << "************************" << endl;
	logFile << endl;

	std::map<std::string, oPType> modelData;

	nDims = InputData->getValue("usedDegreesOfFreedom"); //temp default value

	// Define all sample nodes
	logFile << "***** Setting up Sample Nodes *****" << endl << endl;
	setGridNodes(InputData, logFile);

	logFile << "GridNodes::GridNodes(Data& iData, InputFileData* InputData, ofstream& logFile);"
			":gridGeometry(NULL),MeshlessData(NULL) is not ready yet !!" << endl;
	cout << "GridNodes::GridNodes(Data& iData, InputFileData* InputData, ofstream& logFile):"
			"gridGeometry(NULL),MeshlessData(NULL) is not ready yet !!" << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);

}

/*!****************************************************************************/
/*!****************************************************************************/
GridNodes::GridNodes(InputFileData* InputData, ofstream& logFile)
:gridGeometry(NULL),MeshlessData(NULL) {

//	std::map<std::string,double> modelData;
//	modelData["integrationMethod"] = 1; // Gauss integration scheme.
//	modelData["usedDegreesOfFreedom"] = 3;
//	modelData["usedDimensions"] = 3;
//
//	std::string folderName = "gridNodes/";
//	std::string meshFileName = folderName + "mesh.dat";
//	gridGeometry =
//			new FEMGeometryExt(InputData, modelData, meshFileName, logFile);

	nDims = InputData->getValue("usedDegreesOfFreedom");
	setGridNodes(InputData, logFile);
//	saveNodesToFile(logFile);
}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::addNode(Node myNode) {

	if (nodalList.size() == 0)
		myNode.setId(0);
	else
		myNode.setId(nodalList[nodalList.size() - 1].getId() + 1);

	nodalList.push_back(myNode);

}

/*!****************************************************************************/
/*!****************************************************************************/
int GridNodes::findNode(int val) {

	int location = -1;

	for (int i = 0; i < nodalList.size(); i++) {

		if (nodalList[i].getId() == val) {

			location = i;
			break;
		}
	}
	if (location == -1) {
		cout << endl << "ERROR:Cannot find Data Id[" << val << "] in nodalList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return location;
}

/*!****************************************************************************/
/*!****************************************************************************/
// Read a complete Finite Elements's mesh from the mesh file.
void GridNodes::setGridNodes(InputFileData* InputData, ofstream& logFile) {

	int choice = InputData->getValue("gridNodesType");

	switch (choice) {
	case 1:
		generateGridNodes(InputData,logFile);
		break;

	case 2:
	{

		std::map<std::string, double> modelData;
		modelData["integrationMethod"] = 1; // Gauss integration scheme.
		modelData["usedDegreesOfFreedom"] = 3;
		modelData["usedDimensions"] = 3;

		std::string folderName = "gridNodes/";
		std::string meshFileName = folderName + "mesh.dat";
		std::string inputFileName = folderName + "input.dat";

		gridGeometry = new FEMGeometryExt(InputData, modelData, meshFileName,
				inputFileName, logFile);

		importGridNodes(InputData, logFile);

//		cout << "Num of Gauss Points: " << gridGeometry->getFEMGeoData()->getGaussPointsVec().size() << endl;

//		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}
	default:
		logFile << "In GridNodes::setGridNodes, choice of"
				" gridNodesType is not valid" << endl;
		cout << "In GridNodes::setGridNodes, choice of"
				" gridNodesType is not valid" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::generateGridNodes(InputFileData* InputData, ofstream& logFile){

	// Defining domain (L,B,D)
		oPType maxX, minX, maxY, minY, maxZ, minZ;
		maxX = InputData->getValue("GridNodeMaxX");
		minX = InputData->getValue("GridNodeMinX");

		maxY = InputData->getValue("GridNodeMaxY");
		minY = InputData->getValue("GridNodeMinY");

		maxZ = InputData->getValue("GridNodeMaxZ");
		minZ = InputData->getValue("GridNodeMinZ");


		// Uniformly distributed nodes configuration
		oPType nodalSpacing = InputData->getValue("gridNodesSpacing");

		// Find the number of points along each axis
		int nPointsX = static_cast<int>(round((maxX - minX) / nodalSpacing));
		int nPointsY = static_cast<int>(round((maxY - minY) / nodalSpacing));
		int nPointsZ = static_cast<int>(round((maxZ - minZ) / nodalSpacing));

#ifdef _GridNodesDebugMode_
	logFile << "Setting up grid with the following dimensions:" << endl;
	logFile << "Length:  " << minX << " < L < " << maxX << endl;
	logFile << "Breadth: " << minY << " < B < " << maxY << endl;
	logFile << "Depth:   " << minZ << " < D < " << maxZ << endl;
	logFile << "Num of Points: X[" << nPointsX << "] Y[" << nPointsY << "] Z["
			<< nPointsZ << "]" << endl;
#endif

		dbVector coord;
		resizeArray(coord, nDims);

		int iDCounter = 0;
		// Compute the coordinates of each sample point
		for (int i = 0; i < nPointsX; i++) {
			for (int j = 0; j < nPointsY; j++) {
				for (int k = 0; k < nPointsZ; k++) {
					coord[0] = minX + (i * nodalSpacing);
					coord[1] = minY + (j * nodalSpacing);
					coord[2] = minZ + (k * nodalSpacing);
					Node myNode(coord);
					myNode.setId(iDCounter);
					addNode(myNode);
					iDCounter++;
				}
			}
		}

#ifdef _GridNodesDebugMode_
	logFile << "************ GridNodes::generateGridNodes ************" << endl;
	logFile << "Total Num of Points: "<< nodalList.size() << "]" << endl;
	logFile << "Num of Points: X[" << nPointsX << "] Y[" << nPointsY << "] Z["
			<< nPointsZ << "]" << endl;
	logFile << "The coordinates are:" << endl;
	for (int id = 0; id < nodalList.size(); id++) {
		logFile << "ID[" << nodalList[id].getId() << "]-> Coord("
				<< nodalList[id].getCoords()[0] << ","
				<< nodalList[id].getCoords()[1] << ","
				<< nodalList[id].getCoords()[2] << ")" << endl;
	}
	logFile << endl;
#endif


}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::importGridNodes(InputFileData* InputData, ofstream& logFile){

	using namespace std;

	this->clearNodalList();
	vector<ParticleExt>& ptclsList =  gridGeometry->getNodesVec();

	int iDCounter = 0;
	for(int i = 0; i < ptclsList.size(); i++){
		Node myNode(ptclsList[i].getCoords());
		myNode.setId(iDCounter);
		addNode(myNode);
		iDCounter++;
	}

#ifdef _GridNodesDebugMode_
	logFile << "************ GridNodes::importGridNodes ************" << endl;
	logFile << "Total Num of Points: "<< nodalList.size() << "]" << endl;
	logFile << "The coordinates are:" << endl;
	for (int id = 0; id < nodalList.size(); id++) {
		logFile << "ID[" << nodalList[id].getId() << "]-> Coord("
				<< nodalList[id].getCoords()[0] << ","
				<< nodalList[id].getCoords()[1] << ","
				<< nodalList[id].getCoords()[2] << ")" << endl;
	}
	logFile << endl;
#endif


}

/*!****************************************************************************/
/*!****************************************************************************/
// Read a complete Finite Elements's mesh from the mesh file.
//void GridNodes::calcInterpolants(Data& iData, InputFileData* InputData,
//		ofstream& logFile) {
//
//	std::vector<ParticleExt> ptcls;
//
//	dbMatrix particleCoord = iData.getCoords();
//
//	logFile << "Calculating usedDOF: All DOF[" << iData.getDisplacement().size()
//			<< "]/ All Coord[" << iData.getCoords().size() << "]= " << endl;
//	int usedDOF = iData.getDisplacement().size() / iData.getCoords().size();
//	logFile << usedDOF << endl << endl;
//
//	oPType gamma = 1.8;
//	oPType beta = (particleCoord[particleCoord.size() - 1][2]
//			- particleCoord[particleCoord.size() - 2][2]) / gamma;
//
//	logFile << "Setting up the particles list:" << endl;
//	for (int i = 0; i < particleCoord.size(); i++) {
//		ParticleExt myParticle(usedDOF);
//		myParticle.setID(i);
//		myParticle.setCoords(particleCoord[i][0], particleCoord[i][1],
//				particleCoord[i][2]);
//		//myParticle.setRadii(inflRadius,inflRadius,inflRadius);
//		myParticle.setBeta(beta);
//		ptcls.push_back(myParticle);
//		logFile << "Particle[" << myParticle.getID() << "]->Coord("
//				<< myParticle.getCoords()[0] << "," << myParticle.getCoords()[1]
//				<< "," << myParticle.getCoords()[2] << ") ADDED" << endl;
//	}
//
//#ifdef _GridNodesDebugMode_
//	logFile << "Setting up the particles list:" << endl;
//	for (int id = 0; id < ptcls.size(); id++) {
//		logFile << "Particle[" << ptcls[id].getID() << "]->Coord("
//				<< ptcls[id].getCoords()[0] << "," << ptcls[id].getCoords()[1]
//				<< "," << ptcls[id].getCoords()[2] << ") ADDED" << endl;
//	}
//#endif
//
//	intVector sPtcls;
//	dbVector interpolants;
//	int basisTerm = 1.0;
//	// Find the neighbours and calculate the interpolants of each Node
//	for (int j = 0; j < nodalList.size(); j++) {
//
//		findSuppParticles(nodalList[j].getCoords(), ptcls, sPtcls, logFile);
//
//#ifdef _GridNodesDebugMode_
//		logFile << "Node[" << nodalList[j].getId() << "](";
//		for (int jd = 0; jd < nodalList[j].getCoords().size(); jd++) {
//			logFile << nodalList[j].getCoords()[jd];
//			if (jd < nodalList[j].getCoords().size() - 1)
//				logFile << ",";
//		}
//		logFile << "):";
//
//		for (int jd = 0; jd < sPtcls.size(); jd++) {
//			logFile << sPtcls[jd] << "(";
//			for (int jd_ = 0; jd_ < ptcls.size(); jd_++) {
//				if (ptcls[jd_].getID() == sPtcls[jd]) {
//					for (int jdc = 0; jdc < ptcls[jd_].getCoords().size();
//							jdc++) {
//						logFile << ptcls[jd_].getCoord(jdc);
//						if (jdc < ptcls[jd_].getCoords().size() - 1)
//							logFile << ",";
//					}
//					logFile << ")|";
//					break;
//				}
//
//			}
//		}
//		logFile << endl;
//#endif
//
//		nodalList[j].setSPtcls(sPtcls);
//		int sPtclsSize = sPtcls.size();
//
//		/*EFGShapeFunc::calcShapes(InputData,sPtclsSize,sPtcls,ptcls,
//		 nodalList[j].getCoords()[0],nodalList[j].getCoords()[1],
//		 nodalList[j].getCoords()[2],interpolants,basisTerm,logFile);*/
//
//		/*MaxEntShapeFunc A(InputData,logFile);
//		 A.calcShapes(InputData,sPtcls,ptcls,
//		 nodalList[j].getCoords()[0],nodalList[j].getCoords()[1],
//		 nodalList[j].getCoords()[2],interpolants,logFile);*/
//
//		oPType interpolantsVal = 1.0 / sPtcls.size();
//		for (int k = 0; k < sPtcls.size(); k++) {
//			interpolants.push_back(interpolantsVal);
//		}
//
//		nodalList[j].setInterpolants(interpolants);
//
//		logFile << "The interpolants are:";
//		for (int l = 0; l < interpolants.size(); l++) {
//			logFile << interpolants[l] << ", ";
//		}
//		logFile << endl;
//
//		sPtcls.clear();
//		interpolants.clear();
//
//	}
//	//MPI_Abort(MPI_COMM_WORLD,1);
//
//}

/*!****************************************************************************/
/*!****************************************************************************/
// Read a complete Finite Elements's mesh from the mesh file.
void GridNodes::findSuppParticles(dbVector& coord,
		std::vector<ParticleExt>& ptclList, intVector& sPtcls, ofstream& logFile) {

	dbMatrix coordRange;
	resizeArray(coordRange, nDims, 2);
	int rad = 10;

	for (int i = 0; i < nDims; i++) {
		coordRange[i][0] = coord[i] - rad; // min of range
		coordRange[i][1] = coord[i] + rad; // max of range
		logFile << "Coordinate Range of Dim[" << i << "]: " << coordRange[i][0]
				<< " - " << coordRange[i][1] << endl;
	}

	int counter;
	for (int j = 0; j < ptclList.size(); j++) {
		dbVector ptclCoord = ptclList[j].getCoords();
		counter = 0;
		logFile << "Comparing: ";
		for (int k = 0; k < ptclCoord.size(); k++) {
			logFile << coordRange[k][0] << "<" << ptclCoord[k] << "<"
					<< coordRange[k][1] << " ?";
			if (ptclCoord[k] > coordRange[k][0]
					&& ptclCoord[k] <= coordRange[k][1]) {
				logFile << "YES!";
				counter++;
			}
		}
		logFile << "=>Counter[" << counter << "]";
		if (counter == nDims) {
			logFile << " *SELECTED*" << endl;
			resizeArray(sPtcls, sPtcls.size() + 1);
			sPtcls[sPtcls.size() - 1] = ptclList[j].getID();
		} else
			logFile << endl;
	}
	logFile << endl;
}

/*!****************************************************************************/
/*!****************************************************************************/
//void GridNodes::interpolateDisp(Data& iData, InputFileData* InputData,
//		ofstream& logFile) {
//
//	dbMatrix iDataDisp = iData.getDisplacement();
//
//	int nDof = iData.getDisplacement().size() / iData.getCoords().size();
//
//#ifdef _GridNodesDebugMode_
//	logFile << "Number of DOF: " << nDof << endl << endl;
//#endif
//
//	dbMatrix nodalDisp;
//	resizeArray(nodalDisp, (nodalList.size() * nDof), iDataDisp[0].size());
//	logFile << "Size of nodalDisp: " << nodalDisp.size() << " x "
//			<< nodalDisp[0].size() << endl;
//
//	for (int i = 0; i < nodalList.size(); i++) {
//
//		intVector sPtcls = nodalList[i].getSPtlcs();
//		dbVector interPlnts = nodalList[i].getInterpolants();
//
//#ifdef _GridNodesDebugMode_
//		logFile << "The neighbouring particles along with their interpolants are:"
//				<< endl;
//		if (sPtcls.size() == interPlnts.size()) {
//			for (int m = 0; m < sPtcls.size(); m++) {
//				logFile << " |" << sPtcls[m] << "->" << interPlnts[m];
//			}
//		}
//		logFile << endl << endl;
//#endif
//		logFile << "Interpolating ..." << endl;
//		for (int j = 0; j < iDataDisp[0].size(); j++) {
//			for (int k = 0; k < nDof; k++) {
//				for (int l = 0; l < sPtcls.size(); l++) {
//					nodalDisp[(i * nDof) + k][j] +=
//							iDataDisp[(sPtcls[l] * nDof) + k][j]
//									* interPlnts[l];
//#ifdef _GridNodesDebugMode_
//					logFile << iDataDisp[(sPtcls[l] * nDof) + k][j] << " * "
//							<< interPlnts[l] << " = "
//							<< iDataDisp[(sPtcls[l] * nDof) + k][j]
//									* interPlnts[l] << "("
//							<< nodalDisp[(i * nDof) + k][j] << ") |:|";
//#endif
//				}
//#ifdef _GridNodesDebugMode_
//				logFile << endl;
//#endif
//			}
//		}
//#ifdef _GridNodesDebugMode_
//		logFile << "---------------------------------------------------------\n"
//				<< endl;
//#endif
//	}
//	iData.getDisplacement().clear();
//	iData.setDisplacement(nodalDisp);
//
//}

/*!****************************************************************************/
/*!****************************************************************************/
//dbMatrix GridNodes::calcResultOnGridPoint(Data& iData, InputFileData* InputData,
//		ofstream& logFile) {
//
//
////	setSupportingParticles(iData, InputData, logFile);
//	setSupportingParticles_(iData.getMeshData(),nodalList,InputData,logFile);
//
//	setInterpolantsOnNodes(iData, InputData, logFile);
//
//	interpolateNodalResult(iData, InputData, logFile);
//
//	return assembleNodalResultMatrix(iData, InputData, logFile);
//
//}

/*!****************************************************************************/
/*!****************************************************************************/
//dbMatrix GridNodes::initCalcResultOnParticles(Data& iData,
//		InputFileData* InputData, ofstream& logFile) {
void GridNodes::initCalcResultOnParticles(Data& iData,
			InputFileData* InputData, ofstream& logFile) {

	FEMGeometryExt* gFEMGeometry = gridGeometry;
	gridGeometry = iData.getMeshData();

	importGridNodes(InputData,logFile);

	setSupportingParticles_two(gFEMGeometry,nodalList,InputData,logFile);

	setInterpolantsOnNodes(gFEMGeometry, InputData, logFile);

//	return setInterpolantsOnParticles(iData, InputData, logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
dbMatrix GridNodes::calcResultOnParticles(Data& iData, dbMatrix& interpolantsList,
		InputFileData* InputData, ofstream& logFile) {

	interpolateParticleResult(iData, interpolantsList, InputData, logFile);

	return assemblePtclResultMatrix(iData, InputData, logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::interpolantSetup(Data& iData, InputFileData* InputData,
		ofstream& logFile) {


	setSupportingParticles_two(iData.getMeshData(),nodalList,InputData,logFile);

	setInterpolantsOnNodes(iData.getMeshData(), InputData, logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
dbMatrix GridNodes::interpolateResultOnGridPoint(Data& iData,
		InputFileData* InputData,ofstream& logFile) {

	interpolateNodalResult(iData, InputData, logFile);

	return assembleNodalResultMatrix(iData, InputData, logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::setSupportingParticles(Data& iData, InputFileData* InputData,
		ofstream& logFile) {

	std::map<std::string, oPType> modelData; // To satisfy findPointInGeometry input argument

	int volumeElem = 0;
	intVector nodesOutsideGeo;

	logFile << "Total number of grid nodes: " << nodalList.size() << endl;
	for (int i = 0; i < nodalList.size(); i++) {

		volumeElem = iData.getMeshData()->findPointInGeometry(
				nodalList[i].getCoords(), InputData, modelData, logFile);

		nodalList[i].setVolumeElement(volumeElem);

		if (volumeElem != -1) {

#ifdef _GridNodesDebugMode_
			logFile << "Setting supporting particles ..." << endl;
			logFile << "initial number of sPtcls at node: " << nodalList[i].getSPtlcs().size() << endl;
			logFile << "Num of sPtcls to be set: " << iData.getMeshData()->getNodesElemsVec()[volumeElem].getNodes().size() << endl;
			logFile << "With particles: ";
			for(int b = 0; b < iData.getMeshData()->getNodesElemsVec()[volumeElem].getNodes().size(); b++){
						logFile << iData.getMeshData()->getNodesElemsVec()[volumeElem].getNodes()[b] << " ";
					}
			logFile << endl;
#endif

			nodalList[i].setSPtcls(
					iData.getMeshData()->getNodesElemsVec()[volumeElem].getNodes());

#ifdef _GridNodesDebugMode_
			logFile << "Number of sPtcls set: " << nodalList[i].getSPtlcs().size() << endl;
#endif
		}
		else{
			nodesOutsideGeo.push_back(i);
		}
	}

#ifdef _GridNodesDebugMode_
	for (int j = 0; j < nodalList.size(); j++) {
		logFile << "Node[" << j << "] Supporting ptcls("
				<< nodalList[j].getSPtlcs().size() << "): ";
		if(nodalList[j].getSPtlcs().size() > 0 ){
			for(int k=0; k<nodalList[j].getSPtlcs().size(); k++)
				logFile << nodalList[j].getSPtlcs()[k] << " ";
		}
		else{
			logFile << "NULL" << endl;
		}
		logFile << endl;
	}

	logFile << endl << "Nodes laying outside geometry(Num = " << nodesOutsideGeo.size() <<") : ";
	for(int k=0; k<nodesOutsideGeo.size(); k++){
		logFile << nodesOutsideGeo[k] << ", ";
	}
	logFile << endl;
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::setSupportingParticles_two(FEMGeometryExt* FEMDataExt,
		std::vector<Node>& nodesVec, InputFileData* InputData,
		ofstream& logFile) {

#ifdef _GridNodesDebugMode_
	logFile << "GridNodes::setSupportingParticles_two" << endl;
	FEMDataExt->printVolumePtclsDetails(InputData, logFile);
#endif

	intVector nodesOutsideGeo, nodesInsideGeo;

//	std::string fileName = "mesh_N" +
//			convertIntToString(FEMDataExt->getFEMGeoData()->getNodesVec().size());
//	FEMDataExt->saveToGidMeshFile(fileName,InputData,logFile);

	logFile << "Total number of grid nodes: " << nodalList.size() << endl;
	for (int i = 0; i < nodesVec.size(); i++) {
		dbVector& nodeCoord = nodesVec[i].getCoords();

		intVector sVols;
		intVector sPtcls = FEMDataExt->findSupportingPtcls(nodeCoord, sVols,
				InputData, logFile);

#ifdef _GridNodesDebugMode_
		logFile << "-----------------------------------------" << endl;
		logFile << "For point[" << i << "]: " << nodeCoord[0] << " "
		<< nodeCoord[1] << " " << nodeCoord[2] << " " << endl;

		printVector(sPtcls,"sPtcls",logFile);
		printVector(sVols,"sVols",logFile);
		for(int j = 0; j < sVols.size(); j++) sVols[j]++;
		printVector(sVols,"sVols",logFile);
#endif

		if (sPtcls.size() == 0) {
			nodesOutsideGeo.push_back(i);
//			printVector(nodesOutsideGeo,"nodesOutsideGeo",logFile);
		} else {
			nodesInsideGeo.push_back(i);
//			nodesVec[nodesInsideGeo[i]].getSPtlcs() = sPtcls;
			nodesVec[i].getSPtlcs() = sPtcls;
		}
	}

	int activateFilter = InputData->getValue("filterSupportingParticles");
	if (activateFilter == 1) {
		cout << "Particle Filter activated" << endl;
		logFile << "Particle Filter activated" << endl;
		vector<Particle>& ptclList = FEMDataExt->getFEMGeoData()->getNodesVec();
		filterSupportingParticles(nodesVec, ptclList, InputData, logFile);
	}

#ifdef _GridNodesDebugMode_
	logFile << "********** Supporting List **********" << endl;
	for(int i = 0; i < nodesVec.size(); i++) {
		logFile << "Node[" << i << "]: ";
		intVector& suppPtcls = nodesVec[i].getSPtlcs();
		for(int j = 0; j < suppPtcls.size(); j++) {
			logFile << suppPtcls[j] << " ";
		}
		logFile << endl;
	}

	logFile << "Total nodes outside the geometry: " << nodesOutsideGeo.size() << endl;
	printVector(nodesOutsideGeo,"nodesOutsideGeo",logFile);
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::setSupportingParticles_(FEMGeometryExt* FEMDataExt,
		std::vector<Node>& nodesVec, InputFileData* InputData,
		ofstream& logFile) {

	logFile << "GridNodes::setSupportingParticles_" << endl;

	std::map<std::string, oPType> modelData; // To satisfy findPointInGeometry input argument
	modelData["integrationMethod"] = 1; // Gauss integration scheme.
	modelData["usedDegreesOfFreedom"] = 3;
	modelData["usedDimensions"] = 3;

	int volumeElem = 0;
	intVector nodesOutsideGeo, nodesInsideGeo;

	vector<GaussPoint> gaussPointsVec;
	int counterID = 0;

	FEMGeometry* FEMData = FEMDataExt->getFEMGeoData();

	FEMDataExt->printVolumePtclsDetails(InputData,logFile);

	logFile << "Total number of grid nodes: " << nodalList.size() << endl;
	for (int i = 0; i < nodesVec.size(); i++) {

		dbVector nodeCoord = nodesVec[i].getCoords();
		volumeElem = FEMDataExt->findPointInGeometry(nodeCoord, InputData,
				modelData, logFile);

		nodesVec[i].setVolumeElement(volumeElem);

		if (volumeElem != -1) {

			nodesInsideGeo.resize(nodesInsideGeo.size() + 1);
			nodesInsideGeo[nodesInsideGeo.size() - 1] = i;

			counterID++;

			GaussPoint gPoint;
			gPoint.setGlobalID(counterID);

			gPoint.setCoords(nodeCoord);

			gPoint.elementInfo.resize(1);
			gPoint.elementInfo[0] = volumeElem;

			gaussPointsVec.resize(gaussPointsVec.size() + 1);
			gaussPointsVec[gaussPointsVec.size() - 1] = gPoint;

			intVector volumeNodes =
					FEMDataExt->getNodesElemsVec()[volumeElem].getNodes();

#ifdef _GridNodesDebugMode_
			logFile << "GaussPoint[" << gPoint.getGlobalID()
					<< "](ID:" << i << "|";
			for(int j = 0; j < gPoint.getCoords().size(); j++){
				logFile << gPoint.getCoords()[j] << ",";
			}
			logFile	<< ") is in volume " << gPoint.elementInfo[0]
					<< " with nodes:";
			for (int j = 0; j < volumeNodes.size(); j++) {
				logFile << " " << volumeNodes[j] << ",";
			}
			logFile << endl;
#endif

		} else {
			nodesOutsideGeo.push_back(i);
		}
	}

	saveSelectedNodesToFile(nodesInsideGeo,logFile);

	logFile << "Num of Gauss Points created: " << gaussPointsVec.size() << endl;
	logFile << "Num of Gauss Points to be replaced: " << FEMData->getGaussPointsVec().size() << endl;

	FEMData->getGaussPointsVec() = gaussPointsVec;
	logFile << "GaussPoint list overwritten in FEMData " << endl;

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	// Set up the meshfree approximation data.

	std::map<std::string, double> calcData;
	PetscViewer viewerMPI, viewerSEQ;

	MeshlessData = new MeshlessApproximation(FEMData,
			FEMDataExt->getFEMInputData(), calcData, modelData, logFile,
			viewerMPI, viewerSEQ);

	findSupportingParticles(modelData, FEMDataExt->getFEMInputData(), logFile);

	vector<GaussPoint>& gaussPointsVecMeshless =
										MeshlessData->getGaussPointsVec();

	logFile << "gaussPointsVecMeshless: " << gaussPointsVecMeshless.size()
			<< endl;

//	vector<ParticleExt>& ptclList = FEMDataExt->getNodesVec();
	vector<Particle>& ptclList = MeshlessData->getParticlesVec();

	logFile << "In each particle, elems store:" << endl;
	for(int i = 0; i < ptclList.size(); i++){

		intVector elems = ptclList[i].getElems();
		logFile << i << ")[" << elems.size() << "]:";
		for(int j = 0 ; j < elems.size(); j++){
			logFile << elems[j] << ", ";
		}
		logFile << endl;
	}


//	for (int i = 0; i < nodesInsideGeo.size(); i++) {
//
//		dbVector& nCoord = nodesVec[nodesInsideGeo[i]].getCoords();
//		intVector& supPtcls = gaussPointsVecMeshless[i].getSupportPtcls();
//		printVector(nCoord,"nCoord",logFile);
//		for(int i = 0; i < supPtcls.size(); ){
//			printVector(ptclList[supPtcls[i]].getCoords(),"supPtcls[i]].getCoords()",logFile);
//			bool isSame = true;
//			for (int j = 0; j < nCoord.size(); j++) {
//				if(nCoord[j] != ptclList[supPtcls[i]].getCoords()[j])
//					isSame = false;
//			}
//
//			if (isSame == true) {
//				supPtcls.erase(supPtcls.begin() + i);
//			} else {
//				i++;
//			}
//		}
//
//		nodesVec[nodesInsideGeo[i]].getSPtlcs() = supPtcls;
//	}
#ifdef _GridNodesDebugMode_
	logFile <<"In GridNodes::setSupportingParticles_, particle list is:" << endl;
	for(int i = 0; i < ptclList.size(); i++){
		logFile << "Particle[" << i << "]: ";
		for(int j = 0; j < ptclList[i].getCoords().size(); j++){
			logFile << ptclList[i].getCoords()[j] << ", ";
		}
		logFile << endl;
	}
#endif

	for (int i = 0; i < nodesInsideGeo.size(); i++) {

		dbVector& nCoord = nodesVec[nodesInsideGeo[i]].getCoords();
		intVector& supPtcls = gaussPointsVecMeshless[i].getSupportPtcls();

		// Delete supporting nodes that coincides with gridNodes (Disabled temporarily
		// TODO: Investiage if algorithm below makes a big difference on the results
		//		logFile << " ----------------------------------" << endl;
//		printVector(gaussPointsVecMeshless[i].getSupportPtcls(),
//				"gaussPointsVecMeshless[i].getSupportPtcls()", logFile);
//		printVector(nCoord, "nCoord", logFile);

//		for (int j = 0; j < supPtcls.size();) {
//			logFile << "Considering particle: " << supPtcls[j] << endl;
//			printVector(ptclList[supPtcls[j]].getCoords(),
//					"supPtcls[j]].getCoords()", logFile);
//
//			if (nCoord == ptclList[supPtcls[j]].getCoords()) {
//				logFile << "Particle [" << supPtcls[j] << "|"
//						<< &*(supPtcls.begin() + j) << "] removed" << endl;
//				supPtcls.erase(supPtcls.begin() + j);
//			} else {
//				j++;
//			}
//		}

//		printVector(gaussPointsVecMeshless[i].getSupportPtcls(),
//				"gaussPointsVecMeshless[i].getSupportPtcls()", logFile);

		nodesVec[nodesInsideGeo[i]].getSPtlcs() = supPtcls;
	}

	int activateFilter = InputData->getValue("filterSupportingParticles");
	if(activateFilter == 1){
		cout << "Particle Filter activated" << endl;
		logFile << "Particle Filter activated" << endl;
		filterSupportingParticles(nodesVec,ptclList,InputData,logFile);
	}

#ifdef _GridNodesDebugMode_
	for (int j = 0; j < nodesVec.size(); j++) {
		logFile << "Node[" << j << "] Supporting ptcls("
				<< nodesVec[j].getSPtlcs().size() << "): ";
		if (nodesVec[j].getSPtlcs().size() > 0) {
			for (int k = 0; k < nodesVec[j].getSPtlcs().size(); k++)
				logFile << nodesVec[j].getSPtlcs()[k] << " ";
		} else {
			logFile << "NULL" << endl;
		}
		logFile << endl;
	}
	logFile << endl << "Nodes laying outside geometry(Num = "
			<< nodesOutsideGeo.size() << ") : ";
	for (int k = 0; k < nodesOutsideGeo.size(); k++) {
		logFile << nodesOutsideGeo[k] << ", ";
	}
	logFile << endl;
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::findSupportingParticles(
		std::map<std::string, oPType>& modelData, InputFileData* InputData,
		ofstream& logFile) {

	logFile << "-------- Starting Gauss point support search --------" << endl;

	InputData->setValue("supportComputationMode", 2);

	MeshlessData->BackgroundMesh::setGaussPtcleConn(InputData, modelData,
			logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::filterSupportingParticles(std::vector<Node>& nodesVec,
		vector<Particle>& ptclList, InputFileData* InputData,
		ofstream& logFile) {

	logFile << "-------- GridNodes::filterSupportingParticles --------" << endl;

	int allowedSupPtcls = InputData->getValue("filterSupportingParticlesLimit");

	for (int i = 0; i < nodesVec.size(); i++) {

//		logFile << "----------------------" << endl;
//		logFile << "For node[" << i << "]:" << endl;

		dbVector& nCoords = nodesVec[i].getCoords();
		intVector sPtcls = nodesVec[i].getSPtlcs();
//		printVector(sPtcls,"sPtcls",logFile);
		if (allowedSupPtcls < sPtcls.size()) {
			// find sPtcls distances to nCoords
			dbVector distVec(sPtcls.size(), 0);
			for (int j = 0; j < sPtcls.size(); j++) {
				dbVector& sCoords = ptclList[sPtcls[j]].getCoords();

				double sum = 0;
				for (int k = 0; k < sCoords.size(); k++) {
					sum += pow(nCoords[k] - sCoords[k], 2);
				}
				distVec[j] = sqrt(sum);
			}
//			printVector(distVec,"distVec",logFile);

			while(sPtcls.size() > allowedSupPtcls){
				vector<double>::iterator it = max_element(distVec.begin(),distVec.end());
				int location = it - distVec.begin();

//				logFile << "Max Value found: " << *it << endl;

				distVec.erase(distVec.begin()+location);
				sPtcls.erase(sPtcls.begin()+location);
//				printVector(distVec,"distVec(updated)",logFile);
//				printVector(sPtcls,"sPtcls(updated)",logFile);

			}

			// Select the specified amount of smallest distance
		}
		nodesVec[i].getSPtlcs() = sPtcls;
//		printVector(nodesVec[i].getSPtlcs(),"nodesVec[i].getSPtlcs()",logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::setInterpolantsOnNodes(FEMGeometryExt* FEMData, InputFileData* InputData,
		ofstream& logFile) {

	vector<ParticleExt> ptcls = FEMData->getNodesVec();

	logFile <<"In GridNodes::setInterpolantsOnNodes, particle list is:" << endl;
	for(int i = 0; i < ptcls.size(); i++){
			logFile << "Particle[" << i << "]: ";
			for(int j = 0; j < ptcls[i].getCoords().size(); j++){
				logFile << ptcls[i].getCoords()[j] << ", ";
			}
			logFile << endl;
	}

	for (int i = 0; i < nodalList.size(); i++) {

		intVector supportingPtlcs = nodalList[i].getSPtlcs();

#ifdef _GridNodesDebugMode_
		logFile << "--- Node[" << nodalList[i].getId() << "] ---" << endl;
		printVector(supportingPtlcs, "Supporting particles", logFile);
#endif

		if (supportingPtlcs.size() > 0) {

			dbMatrix sPtclsCoord(nodalList[i].getSPtlcs().size(), dbVector(0));

			for (int j = 0; j < nodalList[i].getSPtlcs().size(); j++) {
//				sPtclsCoord[j] = ptcls[nodalList[i].getSPtlcs()[j]-1].getCoords();
				sPtclsCoord[j] = ptcls[nodalList[i].getSPtlcs()[j]].getCoords();
			}

			// Find radius for interpolants calculation
			dbVector radiusVec(sPtclsCoord[0].size(),0);
			dbVector sPtclsInterpolants;

			oPType maxRad = 0;
			int maxRadId;

			int choice = 2;

			switch (choice) {

			case 1:
				// --Find the farthest particle
				for (int k = 0; k < sPtclsCoord.size(); k++) {

					oPType dist = 0;
					for (int l = 0; l < sPtclsCoord[k].size(); l++) {
						dist += pow(
								sPtclsCoord[k][l] - nodalList[i].getCoords()[l],
								2);
					}

					if (maxRad < dist) {
						maxRad = dist;
						maxRadId = k;
					}
				}

				for (int m = 0; m < sPtclsCoord[maxRadId].size(); m++) {
					radiusVec[m] = (sPtclsCoord[maxRadId][m]
							- nodalList[i].getCoords()[m]) * 1.1;
				}

				break;

			case 2:

				for (int k = 0; k < radiusVec.size(); k++) {
					for(int l = 0; l < sPtclsCoord.size(); l++){

						if(radiusVec[k] < abs(sPtclsCoord[l][k]-
								nodalList[i].getCoords()[k])){
							radiusVec[k] = abs(sPtclsCoord[l][k]-
									nodalList[i].getCoords()[k]) * 1.1;
						}
					}
				}


				break;

			default:
				logFile << "In GridNodes::setInterpolantsOnNodes, choice of"
						" radius determination algorithm is not valid" << endl;
				cout << "In GridNodes::setInterpolantsOnNodes, choice of"
						" radius determination algorithm is not valid" << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}

#ifdef _GridNodesDebugMode_
			logFile << "******* GridNodes::setInterpolantsOnNodes *******"
					<< endl;
			printVector(nodalList[i].getCoords(), "nodalList[i].getCoords()", logFile);
			printMatrix(sPtclsCoord, "sPtclsCoord", logFile);
			printVector(radiusVec, "radiusVec", logFile);
#endif

			int intType = InputData->getValue("gridNodesInterpolationType");

			switch(intType){

			case 1:
				// Calculating interpolants using multi-DIM MLS
				Interpolation(nodalList[i].getCoords(), sPtclsCoord, radiusVec,
						sPtclsInterpolants, InputData, logFile);
				break;

			case 2:
				// Calculating interpolants using Seska's MLS
				calcMatrixFieldMLSApproximants(nodalList[i].getCoords(),
						sPtclsCoord, radiusVec, sPtclsInterpolants, InputData,
						logFile);
				break;

			}

			nodalList[i].setInterpolants(sPtclsInterpolants);
			printVector(sPtclsInterpolants,"sPtclsInterpolants",logFile);
		}
	}
}


/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::calcMatrixFieldMLSApproximants(dbVector& iPoint,
		dbMatrix& sPtclsCoords, dbVector& radiusVec, dbVector& interpolants,
		InputFileData* InputData, ofstream& logFile) {

	dbMatrix MLSMat(1, dbVector(1, 1));

	vector<Point> pointList(1,
			Point((int) InputData->getValue("usedDimensions")));

	pointList[0].getCoords() = iPoint;

	vector<Particle> supportingPtcls(sPtclsCoords.size(),
			Particle((int) InputData->getValue("usedDegreesOfFreedom")));

	for (int k = 0; k < sPtclsCoords.size(); k++) {
		supportingPtcls[k].getCoords() = sPtclsCoords[k];
		supportingPtcls[k].getMLSmat() = MLSMat;
	}

	std::map<std::string, double> modelData;
	modelData["basisPolynomOrder"] = InputData->getValue("basisPolynomOrder");
	modelData["usedDimensions"] = (int) InputData->getValue("usedDimensions");
	modelData["integrationMethod"] = 1; // Gauss integration scheme.

	PetscViewer viewerSEQ;

	// ---------------------------------------------------------------------

	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(logFile.rdbuf()); //redirect std::cout to out.txt!

	MeshlessData->MLSDiscretising::matrixFieldMLSApproximation(supportingPtcls,
			pointList, modelData, logFile, viewerSEQ);

	std::cout.rdbuf(coutbuf); //reset to standard output again

	// ---------------------------------------------------------------------

	logFile << "-------- Checking Points shape functions --------" << endl;
	oPType sum;
	for (int i = 0; i < pointList.size(); i++) {
		logFile << "Point " << i << " (";

		for (int j = 0; j < pointList[i].getCoords().size(); j++) {
			logFile << pointList[i].getCoords()[j] << " ,";
		}
		logFile << ") -> " << endl;

		sum = 0;
		dbVector& shapeFuncs = pointList[i].getShapeFuncs();
		for (int j = 0; j < shapeFuncs.size(); j++) {
			logFile << shapeFuncs[j] << ", ";
			sum += shapeFuncs[j];
		}
		logFile << "-> [" << sum << "]";
		logFile << endl;
	}

	interpolants.clear();
	interpolants = pointList[0].getShapeFuncs();

}

/*!***************************************************************************
!***************************************************************************
void GridNodes::interpolateNodalResult(Data& iData, InputFileData* InputData,
		ofstream& logFile) {

	vector<ParticleExt>& ptcls = iData.getMeshData()->getNodesVec();

	for (int i = 0; i < nodalList.size(); i++) {

		intVector supportingPtcls = nodalList[i].getSPtlcs();

		if (supportingPtcls.size() > 0) {
			dbVector& interpolants = nodalList[i].getInterpolants();

			int numCalcStep =
//			ptcls[supportingPtcls[0] - 1].getStepDOFMat()[0].size();
					ptcls[supportingPtcls[0]].getStepDOFMat()[0].size();

//			nodalList[i].getStepDOFMat().resize(
//					ptcls[supportingPtcls[0] - 1].getStepDOFMat().size(),
//					dbVector(numCalcStep, 0));

			nodalList[i].getStepDOFMat().resize(
					ptcls[supportingPtcls[0]].getStepDOFMat().size(),
					dbVector(numCalcStep, 0));

			for (int j = 0; j < numCalcStep; j++) {
				for (int k = 0; k < supportingPtcls.size(); k++) {
					for (int l = 0; l
//									< ptcls[supportingPtcls[k] - 1].getStepDOFMat().size();
							< ptcls[supportingPtcls[k]].getStepDOFMat().size();
							l++) {
//						nodalList[i].getStepDOFMat()[l][j] +=
//								ptcls[supportingPtcls[k] - 1].getStepDOFMat()[l][j]
//										* interpolants[k];

						nodalList[i].getStepDOFMat()[l][j] +=
								ptcls[supportingPtcls[k]].getStepDOFMat()[l][j]
										* interpolants[k];
					}
				}
			}

		}
	}
#ifdef _GridNodesDebugMode_
	logFile << "******* GridNodes::interpolateNodalResult *******" << endl;

	for (int i = 0; i < nodalList.size(); i++) {
		logFile << "Node: " << i << " ";
		intVector supportingPtcls = nodalList[i].getSPtlcs();
		if (supportingPtcls.size() > 0) {

			logFile << "Interpolants(Supporting Particles): ";
			for (int j = 0; j < supportingPtcls.size(); j++) {
				logFile << nodalList[i].getInterpolants()[j] << "("
						<< supportingPtcls[j] << ") ";
			}
			logFile << endl;

			printMatrix(nodalList[i].getStepDOFMat(),
					"Step-Displacement Matrix", logFile);

		} else {
			logFile << "Supporting Particles: NULL" << endl;
		}
	}
#endif

}*/

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::interpolateNodalResult(Data& iData, InputFileData* InputData,
		ofstream& logFile) {

	logFile << "------------------------------------" << endl;
	logFile << "In GridNodes::interpolateNodalResult" << endl;
	logFile << "------------------------------------" << endl;

	vector<ParticleExt>& ptcls = iData.getMeshData()->getNodesVec();

	for (int i = 0; i < nodalList.size(); i++) {

		logFile << "Node[" << i << "]" << endl;

		intVector& supportingPtcls = nodalList[i].getSPtlcs();

		if (supportingPtcls.size() > 0) {
			dbVector& interpolants = nodalList[i].getInterpolants();

			printVector(supportingPtcls,"supportingPtcls",logFile);
			printVector(interpolants,"interpolants",logFile);


			int numCalcStep =
					ptcls[supportingPtcls[0]].getStepDOFMat()[0].size();

			nodalList[i].getStepDOFMat().resize(
					ptcls[supportingPtcls[0]].getStepDOFMat().size(),
					dbVector(numCalcStep, 0));

			for (int j = 0; j < numCalcStep; j++) {
				for (int k = 0; k < supportingPtcls.size(); k++) {
					for (int l = 0; l
							< ptcls[supportingPtcls[k]].getStepDOFMat().size();
							l++) {

						nodalList[i].getStepDOFMat()[l][j] +=
								ptcls[supportingPtcls[k]].getStepDOFMat()[l][j]
										* interpolants[k];
					}
				}
			}

			printMatrix(nodalList[i].getStepDOFMat(),
					"nodalList[i].getStepDOFMat()",logFile);

		}
	}
#ifdef _GridNodesDebugMode_
	logFile << "******* GridNodes::interpolateNodalResult *******" << endl;

	for (int i = 0; i < nodalList.size(); i++) {
		logFile << "Node: " << i << " ";
		intVector supportingPtcls = nodalList[i].getSPtlcs();
		if (supportingPtcls.size() > 0) {

			logFile << "Interpolants(Supporting Particles): ";
			for (int j = 0; j < supportingPtcls.size(); j++) {
				logFile << nodalList[i].getInterpolants()[j] << "("
						<< supportingPtcls[j] << ") ";
			}
			logFile << endl;

			printMatrix(nodalList[i].getStepDOFMat(),
					"Step-Displacement Matrix", logFile);

		} else {
			logFile << "Supporting Particles: NULL" << endl;
		}
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::setSupportingNodes(Data& iData, InputFileData* InputData,
		ofstream& logFile) {

	intMatrix supportingNodesList(iData.getMeshData()->getNodesNum());

	// Find the neighbouring nodes supporting the particles
	for (int i = 0; i < nodalList.size(); i++) {
		for (int j = 0; j < nodalList[i].getSPtlcs().size(); j++) {
//			supportingNodesList[nodalList[i].getSPtlcs()[j]-1].push_back(i);
			supportingNodesList[nodalList[i].getSPtlcs()[j]].push_back(i);
		}
	}

#ifdef _GridNodesDebugMode_
	logFile << "********** Particles supported by nodes **********" << endl;
	for(int i = 0; i < supportingNodesList.size(); i++){
		logFile << "Ptcle[" << i << "] supported by Node("
				<< supportingNodesList[i].size() << "): ";
		for(int j = 0; j < supportingNodesList[i].size(); j++){
			logFile << supportingNodesList[i][j] << ", ";
		}
		logFile << endl;
	}
#endif

	// Find nodes having insufficient supporting particles
	int minPtclSup = InputData->getValue("minMLSSupport");
	dbVector insufSupportNodalList;
	for(int i=0; i<supportingNodesList.size(); i++) {
		if(supportingNodesList[i].size() < minPtclSup){
			resizeArray(insufSupportNodalList,insufSupportNodalList.size()+1);
			insufSupportNodalList[supportingNodesList.size()] = i;
		}
	}

#ifdef _GridNodesDebugMode_
	logFile << "Nodes having insufficient support(min support = " << minPtclSup
			<< "):" << endl;
	logFile << "------------------------------------------------------" << endl;
	for (int i = 0; i < insufSupportNodalList.size(); i++) {
		logFile << "Node " << insufSupportNodalList[i] << "("
				<< supportingNodesList[insufSupportNodalList[i]].size() << "):";
		for (int j = 0;
				j < supportingNodesList[insufSupportNodalList[i]].size(); j++) {
			logFile << " " << supportingNodesList[insufSupportNodalList[i]][j]
					<< ",";
		}
		logFile << endl;
	}
#endif

	intMatrix nodesInsideVolumeList(
			iData.getMeshData()->getNodesElemsVec().size(),intVector());

	// setup list of volume elements with the associated nodes
	int volElem;
	for(int i = 0; i < nodalList.size(); i++) {

		volElem = nodalList[i].getVolumeElement();

		resizeArray(nodesInsideVolumeList[volElem],
				nodesInsideVolumeList[volElem].size()+1);

		nodesInsideVolumeList[volElem][nodesInsideVolumeList[volElem].size()] = i;
	}

	for(int i = 0; i < insufSupportNodalList.size(); i++){

	}

	iData.setSupportingNodesList(supportingNodesList);

}

/*!****************************************************************************/
/*!****************************************************************************/
dbMatrix GridNodes::setInterpolantsOnParticles(Data& iData,
		InputFileData* InputData, ofstream& logFile) {

	intMatrix supportingNodesList = iData.getSupportingNodesList();
	dbMatrix sNodesInterpolantsList(supportingNodesList.size(), dbVector(0));

	vector<ParticleExt>& ptcls = iData.getMeshData()->getNodesVec();

	for (int i = 0; i < ptcls.size(); i++) {

		intVector supportingNodes = supportingNodesList[i];

		if (supportingNodes.size() > 0) {

			dbMatrix sNodesCoord(supportingNodes.size(), dbVector(0));

			// Extract supporting nodes coordinates
			for (int j = 0; j < supportingNodes.size(); j++) {
				sNodesCoord[j] = nodalList[supportingNodes[j]].getCoords();
			}

			// Find radius for interpolants calculation
			dbVector radiusVec(sNodesCoord[0].size(),0);
			dbVector sNodesInterpolants;

			oPType maxRad = 0;
			int maxRadId;

			dbVector maxRadVec(nDims, 0);
			// --Find the farthest particle
			for (int k = 0; k < sNodesCoord.size(); k++) {

				// A Different radius determination algorithm
				for (int l = 0; l < sNodesCoord[k].size(); l++) {
					if (maxRadVec[l]
							< abs(
									sNodesCoord[k][l]
											- ptcls[i].getCoords()[l])) {
						maxRadVec[l] = abs(
								sNodesCoord[k][l] - ptcls[i].getCoords()[l]);
					}
				}

			}

			for (int m = 0; m < maxRadVec.size(); m++) {
				radiusVec[m] = maxRadVec[m] * 1.1;
			}

#ifdef _GridNodesDebugMode_
			logFile << "******* GridNodes::setInterpolantsOnParticles *******"
					<< endl;
			printVector(ptcls[i].getCoords(), "ptcls[i].getCoords()", logFile);
			printMatrix(sNodesCoord, "sNodesCoord", logFile);
			printVector(radiusVec, "radiusVec", logFile);
#endif

			// Calculating interpolants
			Interpolation(ptcls[i].getCoords(), sNodesCoord, radiusVec,
					sNodesInterpolants, InputData, logFile);

			// nodalList[i].setInterpolants(sNodesInterpolants);

			sNodesInterpolantsList[i] = sNodesInterpolants;

#ifdef _GridNodesDebugMode_
			logFile << "For Particle[" << i << "], Interpolants are" << endl;
			printVector(sNodesInterpolants, "", logFile);
			logFile << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl
					<< endl;
#endif

		}
	}

	return sNodesInterpolantsList;
}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::interpolateParticleResult(Data& iData, dbMatrix interpolantsList,
		InputFileData* InputData, ofstream& logFile) {

	vector<ParticleExt>& ptcls = iData.getMeshData()->getNodesVec();
	intMatrix supportingNodesList = iData.getSupportingNodesList();

	for (int i = 0; i < ptcls.size(); i++) {
		dbVector interpolantsVec = interpolantsList[i];
		intVector sNodes = supportingNodesList[i];

		int numDof =
				nodalList[supportingNodesList[i][0]].getStepDOFMat().size();
		int numSteps =
				nodalList[supportingNodesList[i][0]].getStepDOFMat()[0].size();

		ptcls[i].getStepDOFMat().resize(numDof, dbVector(numSteps, 0));

		for (int j = 0; j < numSteps; j++) {
			for (int k = 0; k < numDof; k++) {
				for (int l = 0; l < interpolantsVec.size(); l++) {

#ifdef _GridNodesDebugMode_
					logFile << "Ptcle[" << i << "] -> "
							<< ptcls[i].getStepDOFMat()[k][j] << " = "
							<< nodalList[supportingNodesList[i][l]].getStepDOFMat()[k][j]
							<< " x " << interpolantsVec[l];
#endif

					ptcls[i].getStepDOFMat()[k][j] +=
							nodalList[supportingNodesList[i][l]].getStepDOFMat()[k][j]
									* interpolantsVec[l];

				}
			}
		}
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::interpolateParticleResult_(Data& iData, dbMatrix interpolantsList,
		InputFileData* InputData, ofstream& logFile) {

	vector<ParticleExt>& ptcls = iData.getMeshData()->getNodesVec();
	intMatrix supportingNodesList = iData.getSupportingNodesList();

	for (int i = 0; i < ptcls.size(); i++) {
		dbVector interpolantsVec = interpolantsList[i];
		intVector sNodes = supportingNodesList[i];

		int numDof =
				nodalList[supportingNodesList[i][0]].getStepDOFMat().size();
		int numSteps =
				nodalList[supportingNodesList[i][0]].getStepDOFMat()[0].size();

		ptcls[i].getStepDOFMat().resize(numDof, dbVector(numSteps, 0));

		for (int j = 0; j < numSteps; j++) {
			for (int k = 0; k < numDof; k++) {
				for (int l = 0; l < interpolantsVec.size(); l++) {

#ifdef _GridNodesDebugMode_
					logFile << "Ptcle[" << i << "] -> "
							<< ptcls[i].getStepDOFMat()[k][j] << " = "
							<< nodalList[supportingNodesList[i][l]].getStepDOFMat()[k][j]
							<< " x " << interpolantsVec[l];
#endif

					ptcls[i].getStepDOFMat()[k][j] +=
							nodalList[supportingNodesList[i][l]].getStepDOFMat()[k][j]
									* interpolantsVec[l];

				}
			}
		}
	}
}


/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::setResultOnNodes(Data& iData, dbMatrix dispMatrix,
		InputFileData* InputData, ofstream& logFile) {

	int numDOF = iData.getNdofs();
	int numSteps = dispMatrix[0].size();

	for (int i = 0; i < nodalList.size(); i++) {

		nodalList[i].getStepDOFMat().resize(numDOF, dbVector(numSteps, 0));

		for (int j = 0; j < numDOF; j++) {
			for (int k = 0; k < numSteps; k++) {
				nodalList[i].getStepDOFMat()[j][k] =
						dispMatrix[(i * numDOF) + j][k];
			}
		}
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::resetNodes(){
	for (int i = 0; i < nodalList.size(); i++) {
		nodalList[i].reset();
	}
}

void GridNodes::resetNodesStepDOFMat(){
	for (int i = 0; i < nodalList.size(); i++) {
		nodalList[i].resetStepDOFMat();
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
dbMatrix GridNodes::assembleNodalResultMatrix(Data& iData, InputFileData* InputData,
		ofstream& logFile){

	int numDof = 0, numSteps = 0;
	for (int i = 0; i < nodalList.size(); i++) {
		if (nodalList[i].getStepDOFMat().size() > 0) {
			numDof = nodalList[i].getStepDOFMat().size();
			numSteps = nodalList[i].getStepDOFMat()[0].size();
			break;
		}
	}

	dbMatrix resultField(nodalList.size()*numDof,dbVector(numSteps,0));

	for (int j = 0; j < nodalList.size(); j++) {
		if (nodalList[j].getStepDOFMat().size() > 0) {
			for (int k = 0; k < numDof; k++) {
				resultField[(j*numDof)+k] = nodalList[j].getStepDOFMat()[k];
			}
		}
	}

	return resultField;
}

/*!****************************************************************************/
/*!****************************************************************************/
dbMatrix GridNodes::assemblePtclResultMatrix(Data& iData,
		InputFileData* InputData, ofstream& logFile) {

	logFile << "***** GridNodes::assemblePtclDispMatrix *****" << endl;

	vector<ParticleExt>& ptcls = iData.getMeshData()->getNodesVec();

	int numDof = 0, numSteps = 0;
	for (int i = 0; i < ptcls.size(); i++) {
		if (ptcls[i].getStepDOFMat().size() > 0) {
			numDof = ptcls[i].getStepDOFMat().size();
			numSteps = ptcls[i].getStepDOFMat()[0].size();
			break;
		}
	}

	logFile << "numDof: " << numDof << endl;
	logFile << "numSteps: " << numSteps << endl;

	dbMatrix dispField(ptcls.size() * numDof, dbVector(numSteps, 0));

	for (int j = 0; j < ptcls.size(); j++) {
		if (ptcls[j].getStepDOFMat().size() > 0) {
			for (int k = 0; k < numDof; k++) {
				dispField[(j * numDof) + k] = ptcls[j].getStepDOFMat()[k];
			}
		}
	}

	logFile << "assembleDispMatrix completed" << endl;

	return dispField;

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::saveNodesToFile(ofstream& logFile){

	string outputFileName = "gridNodes.dat";
	ofstream writeStream(outputFileName.c_str());

	if (writeStream.is_open()) {
		for (int i = 0; i < nodalList.size(); i++) {
			for(int j=0; j<nodalList[i].getCoords().size();j++){
				writeStream <<nodalList[i].getCoords()[j] << " ";
			}
			writeStream << endl;
		}
	}
	writeStream.close();

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::saveSelectedNodesToFile(intVector& nodesVec,ofstream& logFile){

	string outputFileName;

	for(int i = 0; ; i++){
		outputFileName
			  = ((string)"selectedGridNodes" + convertIntToString(i) + ".dat");
		ifstream file_to_check (outputFileName.c_str());
		if(file_to_check == false) break;
	}

	ofstream writeStream(outputFileName.c_str());

	if (writeStream.is_open()) {
		for (int i = 0; i < nodesVec.size(); i++) {
			for(int j=0; j<nodalList[nodesVec[i]].getCoords().size();j++){
				writeStream <<nodalList[nodesVec[i]].getCoords()[j] << " ";
			}
			writeStream << endl;
		}
	}
	writeStream.close();

}

/*!****************************************************************************/
/*!****************************************************************************/
FEMGeometryExt* GridNodes::getGridGeometry(ofstream& logFile){

	if(gridGeometry != NULL)
		return gridGeometry;
	else{
		logFile << "In GridNodes::getGridGeometry, "
				"gridGeometry assigned to NULL." << endl;
		cout << "In GridNodes::getGridGeometry, "
						"gridGeometry assigned to NULL." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::delVolumeElement(intVector& volElemList,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	sortIntVector(volElemList, 0, volElemList.size() - 1);

	// Delete the volume elements
	vector<FEMElementExt>& volFEMList = gridGeometry->getNodesElemsVec();
	for (int i = volElemList.size() - 1; i < 0; i--) {
		volFEMList.erase(volFEMList.begin() + volElemList[i]);
	}

	// Delete the volume elements in the particle's volume element list. If a
	// particle ends up with zero volume element, the particle is deleted.
	vector<ParticleExt>& ptclList = gridGeometry->getNodesVec();
	for (int i = 0; i < ptclList.size();) {
		intVector& ptclVolList = ptclList[i].getElems();
		for (int j = 0; j < ptclVolList.size();) {
			for (int k = 0; k < volElemList.size(); k++) {
				if (ptclVolList[j] == volElemList[k]) {
					ptclVolList.erase(ptclVolList.begin() + j);
				} else {
					j++;
				}
			}
		}

		if (ptclVolList.size() == 0) {
			ptclList.erase(ptclList.begin() + i);
		} else {
			i++;
		}
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodes::clearNodalList(){

	vector<Node>().swap(nodalList);

}
