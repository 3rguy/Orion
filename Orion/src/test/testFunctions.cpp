/*
 * testFunctions.cpp
 *
 *  Created on: Jul 25, 2014
 *      Author: rama
 */

#include "testFunctions.h"


// *****************************************************************************
// *****************************************************************************
void testFunctions(InputFileData* InputData, ofstream& logFile) {

	int choice = InputData->getValue("orionTestMode");
	cout << "TESTING FUNCTION MODE: ";
	logFile << "TESTING FUNCTION MODE: ";
	switch (choice) {
	case 0: {
		//**********************************************************************
		// Point-In-Polygon test
		cout << "Find-Point-In-Polygon test" << endl;
		findPointInPolygonTest(InputData, logFile);
		break;

	}

	case 1: {
		//**********************************************************************
		// snapshot method test
		cout << "POD-Calculation test" << endl;
		PODCalcTest(InputData, logFile);
		break;
	}

	case 2: {
		//**********************************************************************
		// extracting result from FEM.res file
		cout << "read-result-file test" << endl;
		readResultFileTest(InputData, logFile);
		break;
	}

	case 3: {
		//**********************************************************************
		// Find supporting particles of a node in a geometry
		cout << "supporting-particles test" << endl;
		supportingPtclsTest(InputData, logFile);
		break;
	}

	case 4: {
		//**********************************************************************
		// extracting result from FEM.res file
		cout << "Volume-calculation" << endl;
//		calcVolume(InputData, logFile);
		calcVolumeTwo(InputData, logFile);
		break;
	}

	case 5: {
		//**********************************************************************
		// extracting result from FEM.res file
		cout << "MLS Interpolation test" << endl;
//			mlsTest(InputData, logFile);
		break;
	}

	case 6: {
		//**********************************************************************
		// Printing the POMs extracted from result matrix over time
		cout << "Printing snapshot POMs test" << endl;
		printSnapshotPOMs(InputData, logFile);
		break;
	}

	case 7: {
		//**********************************************************************
		// Printing the POMs extracted from result matrix over time
		cout << "Grid Nodes test" << endl;
		GridNodesTest* GNTest = new GridNodesTest(InputData, logFile);
		delete GNTest;
		break;
	}
	default:
		cout << "ERROR: In testFunctions.cpp::testFunctions,"
				" choice is not valid" << endl;
		logFile << "ERROR: In testFunctions.cpp::testFunctions,"
				" choice is not valid" << endl;

	}

}

// *****************************************************************************
// *****************************************************************************
void findPointInPolygonTest(InputFileData* InputData, ofstream& logFile) {

	// Defining parameter vector
	dbVector paramVec(3);
	paramVec[0] = 15000;
	paramVec[1] = 10;
	paramVec[2] = 0.5;

	dbVector anchorPoint(3);
	anchorPoint[0] = 10;
	anchorPoint[1] = 10;
	anchorPoint[2] = 0;

	// Create an instance of Data called myData
	Data myData(paramVec);

	string meshFolderName = "216/";
	myData.setFolderName(meshFolderName);

	myData.readMeshDataFile(InputData, logFile);

	myData.setFolderName(meshFolderName);
	std::map<std::string, double> modelData;

	//==========================================================================
	// Function to find if a point belong to a geometry
	//==========================================================================
	//Point definition
	dbVector pCoord(3);
//	cout << "Enter X-Coord: ";
//	cin >> pCoord[0];
//	cout << "Enter Y-Coord: ";
//	cin >> pCoord[1];
//	cout << "Enter Z-Coord: ";
//	cin >> pCoord[2];
//
//	cout << "Checking if point (" << pCoord[0] << "," << pCoord[1] << ","
//			<< pCoord[2] << ") is in geometry ..." << endl;

	FEMGeometryExt* FEMData = myData.getMeshData();

	// -------------------------------------------------------------------------
	// ----------------- Reset particle coordinates ----------------------------
	vector<ParticleExt>& ptcls = FEMData->getNodesVec();
	for (int i = 0; i < ptcls.size(); i++) {
		for (int j = 0; j < anchorPoint.size(); j++) {
			ptcls[i].getCoord(j) = ptcls[i].getCoord(j) - anchorPoint[j];
		}
	}
	FEMData->writeMeshFile("fem.msh",InputData, logFile);
	// -------------------------------------------------------------------------

	// -------------------------------------------------------------------------
	// --------------------- Read gridNodes.dat --------------------------------
	string coordinatesFile = "gridNodes.dat";
	ifstream coordFile(coordinatesFile.c_str());
	dbMatrix nodesList;
	logFile << "--------------- Reading gridNodes.dat -------------------------"
			<< endl;
	for (;;) {

		if (coordFile.eof()) {
			logFile << "Num of points read: " << nodesList.size() << endl;
			break;
		}

		nodesList.resize(nodesList.size() + 1, dbVector(3));

		cout << "Reading point: " << nodesList.size() << endl;

		//! Read the nodeID, name and the coordinates in each axis
		coordFile >> nodesList[nodesList.size() - 1][0]
				>> nodesList[nodesList.size() - 1][1]
				>> nodesList[nodesList.size() - 1][2];

		logFile << "Nodes[" << nodesList.size() - 1 << "]: "
				<< nodesList[nodesList.size() - 1][0] << " "
				<< nodesList[nodesList.size() - 1][1] << " "
				<< nodesList[nodesList.size() - 1][2];

	}
	// -------------------------------------------------------------------------

	// -------------------------------------------------------------------------
	// --------------------- Read excluded_nodes.dat --------------------------------
	string excludeFile = "excluded_nodes.dat";
	ifstream exFile(excludeFile.c_str());
	intVector excludedNodesList;
	logFile
			<< "--------------- Reading excluded_nodes.dat -------------------------"
			<< endl;
	if (exFile.good()) {
		for (;;) {

			if (exFile.eof()) {
				logFile << "Num of points read: " << excludedNodesList.size() << endl;
				break;
			}

			excludedNodesList.resize(excludedNodesList.size() + 1);

			cout << "Reading point: " << excludedNodesList.size() << endl;

			//! Read the nodeID, name and the coordinates in each axis
			exFile >> excludedNodesList[excludedNodesList.size() - 1];

			logFile << excludedNodesList[excludedNodesList.size() - 1];

		}
	}
	// -------------------------------------------------------------------------

	if (excludedNodesList.size() == 0) {
		for (int i = 0; i < nodesList.size(); i++) {

			cout << "Checking point: " << i << endl;

			pCoord[0] = nodesList[i][0];
			pCoord[1] = nodesList[i][1];
			pCoord[2] = nodesList[i][2];

			int pointElemLocation = FEMData->findPointInGeometry(pCoord,
					InputData, modelData, logFile);

			if (pointElemLocation == -1) {
				logFile << "=====>=>=>=> POINT IS NOT IN GEOMETRY <=<=<=<====="
						<< endl;
				excludedNodesList.push_back(i);
			} else {
				logFile << "Point is located in volume element: "
						<< pointElemLocation << endl;
			}
		}
	}

	logFile << "Nodes laying outside geometry are: " << endl;
	for (int i = 0; i < excludedNodesList.size(); i++) {
		logFile << "Node[" << excludedNodesList[i] << "]: "
				<< nodesList[excludedNodesList[i]][0] << ", "
				<< nodesList[excludedNodesList[i]][1] << ", "
				<< nodesList[excludedNodesList[i]][2] << endl;
	}

	string fileOutput = "podLogFile2";
	ofstream logFile2(fileOutput.c_str(), std::ofstream::out);

	logFile2 << "Re-checking excluded points: " << endl;
	for (int i = 0; i < excludedNodesList.size(); i++) {
		cout << "Re-Checking point: " << excludedNodesList[i] << endl;

		pCoord[0] = nodesList[excludedNodesList[i]][0];
		pCoord[1] = nodesList[excludedNodesList[i]][1];
		pCoord[2] = nodesList[excludedNodesList[i]][2];

		int pointElemLocation = FEMData->findPointInGeometry(pCoord, InputData,
				modelData, logFile2);

		if (pointElemLocation == -1) {
			logFile2 << "=====>=>=>=> POINT IS NOT IN GEOMETRY <=<=<=<====="
					<< endl;
		} else {
			logFile2 << "Point is located in volume element: "
					<< pointElemLocation << endl;
		}
	}

	string excludeCoordFile = "excluded_nodes_coords.dat";
	ofstream exCoordFile(excludeCoordFile.c_str());
	for (int i = 0; i < excludedNodesList.size(); i++) {
		exCoordFile << nodesList[excludedNodesList[i]][0] << " "
					<< nodesList[excludedNodesList[i]][1] << " "
					<< nodesList[excludedNodesList[i]][2] << endl;
	}
	exCoordFile.close();


}

// *****************************************************************************
// *****************************************************************************
void PODCalcTest(InputFileData* InputData, ofstream& logFile) {

	double energy_Level = 95;
	dbMatrix fullMatrix;

	PODCalcTest_readFile(fullMatrix,InputData,logFile);

	printMatrix(fullMatrix,"fullMatrix",logFile);

	// -------------------------------------------------------------------------
	dbMatrix reducedMatrix_one;
	logFile << "Using the SVD Method" << endl;
	logFile << "--------------------" << endl;
	InputData->setValue("PODCalculationType",1);
//	PODCalc* PODCalculation_one
//		= new PODCalc(fullMatrix,reducedMatrix_one,energy_Level,InputData,logFile);

	PODCalc PODCalculation_one(fullMatrix,reducedMatrix_one,energy_Level,InputData,logFile);

	// -------------------------------------------------------------------------
	dbMatrix reducedMatrix_three;
	logFile << "Using the Snapshot Method" << endl;
	logFile << "-------------------------" << endl;
	InputData->setValue("PODCalculationType",3);
//	PODCalc* PODCalculation_three
//		= new PODCalc(fullMatrix,reducedMatrix_three,energy_Level,InputData,logFile);

	PODCalc PODCalculation_three(fullMatrix,reducedMatrix_three,energy_Level,InputData,logFile);


//	printMatrix(PODCalculation_one->getPOMs(),"POMs_one",logFile);
//	printVector(PODCalculation_one->getPOVs(),"POVs_one",logFile);

//	printMatrix(PODCalculation_three->getPOMs(),"POMs_three",logFile);
//	printVector(PODCalculation_three->getPOVs(),"POVs_three",logFile);

	printMatrix(PODCalculation_one.getPOMs(),"POMs_one",logFile);
	printMatrix(PODCalculation_three.getPOMs(),"POMs_three",logFile);

	printVector(PODCalculation_one.getPOVs(),"POVs_one",logFile);
	printVector(PODCalculation_three.getPOVs(),"POVs_three",logFile);


//	delete PODCalculation_one;
//	delete PODCalculation_three;

}

// *****************************************************************************
// *****************************************************************************
void PODCalcTest_readFile(dbMatrix& fullMatrix, InputFileData* InputData,
		ofstream& logFile) {

	string line;
	string fileName = "myMatrix.txt";
	ifstream myfile(fileName.c_str());

	if (myfile.is_open()) {
		while (myfile.good()) {

			if (!getline(myfile, line))
				break;

			fullMatrix.resize(fullMatrix.size()+1,dbVector());

			logFile << "==============================================" << endl;
			logFile << "Reading line: " << line << endl;
			logFile << "==============================================" << endl;

			istringstream ss(line);
			for (;;) {
				string s;
				if (!getline(ss, s, ';'))
					break;
				logFile << "reading:" << atof(s.c_str()) << endl;
				fullMatrix[fullMatrix.size()-1].push_back(atof(s.c_str()));
			}
		}

		myfile.close();
	} else {
		cout << "In main.C, PODCalcTest_readFile is unable to open "
				"or cannot find file:" << fileName << endl;
		logFile << "In main.C, PODCalcTest_readFile is unable to open "
				"or cannot find file:" << fileName << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

// *****************************************************************************
// *****************************************************************************
void readResultFileTest(InputFileData* InputData,ofstream& logFile){

	std::chrono::time_point<std::chrono::system_clock> start, end;

	// -------------------------------------------------------------------------
	start = std::chrono::system_clock::now();

	readResultFileTest_M2(InputData, logFile);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "readResultFileTest_M2 calculation time: " << elapsed_seconds.count() << " sec"
			<< endl;
	logFile << "readResultFileTest_M2 calculation time: " << elapsed_seconds.count() << " sec"
			<< endl;

	// -------------------------------------------------------------------------
	start = std::chrono::system_clock::now();

	readResultFileTest_M1(InputData, logFile);

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;

	cout << "readResultFileTest_M1 calculation time: " << elapsed_seconds.count() << " sec"
	     << endl;
	logFile << "readResultFileTest_M1 calculation time: " << elapsed_seconds.count() << " sec"
			<< endl;

}

// *****************************************************************************
// *****************************************************************************
void readResultFileTest_M1(InputFileData* InputData,ofstream& logFile){
	Data myData;

	vector<vector<string> > resultNameList(11,vector<string>(2));
	resultNameList[0][0] = "displacement"; 		resultNameList[0][1] = "Vector";
	resultNameList[1][0] = "strain";			resultNameList[1][1] = "Matrix";
	resultNameList[2][0] = "effective strain";	resultNameList[2][1] = "Scalar";
	resultNameList[3][0] = "fibre strain E_11";	resultNameList[3][1] = "Scalar";
	resultNameList[4][0] = "fibre strain E_22";	resultNameList[4][1] = "Scalar";
	resultNameList[5][0] = "fibre strain E_33";	resultNameList[5][1] = "Scalar";
	resultNameList[6][0] = "stress";			resultNameList[6][1] = "Matrix";
	resultNameList[7][0] = "effective stress";	resultNameList[7][1] = "Scalar";
	resultNameList[8][0] = "fibre stress 1";	resultNameList[8][1] = "Scalar";
	resultNameList[9][0] = "fibre stress 2";	resultNameList[9][1] = "Scalar";
	resultNameList[10][0] = "fibre stress 3";	resultNameList[10][1] = "Scalar";

	string inputFileName = "fem.res";


	for(int i=0; i<resultNameList.size();i++){

		dbMatrix resultMatrix;
		dbVector resultStepVector;
		int numDofs ;

		cout << "Reading result: " << resultNameList[i][0] << endl;
		logFile << "Reading result: " << resultNameList[i][0] << endl;

		bool isResultFound = myData.readResultFile_resFormat_specificResult(inputFileName,
				resultNameList[i][0], resultMatrix, numDofs, resultStepVector,
				logFile);

		if (isResultFound == true)
			readResultFileTest_writeResultToFile(resultNameList[i][0],
					resultStepVector, resultMatrix, numDofs, InputData,
					logFile);

		printMatrix(resultMatrix,resultNameList[i][0].c_str(),logFile);

	}

}

// *****************************************************************************
// *****************************************************************************
void readResultFileTest_M2(InputFileData* InputData,ofstream& logFile){
	Data myData;

	vector<string> resultNameList(11);
	resultNameList[0] = "displacement";
	resultNameList[1] = "strain";
	resultNameList[2] = "effective strain";
	resultNameList[3] = "fibre strain E_11";
	resultNameList[4] = "fibre strain E_22";
	resultNameList[5] = "fibre strain E_33";
	resultNameList[6] = "stress";
	resultNameList[7] = "effective stress";
	resultNameList[8] = "fibre stress 1";
	resultNameList[9] = "fibre stress 2";
	resultNameList[10] = "fibre stress 3";

	string inputFileName = "fem.res";

//	dbMatrix resultMatrix;
//	dbVector resultStepVector;
//	map<double, map<string, dbMatrix> > resultsData;
//
//	myData.readResultFile_resFormat_allResult(inputFileName, resultNameList,
//			resultsData, logFile);
//
//	logFile << "Num of steps: " << resultsData.size() << endl;
//
//	map<double, map<string, dbMatrix> >::iterator it_stepVal,it_stepVal_beg,it_stepVal_end;
//	it_stepVal_beg = resultsData.begin();
//	it_stepVal_end = resultsData.end();
//	map<string, dbMatrix>::iterator it_result,it_result_beg,it_result_end;
//
//	for(it_stepVal = it_stepVal_beg; it_stepVal != it_stepVal_end; it_stepVal++){
//		logFile << "--------------------------------------------------" << endl;
//		logFile << "Step: " << it_stepVal->first << endl;
//		logFile << "--------------------------------------------------" << endl;
//
//		it_result_beg = it_stepVal->second.begin();
//		it_result_end = it_stepVal->second.end();
//
//		for(it_result = it_result_beg; it_result != it_result_end; it_result++){
//			dbMatrix& resultMat = it_result->second;
//
//			logFile << "Result: " << it_result->first.c_str() << endl;
//
//			for(int i=0; i < resultMat.size(); i++){
//				logFile << i << " " ;
//				for(int j=0; j < resultMat[i].size(); j++){
//					logFile << resultMat[i][j] << " ";
//				}
//				logFile << endl;
//			}
//		}
//
//	}

	logFile << "readResultFile_resFormat_2" << endl;
	myData.readResultFile_resFormat(inputFileName, logFile);

//	printVector(myData.getStepValueVec(),"StepValueVec",logFile);
//
//	vector<string> resultNameVec = myData.getResultNameList();
//
//	for(int i=0; i<resultNameVec.size(); i++){
//		logFile << "--------------------------------------------------" << endl;
//		logFile << "Result: " << resultNameVec[i] << "\t nDOF: "
//				<< myData.getResultDOF(resultNameVec[i].c_str()) << endl;
//		logFile << "--------------------------------------------------" << endl;
//
//		printMatrix(myData.getResult(resultNameVec[i].c_str()),"",logFile);
//	}

}

// *****************************************************************************
// *****************************************************************************
void readResultFileTest_writeResultToFile(string& resultName,dbVector& resultStepVector,
		dbMatrix& resultMatrix,int& numDofs,InputFileData* InputData,ofstream& logFile){

		ofstream writeToFemRes(resultName.c_str());
		writeToFemRes.precision(12);
		writeToFemRes.setf(ios_base::scientific,ios_base::floatfield);

		// Result details
		string analysis_name = " ";
		string my_type = "Vector";
		string my_location = "OnNodes";

		if (writeToFemRes.is_open()) {

			for (int j = 0; j < resultStepVector.size(); j++) {

				// Setting up the Result line
				writeToFemRes << "Result " << "\"" << resultName << "\" \""
						<< analysis_name << "\" " << resultStepVector[j] << " "
						<< my_type << " " << my_location << endl;

				writeToFemRes << "Values" << endl;

				int nodeNum = 1,counter = 0;
				for (int r = 0; r < resultMatrix.size(); r++) {

					if(counter == 0){
						writeToFemRes << nodeNum;
					}

					writeToFemRes << " " << resultMatrix[r][j];
					counter++;

					if(counter == numDofs){
						writeToFemRes << endl;
						counter = 0;
						nodeNum++;
					}
				}

				writeToFemRes << "End Values" << endl;
			}
			writeToFemRes.close();
		}

}

// *****************************************************************************
// *****************************************************************************
void supportingPtclsTest(InputFileData* InputData, ofstream& logFile) {

	string folderName = "";
	string meshFileName = folderName + "mesh.dat";
	string inputFileName = folderName + "input.dat";

	std::map<std::string, double> modelData;
	FEMGeometryExt* meshData = new FEMGeometryExt(InputData, modelData,
			meshFileName, inputFileName, logFile);

	while (1) {
		dbVector pCoords(3, 0);
		cout << "Enter a point's coordinates to find its supporting particles"
				<< endl;

		cout << "X: ";
		cin >> pCoords[0];

		if (pCoords[0] == 999) break;

		cout << "Y: ";
		cin >> pCoords[1];
		cout << "Z: ";
		cin >> pCoords[2];



		intVector sVols;
		intVector sPtcls = meshData->findSupportingPtcls(pCoords,sVols,InputData,
				logFile);

		if (sPtcls.size() == 0) {
			cout << "No supporting particles was found. Point may probably \n";
			cout << "lie outside the geometry" << endl;

			logFile << "No supporting particles was found. Point may probably \n";
			logFile << "lie outside the geometry" << endl;

		} else {
			cout << "The supporting volume elements are[Orion Format]:" << endl;
			logFile << "The supporting volume elements are[Orion Format]:" << endl;
			for (int i = 0; i < sVols.size(); i++) {
				cout << sVols[i] << " ";
				logFile << sVols[i] << " ";
			}
			cout << endl;
			logFile << endl;

			cout << "The supporting volume elements are[GID Format]:" << endl;
			logFile << "The supporting volume elements are[GID Format]:" << endl;
			for (int i = 0; i < sVols.size(); i++) {
				cout << sVols[i]+1 << " ";
				logFile << sVols[i]+1 << " ";
			}
			cout << endl;
			logFile << endl;

			cout << "The supporting particles are:" << endl;
			logFile << "The supporting particles are:" << endl;
			for (int i = 0; i < sPtcls.size(); i++) {
				cout << sPtcls[i] << " ";
				logFile << sPtcls[i] << " ";
			}
			cout << endl;
			logFile << endl;
		}
	}

}

// *****************************************************************************
// *****************************************************************************
void calcVolume(InputFileData* InputData,ofstream& logFile){

//	dbMatrix& loadMat = InputData->getSurfacePressureLoads();
//	printMatrix(loadMat,"Load Matrix",logFile);

#ifdef _InputFileType_
	dbMatrix& allLoadMat = InputData->getSurfacePressureLoads();
#else
	dbMatrix allLoadMat;
#endif
	printMatrix(allLoadMat,"allLoadMat",logFile);

	dbMatrix loadMat;

	for(int i = 0 ; i < allLoadMat.size(); i++){
		if(allLoadMat[i][4] == (double)1.0){
			loadMat.push_back(allLoadMat[i]);
		}
	}

	printMatrix(loadMat,"loadMat",logFile);

	int numOfFaces = loadMat.size();

	string folderName = "";
	string meshFileName = folderName + "mesh.dat";
	string inputFileName = folderName + "input.dat";


	InputData->setValue("usedDimensions",3);
	InputData->setValue("FEMGeometrySetupType",2);
	std::map<std::string, double> modelData;
	FEMGeometryExt* meshData = new FEMGeometryExt(InputData, modelData,
				meshFileName, inputFileName, logFile);

	vector<Particle> & ptcls = meshData->getFEMGeoData()->getNodesVec();

	double volume = 0 ;

	int x = 0;
	int y = 1;
	int z = 2;

	for(int i = 0 ; i < numOfFaces; i++){

		dbVector& v1 = ptcls[loadMat[i][1]-1].getCoords();
		dbVector& v2 = ptcls[loadMat[i][2]-1].getCoords();
		dbVector& v3 = ptcls[loadMat[i][3]-1].getCoords();

		volume += ((v2[y] - v1[y])*(v3[z] - v1[z]) -
				   (v2[z] - v1[z])*(v3[y] - v1[y]))*
				   (v1[x] + v2[x] + v3[x] );
	}

	volume = volume / 6;

	cout << "Geometry volume =" << volume << endl;
}

// *****************************************************************************
// *****************************************************************************
void calcVolumeTwo(InputFileData* InputData,ofstream& logFile){

	InputData->setValue("usedDimensions",3);
	InputData->setValue("FEMGeometrySetupType",2);

	dbVector param(3);
	param[0] = 0;
	param[1] = 1;
	param[2] = 2;

	Data myData(param);
	myData.getFolderName() = "./";

	myData.readMeshDataFile(InputData,logFile);
	myData.readResultFile(InputData,logFile);
	myData.readVentriclesPVGraphResultFile(InputData,logFile);

	dbMatrix& resultMatrix = myData.getResult("displacement");
	printMatrix(resultMatrix,"displacement Matrix",logFile);

	vector<dbMatrix> resultMatList(resultMatrix[0].size(),dbMatrix(resultMatrix.size()/3, dbVector(3)));


	for(int i = 0 ; i < resultMatrix[0].size(); i++){
		dbMatrix& mat = resultMatList[i];

		int row = 0; int col = 0;
		for(int j = 0 ; j < resultMatrix.size() ; j++){

			mat[row][col] = resultMatrix[j][i];
			col++;

			if(col>2){
				col = 0;
				row++;
			}
		}

//		logFile << " -------------------------------------------" << endl;
//		logFile << " STEP: " << i << endl;
//		printMatrix(mat,"mat",logFile);

	}


	FEMGeometry* myData_MeshData =  myData.getMeshData()->getFEMGeoData();
	InputFileData* myData_InputData =  myData.getMeshData()->getFEMInputData();

#ifdef _InputFileType_
	dbMatrix& allLoadMat = myData_InputData->getSurfacePressureLoads();
#else
	dbMatrix allLoadMat;
#endif

	printMatrix(allLoadMat,"allLoadMat",logFile);

	dbMatrix loadMat;

	int choice = 0;
	if (choice == 0) {

		// Left ventricle cavity volume calculation
		for (int i = 0; i < allLoadMat.size(); i++) {
			if (allLoadMat[i][4] == (double) 1.066) {
				loadMat.push_back(allLoadMat[i]);
			}
		}
	} else if (choice == 1) {

		// Right ventricle cavity volume calculation
		for (int i = 0; i < allLoadMat.size(); i++) {
			if (allLoadMat[i][4] == (double) 0.533) {
				loadMat.push_back(allLoadMat[i]);
			}
		}

	} else {
		loadMat = allLoadMat;
	}

	printMatrix(loadMat,"loadMat",logFile);

#ifdef _InputFileType_
	dbMatrix& lineLoadMat = myData_InputData->getLineDispBoundConds();
#else
	dbMatrix lineLoadMat;
#endif

	printMatrix(lineLoadMat,"lineLoadMat",logFile);

	int numOfFaces = loadMat.size();

	int numOfSteps = resultMatList.size();
	dbVector volumeVec(numOfSteps,0);

	logFile << "NumOfSteps: " << numOfSteps << endl;

	for(int i = 0 ; i < numOfSteps; i++){

		dbMatrix& rMatrix = resultMatList[i];

		vector<Particle> ptcls = myData_MeshData->getNodesVec();
		// Update the coordinates
		for(int j = 0 ; j < ptcls.size(); j++){
			for(int k = 0 ; k < 3; k++){
				ptcls[j].getCoord(k) += rMatrix[j][k];
			}
		}


		double volume = 0 ;

		int x = 0;
		int y = 1;
		int z = 2;

		for(int m = 0 ; m < numOfFaces; m++){

			dbVector& v1 = ptcls[loadMat[m][1]-1].getCoords();
			dbVector& v3 = ptcls[loadMat[m][2]-1].getCoords();
			dbVector& v2 = ptcls[loadMat[m][3]-1].getCoords();

			volume += ((v2[y] - v1[y])*(v3[z] - v1[z]) -
					   (v2[z] - v1[z])*(v3[y] - v1[y]))*
					   (v1[x] + v2[x] + v3[x] );
		}

		volumeVec[i] = volume / 6.0;

		cout << "[" << i << "] Geometry volume = " << volumeVec[i] << endl;
	}

	dbVector& pressureVec = myData.getLeftCavityPressures();

	ofstream writeGraph("defVolLoad_myAlgo.grf");
	writeGraph.precision(15);

	if(writeGraph.good()){
		for(int i = 0 ; i < pressureVec.size(); i++){
			writeGraph << volumeVec[i] << " " << pressureVec[i] << endl;
		}
	}

	writeGraph.close();
}

// *****************************************************************************
// *****************************************************************************
dbMatrix decompPolygonToTriangles(dbVector polyNodeList,
		InputFileData* InputData,ofstream& logFile){

	if(polyNodeList.size() < 3){
		logFile << "In decompPolygonToTriangles, polyNodeList.size() < 3" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int numFaces = polyNodeList.size() - 1 ;
	int numTriangles = numFaces - 2 ;

	dbMatrix triangles(numTriangles,dbVector(3,0));

	for(int i = 0 ; i < numTriangles; i++){
		triangles[i][0] = polyNodeList[0];
	}

	int triangleID = 0;
	int count = 1 ;
	for(int i = 1 ; i < polyNodeList.size(); i++){ // Skip the first entry
		triangles[triangleID][count] = polyNodeList[i];
		count++;

		if(count > 2){
			i--;
			count = 1;
			triangleID++;
		}
	}

	printMatrix(triangles,"triangles formed:",logFile);

	return triangles;

}

// *****************************************************************************
// *****************************************************************************
void printSnapshotPOMs(InputFileData* InputData,ofstream& logFile){

	Data myData;

	vector<string> resultNameList(11);
	resultNameList[0] = "displacement";
	resultNameList[1] = "strain";
	resultNameList[2] = "effective strain";
	resultNameList[3] = "fibre strain E_11";
	resultNameList[4] = "fibre strain E_22";
	resultNameList[5] = "fibre strain E_33";
	resultNameList[6] = "stress";
	resultNameList[7] = "effective stress";
	resultNameList[8] = "fibre stress 1";
	resultNameList[9] = "fibre stress 2";
	resultNameList[10] = "fibre stress 3";

	string inputFileName = "fem.res";

	logFile << "read Result File: " << inputFileName << endl;
	cout << "read Result File: " << inputFileName << endl;
	myData.readResultFile_resFormat(inputFileName, logFile);
	logFile << "read Result File completed" << endl;
	cout << "read Result File completed" << endl;

	vector<string> resultAvailableList = myData.getResultNameList();


	InputData->setValue("interpolantionType",
			InputData->getValue("dbInterpolantionType"));
	InputData->setValue("MLSCalculationType",
			InputData->getValue("dbMLSCalculationType"));
	InputData->setValue("MLSPolynomialDegree",
			InputData->getValue("dbMLSPolynomialDegree"));
	InputData->setValue("MLSWeightFunc",
			InputData->getValue("dbMLSWeightFunc"));
	InputData->setValue("parameterPolynomialDegree",
			InputData->getValue("dbparameterPolynomialDegree"));

	DataContainer* problemData = new DataContainer();

	problemData->setValue("PODIStepValueVec",myData.getStepValueVec());


	// ------------------------------------------------------------------------
	// Set random myParameters
	dbVector& stepValues = myData.getStepValueVec();
	double minVal, maxVal;
	minVal = *stepValues.begin();
	maxVal = *stepValues.end();

	dbVector myParameters(1);

	srand(time(NULL));
	double randNum = ((double) rand() / (RAND_MAX));

	myParameters[0] = minVal + ((maxVal-minVal)*0.5);

	problemData->setValue("myParameters",myParameters);


	// ------------------------------------------------------------------------
	dbVector parameterRadii(1);
	parameterRadii[0] = (stepValues[1] - stepValues[0])*2;

	problemData->setValue("parameterRadii",parameterRadii);

	// ------------------------------------------------------------------------
	intVector supportDataID;
	for(int i = 0 ; i < stepValues.size(); i++){
		supportDataID.push_back(i);
	}

	problemData->setValue("supportDataID",supportDataID);

	// ------------------------------------------------------------------------
	dbMatrix dataParametersList(stepValues.size(),dbVector(1));
	for(int i=0; i<stepValues.size(); i++){
		dataParametersList[i][0] = stepValues[i];
	}

	problemData->setValue("dataParametersList",dataParametersList);



	for (int i = 0; i < resultAvailableList.size(); i++) {

		logFile << "Processing result: " << resultAvailableList[i] << endl;
		cout << "Processing result: " << resultAvailableList[i] << endl;

		vector<dbMatrix> resultList(1,dbMatrix());
		resultList[0] = myData.getResult(resultAvailableList[i].c_str());
		problemData->setValue("resultList",resultList);

		problemData->setValue("PODIResultName", resultAvailableList[i]);
		problemData->setValue("PODIDofPerNode",
				myData.getResultDOF(resultAvailableList[i].c_str()));

		problemData->setValue("calcResultList",dbMatrix());

		PODICalc* Podi = new PODICalc(problemData->getDbVector("myParameters"),
				problemData->getDbVector("parameterRadii"),
				problemData->getIntVector("supportDataID"),
				problemData->getDbMatrix("dataParametersList"),
				problemData->getDbMatrixVec("resultList"),
				problemData->getDbMatrix("calcResultList"), problemData,
				InputData, logFile);

		delete Podi;

//		myData.setResult(resultAvailableList[i].c_str(),
//						problemData->getDbMatrix("calcResultList"));
//		dbMatrix().swap(problemData->getDbMatrix("calcResultList"));

		problemData->deleteDbMatrix("calcResultList");
		problemData->deleteDbMatrixVec("resultList");

		problemData->deleteString("PODIResultName");
		problemData->deleteInt("PODIDofPerNode");

	}

}


// *****************************************************************************
// *****************************************************************************
//void pcl_test(InputFileData* InputData,ofstream& logFile){
//
//	  pcl::PointCloud<pcl::PointXYZ> cloud;
//
//	  // Fill in the cloud data
//	  cloud.width    = 5;
//	  cloud.height   = 1;
//	  cloud.is_dense = false;
//	  cloud.points.resize (cloud.width * cloud.height);
//
//	  for (size_t i = 0; i < cloud.points.size (); ++i)
//	  {
//	    cloud.points[i].x = 1024 * rand () / (RAND_MAX + 1.0f);
//	    cloud.points[i].y = 1024 * rand () / (RAND_MAX + 1.0f);
//	    cloud.points[i].z = 1024 * rand () / (RAND_MAX + 1.0f);
//	  }
//
//	  pcl::io::savePCDFileASCII ("test_pcd.pcd", cloud);
//	  std::cerr << "Saved " << cloud.points.size () << " data points to test_pcd.pcd." << std::endl;
//
//	  for (size_t i = 0; i < cloud.points.size (); ++i)
//	    std::cerr << "    " << cloud.points[i].x << " " << cloud.points[i].y << " " << cloud.points[i].z << std::endl;
//
//
//}




//=============================================================================
// Function to read fem.res file
//=============================================================================
/*
 string fileName = "fem_mod.res";
 ifstream readFemRes(fileName.c_str());
 string line;

 vector<dbMatrix> listOfDisp;
 dbMatrix disp;
 dbVector step_value_vec;

 bool is_name_found, is_value_found, is_endValue_found;
 is_name_found = false;
 is_value_found = false;
 is_endValue_found = false;

 //string result_name = "displacement";
 string result_name = "displacement";
 string analysis_name;
 double step_value = 0;
 string my_result_type, my_location;

 int count = 0;

 // Look for file
 if (!readFemRes) {
 logFile << "File " << fileName << "does not exist or is empty" << endl;
 MPI_Abort(MPI_COMM_WORLD, 1);
 }

 // Read file
 while (readFemRes.is_open()) {
 getline(readFemRes, line);
 logFile << "Reading line => " << line << endl;
 if (readFemRes.good()) {

 //logFile << line << endl;
 istringstream is_line(line);
 string is_line_str;
 getline(is_line, is_line_str, '"');
 logFile << "is_line_str: " << is_line_str << endl;
 if (is_line_str == "Result ") {
 getline(is_line, is_line_str, '"');
 logFile << "***Reading Line: "<< is_line_str <<endl;
 if (is_line_str == result_name) {
 is_name_found = true;

 // Extract analysis_name
 getline(is_line, analysis_name, '"');
 logFile << "***analysis_name: "<< analysis_name <<endl;

 // Skip blank space
 getline(is_line, is_line_str, '"');
 logFile << "***Reading Line: "<< is_line_str <<endl;

 // Extract line containing step_value, my_result_type
 // and my_location;
 getline(is_line, is_line_str, '"');
 logFile << "***Reading Line: "<< is_line_str <<endl;
 istringstream is_line_ws(is_line_str);
 is_line_ws >> step_value >> my_result_type >> my_location;

 logFile << "step_value["<<step_value<<"] my_result_type["
 <<my_result_type<<"] my_location["<<my_location<<"]"
 <<endl;

 // Record the step values
 step_value_vec.push_back(step_value);
 }
 }

 if (is_name_found == true && is_value_found == true) {

 if (line == "End Values") {
 is_name_found = false;
 logFile << "===============> is_name_found = false" << endl;
 is_value_found = false;
 logFile << "===============> is_value_found = false"
 << endl;
 count++;
 listOfDisp.push_back(disp);
 disp.clear();

 } else {

 logFile
 << "**************************************************"
 << endl;
 logFile << line << endl;
 logFile
 << "**************************************************"
 << endl;

 istringstream is(line);
 double is_n;
 vector<double> temp;
 while (is >> is_n) {
 temp.push_back(is_n);
 logFile << "reading:" << is_n << endl;
 }

 disp.push_back(temp);
 logFile << "Size of: temp(" << temp.size() << ") & disp("
 << disp.size() << endl;
 }

 }

 // Locate the values and indicate if found
 if (line == "Values" && is_name_found == true) {
 is_value_found = true;
 }

 } else {
 break;
 }
 }

 // Close .res file
 readFemRes.close();

 // Remove all Nodal ID from the displacement vectors
 if (result_name == "displacement") {
 for (int i = 0; i < listOfDisp.size(); i++) {
 for (int j = 0; j < listOfDisp[i].size(); j++) {
 listOfDisp[i][j].erase(listOfDisp[i][j].begin());
 }
 }
 }

 #ifdef _PODDebugMode_
 logFile << "Printing the results read from fem.res" << endl;
 logFile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
 for (int k = 0; k < listOfDisp.size(); k++) {
 logFile << "Step Value" << step_value_vec[k]<<endl;
 printMatrix(listOfDisp[k], logFile);
 logFile << "=========================================================="
 << endl;
 }
 #endif
 */

//==========================================================================
// Function to find if a point belong to a geometry
//==========================================================================
//Point definition
/*dbVector pCoord(3);
 cout << "Enter X-Coord: ";
 cin >> pCoord[0];
 cout << "Enter Y-Coord: ";
 cin >> pCoord[1];
 cout << "Enter Z-Coord: ";
 cin >> pCoord[2];

 cout << "Checking if point (" << pCoord[0] << "," << pCoord[1] << ","
 << pCoord[2] << ") is in geometry ..." << endl;

 int pointElemLocation = FEMData->findPointInGeometry(pCoord, InputData,
 modelData, logFile);

 if (pointElemLocation == -1) {
 cout << "=====>=>=>=> POINT IS NOT IN GEOMETRY <=<=<=<=====" << endl;
 } else {
 cout << "Point is located in volume element: " << pointElemLocation
 << endl;
 }*/

//==========================================================================
// Function to read mesh.dat file using FEMGeometry class file from Seska
//==========================================================================
/* string meshFileName = "mesh.dat";

 myData.setFileName(meshFileName);
 std::map<std::string, double> modelData;

 myData.readMeshFile(InputData, modelData, logFile);

 //FEMGeometry* FEMData = new FEMGeometry(InputData,modelData,meshFileName,logFile);
 FEMGeometry* FEMData = myData.meshDataAccess();
 cout << "Reading and Building of FEMGeometry completed" << endl;

 vector<Particle> myPtclList = FEMData->getNodesVec();
 vector<FEMElement> myElemList = FEMData->getNodesElemsVec();

 logFile << "PARTICLE DETAILS" << endl;
 logFile << "Number of Particles: " << myPtclList.size() << endl;
 for (int i = 0; i < myPtclList.size(); i++) {
 logFile << "Particle[" << myPtclList[i].getID() << "]->("
 << myPtclList[i].getCoord(0) << "," << myPtclList[i].getCoord(1)
 << "," << myPtclList[i].getCoord(2) << ") Elements:";

 intVector ptclElem = myPtclList[i].getElems();

 for (int j = 0; j < ptclElem.size(); j++) {
 logFile << ptclElem[j] << ",";
 }
 logFile << endl;
 }

 logFile << "============================================================="
 << endl;
 logFile << "ELEMENT DETAILS" << endl;
 logFile << "Number of Elements: " << myElemList.size() << endl;
 for (int k = 0; k < myElemList.size(); k++) {
 logFile << "Element[" << myElemList[k].getGlobalID() << "]->Nodes: ";
 for (int l = 0; l < myElemList[k].getNodes().size(); l++) {
 logFile << myElemList[k].getNodes()[l] << ",";
 }
 logFile << endl;
 }

 logFile << "============================================================="
 << endl;
 logFile << "SURFACE ELEMENT DETAILS" << endl;
 logFile << "Number of Particles: " << myPtclList.size() << endl;
 for (int k = 0; k < myPtclList.size(); k++) {
 logFile << "Particle[" << myPtclList[k].getID() << "]->Normals["
 << myPtclList[k].getAllSurfaceNormals().size() << "]: ";
 for (int l = 0; l < myPtclList[k].getAllSurfaceNormals().size(); l++) {
 logFile << l << ": ";
 for (int m = 0; m < myPtclList[k].getAllSurfaceNormals()[l].size();
 m++) {
 logFile << myPtclList[k].getAllSurfaceNormals()[l][m] << ",";
 }
 logFile << endl;
 }
 logFile << endl;
 }
 logFile << "============================================================="
 << endl;
 */

//==========================================================================
// Function to generate interpolants
//==========================================================================
//Interpolation* interpolants = new Interpolation(iPoint,coords,radiusVec,
//		InputData,logFile);

//==========================================================================

//dbMatrix readMatrix(){
//	string fileInput = "matrixFile.csv";
//		ifstream readFile(fileInput.c_str(), std::ofstream::out);
//
//		int row,col,temp;
//		dbMatrix myMatrix;
//		string line;
//		if (readFile.is_open()) {
//				readFile >> row >> col ;
//				cout << "Rows:" << row << endl;
//				cout << "Columns:" << col << endl;
//
//				getline(readFile, line); // Get rid of ghost line
//
//				resizeArray(myMatrix, row, col);
//
//				for (int i = 0; i < row; i++) {
//					getline(readFile, line);
//					istringstream ss(line);
//					for (int j = 0; j < row; j++) {
//						string s;
//						if (!getline(ss, s, ' '))
//							break;
//						cout << "reading:" << atof(s.c_str()) << endl;
//						myMatrix[i][j] = atof(s.c_str());
//					}
//				}
//		}
//		return myMatrix;
//}

//dbMatrix myInverseCalc(dbMatrix& myMatrix, ofstream& logFile){
//
//	string myS;
//
//	// CHECK: Matrix dimensions
//		for(int i=1;i<myMatrix.size();i++){
//			if(myMatrix[i-1].size() != myMatrix[i].size()){
//				cout << "ERROR: Matrix is missing some elements" << endl;
//				logFile << "ERROR: Matrix is missing some elements" << endl;
//				MPI_Abort(MPI_COMM_WORLD, 1);
//			}
//		}
//
//		// LU factorisation of matrix
//		// Define variables
//		int M = myMatrix.size(); 	// Num of rows of matrix
//		int N = myMatrix[0].size(); // Num of columns of matrix
//		dbVector A(M*N);			// Array of matrix
//		int LDA = max(1,M);			// Leading Dimensional Array of A
//		intVector IPIV(min(M,N));	// Pivot indices
//		int INFO;					// Outcome of the LU factorisation
//
//		// Define Array vector A
//		// The columns of matrix are laid out one after each other
//		for(int i=0; i<M; i++){
//			for(int j=0; j<N; j++){
//				A[i*LDA + j] = myMatrix[j][i];
//			}
//		}
//
//		// Debug outputs
//		myS = "A-array:"; printVector(A,myS,logFile);
//
//		// Compute LU factorisation
//		dgetrf_(M,N,&A[0],LDA,&IPIV[0],INFO);
//
//		if(INFO < 0){
//			cout << "ERROR: The "<< INFO <<"-th argument had an illegal value" << endl;
//			logFile << "ERROR: The "<< INFO <<"-th argument had an illegal value" << endl;
//			MPI_Abort(MPI_COMM_WORLD, 1);
//		}
//		else if (INFO > 0) {
//			cout <<"WARNING: U("<< INFO <<","<< INFO <<") is exactly zero. The factorization "
//					"has been completed, but the factor U is exactly singular, and "
//					"division by zero will occur if it is used to solve a system of "
//					"equations." <<endl;
//			logFile <<"WARNING: U("<< INFO <<","<< INFO <<") is exactly zero. The factorization "
//							"has been completed, but the factor U is exactly singular, and "
//							"division by zero will occur if it is used to solve a system of "
//							"equations." <<endl;
//		}
//
//		// Debug outputs
//		myS = "A-array results:"; printVector(A,myS,logFile);
//		myS = "IPIV:"; printVector(IPIV,myS,logFile);
//
//		// Define additional variables for inverse calculation
//		int LWORK = max(1,N);
//		dbVector WORK(LWORK);
//
//		// Carry out the inverse calculation
//		dgetri_(N,&A[0],LDA,&IPIV[0],&WORK[0],LWORK,INFO);
//		myS = "A inverse:"; printVector(A,myS,logFile);
//
//		if(INFO < 0){
//			cout <<"WARNING: The "<< INFO <<"-th argument had an illegal value"<< endl;
//			logFile <<"WARNING: The "<< INFO <<"-th argument had an illegal value"<< endl;
//			MPI_Abort(MPI_COMM_WORLD, 1);
//		}
//		else if (INFO > 0){
//			cout <<"WARNING: U(<<"<< INFO <<","<< INFO <<") is exactly zero; the matrix is "
//					"singular and its inverse could not be computed." << endl;
//			logFile <<"WARNING: U(<<"<< INFO <<","<< INFO <<") is exactly zero; the matrix is "
//					"singular and its inverse could not be computed." << endl;
//			MPI_Abort(MPI_COMM_WORLD, 1);
//		}
//
//		dbMatrix myMatrixInverse(M,dbVector(N));
//		for(int i=0; i<M; i++){
//			for(int j=0; j<N; j++){
//				myMatrixInverse[j][i] = A[i*LDA + j];
//			}
//		}
//
//		clock_t myTime = clock();
//		clock_t diffTime;
//		for(;;){
//			diffTime = clock() - myTime;
//			if((((float)diffTime)/CLOCKS_PER_SEC) > 2)
//				break;
//		}
//
//
//		return myMatrixInverse;
//}

//==========================================================================
// Finding limits of each variable type
//==========================================================================
//cout << "int: " << std::numeric_limits<int>::max() << endl;
//cout << "unsigned int: " << std::numeric_limits<unsigned int>::max() << endl;
//cout << "long: " << std::numeric_limits<long>::max() << endl;
//cout << "long int: " << std::numeric_limits<long int>::max() << endl;
//cout << "long unsigned int: " << std::numeric_limits<long unsigned int>::max() << endl;

//==========================================================================
// Function to test inverse lapack functions
//==========================================================================

//// Read matrix from files
//dbMatrix myMatrix = readMatrix();
//myS = "myMatrix:"; printMatrix(myMatrix,myS,logFile);
//
//dbMatrix seskaMatrixInverse;
//clock_t sTimer = clock();
//calcInvDouble(myMatrix,seskaMatrixInverse,logFile);
//sTimer = clock() - sTimer;
//myS = "seskaMatrixInverse:"; printMatrix(seskaMatrixInverse,myS,logFile);
//
//clock_t mTimer = clock();
//dbMatrix myMatrixInverse = myInverseCalc(myMatrix,logFile);
//mTimer = clock() - mTimer;
//myS = "myMatrixInverse:"; printMatrix(myMatrixInverse,myS,logFile);
//
//cout << "Duration: seska["<< sTimer
//		<<" clicks] v/s myMatrix["<< mTimer <<" clicks]" << endl;
