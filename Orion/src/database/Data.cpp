#include "Data.h"

Data::Data() {
	id = -1;
	nDims = 3;
	nDofs = 3;
	nNodes = 0;
	meshData = NULL;
	folderName = "";
	anchorPoint = -1;


	allResultsNameList.resize(11, vector<string>(2));

	allResultsNameList[0][0] = "displacement";
	allResultsNameList[0][1] = "Vector";

	allResultsNameList[1][0] = "strain";
	allResultsNameList[1][1] = "Matrix";

	allResultsNameList[2][0] = "effective strain";
	allResultsNameList[2][1] = "Scalar";

	allResultsNameList[3][0] = "fibre strain E_ff";
	allResultsNameList[3][1] = "Scalar";

	allResultsNameList[4][0] = "fibre strain E_tt";
	allResultsNameList[4][1] = "Scalar";

	allResultsNameList[5][0] = "fibre strain E_nn";
	allResultsNameList[5][1] = "Scalar";

	allResultsNameList[6][0] = "stress";
	allResultsNameList[6][1] = "Matrix";

	allResultsNameList[7][0] = "effective stress";
	allResultsNameList[7][1] = "Scalar";

	allResultsNameList[8][0] = "fibre stress S_ff";
	allResultsNameList[8][1] = "Scalar";

	allResultsNameList[9][0] = "fibre stress S_tt";
	allResultsNameList[9][1] = "Scalar";

	allResultsNameList[10][0] = "fibre stress S_nn";
	allResultsNameList[10][1] = "Scalar";




//	allResultsNameList.resize(5, vector<string>(2));
//
//	allResultsNameList[0][0] = "displacement";
//	allResultsNameList[0][1] = "Vector";
//
//	allResultsNameList[1][0] = "strain";
//	allResultsNameList[1][1] = "Matrix";
//
//	allResultsNameList[2][0] = "effective strain";
//	allResultsNameList[2][1] = "Scalar";
//
//	allResultsNameList[3][0] = "stress";
//	allResultsNameList[3][1] = "Matrix";
//
//	allResultsNameList[4][0] = "effective stress";
//	allResultsNameList[4][1] = "Scalar";
}

Data::Data(dbVector paramVec) {

	id = -1;
	nDims = 3;
	nDofs = 3;
	nNodes = 0;
	meshData = NULL;
	folderName = "";
	paramValues = paramVec;
	anchorPoint = -1;

	allResultsNameList.resize(11, vector<string>(2));

	allResultsNameList[0][0] = "displacement";
	allResultsNameList[0][1] = "Vector";

	allResultsNameList[1][0] = "strain";
	allResultsNameList[1][1] = "Matrix";

	allResultsNameList[2][0] = "effective strain";
	allResultsNameList[2][1] = "Scalar";

	allResultsNameList[3][0] = "fibre strain E_ff";
	allResultsNameList[3][1] = "Scalar";

	allResultsNameList[4][0] = "fibre strain E_tt";
	allResultsNameList[4][1] = "Scalar";

	allResultsNameList[5][0] = "fibre strain E_nn";
	allResultsNameList[5][1] = "Scalar";

	allResultsNameList[6][0] = "stress";
	allResultsNameList[6][1] = "Matrix";

	allResultsNameList[7][0] = "effective stress";
	allResultsNameList[7][1] = "Scalar";

	allResultsNameList[8][0] = "fibre stress S_ff";
	allResultsNameList[8][1] = "Scalar";

	allResultsNameList[9][0] = "fibre stress S_tt";
	allResultsNameList[9][1] = "Scalar";

	allResultsNameList[10][0] = "fibre stress S_nn";
	allResultsNameList[10][1] = "Scalar";

//	allResultsNameList.resize(5, vector<string>(2));
//
//		allResultsNameList[0][0] = "displacement";
//		allResultsNameList[0][1] = "Vector";
//
//		allResultsNameList[1][0] = "strain";
//		allResultsNameList[1][1] = "Matrix";
//
//		allResultsNameList[2][0] = "effective strain";
//		allResultsNameList[2][1] = "Scalar";
//
//		allResultsNameList[3][0] = "stress";
//		allResultsNameList[3][1] = "Matrix";
//
//		allResultsNameList[4][0] = "effective stress";
//		allResultsNameList[4][1] = "Scalar";
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read the FEM.res file using the FEMGeometry class
void Data::readMeshDataFile(InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	std::map<std::string, double> modelData;

	string meshFileName = folderName + "mesh.dat";
	string inputFileName = folderName + "input.dat";

	meshData = new FEMGeometryExt(InputData, modelData, meshFileName,
														inputFileName,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read the FEM.res file using the FEMGeometry class
void Data::readMeshDataFile(string& meshFileName,string& inputFileName,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	std::map<std::string, double> modelData;

	meshData = new FEMGeometryExt(InputData, modelData, meshFileName,
														inputFileName,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::delMeshData(){

	using namespace std;

	if(meshData != NULL){
		delete meshData;
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read the FEM.res file using the FEMGeometry class
void Data::readTransformedMeshDataFile(InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	std::map<std::string, double> modelData;

	string meshFileName = folderName + "mesh_map.dat";
	string inputFileName = folderName + "input.dat";

	meshData = new FEMGeometryExt(InputData, modelData, meshFileName,
														inputFileName,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Write to FEM.msh file using data from FEMGeometry class
void Data::writeMeshToMSHFile(string& outputFileName,InputFileData* InputData,
		ofstream& logFile){

	ofstream writeGraphRes(outputFileName);

	vector<ParticleExt>& ptcls = meshData->getNodesVec();
	vector<FEMElementExt>& elems = meshData->getNodesElemsVec();

	int dim = ptcls[0].getCoords().size();
	int nodesPerElem = elems[0].getNodes().size();

	string elemType;
	if(nodesPerElem == 3)
		elemType = "Tetrahedra";
	else if (nodesPerElem == 4)
		elemType = "Hexahedra";
	else{
		logFile << "In Data::writeMeshToMSHFile, nodesPerElem = "
				<< nodesPerElem << " and is not supported." << endl;
		cout << "In Data::writeMeshToMSHFile, nodesPerElem = "
				<< nodesPerElem << " and is not supported." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	string header_line = "MESH  dimension " + to_string(dim) + " ElemType "
			+ elemType + " Nnode " + to_string(nodesPerElem);

	writeGraphRes << header_line << endl;

	// Writing coordinates
	writeGraphRes << "Coordinates" << endl;
	for(int i = 0 ; i < ptcls.size(); i++){
		writeGraphRes << i << " ";

		dbVector coords = ptcls[i].getCoords();
		for(int j = 0 ; j < coords.size(); j++){
			writeGraphRes << coords[j] << " ";
		}

		writeGraphRes << endl;
	}
	writeGraphRes << "end coordinates" << endl;


	// Writing elements
	writeGraphRes << "Elements" << endl;
	for(int i = 0 ; i < elems.size(); i++){

		writeGraphRes << i << " ";

		intVector elemNodes = elems[i].getNodes();
		for(int j=0; j< elemNodes.size(); j++){
			writeGraphRes << elemNodes[j] << " ";
		}
		writeGraphRes << "0" << endl;

	}
	writeGraphRes << "end elements" << endl;



}

/*!****************************************************************************/
/*!****************************************************************************/
//! Print the displacement matrix
//void Data::printDispMatrix(ofstream& logFile) {
//	cout << "Data::printDispMatrix" << endl;
//
//	logFile << "Displacement matrix: " << endl;
//	for (int i = 0; i < dispMatrixList.size(); i++) {
//		for (int j = 0; j < dispMatrixList[i].size(); j++) {
//			logFile << dispMatrixList[i][j] << ", ";
//		}
//		logFile << endl;
//	}
//}

/*!****************************************************************************/
/*!****************************************************************************/
//! Plot the displacement field of a particular node
//void Data::plotNode(int dof, ofstream& logFile) {
//
//	logFile << "Displacement Node[" << dof << "] of Data[" << getId() << "]:"
//			<< endl;
//		for (int j = 0; j < dispMatrixList[dof - 1].size(); j++) {
//			logFile << dispMatrixList[dof - 1][j] << "\t";
//			logFile << j << "\t" << dispMatrixList[dof - 1][j] << endl;
//		}
//
//	logFile << endl;
//}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read the displacement data file form different file formats
void Data::readResultFile(InputFileData* InputData, ofstream& logFile) {

	int choice = InputData->getValue("resultFileType");

	switch (choice) {
	case 1:
		readResultFile_txtFormat(logFile);
		break;

	case 2: {

		string inputFileName = folderName + "fem.res";
		readResultFile_resFormat(inputFileName, logFile);

		break;
	}
	default:
		cout << "ERROR: Specified file format cannot be read" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readResultFile_resFormat_resultwise(string& inputFileName, ofstream& logFile) {

	int numResultFound = 0;

	for (int i = 0; i < allResultsNameList.size(); i++) {

		dbMatrix resultMatrix;
		dbVector resultStepVector;
		int numDofs;

		cout << "Reading result: " << allResultsNameList[i][0] << " .... ";
		logFile << "Reading result: " << allResultsNameList[i][0] << " .... ";

		bool isResultFound = readResultFile_resFormat_specificResult(
				inputFileName, allResultsNameList[i][0], resultMatrix, numDofs,
				resultStepVector, logFile);

		if (isResultFound == true) {

			bool isStepValueVecSimilar = true;

			// Record the first series of step values
			if (step_value_vec.size() == 0) {
				step_value_vec = resultStepVector;
			}
			// Else check if each step value are the same wrt to previously
			// recorded ones
			else {

				// Carrying out check
				if (step_value_vec.size() == resultStepVector.size()) {
					for (int i = 0; i < step_value_vec.size(); i++) {
						logFile << "Step[" << i <<"]: Comparing (" << resultStepVector[i] << ") Default(" << step_value_vec[i] << ")" << endl;
						if (fabs(step_value_vec[i] - resultStepVector[i]) > 1e-10) {
							isStepValueVecSimilar = false;
							break;
						}
					}
				} else {
					isStepValueVecSimilar = false;
					logFile << "WARNING: Default step size is not the same with the actual one" << endl;
				}

				// Discard result if step values are incompatible
				if (isStepValueVecSimilar == false) {

					cout << "Skipped(Check logFile for details)" << endl;
					logFile << "Skipped(Check logFile for details)" << endl;

					logFile << "WARNING: The step values of '"
							<< allResultsNameList[i][0] << "' in '"
							<< inputFileName
							<< "' are not the same compared to the other"
									" result types" << endl << endl ;

				}

			}

			// Record result if step values are compatible
		if(isStepValueVecSimilar == true) {

			numResultFound++;

			cout << "Found" << endl;
			logFile << "Found" << endl;

			resultNameList.push_back(allResultsNameList[i][0]);
			setResult(allResultsNameList[i][0].c_str(), resultMatrix);
			setResultDOF(allResultsNameList[i][0].c_str(), numDofs);

		}

#ifdef _DataDebugMode_
		printMatrix(resultMatrix, allResultsNameList[i][0].c_str(), logFile);
#endif

		}
		else{
			cout << "Not found" << endl;
			logFile << "Not found" << endl;
		}

	}

	if(numResultFound == 0){
		cout << "ERROR: No result has been found" << endl;
		logFile << "ERROR: No result has been found" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readResultFile_resFormat(string& inputFileName, ofstream& logFile) {

	dbMatrix resultMatrix;
	dbVector resultStepVector;
	map<double, map<string, dbMatrix> > resultsData;

	readResultFile_resFormat_allResult(inputFileName, resultNameList,
			resultsData, logFile);

	intVector totalDOFVec(resultNameList.size());

	// extract the number of degrees of freedom
	map<string, dbMatrix>::iterator it_result;
	for (int i = 0; i < resultNameList.size(); i++) {

		it_result = resultsData.begin()->second.find(resultNameList[i]);

		int nCol = it_result->second[0].size();
		int nRow = it_result->second.size();

		// Record num of DOFs per node
		setResultDOF(resultNameList[i].c_str(), nCol);

		// Record total number of DOFs
		totalDOFVec[i] = nCol * nRow;

	}

	int numSteps = resultsData.size();
	map<double, map<string, dbMatrix> >::iterator it_step, it_step_begin,
			it_step_end;
	it_step_begin = resultsData.begin();
	it_step_end = resultsData.end();

	// Save the result matrices
	for (int i = 0; i < resultNameList.size(); i++) {
		dbMatrix resultMatrix(totalDOFVec[i], dbVector(numSteps));

		int stepCounter = 0;
		for (it_step = it_step_begin; it_step != it_step_end; it_step++) {

			it_result = it_step->second.find(resultNameList[i]);
			dbMatrix& resultMatTemp = it_result->second;

			int dofCounter = 0;
			for (int r = 0; r < resultMatTemp.size(); r++) {
				for (int c = 0; c < resultMatTemp[r].size(); c++) {
					resultMatrix[dofCounter][stepCounter] = resultMatTemp[r][c];
					dofCounter++;
				}
			}

			stepCounter++;

			it_step->second.erase(it_result);

		}

		setResult(resultNameList[i].c_str(), resultMatrix);
	}

	// Save the step values
	dbVector stepVector;
	for (it_step = it_step_begin; it_step != it_step_end; it_step++) {
		stepVector.push_back(it_step->first);
	}
	setStepValueVec(stepVector);

//	printVector(this->getStepValueVec(), "StepValueVec", logFile);
//
//	vector<string> resultNameVec = this->getResultNameList();
//
//	for (int i = 0; i < resultNameVec.size(); i++) {
//		logFile << "--------------------------------------------------" << endl;
//		logFile << "Result: " << resultNameVec[i] << "\t nDOF: "
//				<< this->getResultDOF(resultNameVec[i].c_str()) << endl;
//		logFile << "--------------------------------------------------" << endl;
//
//		printMatrix(this->getResult(resultNameVec[i].c_str()), "", logFile);
//	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readResultFile_resFormat_allResult(string& inputFileName,
		vector<string>& resultNameList,
		map<double, map<string, dbMatrix> >& resultsData, ofstream& logFile) {

	resultNameList.clear();

	ifstream readFemRes(inputFileName.c_str());
	string line, result_name, analysis_name, my_result_type, my_location;
	double step_value;

	bool resultNameFound = false;

	if (!readFemRes) {
		logFile << "File " << inputFileName << " does not exist or is empty"
				<< endl;
		cout << "File " << inputFileName << " does not exist or is empty"
						<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Read file
	getline(readFemRes, line);
	while (readFemRes.good()) {

#ifdef _DataDebugMode_
		logFile << "Reading line => " << line << endl;
#endif

		//logFile << line << endl;
		istringstream is_line(line);
		string is_line_str;
		getline(is_line, is_line_str, '"');

#ifdef _DataDebugMode_
		logFile << "is_line_str: " << is_line_str << endl;
#endif

		if (is_line_str == "Result ") {
			getline(is_line, result_name, '"');

			// Extract analysis_name
			getline(is_line, analysis_name, '"');

			// Skip blank space
			getline(is_line, is_line_str, '"');

			// Extract line containing step_value, my_result_type
			// and my_location;
			getline(is_line, is_line_str, '"');
			istringstream is_line_ws(is_line_str);
			is_line_ws >> step_value >> my_result_type >> my_location;

#ifdef _DataDebugMode_
			logFile << "result_name: " << result_name << endl;
			logFile << "analysis_name: " << analysis_name << endl;
			logFile << "step_value: " << step_value << endl;
			logFile << "my_result_type: " << my_result_type << endl;
			logFile << "my_location: " << my_location << endl;
#endif

			// Find result name in resultNameList
			for (int i = 0; i < allResultsNameList.size(); i++) {
				if (result_name == allResultsNameList[i][0]) {

//					if (resultsData.begin() == resultsData.end()) {
//						resultNameList.push_back(result_name);
//						resultNameFound = true;
//						break;
//					} else {
//						for (int j = 0; j < resultNameList.size(); j++) {
//							if (result_name == resultNameList[j]) {
								resultNameFound = true;
//								break;
//							}
//						}
//					}

				}
			}

			// If the result name has been found in list, read the result field
			if (resultNameFound == true) {

				// --------------------------------------------------------
				// Read result field
				getline(readFemRes, line); // skip component names
				getline(readFemRes, line); // skip "values"

				// Read result field
				dbMatrix resultFieldMatrix;
				getline(readFemRes, line);
				while (line != "End Values") {

#ifdef _DataDebugMode_
					logFile << "Reading: " << line << " ";
#endif

					istringstream is(line);

					double nodeId;
					is >> nodeId;

					double result_field_value;
					dbVector resultFieldVector;

#ifdef _DataDebugMode_
					logFile << "Recording: " << nodeId << " ";
#endif
					while (is >> result_field_value) {
						resultFieldVector.push_back(result_field_value);
#ifdef _DataDebugMode_
						logFile << result_field_value << " ";
#endif
					}
#ifdef _DataDebugMode_
					logFile << endl;
#endif

					resultFieldMatrix.push_back(resultFieldVector);

					getline(readFemRes, line);

				}

				// --------------------------------------------------------
				// Save result field + stepvalue + result Name
				map<double, map<string, dbMatrix> >::iterator it_step;
				it_step = resultsData.find(step_value);

#ifdef _DataDebugMode_
				logFile << "Saving: " << result_name << " Timestep: " << step_value
						<< endl;
#endif

				if (it_step == resultsData.end()) {
					// add to resultsData
					map<string, dbMatrix> temp;
					temp[result_name] = resultFieldMatrix;
					resultsData[step_value] = temp;
				} else {
					// Look for result name
					map<string, dbMatrix>::iterator it_resultName;
					it_resultName = it_step->second.find(result_name);

					if (it_resultName == it_step->second.end()) {
						it_step->second[result_name] = resultFieldMatrix;
					} else {
						cout << "In " << inputFileName << ", result: '"
								<< result_name << "' has two result field"
										" of the same step value: '"
								<< step_value << endl;
						logFile << "In " << inputFileName << ", result: '"
								<< result_name << "' has two result field"
										" of the same step value: '"
								<< step_value << endl;
						MPI_Abort(MPI_COMM_WORLD, 1);
					}
				}

			}
			resultNameFound = false;

		}
		getline(readFemRes, line);
	}

	//

	map<double, map<string, dbMatrix> >::iterator it_stepVal,it_stepVal_beg,it_stepVal_end;
	it_stepVal_beg = resultsData.begin();
	it_stepVal_end = resultsData.end();
	map<string, dbMatrix>::iterator it_result,it_result_beg,it_result_end;

	vector<string> nameToDeleteList;
	intVector nameToDeleteListID;
	for (it_stepVal = it_stepVal_beg; it_stepVal != it_stepVal_end;
			it_stepVal++) {

		if (it_stepVal == it_stepVal_beg) {
			it_result_beg = it_stepVal->second.begin();
			it_result_end = it_stepVal->second.end();
#ifdef _DataDebugMode_
			logFile << "------------------------------------------" << endl;
			logFile << "Step: " << it_stepVal->first;
#endif
			for (it_result = it_result_beg; it_result != it_result_end;
					it_result++) {
#ifdef _DataDebugMode_
				logFile << it_result->first << endl;
#endif
				resultNameList.push_back(it_result->first);
			}
		} else {
#ifdef _DataDebugMode_
			logFile << "------------------------------------------" << endl;
			logFile << "Step: " << it_stepVal->first << endl;
			for (it_result = it_result_beg; it_result != it_result_end;
								it_result++) {
				logFile << it_result->first << endl;
			}
#endif
			for(int i=0; i < resultNameList.size(); i++){
				it_result_beg = it_stepVal->second.begin();
				it_result_end = it_stepVal->second.end();

#ifdef _DataDebugMode_
				logFile << "looking for: " << resultNameList[i] << endl;
#endif
				if(it_stepVal->second.find(resultNameList[i])
						== it_stepVal->second.end()){
#ifdef _DataDebugMode_
					logFile << "Not found !!" << endl;
#endif
					nameToDeleteList.push_back(resultNameList[i]);
					nameToDeleteListID.push_back(i);
				}
			}

		}
	}

	cout << "deleting incompatible result names" << endl;

	// Delete incompatible result names
#ifdef _DataDebugMode_
	logFile << "size of nameToDeleteList: " << nameToDeleteList.size() << endl;
#endif

	for(int i = nameToDeleteListID.size() ; i != 0; i--){

#ifdef _DataDebugMode_
		logFile << "resultNameList" << endl;
		for(int j = 0; j < resultNameList.size(); j++){
			logFile << "[" << j << "]" << resultNameList[j] << endl;
		}
		logFile << "deleting entry: " << nameToDeleteListID[i] << endl;
#endif

		resultNameList.erase(resultNameList.begin()+nameToDeleteListID[i]);

#ifdef _DataDebugMode_
		logFile << "resultNameList" << endl;
		for(int j = 0; j < resultNameList.size(); j++){
			logFile << "[" << j << "]" << resultNameList[j] << endl;
		}
#endif
	}

#ifdef _DataDebugMode_
	logFile << "Deleting elements in result map " << endl;
#endif
	for (it_stepVal = it_stepVal_beg; it_stepVal != it_stepVal_end;
				it_stepVal++) {

		it_result_beg = it_stepVal->second.begin();
		it_result_end = it_stepVal->second.end();

		for(int i=0; i < nameToDeleteList.size(); i++){

#ifdef _DataDebugMode_
			logFile << "Here" << endl;
#endif

			it_result = it_stepVal->second.find(nameToDeleteList[i]);

			if(it_result != it_result_end){
#ifdef _DataDebugMode_
				logFile << "Result excluded: " << it_result->first;
				cout << "Result excluded: " << it_result->first;
#endif
				it_stepVal->second.erase(it_result);
			}
		}


		it_result = it_result_end;
		it_result--;
		for(; it_result != it_result_beg; it_result--){
			bool isFound = false;
			for(int i = 0 ; i < resultNameList.size(); i++){
#ifdef _DataDebugMode_
				logFile << "Comparing resultNameList[i]: " << resultNameList[i]
				        << " with it_result->first: " << it_result->first << endl;
#endif
				if(resultNameList[i] == it_result->first){
#ifdef _DataDebugMode_
					logFile << "Found !!" << endl;
#endif
					isFound = true;
					break;
				}
			}

			if(isFound == false){
#ifdef _DataDebugMode_
				logFile << "Step[" << it_stepVal->first << "] Deleting result:" << it_result->first << endl;
#endif
				it_stepVal->second.erase(it_result);
			}
		}
	}

	// Close .res file
	readFemRes.close();
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readResultFile_txtFormat(ofstream& logFile) {
	double youngModulus, length, timeSteps;
	string temp;

	dbMatrix displacementMatrix;

	string line;
	string fileName = folderName + "displacement.txt";
	ifstream myfile(fileName.c_str());
	if (myfile.is_open()) {
		myfile >> temp >> youngModulus;
		//        cout<<"E:"<<youngModulus<<endl;
		myfile >> temp >> length;
		//        cout<<"Length:"<<length<<endl;
		myfile >> temp >> nDofs >> timeSteps;
		cout << "Dofs:" << nDofs << endl;
		cout << "timeSteps:" << timeSteps << endl;

		getline(myfile, line); // Get rid of ghost line

		resizeArray(displacementMatrix, nDofs, timeSteps);

		for (int i = 0; i < nDofs; i++) {
			getline(myfile, line);
			istringstream ss(line);
			for (int j = 0; j < timeSteps; j++) {
				string s;
				if (!getline(ss, s, ';'))
					break;
				cout << "reading:" << atof(s.c_str()) << endl;
				displacementMatrix[i][j] = atof(s.c_str());
			}
		}

		myfile >> temp >> nNodes >> nDims;
		cout << "nNodes:" << nNodes << endl;
		cout << "nDims:" << nDims << endl;

		getline(myfile, line); // Get rid of ghost line

		resizeArray(coordsList, nNodes, nDims);

		for (int i = 0; i < nNodes; i++) {
			getline(myfile, line);
			istringstream ss(line);
			for (int j = 0; j < nDims; j++) {
				string s;
				if (!getline(ss, s, ';'))
					break;
				cout << "reading:" << atof(s.c_str()) << endl;
				coordsList[i][j] = atof(s.c_str());
			}
		}

		myfile.close();
	} else {
		cout << "Unable to open Data file:" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	cout << "Data::readDispFile" << endl;

	cout << "~~~ Displacement matrix ~~~" << endl;
	for (int i = 0; i < displacementMatrix.size(); i++) {
		for (int j = 0; j < displacementMatrix[i].size(); j++) {
			cout << displacementMatrix[i][j] << ";";
		}
		cout << endl;
	}
	cout << endl;

	cout << "~~~ Coordinate List ~~~" << endl;
	for (int i = 0; i < coordsList.size(); i++) {
		for (int j = 0; j < coordsList[i].size(); j++) {
			cout << coordsList[i][j] << ";";
		}
		cout << endl;
	}
	cout << endl;

//	dispMatrixList.push_back(displacementMatrix);
	logFile << "Displacement not assigned to Data" << endl;
	cout << "Displacement not assigned to Data" << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readGraphFile_grfFormat(std::string& fileName, dbMatrix& grfMatrix,
		ofstream& logFile) {

	ifstream readGraphFile(fileName.c_str());
	string line;

	dbMatrix resMat;
	string blank_line = "";

	int counter = 0;
	if (readGraphFile.is_open()) {
		while (readGraphFile.good()) {
			getline(readGraphFile,line);
			string trimLine = line;
			trimLine.erase(remove(trimLine.begin(), trimLine.end(), '\t'), trimLine.end());
			trimLine.erase(remove(trimLine.begin(), trimLine.end(), ' '), trimLine.end());
			if (trimLine.compare(0,1,"#") != 0 && trimLine.compare(blank_line) != 1) {
//				logFile << "----------------------------------" << endl;
//				logFile << "Line: " << line << endl;
				istringstream ss(line);
				dbVector resultLine;
				for (;;) {
					string s;
					if (!getline(ss, s, ' '))
						break;
					else {
//						cout << "reading:" << atof(s.c_str()) << endl;
						resultLine.push_back(atof(s.c_str()));
					}
				}
				if(resultLine.size() > 0)
					resMat.push_back(resultLine);

			}
		}
		readGraphFile.close();
	} else {
		cout << "Unable to open " << fileName << " file:" << endl;
		logFile << "Unable to open " << fileName << " file:" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

//	logFile << "Data::readGraphFile_grfFormat" << endl;
//	printMatrix(resMat,"",logFile);

	// Transpose resMat
	grfMatrix.resize(resMat[0].size(),dbVector(resMat.size(),0));
	for(int i = 0; i < resMat.size(); i++){
		for(int j = 0; j < resMat[0].size(); j++){
			grfMatrix[j][i] = resMat[i][j];
		}
	}

#ifdef _DataDebugMode_
		logFile << "Data::readGraphFile_grfFormat" << endl;
		printMatrix(grfMatrix,"",logFile);
#endif


}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readVentriclesPVGraphResultFile(InputFileData* InputData,
		ofstream& logFile){

	readLeftVentriclePVGraphResultFile(InputData,logFile);
	readRightVentriclePVGraphResultFile(InputData,logFile);
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readLeftVentriclePVGraphResultFile(InputFileData* InputData,
		ofstream& logFile) {

	std::string vFileName = folderName + "volume_time-dirichletID_0.grf";
	dbMatrix vGraphResult;
	readGraphFile_grfFormat(vFileName, vGraphResult, logFile);
	setGraph("LVTimeSteps",vGraphResult[0]);
	leftCavityVolumes = vGraphResult[1];


	std::string pFileName = folderName + "load_time-loadID_0.grf";
	dbMatrix pGraphResult;
	readGraphFile_grfFormat(pFileName, pGraphResult, logFile);
	setGraph("LPTimeSteps",pGraphResult[0]);
	leftCavityPressures = pGraphResult[1];

#ifdef _DataDebugMode_
	printVector(leftCavityVolumes, "leftCavityVolumes", logFile);
	printVector(leftCavityPressures, "leftCavityPressures", logFile);
#endif
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from a simple text file format
void Data::readRightVentriclePVGraphResultFile(InputFileData* InputData,
		ofstream& logFile){

	std::string vFileName = folderName + "volume_time-dirichletID_1.grf";
	dbMatrix vGraphResult;
	readGraphFile_grfFormat(vFileName, vGraphResult, logFile);
	printMatrix(vGraphResult, "vGraphResult", logFile);
	setGraph("RVTimeSteps", vGraphResult[0]);
	rightCavityVolumes = vGraphResult[1];

	std::string pFileName = folderName + "load_time-loadID_1.grf";
	dbMatrix pGraphResult;
	readGraphFile_grfFormat(pFileName, pGraphResult, logFile);
	printMatrix(pGraphResult, "pGraphResult", logFile);
	setGraph("RPTimeSteps", pGraphResult[0]);
	rightCavityPressures = pGraphResult[1];

	printVector(rightCavityVolumes,"rightCavityVolumes",logFile);
	printVector(rightCavityPressures,"rightCavityPressures",logFile);
}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::saveGraphResultsToFile_grf_format(std::string outputFileName,
		dbVector& graphVecOne,dbVector& graphVecTwo,ofstream& logFile){

	printVector(graphVecOne,"graphVecOne",logFile);
	printVector(graphVecTwo,"graphVecTwo",logFile);

	if(graphVecOne.size() != graphVecTwo.size()){
		cout << "In Data::saveGraphResultsToFile_grf_format, graphVecOne and"
						" graphVecTwo are not of the same size" << endl;
		logFile << "In Data::saveGraphResultsToFile_grf_format, graphVecOne and"
				" graphVecTwo are not of the same size" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ofstream writeGraphRes(outputFileName);

	for(int i = 0 ; i < graphVecOne.size(); i++){
			writeGraphRes << graphVecOne[i] << " " << graphVecTwo[i] << endl;
	}

	writeGraphRes.close();

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from the GiD "res" file format
//void Data::readResultFile_resFormat(ofstream& logFile) {
////	cout<< "void Data::readDispFile_GiD_res_format function not yet implemented"
////			<< endl;
////	MPI_Abort(MPI_COMM_WORLD, 1);
//
//	//string inputFileName = "fem.res";
//	string inputFileName = folderName + "fem.res";
//	ifstream readFemRes(inputFileName.c_str());
//	string line;
//
//	dbMatrix listOfDisp;
//	dbVector dispVec;
//
//	bool is_name_found, is_value_found, is_endValue_found, nDofs_recorded,
//	nNodes_recorded;
//	is_name_found = false;
//	is_value_found = false;
//	is_endValue_found = false;
//	nDofs_recorded = false;
//	nNodes_recorded = false;
//
//	string result_name = "displacement";
//	string analysis_name;
//	double step_value = 0;
//	string my_result_type, my_location;
//
//	int numNodes, dofCounter, count = 0;
//
//	// Look for file
//	if (!readFemRes) {
//		logFile << "File " << inputFileName << " does not exist or is empty"
//				<< endl;
//		MPI_Abort(MPI_COMM_WORLD, 1);
//	}
//
//	// Read file
//	while (readFemRes.is_open()) {
//		getline(readFemRes, line);
//
//#ifdef _DataDebugMode_
//		logFile << "Reading line => " << line << endl;
//#endif
//
//		if (readFemRes.good()) {
//
//			//logFile << line << endl;
//			istringstream is_line(line);
//			string is_line_str;
//			getline(is_line, is_line_str, '"');
//
//#ifdef _DataDebugMode_
//			logFile << "is_line_str: " << is_line_str << endl;
//#endif
//
//			if (is_line_str == "Result ") {
//				getline(is_line, is_line_str, '"');
//				if (is_line_str == result_name) {
//					is_name_found = true;
//
//					// Extract analysis_name
//					getline(is_line, analysis_name, '"');
//
//					// Skip blank space
//					getline(is_line, is_line_str, '"');
//
//					// Extract line containing step_value, my_result_type
//					// and my_location;
//					getline(is_line, is_line_str, '"');
//					istringstream is_line_ws(is_line_str);
//					is_line_ws >> step_value >> my_result_type >> my_location;
//
//#ifdef _DataDebugMode_
//					logFile << "step_value: " << step_value << endl;
//					logFile << "my_result_type: " << my_result_type << endl;
//					logFile << "my_location: " << my_location << endl;
//#endif
//
//					// Record the step values
//					step_value_vec.push_back(step_value);
//				}
//			}
//
//			// If line containing step_value, my_result_type and my_location
//			// has been found, proceed with reading the degrees of freedom
//			if (is_name_found == true && is_value_found == true) {
//
//				// If end of degrees of freedom has been reached
//				if (line == "End Values") {
//
//					is_name_found = false;
//					is_value_found = false;
//
//					if (nNodes_recorded == false){
//						nNodes = numNodes;
//						nNodes_recorded = true;
//					}
//					else if (numNodes != nNodes) {
//						logFile << "nDofs: " << nDofs << " <=> " << "numNodes: "<< numNodes << endl;
//						cout << " Number of nodes across steps values "
//								"in " << inputFileName << " does not match"
//								<< endl;
//						logFile << " Number of nodes across step values "
//								"in " << inputFileName << " does not match"
//								<< endl;
//						MPI_Abort(MPI_COMM_WORLD, 1);
//					}
//
//#ifdef _DataDebugMode_
//					logFile << "===> is_name_found = false" << endl;
//					logFile << "===> is_value_found = false" << endl;
//#endif
//
//					count++;
//					listOfDisp.push_back(dispVec);
//					dispVec.clear();
//
//				} else {	string inputFileName = folderName + "fem.res";
//
//#ifdef _DataDebugMode_
//					logFile << "***********************************************"
//							"***" << endl;
//					logFile << line << endl;
//					logFile << "***********************************************"
//							"***" << endl;
//#endif
//
//					istringstream is(line);
//					double is_n;
//					is >> numNodes;		// Record the number of nodes
//					// vector<double> temp;
//
//					dofCounter = 0;
//					while (is >> is_n) {
//						dispVec.push_back(is_n);
//						dofCounter++;
//#ifdef _DataDebugMode_
//						logFile << "reading:" << is_n << endl;
//#endif
//					}
//
//					if(nDofs_recorded == false){
//						nDofs = dofCounter;
//						nDofs_recorded = true;
//					}
//					else if(dofCounter != nDofs){
//						cout << " Number of Degrees of freedoms for each "
//								"node in " << inputFileName << " does not match"
//								<< endl;
//						logFile << " Number of Degrees of freedoms for each "
//								"node in " << inputFileName << " does not match"
//								<< endl;
//						MPI_Abort(MPI_COMM_WORLD, 1);
//					}
//
//					//dispVec.push_back(temp);
//#ifdef _DataDebugMode_
//					logFile << "Size of dispVec: " << dispVec.size() << endl;
//#endif
//				}
//			}
//
//			// Indicate if the DOFs values have been found
//			if (line == "Values" && is_name_found == true) {
//				is_value_found = true;
//				numNodes = 0;
//			}
//
//		} else {
//			break;
//		}
//	}
//
//	// Close .res file
//	readFemRes.close();
//
//	if(listOfDisp.size() == 0){
//		logFile << "ERROR: " << result_name << " has not been found in "
//				<< inputFileName << endl;
//		cout << "ERROR: " << result_name << " has not been found in "
//				<< inputFileName << endl;
//		MPI_Abort(MPI_COMM_WORLD, 1);
//	}
//
//	// Set the step as columns
//	dispMatrixList.resize(listOfDisp[0].size(),dbVector(listOfDisp.size(),0));
//
//	for(int i = 0 ; i < listOfDisp.size(); i++){
//		for(int j = 0 ; j < listOfDisp[i].size(); j++){
//			dispMatrixList[j][i] = listOfDisp[i][j];
//		}
//	}
//
////	displacementMatrix = listOfDisp[0];
//
//#ifdef _DataDebugMode_
//	logFile << "Printing the results read from fem.res" << endl;
//	logFile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//	int counter = 0 ;
//	for (int i = 0; i < dispMatrixList[0].size(); i++) {
//		logFile << "Step Value :" << step_value_vec[i] << endl;
//		for (int j = 0; j < dispMatrixList.size(); j++) {
//			logFile << dispMatrixList[j][i] << "\t";
//			counter++;
//			if(nDofs == counter){
//				logFile << endl;
//				counter = 0 ;
//			}
//		}
//		logFile << endl
//				<< "--------------------------------------------------"
//				<< endl;
//	}
//#endif
//}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from the GiD "res" file format
bool Data::readResultFile_resFormat_specificResult(string& inputFileName,
		string& result_name,dbMatrix& resultMatrix,int& numDofs,
		dbVector& resultStepVector, ofstream& logFile) {

	ifstream readFemRes(inputFileName.c_str());
	string line;

	dbMatrix listOfDisp;
	dbVector dispVec;

	bool is_name_found, is_value_found, is_endValue_found, nDofs_recorded,
	nNodes_recorded;
	is_name_found = false;
	is_value_found = false;
	is_endValue_found = false;
	nDofs_recorded = false;
	nNodes_recorded = false;

	string analysis_name;
	double step_value = 0;
	string my_result_type, my_location;

	int numNodes, dofCounter, count = 0;

	numDofs = 0;

	// Look for file
	if (!readFemRes) {
		logFile << "File " << inputFileName << " does not exist or is empty"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Read file
	while (readFemRes.is_open()) {
		getline(readFemRes, line);

#ifdef _DataDebugMode_
		logFile << "Reading line => " << line << endl;
#endif

		if (readFemRes.good()) {

			//logFile << line << endl;
			istringstream is_line(line);
			string is_line_str;
			getline(is_line, is_line_str, '"');

#ifdef _DataDebugMode_
			logFile << "is_line_str: " << is_line_str << endl;
#endif

			if (is_line_str == "Result ") {
				getline(is_line, is_line_str, '"');
				if (is_line_str == result_name) {
					is_name_found = true;

					// Extract analysis_name
					getline(is_line, analysis_name, '"');

					// Skip blank space
					getline(is_line, is_line_str, '"');

					// Extract line containing step_value, my_result_type
					// and my_location;
					getline(is_line, is_line_str, '"');
					istringstream is_line_ws(is_line_str);
					is_line_ws >> step_value >> my_result_type >> my_location;

#ifdef _DataDebugMode_
					logFile << "step_value: " << step_value << endl;
					logFile << "my_result_type: " << my_result_type << endl;
					logFile << "my_location: " << my_location << endl;
#endif

					// Record the step values
					resultStepVector.push_back(step_value);
				}
			}

			// If line containing step_value, my_result_type and my_location
			// has been found, proceed with reading the degrees of freedom
			if (is_name_found == true && is_value_found == true) {

				// If end of degrees of freedom has been reached
				if (line == "End Values") {

					is_name_found = false;
					is_value_found = false;

					if (nNodes_recorded == false){
						nNodes = numNodes;
						nNodes_recorded = true;
					}
					else if (numNodes != nNodes) {
						logFile << "nDofs: " << numDofs << " <=> " << "numNodes: "<< numNodes << endl;
						cout << " Number of nodes across steps values "
								"in " << inputFileName << " does not match"
								<< endl;
						logFile << " Number of nodes across step values "
								"in " << inputFileName << " does not match"
								<< endl;
						MPI_Abort(MPI_COMM_WORLD, 1);
					}

#ifdef _DataDebugMode_
					logFile << "===> is_name_found = false" << endl;
					logFile << "===> is_value_found = false" << endl;
#endif

					count++;
					listOfDisp.push_back(dispVec);
					dispVec.clear();

				} else {

#ifdef _DataDebugMode_
					logFile << "***********************************************"
							"***" << endl;
					logFile << line << endl;
					logFile << "***********************************************"
							"***" << endl;
#endif

					istringstream is(line);
					double is_n;
					is >> numNodes;		// Record the number of nodes
					// vector<double> temp;

					dofCounter = 0;
					while (is >> is_n) {
						dispVec.push_back(is_n);
						dofCounter++;
#ifdef _DataDebugMode_
						logFile << "reading:" << is_n << endl;
#endif
					}

					if(nDofs_recorded == false){
						numDofs = dofCounter;
						nDofs_recorded = true;
					}
					else if(dofCounter != numDofs){
						cout << " Number of Degrees of freedoms for each "
								"node in " << inputFileName << " does not match"
								<< endl;
						logFile << " Number of Degrees of freedoms for each "
								"node in " << inputFileName << " does not match"
								<< endl;
						MPI_Abort(MPI_COMM_WORLD, 1);
					}

					//dispVec.push_back(temp);
#ifdef _DataDebugMode_
					logFile << "Size of dispVec: " << dispVec.size() << endl;
#endif
				}
			}

			// Indicate if the DOFs values have been found
			if (line == "Values" && is_name_found == true) {
				is_value_found = true;
				numNodes = 0;
			}

		} else {
			break;
		}
	}

	// Close .res file
	readFemRes.close();

	bool isResultFound = true;

	if (listOfDisp.size() == 0) {
		logFile << "WARNING: " << result_name << " has not been found in "
				<< inputFileName << endl;
		isResultFound = false;
	} else {

		// Set the step as columns
		resultMatrix.resize(listOfDisp[0].size(),
				dbVector(listOfDisp.size(), 0));

		for (int i = 0; i < listOfDisp.size(); i++) {
			for (int j = 0; j < listOfDisp[i].size(); j++) {
				resultMatrix[j][i] = listOfDisp[i][j];
			}
		}

	}

	return isResultFound;



#ifdef _DataDebugMode_
	if (isResultFound == true) {
		logFile << "Printing the results read from fem.res" << endl;
		logFile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		int counter = 0;
		for (int i = 0; i < resultMatrix[0].size(); i++) {
			logFile << "Step Value :" << resultStepVector[i] << endl;
			for (int j = 0; j < resultMatrix.size(); j++) {
				logFile << resultMatrix[j][i] << "\t";
				counter++;
				if (nDofs == counter) {
					logFile << endl;
					counter = 0;
				}
			}
			logFile << endl
					<< "--------------------------------------------------"
					<< endl;
		}
	}
#endif
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Read and Extract results from the GiD "res" file format
vector<vector<string> > Data::readResultFile_resFormat_HeadersOnly(
		string& inputFileName, ofstream& logFile) {

	// string inputFileName = folderName + "fem.res";
	ifstream readFemRes(inputFileName.c_str());
	string line;

	string resultName;

	vector<vector<string> > resultNameList;

	bool repeat_block = false;
//	bool is_name_found, is_value_found, is_endValue_found, nDofs_recorded,
//			nNodes_recorded;
//	is_name_found = false;
//	is_value_found = false;
//	is_endValue_found = false;
//	nDofs_recorded = false;
//	nNodes_recorded = false;

	string analysis_name;
	double step_value = 0;
	string my_result_type, my_location;

//	int numNodes, dofCounter, count = 0;

	// Look for file
	if (!readFemRes) {
		logFile << "File " << inputFileName << " does not exist or is empty"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Read file
	while (readFemRes.is_open()) {
		getline(readFemRes, line);

#ifdef _DataDebugMode_
		logFile << "Reading line => " << line << endl;
#endif

		if (readFemRes.good()) {

			//logFile << line << endl;
			istringstream is_line(line);
			string is_line_str;
			getline(is_line, is_line_str, '"');

#ifdef _DataDebugMode_
			logFile << "is_line_str: " << is_line_str << endl;
#endif

			if (is_line_str == "Result ") {
				getline(is_line, resultName, '"');

				// Extract analysis_name
				getline(is_line, analysis_name, '"');

				// Skip blank space
				getline(is_line, is_line_str, '"');

				// Extract line containing step_value, my_result_type
				// and my_location;
				getline(is_line, is_line_str, '"');
				istringstream is_line_ws(is_line_str);
				is_line_ws >> step_value >> my_result_type >> my_location;

#ifdef _DataDebugMode_
				logFile << "step_value: " << step_value << endl;
				logFile << "my_result_type: " << my_result_type << endl;
				logFile << "my_location: " << my_location << endl;
#endif

				if (my_location == "OnNodes") {

					bool isInList = false;
					for (int i = 0; resultNameList.size(); i++) {
						if (resultNameList[i][0] == my_location) {
							isInList = true;
							repeat_block = true;
						}

					}

					if (isInList == false) {
						resultNameList.push_back(vector<string>(2));
						resultNameList[resultNameList.size() - 1][1] =
								my_location;
						resultNameList[resultNameList.size() - 1][1] =
								my_result_type;
					}

				}
			}
		}

		// Repeating result block found, stop reading file!
		if(repeat_block == true){
			break;
		}
	}

	// Close .res file
	readFemRes.close();


	return resultNameList;

}




/*!****************************************************************************/
/*!****************************************************************************/
//! Save the displacement matrix to the a specific file format
//void Data::saveDispFile(ofstream& logFile) {
//
//	int choice = 1;
//
//	switch (choice) {
//	case 1:
//		//Save matrix to .res file format
//		saveDispFile_res_format("fem_orion.res",logFile);
//		break;
//
//	default:
//		logFile << "ERROR: Specified file format cannot be written to" << endl;
//		MPI_Abort(MPI_COMM_WORLD, 1);
//		break;
//	}
//
//}

/*!****************************************************************************/
/*!****************************************************************************/
// Writing generated displacement matrix to files
//void Data::saveDispFile_res_format(const char* outputFileName,ofstream& logFile) {
//
//	system("cp -f fem.msh fem_orion.msh");
//
//	ofstream readFemRes(outputFileName);
//
//	string headerLine = "GiD Post Results File 1.0";
//
//	// Result details
//	string result_name = "displacement";
//	string analysis_name = " ";
//	string my_type = "Vector";
//	string my_location = "OnNodes";
//
//	// Setting up the result type header line
//	string my_type_header = "ComponentNames \"DOF_1\",";
//	if (nDofs > 1) {
//		for (int i = 2; i < nDofs + 1; i++) {
//			stringstream i_ss;
//			i_ss << i;
//			my_type_header += " \"DOF_" + i_ss.str() + "\",";
//		}
//	}
//
//	if (readFemRes.is_open()) {
//
//		readFemRes << headerLine << endl;
//
//#ifdef _DataDebugMode_
//		logFile << "size of stepValue vector: " << step_value_vec.size() << endl;
//		logFile << "nDofs: " << nDofs << endl;
//		logFile << "nNodes: " << nNodes << endl;
//		logFile << "dispMatrixList: " << dispMatrixList.size() << " x "
//				<< dispMatrixList[0].size() << endl;
//#endif
//
//		for (int j = 0; j < step_value_vec.size(); j++) {
//
//			// Setting up the Result line
//			readFemRes << "Result " << "\"" << result_name << "\" \""
//					<< analysis_name << "\" " << step_value_vec[j] << " "
//					<< my_type << " " << my_location << endl;
//
//			readFemRes << my_type_header << endl;
//			readFemRes << "Values" << endl;
//
//			int nodeNum = 1,counter = 0;
//			for (int r = 0; r < dispMatrixList.size(); r++) {
//
//				if(counter == 0){
//					readFemRes << nodeNum;
//				}
//
//				readFemRes << " " << dispMatrixList[r][j];
//				counter++;
//
//				if(counter == nDofs){
//					readFemRes << endl;
//					counter = 0;
//					nodeNum++;
//				}
//			}
//
//			readFemRes << "End Values" << endl;
//		}
//		readFemRes.close();
//	}
//}

/*!****************************************************************************/
/*!****************************************************************************/
//! Save the displacement matrix to the a specific file format
void Data::saveResultsToFile(InputFileData* InputData,ofstream& logFile) {

	int choice = 1;

	switch (choice) {
	case 1:
	{
		//Save matrix to .res file format
		string str = folderName + "fem_orion.res";
		saveAllResultsToFile_res_format(str.c_str(),InputData,logFile);
		break;
	}
	default:
		logFile << "ERROR: Specified file format cannot be written to" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

}

void Data::saveResultsToFile_step(InputFileData* InputData,ofstream& logFile) {

	int choice = 1;

	switch (choice) {
	case 1:
	{
		//Save matrix to .res file format
		string str = folderName + "fem_orion_step.res";
		saveAllResultsToFile_res_format(str.c_str(),InputData,logFile);
		break;
	}
	default:
		logFile << "ERROR: Specified file format cannot be written to" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// Writing generated displacement matrix to files
void Data::saveAllResultsToFile_res_format(const char* outputFileName,
		InputFileData* InputData,ofstream& logFile) {

	std::cout << "Saving fem_orion.res to: " << outputFileName << endl;
	saveResultsToFile_res_format(outputFileName,resultNameList,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
// Writing generated displacement matrix to files
void Data::saveResultsToFile_res_format(const char* outputFileName,
		vector<string>& saveResultNameList, InputFileData* InputData,
		ofstream& logFile) {

#ifdef _DataDebugMode_
	logFile << "****** All Results available ******" << endl;

	map<string, dbMatrix>::iterator it_dbMatrix;
	map<string, int>::iterator it_int;

for(it_dbMatrix = resultList.begin();it_dbMatrix != resultList.end();it_dbMatrix++){
	logFile << "Number of DOFs: "
			<< resultDofList.find(it_dbMatrix->first)->second << endl;
	printMatrix(it_dbMatrix->second,it_dbMatrix->first.c_str(),logFile);

	logFile << "------------------------------------------------------" << endl;
}

#endif


	vector<vector<string> > ResultNameAndTypeToBeSavedList;

	for (int i = 0; i < saveResultNameList.size(); i++) {
		for (int j = 0; j < allResultsNameList.size(); j++) {
			if (saveResultNameList[i] == allResultsNameList[j][0]) {
				ResultNameAndTypeToBeSavedList.resize(
						ResultNameAndTypeToBeSavedList.size() + 1);
				ResultNameAndTypeToBeSavedList[ResultNameAndTypeToBeSavedList.size()
						- 1] = allResultsNameList[j];
				break;
			}
		}
	}

	ofstream readFemRes(outputFileName);

	string headerLine = "GiD Post Results File 1.0";
	string my_location = "OnNodes";

	logFile << "step_value_vec: " << step_value_vec.size() << endl;

	if (readFemRes.is_open()) {

		readFemRes << headerLine << endl;

		int startIndex = 0;
		if(InputData->getValue("insertZeroResultFields") == 1)
			startIndex = 1;

		for (int j = startIndex; j < step_value_vec.size(); j++) {
			for (int k = 0; k < ResultNameAndTypeToBeSavedList.size(); k++) {

				// Result details
				string analysis_name = " ";
				int& numDofs = this->getResultDOF(
						ResultNameAndTypeToBeSavedList[k][0].c_str());

				// Setting up the result type header line
				string my_type_header = "ComponentNames \"DOF_1\",";
				if (numDofs > 1) {
					for (int i = 2; i < numDofs + 1; i++) {
						stringstream i_ss;
						i_ss << i;
						my_type_header += " \"DOF_" + i_ss.str() + "\",";
					}
				}

				// Setting up the Result line
				readFemRes << "Result " << "\""
						<< ResultNameAndTypeToBeSavedList[k][0] << "\" \""
						<< analysis_name << "\" " << step_value_vec[j] << " "
						<< ResultNameAndTypeToBeSavedList[k][1] << " "
						<< my_location << endl;

				readFemRes << my_type_header << endl;

				dbMatrix& resultMatrix =
						this->getResult(ResultNameAndTypeToBeSavedList[k][0].c_str());

#ifdef _DataDebugMode_
				logFile << "size of stepValue vector: " << step_value_vec.size()
						<< endl;
				logFile << "numDofs: " << numDofs << endl;
				logFile << "nNodes: " << nNodes << endl;
				logFile << "resultMatrix: " << resultMatrix.size() << " x "
						<< resultMatrix[0].size() << endl;
#endif


				readFemRes << "Values" << endl;

				int nodeNum = 1, counter = 0;
				for (int r = 0; r < resultMatrix.size(); r++) {

					if (counter == 0) {
						readFemRes << nodeNum;
					}

					readFemRes << " " << resultMatrix[r][j];
					counter++;

					if (counter == numDofs) {
						readFemRes << endl;
						counter = 0;
						nodeNum++;
					}
				}

				readFemRes << "End Values" << endl;
			}
		}
		readFemRes.close();
	}
}


/*!****************************************************************************/
/*!****************************************************************************/
//void Data::assignDispToParticles(InputFileData* InputData, ofstream& logFile){
//
//#ifdef _DataDebugMode_
//	logFile << "Before assigning disp to particles" << endl;
//		printMatrix(dispMatrixList,"",logFile);
//#endif
//
//	vector<ParticleExt>& particleList = meshData->getNodesVec();
//
//	if(dispMatrixList.size() > 0 && particleList.size()*nDofs == dispMatrixList.size()){
//
//		for(int i=0; i< particleList.size(); i++){
//			particleList[i].getStepDOFMat().resize(nDofs,dbVector(step_value_vec.size()));
//
//			for(int j=0; j<step_value_vec.size(); j++){
//				for(int k=0; k<nDofs; k++){
//					particleList[i].getStepDOFMat()[k][j] = dispMatrixList[(i*nDofs)+k][j];
//				}
//			}
//
//#ifdef _DataDebugMode_
//
//			logFile << "******** Particles stepDOFs ********" << endl;
//			for(int l=0; l<particleList.size(); l++){
//				logFile << " For particle[" << i << "]" << endl;
//				printMatrix(particleList[i].getStepDOFMat(),"",logFile);
//				logFile << endl;
//			}
//
//#endif
//
//		}
//	}
//
//}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::assignResultToParticles(const char* resultName,
		InputFileData* InputData, ofstream& logFile) {

	vector<ParticleExt>& particleList = meshData->getNodesVec();

	dbMatrix& resultMatrix = this->getResult(resultName);
	int numDofs = this->getResultDOF(resultName);
	int numSteps = step_value_vec.size();

#ifdef _DataDebugMode_
	logFile << "Before assigning '" << resultName << "' to particles" << endl;
	printMatrix(resultMatrix, "", logFile);
#endif

	if (resultMatrix.size() > 0
			&& particleList.size() * numDofs == resultMatrix.size()) {

		for (int i = 0; i < particleList.size(); i++) {

			particleList[i].getStepDOFMat().resize(numDofs,
					dbVector(numSteps, 0));

			for (int j = 0; j < step_value_vec.size(); j++) {
				for (int k = 0; k < numDofs; k++) {
					particleList[i].getStepDOFMat()[k][j] = resultMatrix[(i
							* numDofs) + k][j];
				}
			}
		}

#ifdef _DataDebugMode_

		logFile << "******** Particles stepDOFs ********" << endl;
		for(int l=0; l<particleList.size(); l++) {
			std::string msg = "For particle[";
			ostringstream convert;
			convert << "For particle[" << l << "]";
			msg = convert.str();
			printMatrix(particleList[l].getStepDOFMat(),msg.c_str(),logFile);
			logFile << endl;
		}

#endif

	} else {
		logFile << "Result matrix is empty or the size of resultMatrix["
				<< resultMatrix.size() << "] and particleList.size()*numDofs["
				<< particleList.size() * numDofs << "] does not correspond."
				<< endl;
		cout << "Result matrix is empty or the size of resultMatrix["
				<< resultMatrix.size() << "] and particleList.size()*numDofs["
				<< particleList.size() * numDofs << "] does not correspond."
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::setResult(const char* resultName, dbMatrix result) {

	string str(resultName);

	map<string, dbMatrix>::iterator it = resultList.find(str);

	if (it == resultList.end())
		resultList[str] = result ;
	else{
		cout << "ERROR: In Data::setResult, '" << resultName
				<< "' already exists in resultList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
dbMatrix& Data::getResult(const char* resultName) {

	string str(resultName);

	map<string, dbMatrix>::iterator it = resultList.find(str);

	if (it != resultList.end())
		return it->second;
	else {
		cout << "ERROR: In Data::getResult, '" << resultName
				<< "' does not exist in resultList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
dbVector& Data::getResultAcrossSteps(const char* resultName,int dofID) {

	string str(resultName);

	map<string, dbMatrix>::iterator it = resultList.find(str);

	if (it == resultList.end()){
		cout << "ERROR: In Data::getResult, '" << resultName
				<< "' does not exist in resultList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	dbMatrix& result = it->second;
	return result[dofID];
}

/*!****************************************************************************/
/*!****************************************************************************/
dbVector& Data::getResultAcrossDOFs(const char* resultName,int step) {

	string str(resultName);

	map<string, dbMatrix>::iterator it = resultList.find(str);

	if (it == resultList.end()){
		cout << "ERROR: In Data::getResult, '" << resultName
				<< "' does not exist in resultList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	dbMatrix& result = it->second;
	dbVector resultVector(result.size(),0);

	for(int i=0; i < result.size() ; i++){
		resultVector[i] = result[i][step];
	}

	return resultVector;
}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::deleteResult(const char* resultName){

	string str(resultName);

	map<string, dbMatrix>::iterator it = resultList.find(str);

	if (it != resultList.end())
		resultList.erase(it);
	else {
		cout << "ERROR: In Data::deleteResult, '" << resultName
				<< "' does not exist in resultList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::setResultDOF(const char* resultName, int& dof) {

	string str(resultName);

	map<string, int>::iterator it = resultDofList.find(str);

	if (it == resultDofList.end())
		resultDofList[str] = dof;
	else {
		cout << "ERROR: In Data::setResultDOF, '" << resultName
				<< "' already exists in resultDofList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
int& Data::getResultDOF(const char* resultName) {

	string str(resultName);

	map<string, int>::iterator it = resultDofList.find(str);

	if (it != resultDofList.end())
		return it->second;
	else {
		cout << "ERROR: In Data::getResultDOF, '" << resultName
				<< "' does not exist in resultDofList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::deleteResultDOF(const char* resultName){

	string str(resultName);

	map<string, int>::iterator it = resultDofList.find(str);

	if (it != resultDofList.end())
		resultDofList.erase(it);
	else {
		cout << "ERROR: In Data::deleteResultDOF, '" << resultName
				<< "' does not exist in resultDofList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// DO NOT USE: Data::calcCavityVolumes is based on the old input file. Hence,
// Data::calcCavityVolumes is superseeded by Data::calcLeftAndRightCavityVolumes
void Data::calcCavityVolumes(InputFileData* InputData,ofstream& logFile){

	using namespace std;

	cout << "Calculating the cavities volumes" << endl;
	logFile << "Calculating the cavities volumes" << endl;

	dbMatrix& resultMatrix = this->getResult("displacement");
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

	}


	FEMGeometry* myData_MeshData =  this->getMeshData()->getFEMGeoData();
	InputFileData* myData_InputData =  this->getMeshData()->getFEMInputData();

	logFile << "In Data::calcCavityVolumes, InputData->getLineDispBoundConds()"
				"problems have not been resolved yet." << endl;
	cout << "In Data::calcCavityVolumes, InputData->getLineDispBoundConds()"
				"problems have not been resolved yet." << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);

	dbMatrix allLoadMat;
	printMatrix(allLoadMat,"allLoadMat",logFile);

	dbMatrix loadMat;

	int choice = 1;
	if (choice == 0) {
		for (int i = 0; i < allLoadMat.size(); i++) {
			if (allLoadMat[i][4] == (double) 1.0) {
				loadMat.push_back(allLoadMat[i]);
			}
		}
	} else {
		loadMat = allLoadMat;
	}

	printMatrix(loadMat,"loadMat",logFile);

	logFile << "In Data::calcCavityVolumes, InputData->getLineDispBoundConds()"
			"problems have not been resolved yet." << endl;
	cout << "In Data::calcCavityVolumes, InputData->getLineDispBoundConds()"
			"problems have not been resolved yet." << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);

	dbMatrix lineLoadMat;
	printMatrix(lineLoadMat,"lineLoadMat",logFile);

	int numOfFaces = loadMat.size();

	int numOfSteps = resultMatList.size();
	dbVector volumeVec(numOfSteps,0);

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
			dbVector& v2 = ptcls[loadMat[m][2]-1].getCoords();
			dbVector& v3 = ptcls[loadMat[m][3]-1].getCoords();

			volume += ((v2[y] - v1[y])*(v3[z] - v1[z]) -
					   (v2[z] - v1[z])*(v3[y] - v1[y]))*
					   (v1[x] + v2[x] + v3[x] );
		}

		volumeVec[i] = abs(volume / 6.0);

		logFile << "[" << i << "] Geometry volume = " << volumeVec[i] << endl;
	}

//	dbMatrix& graphResult = this->getGraphResultList();
//	printMatrix(graphResult,"graphResult", logFile);

	string graphFileName = folderName + "defVolLoad_myAlgo.grf";

	ofstream writeGraph(graphFileName);
	writeGraph.precision(15);

	if(writeGraph.good()){
		//for(int i = 0 ; i < graphResult[0].size(); i++){
			//writeGraph << volumeVec[i] << " " << graphResult[0][i] << endl;
		//}

		for(int i = 0 ; i < step_value_vec.size(); i++){
			writeGraph << volumeVec[i] << " " << step_value_vec[i] << endl;
		}
	}

	writeGraph.close();
}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::calcLeftAndRightCavityVolumes(InputFileData* InputData, ofstream& logFile) {

	cout << "Calculating cavity volumes" << endl;

	calcLeftCavityVolumes(InputData,logFile);
	calcRightCavityVolumes(InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::calcLeftCavityVolumes(InputFileData* InputData, ofstream& logFile) {

	using namespace std;

//	logFile << "Calculating left cavity volumes" << endl;

	dbMatrix& resultMatrix = this->getResult("displacement");
	//printMatrix(resultMatrix,"displacement Matrix",logFile);

	// Convert each result at each step from a vector to a matrix
	vector<dbMatrix> resultMatList(resultMatrix[0].size(),
				dbMatrix(resultMatrix.size() / 3, dbVector(3,0)));
	for (int i = 1; i < resultMatrix[0].size(); i++) {
		dbMatrix& mat = resultMatList[i];

		int row = 0;
		int col = 0;
		for (int j = 0; j < resultMatrix.size(); j++) {

			mat[row][col] = resultMatrix[j][i];
			col++;

			if (col > 2) {
				col = 0;
				row++;
			}
		}
	}

	// Extract the surface of the left ventricle
	intVector surfaceIDList;
	intMatrix surfaceNodes;
	meshData->getSpecificSurfaceNodes("Endocaridum_LV_Surface", surfaceIDList,
			surfaceNodes, InputData, logFile);

//	logFile << "List of surfaces" << endl;
//	for(int i=0; i<surfaceIDList.size(); i++){
//		logFile << "surfaceIDList: " << surfaceNodes[i][0] << ", " << surfaceNodes[i][1] << ", " << surfaceNodes[i][2] << endl;
//
//	}

	int numOfSteps = resultMatList.size();
	dbVector volumeVec(numOfSteps, 0);


	for (int i = 0; i < numOfSteps; i++) {

		dbMatrix& rMatrix = resultMatList[i];
		vector<ParticleExt> ptcls = meshData->getNodesVec();

		// Update the coordinates
		for (int j = 0; j < ptcls.size(); j++) {
			for (int k = 0; k < 3; k++) {
				ptcls[j].getCoord(k) += rMatrix[j][k];
			}
		}

		// Compute the volume of the cavity
		double volume = 0;

		int x = 0;
		int y = 1;
		int z = 2;

		int numOfsurfaces = surfaceIDList.size();
		for (int m = 0; m < numOfsurfaces; m++) {

			dbVector& v1 = ptcls[surfaceNodes[m][0] - 1].getCoords();
			dbVector& v2 = ptcls[surfaceNodes[m][1] - 1].getCoords();
			dbVector& v3 = ptcls[surfaceNodes[m][2] - 1].getCoords();

			volume += ((v2[y] - v1[y]) * (v3[z] - v1[z])
					- (v2[z] - v1[z]) * (v3[y] - v1[y]))
					* (v1[x] + v2[x] + v3[x]);
		}

		volumeVec[i] = abs(volume / 6.0);
	}

	string timeGraphFileName = folderName + "time_volume-LV.grf";
	cout << "Writing to: " << timeGraphFileName << endl;

	ofstream writeTimeGraph(timeGraphFileName);
	writeTimeGraph.precision(15);
	writeTimeGraph.setf( std::ios::scientific, std:: ios::floatfield );

	string pressureGraphFileName = folderName + "pressure_volume-LV.grf";
	cout << "Writing to: " << pressureGraphFileName << endl;

	ofstream writePressureGraph(pressureGraphFileName);
	writePressureGraph.precision(15);
	writePressureGraph.setf( std::ios::scientific, std:: ios::floatfield );


	if (writeTimeGraph.good() && writePressureGraph.good()) {
		for (int i = 0; i < step_value_vec.size(); i++) {
			writeTimeGraph << step_value_vec[i] << " " << volumeVec[i] << endl;
			writePressureGraph << volumeVec[i] << " " << leftCavityPressures[i] << endl;
		}
	}

	writeTimeGraph.close();
	writePressureGraph.close();

	leftCavityVolumes = volumeVec;

}

void Data::calcLeftCavityVolumes_step(InputFileData* InputData, ofstream& logFile) {

	using namespace std;

//	logFile << "Calculating left cavity volumes" << endl;

	dbMatrix& resultMatrix = this->getResult("displacement");

	// Convert each result at each step from a vector to a matrix
	vector<dbMatrix> resultMatList(resultMatrix[0].size(),
				dbMatrix(resultMatrix.size() / 3, dbVector(3,0)));
	for (int i = 1; i < resultMatrix[0].size(); i++) {
		dbMatrix& mat = resultMatList[i];

		int row = 0;
		int col = 0;
		for (int j = 0; j < resultMatrix.size(); j++) {

			mat[row][col] = resultMatrix[j][i];
			col++;

			if (col > 2) {
				col = 0;
				row++;
			}
		}
	}

	// Extract the surface of the left ventricle
	intVector surfaceIDList;
	intMatrix surfaceNodes;
	meshData->getSpecificSurfaceNodes("Endocaridum_LV_Surface", surfaceIDList,
			surfaceNodes, InputData, logFile);

	int numOfSteps = resultMatList.size();
	dbVector volumeVec(numOfSteps, 0);


	for (int i = 0; i < numOfSteps; i++) {

		dbMatrix& rMatrix = resultMatList[i];
		vector<ParticleExt> ptcls = meshData->getNodesVec();

		// Update the coordinates
		for (int j = 0; j < ptcls.size(); j++) {
			for (int k = 0; k < 3; k++) {
				ptcls[j].getCoord(k) += rMatrix[j][k];
			}
		}

		// Compute the volume of the cavity
		double volume = 0;

		int x = 0;
		int y = 1;
		int z = 2;

		int numOfsurfaces = surfaceIDList.size();
		for (int m = 0; m < numOfsurfaces; m++) {

			dbVector& v1 = ptcls[surfaceNodes[m][0] - 1].getCoords();
			dbVector& v2 = ptcls[surfaceNodes[m][1] - 1].getCoords();
			dbVector& v3 = ptcls[surfaceNodes[m][2] - 1].getCoords();

			volume += ((v2[y] - v1[y]) * (v3[z] - v1[z])
					- (v2[z] - v1[z]) * (v3[y] - v1[y]))
					* (v1[x] + v2[x] + v3[x]);
		}

		volumeVec[i] = abs(volume / 6.0);
	}

	string timeGraphFileName = folderName + "time_volume-LV_step.grf";
	logFile << "Writing to: " << timeGraphFileName << endl;
	cout << "Writing to: " << timeGraphFileName << endl;

	ofstream writeTimeGraph(timeGraphFileName);
	writeTimeGraph.precision(15);
	writeTimeGraph.setf( std::ios::scientific, std:: ios::floatfield );

	string pressureGraphFileName = folderName + "pressure_volume-LV_step.grf";
	logFile << "Writing to: " << pressureGraphFileName << endl;
	cout << "Writing to: " << pressureGraphFileName << endl;

	ofstream writePressureGraph(pressureGraphFileName);
	writePressureGraph.precision(15);
	writePressureGraph.setf( std::ios::scientific, std:: ios::floatfield );


	if (writeTimeGraph.good() && writePressureGraph.good()) {
		for (int i = 0; i < step_value_vec.size(); i++) {
			writeTimeGraph << step_value_vec[i] << " " << volumeVec[i] << endl;
			writePressureGraph << volumeVec[i] << " " << leftCavityPressures[i] << endl;
		}
	}

	writeTimeGraph.close();
	writePressureGraph.close();

	leftCavityVolumes = volumeVec;

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::calcRightCavityVolumes(InputFileData* InputData, ofstream& logFile) {

	using namespace std;

//	logFile << "Calculating right cavity volumes" << endl;

	dbMatrix& resultMatrix = this->getResult("displacement");
//	printMatrix(resultMatrix,"displacement Matrix",logFile);
//	MPI_Abort(MPI_COMM_WORLD, 1);

	// Convert each result at each step from a vector to a matrix
	vector<dbMatrix> resultMatList(resultMatrix[0].size(),
			dbMatrix(resultMatrix.size() / 3, dbVector(3,0)));
	for (int i = 0; i < resultMatrix[0].size(); i++) {
		dbMatrix& mat = resultMatList[i];

		int row = 0;
		int col = 0;
		for (int j = 0; j < resultMatrix.size(); j++) {

			mat[row][col] = resultMatrix[j][i];
			col++;

			if (col > 2) {
				col = 0;
				row++;
			}
		}
	}

	// Extract the surface of the left ventricle
	intVector surfaceIDList_EndoRV;
	intMatrix surfaceNodes_EndoRV;
	meshData->getSpecificSurfaceNodes("Endocaridum_RV_Surface",
			surfaceIDList_EndoRV,surfaceNodes_EndoRV, InputData, logFile);

	intVector surfaceIDList_FreeWall;
	intMatrix surfaceNodes_FreeWall;
	meshData->getSpecificSurfaceNodes("Endocaridum_Freewall_Surface",
			surfaceIDList_FreeWall,surfaceNodes_FreeWall, InputData, logFile);

	// Merge Endocaridum_RV_Surface and Endocaridum_Freewall_Surface data
	intVector surfaceIDList;
	intMatrix surfaceNodes;

	surfaceIDList.reserve( surfaceIDList_EndoRV.size() + surfaceIDList_FreeWall.size() );
	surfaceIDList.insert( surfaceIDList.end(), surfaceIDList_EndoRV.begin(), surfaceIDList_EndoRV.end() );
	surfaceIDList.insert( surfaceIDList.end(), surfaceIDList_FreeWall.begin(), surfaceIDList_FreeWall.end() );

	surfaceNodes.reserve( surfaceNodes_EndoRV.size() + surfaceNodes_FreeWall.size() ); // preallocate memory
	surfaceNodes.insert( surfaceNodes.end(), surfaceNodes_EndoRV.begin(), surfaceNodes_EndoRV.end() );
	surfaceNodes.insert( surfaceNodes.end(), surfaceNodes_FreeWall.begin(), surfaceNodes_FreeWall.end() );


	int numOfSteps = resultMatList.size();
	dbVector volumeVec(numOfSteps, 0);


	for (int i = 0; i < numOfSteps; i++) {

		dbMatrix& rMatrix = resultMatList[i];

		vector<ParticleExt> ptcls = meshData->getNodesVec();

		// Update the coordinates
		for (int j = 0; j < ptcls.size(); j++) {
			for (int k = 0; k < 3; k++) {
				ptcls[j].getCoord(k) += rMatrix[j][k];
			}
		}

		// Compute the volume of the cavity
		double volume = 0;

		int x = 0;
		int y = 1;
		int z = 2;

		int numOfsurfaces = surfaceIDList.size();
		for (int m = 0; m < numOfsurfaces; m++) {

			dbVector& v1 = ptcls[surfaceNodes[m][0] - 1].getCoords();
			dbVector& v2 = ptcls[surfaceNodes[m][1] - 1].getCoords();
			dbVector& v3 = ptcls[surfaceNodes[m][2] - 1].getCoords();

			volume += ((v2[y] - v1[y]) * (v3[z] - v1[z])
					- (v2[z] - v1[z]) * (v3[y] - v1[y]))
					* (v1[x] + v2[x] + v3[x]);
		}

		volumeVec[i] = abs(volume / 6.0);

//		logFile << "[" << i << "] Geometry volume = " << volumeVec[i] << endl;
	}

	string timeGraphFileName = folderName + "time_volume-RV.grf";
	logFile << "Writing to: " << timeGraphFileName << endl;
	cout << "Writing to: " << timeGraphFileName << endl;

	ofstream writeTimeGraph(timeGraphFileName);
	writeTimeGraph.precision(15);
	writeTimeGraph.setf(std::ios::scientific, std::ios::floatfield);


	string pressureGraphFileName = folderName + "pressure_volume-RV.grf";
	logFile << "Writing to: " << pressureGraphFileName << endl;
	cout << "Writing to: " << pressureGraphFileName << endl;

	ofstream writePressureGraph(pressureGraphFileName);
	writePressureGraph.precision(15);
	writePressureGraph.setf( std::ios::scientific, std:: ios::floatfield );


	if (writeTimeGraph.good() && writePressureGraph.good()) {
		for (int i = 0; i < step_value_vec.size(); i++) {
			writeTimeGraph << step_value_vec[i] << " " << volumeVec[i] << endl;
			writePressureGraph << volumeVec[i] << " " << rightCavityPressures[i] << endl;
		}
	}

	writeTimeGraph.close();
	writePressureGraph.close();

	rightCavityVolumes = volumeVec;

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::calcRightCavityVolumes_step(InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	logFile << "Calculating right cavity volumes" << endl;

	dbMatrix& resultMatrix = this->getResult("displacement");
//	printMatrix(resultMatrix,"displacement Matrix",logFile);
//	MPI_Abort(MPI_COMM_WORLD, 1);

	// Convert each result at each step from a vector to a matrix
	vector<dbMatrix> resultMatList(resultMatrix[0].size(),
			dbMatrix(resultMatrix.size() / 3, dbVector(3,0)));
	for (int i = 0; i < resultMatrix[0].size(); i++) {
		dbMatrix& mat = resultMatList[i];

		int row = 0;
		int col = 0;
		for (int j = 0; j < resultMatrix.size(); j++) {

			mat[row][col] = resultMatrix[j][i];
			col++;

			if (col > 2) {
				col = 0;
				row++;
			}
		}
	}

	// Extract the surface of the left ventricle
	intVector surfaceIDList_EndoRV;
	intMatrix surfaceNodes_EndoRV;
	meshData->getSpecificSurfaceNodes("Endocaridum_RV_Surface",
			surfaceIDList_EndoRV,surfaceNodes_EndoRV, InputData, logFile);

	intVector surfaceIDList_FreeWall;
	intMatrix surfaceNodes_FreeWall;
	meshData->getSpecificSurfaceNodes("Endocaridum_Freewall_Surface",
			surfaceIDList_FreeWall,surfaceNodes_FreeWall, InputData, logFile);

	// Merge Endocaridum_RV_Surface and Endocaridum_Freewall_Surface data
	intVector surfaceIDList;
	intMatrix surfaceNodes;

	surfaceIDList.reserve( surfaceIDList_EndoRV.size() + surfaceIDList_FreeWall.size() );
	surfaceIDList.insert( surfaceIDList.end(), surfaceIDList_EndoRV.begin(), surfaceIDList_EndoRV.end() );
	surfaceIDList.insert( surfaceIDList.end(), surfaceIDList_FreeWall.begin(), surfaceIDList_FreeWall.end() );

	surfaceNodes.reserve( surfaceNodes_EndoRV.size() + surfaceNodes_FreeWall.size() ); // preallocate memory
	surfaceNodes.insert( surfaceNodes.end(), surfaceNodes_EndoRV.begin(), surfaceNodes_EndoRV.end() );
	surfaceNodes.insert( surfaceNodes.end(), surfaceNodes_FreeWall.begin(), surfaceNodes_FreeWall.end() );


	int numOfSteps = resultMatList.size();
	dbVector volumeVec(numOfSteps, 0);


	for (int i = 0; i < numOfSteps; i++) {

		dbMatrix& rMatrix = resultMatList[i];

		vector<ParticleExt> ptcls = meshData->getNodesVec();

		// Update the coordinates
		for (int j = 0; j < ptcls.size(); j++) {
			for (int k = 0; k < 3; k++) {
				ptcls[j].getCoord(k) += rMatrix[j][k];
			}
		}

		// Compute the volume of the cavity
		double volume = 0;

		int x = 0;
		int y = 1;
		int z = 2;

		int numOfsurfaces = surfaceIDList.size();
		for (int m = 0; m < numOfsurfaces; m++) {

			dbVector& v1 = ptcls[surfaceNodes[m][0] - 1].getCoords();
			dbVector& v2 = ptcls[surfaceNodes[m][1] - 1].getCoords();
			dbVector& v3 = ptcls[surfaceNodes[m][2] - 1].getCoords();

			volume += ((v2[y] - v1[y]) * (v3[z] - v1[z])
					- (v2[z] - v1[z]) * (v3[y] - v1[y]))
					* (v1[x] + v2[x] + v3[x]);
		}

		volumeVec[i] = abs(volume / 6.0);

//		logFile << "[" << i << "] Geometry volume = " << volumeVec[i] << endl;
	}

	string timeGraphFileName = folderName + "time_volume-RV_step.grf";
	cout << "Writing to: " << timeGraphFileName << endl;

	ofstream writeTimeGraph(timeGraphFileName);
	writeTimeGraph.precision(15);
	writeTimeGraph.setf(std::ios::scientific, std::ios::floatfield);


	string pressureGraphFileName = folderName + "pressure_volume-RV_step.grf";
	cout << "Writing to: " << pressureGraphFileName << endl;

	ofstream writePressureGraph(pressureGraphFileName);
	writePressureGraph.precision(15);
	writePressureGraph.setf( std::ios::scientific, std:: ios::floatfield );


	if (writeTimeGraph.good() && writePressureGraph.good()) {
		for (int i = 0; i < step_value_vec.size(); i++) {
			writeTimeGraph << step_value_vec[i] << " " << volumeVec[i] << endl;
			writePressureGraph << volumeVec[i] << " " << rightCavityPressures[i] << endl;
		}
	}

	writeTimeGraph.close();
	writePressureGraph.close();

	rightCavityVolumes = volumeVec;

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::setGraph(const char* graphName,dbVector graph) {

	string str(graphName);

	map<string, dbVector>::iterator it = graphList.find(str);

	if (it == graphList.end())
		graphList[str] = graph ;
	else{
		cout << "ERROR: In Data::setGraph, '" << graphName
				<< "' already exists in graphList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
dbVector& Data::getGraph(const char* graphName) {

	string str(graphName);

	map<string, dbVector>::iterator it = graphList.find(str);

	if (it != graphList.end())
		return it->second;
	else {
		cout << "ERROR: In Data::getGraph, '" << graphName
				<< "' does not exist in graphList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::deleteGraph(const char* graphName){

	string str(graphName);

	map<string, dbVector>::iterator it = graphList.find(str);

	if (it != graphList.end())
		graphList.erase(it);
	else {
		cout << "ERROR: In Data::deleteGraph, '" << graphName
				<< "' does not exist in graphList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::syncCardiacTimeStepsAndResults(InputFileData* InputData, ofstream& logFile){

	logFile << "Synchronising result timesteps and graph results." << endl;
	cout    << "Synchronising result timesteps and graph results." << endl;

#ifdef _DataDebugMode_
	printVector(leftCavityVolumes,"leftCavityVolumes(Before)",logFile);
	printVector(leftCavityPressures,"leftCavityPressures(Before)",logFile);

	printVector(rightCavityVolumes,"rightCavityVolumes(Before)",logFile);
	printVector(rightCavityPressures,"rightCavityPressures(Before)",logFile);
#endif


	dbVector& LVTimeSteps = getGraph("LVTimeSteps");
	intVector selectedLVIndex;

	logFile << "LVTimeSteps.size() = " << LVTimeSteps.size() << endl;
	logFile << "step_value_vec.size() = " << step_value_vec.size() << endl;

	if(LVTimeSteps.size()>1){
		for(int i=0; i<step_value_vec.size();i++){
			for(int j=0; j<LVTimeSteps.size();j++){
				if(step_value_vec[i] == LVTimeSteps[j]){
					selectedLVIndex.push_back(j);
				}
			}
		}

		if(selectedLVIndex.size() != step_value_vec.size()){
			logFile << "In Data::syncCardiacTimeStepsAndResults, not all the graph "
					"timesteps were found in the step_value_vec vector." << endl;
			cout << "In Data::syncCardiacTimeStepsAndResults, not all the graph "
					"timesteps were found in the step_value_vec vector." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		dbVector tempVolVec(selectedLVIndex.size(),0), tempPreVec(selectedLVIndex.size(),0);
		for(int i=0;i<selectedLVIndex.size();i++){
			int index = selectedLVIndex[i];

			tempVolVec[i] = leftCavityVolumes[index];
			tempPreVec[i] = leftCavityPressures[index];
		}

		leftCavityVolumes = tempVolVec;
		leftCavityPressures = tempPreVec;

		deleteGraph("LVTimeSteps");

	}
	else{
		logFile << "In Data::syncCardiacTimeStepsAndResults, LV timesteps is empty." << endl;
	}





//	dbVector& RVTimeSteps = getGraph("RVTimeSteps");
//	intVector selectedRVIndex;
//
//	logFile << "RVTimeSteps.size() = " << RVTimeSteps.size() << endl;
//	logFile << "step_value_vec.size() = " << step_value_vec.size() << endl;
//
//	if(RVTimeSteps.size()>1){
//		for(int i=0; i<step_value_vec.size();i++){
//			for(int j=0; j<RVTimeSteps.size();j++){
//				if(step_value_vec[i] == RVTimeSteps[j]){
//					selectedRVIndex.push_back(j);
//				}
//			}
//		}
//
//		if(selectedRVIndex.size() != step_value_vec.size()){
//			logFile << "In Data::syncCardiacTimeStepsAndResults, not all the graph "
//					"timesteps were found in the step_value_vec vector." << endl;
//			cout << "In Data::syncCardiacTimeStepsAndResults, not all the graph "
//					"timesteps were found in the step_value_vec vector." << endl;
//			MPI_Abort(MPI_COMM_WORLD, 1);
//		}
//
//		dbVector tempVolVec(selectedRVIndex.size(),0), tempPreVec(selectedRVIndex.size(),0);
//		for(int i=0;i<selectedRVIndex.size();i++){
//			int index = selectedRVIndex[i];
//
//			tempVolVec[i] = rightCavityVolumes[index];
//			tempPreVec[i] = rightCavityPressures[index];
//		}
//
//		rightCavityVolumes = tempVolVec;
//		rightCavityPressures = tempPreVec;
//
//		deleteGraph("RVTimeSteps");
//
//	}
//	else{
//		logFile << "In Data::syncCardiacTimeStepsAndResults, RV timesteps is empty." << endl;
//	}

#ifdef _DataDebugMode_
	printVector(leftCavityVolumes,"leftCavityVolumes(after)",logFile);
	printVector(leftCavityPressures,"leftCavityPressures(after)",logFile);
#endif

//	printVector(rightCavityVolumes,"rightCavityVolumes(after)",logFile);
//	printVector(rightCavityPressures,"rightCavityPressures(after)",logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::insertZeroResultFields(InputFileData* InputData, ofstream& logFile){

	for(int i=0;i<resultNameList.size();i++){
		dbMatrix& result = this->getResult(resultNameList[i].c_str());

		for(int j=0; j<result.size(); j++){
			result[j].insert(result[j].begin(),0);
		}
	}

	step_value_vec.insert(step_value_vec.begin(),0);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::removeZeroResultFields(InputFileData* InputData, ofstream& logFile){

	for(int i=0;i<resultNameList.size();i++){
		dbMatrix& result = this->getResult(resultNameList[i].c_str());

		for(int j=0; j<result.size(); j++){
			result[j].insert(result[j].begin(),0);
		}
	}

	step_value_vec.insert(step_value_vec.begin(),0);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Data::plotPostProcessGraph(InputFileData* InputData, ofstream& logFile){

	// Read influence factors from input file
	string graphX_start = "graphX";
	string graphY_start = "graphY";
	string graphNode_start = "graphNode";
	string graphDOF_start = "graphDOF";

	int nGraphs = InputData->getValue("nGraphs");

	for(int i=0; i < nGraphs ; i ++){

		ostringstream convert;   // stream used for the conversion
		convert << i;

		string graphX_name = graphX_start + convert.str();
		string graphY_name = graphY_start + convert.str();
		string graphNode_name = graphNode_start + convert.str();
		string graphDOF_name = graphDOF_start + convert.str();

		int graphX = InputData->getValue(graphX_name.c_str());
		int graphY = InputData->getValue(graphY_name.c_str());
		int graphNode =InputData->getValue(graphNode_name.c_str());
		int graphDOF = InputData->getValue(graphDOF_name.c_str());

		string resultNameX, resultNameY;
		dbVector graphX_data = getGraphData(graphX,graphNode,graphDOF,resultNameX,InputData,logFile);
		dbVector graphY_data = getGraphData(graphY,graphNode,graphDOF,resultNameY,InputData,logFile);

		ostringstream nodeConvert;   // stream used for the conversion
		nodeConvert << graphNode;

		ostringstream DOFConvert;   // stream used for the conversion
		DOFConvert << graphDOF;

		string graphFileName = resultNameX + "_" + resultNameY + "-node" + nodeConvert.str() + "-DOF" + DOFConvert.str();
		saveGraphResultsToFile_grf_format(graphFileName, graphX_data,
				graphY_data, logFile);

	}

}

/*!****************************************************************************/
/*!****************************************************************************/
dbVector& Data::getGraphData(int graphType, int node, int DOF, string& resultName,
		InputFileData* InputData, ofstream& logFile){

	if (graphType == -2){
		resultName = "load";
		return this->getGraph("surfaceLoad");
	}
	else if (graphType == -1){
		resultName = "time";
		return step_value_vec;
	}
	else if (graphType > -1 && graphType < allResultsNameList.size()){
		resultName = allResultsNameList[graphType][0];

		int numDofs = this->getResultDOF(resultName.c_str());
		int DOFIndex = ((node - 1)*numDofs) + (DOF - 1);

		return this->getResultAcrossSteps(resultName.c_str(),DOFIndex);

	}
	else {
		logFile << "In Data::plotPostProcessGraph, graphType = " << graphType
				<< " and is not supported." << endl;
		cout << "In Data::plotPostProcessGraph, graphType = " << graphType
				<< " and is not supported." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}
