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

	allResultsNameList[3][0] = "fibre strain E_11";
	allResultsNameList[3][1] = "Scalar";

	allResultsNameList[4][0] = "fibre strain E_22";
	allResultsNameList[4][1] = "Scalar";

	allResultsNameList[5][0] = "fibre strain E_33";
	allResultsNameList[5][1] = "Scalar";

	allResultsNameList[6][0] = "stress";
	allResultsNameList[6][1] = "Matrix";

	allResultsNameList[7][0] = "effective stress";
	allResultsNameList[7][1] = "Scalar";

	allResultsNameList[8][0] = "fibre stress 1";
	allResultsNameList[8][1] = "Scalar";

	allResultsNameList[9][0] = "fibre stress 2";
	allResultsNameList[9][1] = "Scalar";

	allResultsNameList[10][0] = "fibre stress 3";
	allResultsNameList[10][1] = "Scalar";

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

	allResultsNameList[3][0] = "fibre strain E_11";
	allResultsNameList[3][1] = "Scalar";

	allResultsNameList[4][0] = "fibre strain E_22";
	allResultsNameList[4][1] = "Scalar";

	allResultsNameList[5][0] = "fibre strain E_33";
	allResultsNameList[5][1] = "Scalar";

	allResultsNameList[6][0] = "stress";
	allResultsNameList[6][1] = "Matrix";

	allResultsNameList[7][0] = "effective stress";
	allResultsNameList[7][1] = "Scalar";

	allResultsNameList[8][0] = "fibre stress 1";
	allResultsNameList[8][1] = "Scalar";

	allResultsNameList[9][0] = "fibre stress 2";
	allResultsNameList[9][1] = "Scalar";

	allResultsNameList[10][0] = "fibre stress 3";
	allResultsNameList[10][1] = "Scalar";
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
void Data::readResultFile_resFormat(string& inputFileName, ofstream& logFile) {

	for (int i = 0; i < allResultsNameList.size(); i++) {

		dbMatrix resultMatrix;
		dbVector resultStepVector;
		int numDofs;

		cout << "Reading result: " << allResultsNameList[i][0] << endl;
		logFile << "Reading result: " << allResultsNameList[i][0] << endl;

		bool isResultFound = readResultFile_resFormat_specificResult(
				inputFileName, allResultsNameList[i][0], resultMatrix, numDofs,
				resultStepVector, logFile);

		if (isResultFound == true) {
			resultNameList.push_back(allResultsNameList[i][0]);
			setResult(allResultsNameList[i][0].c_str(), resultMatrix);
			setResultDOF(allResultsNameList[i][0].c_str(), numDofs);

			// Record the first series of step values
			if (step_value_vec.size() == 0) {
				step_value_vec = resultStepVector;
			}
			// Else check if each step value are the same wrt to previously
			// recorded ones
			else {

				bool isStepValueVecSimilar = true;
				if (step_value_vec.size() == resultStepVector.size()) {
					for (int i = 0; i < step_value_vec.size(); i++) {
						if (step_value_vec[i] != resultStepVector[i]) {
							isStepValueVecSimilar = false;
							break;
						}
					}
				} else {
					isStepValueVecSimilar = false;
				}

				if (isStepValueVecSimilar == false) {
					cout << "ERROR: The step values of '"
							<< allResultsNameList[i][0] << "' in '"
							<< inputFileName
							<< "' is not the same compared to the other"
									" result types" << endl;

					logFile << "ERROR: The step values of '"
							<< allResultsNameList[i][0] << "' in '"
							<< inputFileName
							<< "' are not the same compared to the other"
									" result types" << endl;

					MPI_Abort(MPI_COMM_WORLD, 1);
				}

			}

#ifdef _DataDebugMode_
		printMatrix(resultMatrix, allResultsNameList[i][0].c_str(), logFile);
#endif

		}

	}

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
void Data::saveResultsToFile(ofstream& logFile) {

	int choice = 1;

	switch (choice) {
	case 1:
		//Save matrix to .res file format
		saveAllResultsToFile_res_format("fem_orion.res",logFile);
		break;

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
		ofstream& logFile) {

	saveResultsToFile_res_format(outputFileName,resultNameList,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
// Writing generated displacement matrix to files
void Data::saveResultsToFile_res_format(const char* outputFileName,
		vector<string>& saveResultNameList, ofstream& logFile) {

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

		for (int j = 0; j < step_value_vec.size(); j++) {
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

//#ifdef _DataDebugMode_
	logFile << "Before assigning '" << resultName << "' to particles" << endl;
	printMatrix(resultMatrix, "", logFile);
//#endif

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

//#ifdef _DataDebugMode_

		logFile << "******** Particles stepDOFs ********" << endl;
		for(int l=0; l<particleList.size(); l++) {
			std::string msg = "For particle[";
			ostringstream convert;
			convert << "For particle[" << l << "]";
			msg = convert.str();
			printMatrix(particleList[l].getStepDOFMat(),msg.c_str(),logFile);
			logFile << endl;
		}

//#endif

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
				<< "' already exists in resultList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

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
				<< "' already exists in resultDofList" << endl;
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

