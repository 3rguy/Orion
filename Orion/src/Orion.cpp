//============================================================================
// Name        : Orion.cpp
// Author      : Ritesh Rao Rama
// Version     :
// Copyright   : 
// Description : POD-based real-time calculation
//============================================================================
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include <time.h>
#include <chrono>
#include <ctime>

//#include "clapack.h"
#define PI 3.14159265

//#include "commonTypedefs.h"
//#include "commonFunctions.h"
//#include "fortranFunctions.h"

#include <limits>

//#include "Database.h"
#include "Data.h"
//#include "PODCalc.h"
//#include "FEMGeometryX.h"
//#include "ParticleX.h"
//#include "FEMElementX.h"

#include "defs.h"
#include "ROMCalc.h"
#include "Interpolation.h"

using namespace std;

// =============================================================================
// Function declaration
void testFunctions(InputFileData* InputData,ofstream& logFile);
void findPointInPolygonTest(InputFileData* InputData, ofstream& logFile);
void PODCalcTest_readFile(dbMatrix& fullMatrix, InputFileData* InputData,
		ofstream& logFile);
void PODCalcTest(InputFileData* InputData,ofstream& logFile);
void readResultFileTest(InputFileData* InputData,ofstream& logFile);
void readResultFileTest_writeResultToFile(string& resultName,dbVector& resultStepVector,
		dbMatrix& resultMatrix,int& numDofs,InputFileData* InputData,ofstream& logFile);
// =============================================================================


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

int main(int argc, char* argv[]) {

	int choice = 0;

	PetscInitialize(&argc, &argv, (char*) 0, 0);

	// *************************************************************************
	// Log writer
	string fileOutput = "oLog.dat";
	ofstream logFile(fileOutput.c_str(), std::ofstream::out);
	logFile.precision(15);

	time_t tim;  //create variable of time_t
	time(&tim); //pass variable tim to time function

	cout << "================================" << endl;
	cout << "ORION: " << ctime(&tim);
	cout << "--------------------------------" << endl;

	logFile << "================================" << endl;
	logFile << "ORION: " << ctime(&tim);
	logFile << "--------------------------------" << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	//*************************************************************************
	// Read input.dat file
	string filename = "input.dat";
	InputFileData* InputData = new InputFileData(filename,logFile);

	switch (choice) {
	case 0: {
		//*************************************************************************
		// Reduced Order Method Calculation Calculation
		ROMCalc* reducedOrderCalculation = new ROMCalc(InputData, logFile);
		delete reducedOrderCalculation;
		break;
	}

	case 1: {
		//*************************************************************************
		// Testing specific algorithm in Orion
		testFunctions(InputData, logFile);
		break;
	}

	default:
		cout << "ERROR: In main(), choice is not valid" << endl;
		logFile << "ERROR: In main(), choice is not valid" << endl;
	}

	delete InputData;

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Total calculation time: " << elapsed_seconds.count() << " sec"
			<< endl;
	logFile << "Total calculation time: " << elapsed_seconds.count() << " sec"
			<< endl;

	cout << "--------------------------------" << endl;
	logFile << "--------------------------------" << endl;

	return 0;

	//	cout<<"*** THE END ***"<<endl;
	//	MPI_Abort(MPI_COMM_WORLD, 1);
}


void testFunctions(InputFileData* InputData, ofstream& logFile) {

	int choice = 1;

	switch (choice) {
	case 0: {
		//*************************************************************************
		// Point-In-Polygon test
		findPointInPolygonTest(InputData, logFile);
		break;

	}

	case 1: {
		//*************************************************************************
		// snapshot method test
		PODCalcTest(InputData, logFile);
		break;
	}

	case 2: {
		//*************************************************************************
		// extracting result from FEM.res file
		readResultFileTest(InputData, logFile);
		break;
	}
	default:
		cout << "ERROR: In main()-testFunctions, choice is not valid" << endl;
		logFile << "ERROR: In main()-testFunctions, choice is not valid"
				<< endl;

	}

}


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

	string meshFileName = "216/";
	myData.setFileName(meshFileName);

	myData.readMeshDataFile(InputData, logFile);

	myData.setFileName(meshFileName);
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
	FEMData->writeMeshFile(InputData, logFile);
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

void readResultFileTest(InputFileData* InputData,ofstream& logFile){

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

