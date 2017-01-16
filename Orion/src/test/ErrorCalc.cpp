/*
 * ErrorCalc.cpp
 *
 *  Created on: Oct 21, 2014
 *      Author: rama
 */

#include <ErrorCalc.h>

//ErrorCalc::ErrorCalc(InputFileData* InputData, ofstream& logFile) {
//
//	using namespace std;
//
//	cout << "ERROR CALCULATION MODE" << endl;
//	logFile << "ERROR CALCULATION MODE" << endl;
//
//	string errorFileName = "oError.dat";
//	vector<string> fileList;
//
//	ifstream inputFile(errorFileName.c_str());
//
//	if (inputFile.eof()) {
//		logFile << "Input file: " << errorFileName << " contains no data!"
//				<< endl;
//		MPI_Abort(MPI_COMM_WORLD, 1);
//	}
//
//	if (inputFile) {
//
//		logFile << "Reading: " << errorFileName << endl;
//		cout << "Reading: " << errorFileName << endl;
//
//		string name;
//
//		do {
//
//			fileList.push_back(name);
//
//			logFile << "Error file read: " << name << endl;
//
//			inputFile >> name;
//
//		} while (inputFile.good());
//
//	}
//
//	inputFile.close();
//
//	cout << "Reading: " << fileList[0] << endl;
//	logFile << "Reading: " << fileList[0] << endl;
//	Data* exactData = new Data();
//	exactData->readResultFile_resFormat(fileList[0], logFile);
//
//	for (int i = 1; i < fileList.size(); i++) {
//		cout << "Reading: " << fileList[i] << endl;
//		logFile << "Reading: " << fileList[i] << endl;
//		Data* approxData = new Data();
//		approxData->readResultFile_resFormat(fileList[i], logFile);
//
////		compareTimeSteps(exactData, approxData, InputData, logFile);
//		vector<string> commonNameList = findCommonDataResults(exactData,
//				approxData, InputData, logFile);
//
//		calculateErrors(exactData, approxData, commonNameList, InputData,
//				logFile);
//
//		delete approxData;
//	}
//
//}

ErrorCalc::ErrorCalc(InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	cout << "ERROR CALCULATION MODE" << endl;
	logFile << "ERROR CALCULATION MODE" << endl;

	DataContainer* problemData = new DataContainer();

	readingErrorFile(problemData, InputData, logFile);

	resErrorCalc(problemData, InputData, logFile);

	grfErrorCalc(problemData, InputData, logFile);

}

// *****************************************************************************
// *****************************************************************************
void ErrorCalc::readingErrorFile(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	using namespace std;

	string errorFileName = "oError.dat";
	vector<string> fileList;

	ifstream inputFile(errorFileName.c_str());

	if (inputFile.eof()) {
		logFile << "Input file: " << errorFileName << " contains no data!"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else if(inputFile.good() == 0){
		logFile << "Input file: " << errorFileName << " is missing!"
						<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (inputFile) {

		logFile << "Reading: " << errorFileName << endl;
		cout << "Reading: " << errorFileName << endl;

		string name;

		inputFile >> name;

		while (inputFile.good()) {

			fileList.push_back(name);

			logFile << "Error file read: " << name << endl;

			inputFile >> name;

		}

	}

	inputFile.close();

	logFile << "In fileList:" << endl;
	for(int i = 0; i < fileList.size(); i++){
		logFile << fileList[i] << endl;
	}

	// -------------------------------------------------------------------------
	// Sorting errorFiles

	vector<string> resFileNameList;
	vector<string> grfFileNameList;

	for(int i = 0; i < fileList.size(); i++){
		if(fileList[i] == "RES_FILES"){
			for(;;){
				i++;
				if(fileList[i] == "END")
					break;
				else
					resFileNameList.push_back(fileList[i]);

			}
		}
		else if(fileList[i] == "GRF_FILES"){
			for(;;){
				i++;
				if(fileList[i] == "END")
					break;
				else
					grfFileNameList.push_back(fileList[i]);

			}
		}
	}

	problemData->setValue("resFileNameList",resFileNameList);
	problemData->setValue("grfFileNameList",grfFileNameList);

}

// *****************************************************************************
// *****************************************************************************
void ErrorCalc::resErrorCalc(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	cout << endl << "########## RES Files error calculation ##########" << endl;
	logFile << endl << "########## RES Files error calculation ##########" << endl;

	vector<string>& fileList = problemData->getStringVec("resFileNameList");

	if (fileList.size() == 0) {
		cout << "No res file found" << endl;
		logFile << "No res file found" << endl;

		return;
	}

	cout << "Reading: " << fileList[0] << endl;
	logFile << "Reading: " << fileList[0] << endl;
	Data* exactData = new Data();
	exactData->readResultFile_resFormat(fileList[0], logFile);

	for (int i = 1; i < fileList.size(); i++) {

		cout << "*************************************************" << endl;
		logFile << "*************************************************" << endl;

		cout << "Reading: " << fileList[i] << endl;
		logFile << "Reading: " << fileList[i] << endl;
		Data* approxData = new Data();
		approxData->readResultFile_resFormat(fileList[i], logFile);

		//		compareTimeSteps(exactData, approxData, InputData, logFile);
		vector<string> commonNameList = findCommonDataResults(exactData,
				approxData, InputData, logFile);

		calculateErrors_res_Format(exactData, approxData, commonNameList,
				InputData, logFile);

		delete approxData;
	}

}

// *****************************************************************************
// *****************************************************************************
void ErrorCalc::grfErrorCalc(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	cout << endl << "########## GRF Files error calculation ##########" << endl;
	logFile << endl << "########## GRF Files error calculation ##########" << endl;

	vector<string>& fileList = problemData->getStringVec("grfFileNameList");

	if(fileList.size() == 0){
		cout << "No grf file found" << endl;
		logFile << "No grf file found" << endl;

		return;
	}

	if(fileList.size() % 2 != 0){
		cout << "ERROR:Number of GRF files should be even." << endl;
		logFile << "ERROR:Number of GRF files should be even." << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	int nGRFPairs = fileList.size() / 2 ;

	Data* grfReadData = new Data();

	for (int i = 0; i < nGRFPairs; i++) {

		cout << "*************************************************" << endl;
		logFile << "*************************************************" << endl;

		dbMatrix exactGrfMatrix, approxGrfMatrix;

		cout << "Reading: " << fileList[(i*2)] << endl;
		logFile << "Reading: " << fileList[(i*2)] << endl;
		grfReadData->readGraphFile_grfFormat(fileList[(i*2)], exactGrfMatrix, logFile);

		cout << "Reading: " << fileList[(i*2)+1] << endl;
		logFile << "Reading: " << fileList[(i*2)+1] << endl;
		grfReadData->readGraphFile_grfFormat(fileList[(i*2)+1], approxGrfMatrix, logFile);

		calculateErrors_grf_Format(exactGrfMatrix, approxGrfMatrix,
						InputData, logFile);

	}

}


// *****************************************************************************
// *****************************************************************************
void ErrorCalc::compareTimeSteps(Data* exactData, Data* approxData,
		InputFileData* InputData, ofstream& logFile){

	dbVector& exactStepVec = exactData->getStepValueVec();
	dbVector& approxStepVec = approxData->getStepValueVec();

	if(exactStepVec.size() == approxStepVec.size()){
		for(int i=0; i < exactStepVec.size(); i++){
			if(fabs(exactStepVec[i] - approxStepVec[i]) > 1e-8){
				logFile <<"In ErrorCalc::compareTimeSteps, exacStepVec("
						<< exactStepVec[i] << " != approxData ("
						<< approxStepVec[i] <<")" << endl;
				cout <<"In ErrorCalc::compareTimeSteps, exacStepVec("
						<< exactStepVec[i] << " != approxData ("
						<< approxStepVec[i] <<")" << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
	}
	else{
		logFile << "In ErrorCalc::compareTimeSteps, size of exacStepVec("
				<< exactStepVec.size() << " != approxData ("
				<< approxStepVec.size() << ")" << endl;
		cout << "In ErrorCalc::compareTimeSteps, size of exacStepVec("
				<< exactStepVec.size() << " != approxData ("
				<< approxStepVec.size() << ")" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

// *****************************************************************************
// *****************************************************************************
vector<string> ErrorCalc::findCommonDataResults(Data* exactData,
		Data* approxData, InputFileData* InputData, ofstream& logFile){

	vector<string>& exactNameList = exactData->getResultNameList();
	vector<string>& approxNameList = approxData->getResultNameList();

	vector<string> commonNameList;

	for(int i=0; i<exactNameList.size(); i++){
		for(int j=0; j<approxNameList.size(); j++){
			if(exactNameList[i] == approxNameList[j]){
				commonNameList.push_back(exactNameList[i]);
			}
		}
	}

	logFile << "commonNameList:" << endl;
	for(int i=0; i < commonNameList.size(); i++ ){
		logFile << "(" << i << ") " << commonNameList[i] << endl;
	}

	return commonNameList;
}

// *****************************************************************************
// *****************************************************************************
void ErrorCalc::calculateErrors_res_Format(Data* exactData, Data* approxData,
		vector<string>& commonNameList, InputFileData* InputData,
		ofstream& logFile){

	for(int i=0; i < commonNameList.size(); i++){

		// *********************************************************************
		double rel2Norm =
				relativeL2norm(exactData->getResult(commonNameList[i].c_str()),
				approxData->getResult(commonNameList[i].c_str()),
				InputData,logFile);

		cout << "Relative L2 Norm of " << commonNameList[i] << " : " << rel2Norm << endl;
		logFile << "Relative L2 Norm of " << commonNameList[i] << " : " << rel2Norm << endl;

	}
}

// *****************************************************************************
// *****************************************************************************
void ErrorCalc::calculateErrors_grf_Format(dbMatrix& exactGrfMatrix,
		dbMatrix& approxGrfMatrix,InputFileData* InputData, ofstream& logFile){

	dbMatrix exactGrfMatrix_volOnly(1,dbVector());
	exactGrfMatrix_volOnly[0] = exactGrfMatrix[0];

	dbMatrix approxGrfMatrix_volOnly(1,dbVector());
	approxGrfMatrix_volOnly[0] = approxGrfMatrix[0];

#ifdef _errorCalcMode_
	printMatrix(exactGrfMatrix_volOnly,"exactGrfMatrix_volOnly",logFile);
	printMatrix(approxGrfMatrix_volOnly,"approxGrfMatrix_volOnly",logFile);
#endif

	// *************************************************************************
	double rel2Norm_first =
		relativeL2norm(exactGrfMatrix_volOnly,approxGrfMatrix_volOnly,
				InputData,logFile);

	cout << "First graph column" << endl;
	cout << "Relative L2 Norm: " << rel2Norm_first << endl;
	logFile << "First graph column" << endl;
	logFile << "Relative L2 Norm: " << rel2Norm_first << endl;
	// *************************************************************************


	exactGrfMatrix_volOnly[0] = exactGrfMatrix[1];
	approxGrfMatrix_volOnly[0] = approxGrfMatrix[1];

#ifdef _errorCalcMode_
	printMatrix(exactGrfMatrix_volOnly, "exactGrfMatrix_volOnly", logFile);
	printMatrix(approxGrfMatrix_volOnly, "approxGrfMatrix_volOnly", logFile);
#endif

	// *********************************************************************
	double rel2Norm_second = relativeL2norm(exactGrfMatrix_volOnly,
			approxGrfMatrix_volOnly, InputData, logFile);

	cout << "Second graph column" << endl;
	cout << "Relative L2 Norm: " << rel2Norm_second << endl;
	logFile << "Second graph column" << endl;
	logFile << "Relative L2 Norm: " << rel2Norm_second << endl;
	// *************************************************************************

	// Combined relative L2Norm
	double rel2Norm = sqrt(pow(rel2Norm_first,2)+pow(rel2Norm_second,2));
	cout << "Combined graph column" << endl;
	cout << "Relative L2 Norm: " << rel2Norm << endl;
	logFile << "Second graph column" << endl;
	logFile << "Relative L2 Norm: " << rel2Norm << endl;

}

// *****************************************************************************
// *****************************************************************************
double ErrorCalc::relativeL2norm(dbMatrix& exactMat, dbMatrix& approxMat,
		InputFileData* InputData, ofstream& logFile){

	double sum_sqr_diff = 0;
	double sum_sqr_exact = 0;

//	cout << "Size of exactMat: " << exactMat.size() << " x " << exactMat[0].size() << endl;
//	cout << "Size of approxMat: " << approxMat.size() << " x " << approxMat[0].size() << endl;

	for(int i=0;i<exactMat.size();i++){
		for(int j=0; j<exactMat[i].size();j++){
			double diff = approxMat[i][j]-exactMat[i][j];
			sum_sqr_diff += pow(diff,2);
			sum_sqr_exact += pow(exactMat[i][j],2);
		}
	}

	double relL2Norm = sqrt(sum_sqr_diff)/sqrt(sum_sqr_exact);

	return relL2Norm;

}
