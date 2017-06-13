/*
 * MicroRVE.cpp
 *
 *  Created on: 14 Mar 2017
 *      Author: rama
 */

#include <MicroRVE.h>

MicroRVE::MicroRVE(Database& myDB, DataContainer* problemData,
		 InputFileData* InputData, ofstream& logFile) {

	cout << "Model Calculation: MicroRVE" << endl;
	logFile << "Model Calculation: MicroRVE" << endl;

	myDatabase = myDB;

}


/*!****************************************************************************/
/*!****************************************************************************/
void MicroRVE::preProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	initPreProcessing(myDatabase,problemData,InputData,logFile);

	defaultPreProcessing(myDatabase,problemData,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void MicroRVE::defaultPreProcessing(Database& myDatabase,
	   DataContainer* problemData, InputFileData* InputData, ofstream& logFile){
	// Search database and extract the neighbouring nodes

	searchDatabase(problemData, myDatabase, InputData, logFile);

	// Carry out a preliminary check of the quality of MLS interpolants
	// based on the selected datasets
	preliminaryMLSInterpolantsCalc(myDatabase, problemData, InputData, logFile);

	dbVector& myParameters = problemData->getDbVector("myParameters");
	myData.setParamValuesVec(myParameters);

	loadSelectedData(problemData, myDatabase, InputData, logFile);

	assembleResultsOnly(problemData, InputData, logFile);

	if(InputData->getValue("savePreProcessedData") == 1){
		savePreProcessedData(problemData,InputData,logFile);
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
void MicroRVE::loadSelectedData(DataContainer* problemData,
		Database& myDatabase, InputFileData* InputData, ofstream& logFile) {

	cout << "Cardiac::loadSelectedData" << endl;

	int isRV = InputData->getValue("isRightVentriclePresent");

	intVector& supportDataID = problemData->getIntVector("supportDataID");

	dbMatrix& dataParametersList = problemData->getDbMatrix(
			"dataParametersList");
	dataParametersList.resize(supportDataID.size());

	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");
	stepHistoryList.resize(supportDataID.size(), dbVector());

	dbMatrix& leftCavityVolumesList = problemData->getDbMatrix(
			"leftCavityVolumesList");
	leftCavityVolumesList.resize(supportDataID.size(), dbVector());

	dbMatrix& leftCavityPressuresList = problemData->getDbMatrix(
			"leftCavityPressuresList");
	leftCavityPressuresList.resize(supportDataID.size(), dbVector());

	dbMatrix& rightCavityVolumesList = problemData->getDbMatrix(
			"rightCavityVolumesList");
	rightCavityVolumesList.resize(supportDataID.size(), dbVector());

	dbMatrix& rightCavityPressuresList = problemData->getDbMatrix(
			"rightCavityPressuresList");
	rightCavityPressuresList.resize(supportDataID.size(), dbVector());

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "-----------------------------" << endl;
		logFile << "[" << i << "]Loading Data folder: "
				<< mainData.getFolderName() << endl;
		logFile << "-----------------------------" << endl;

		cout << "-----------------------------" << endl;
		cout << "[" << i << "]Loading Data folder: " << mainData.getFolderName()
				<< endl;
		cout << "-----------------------------" << endl;

		dataParametersList[i] = mainData.getParamValuesVec();

		// Read result from file
		mainData.readResultFile(InputData, logFile);
	}

	cout << "-----------------------------" << endl;

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Reading selected datasets completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "Reading selected datasets completed in: "
			<< elapsed_seconds.count() << " sec" << endl;

}


/*!****************************************************************************/
/*!****************************************************************************/
void MicroRVE::reducedCalculation(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	calculation(myDatabase,problemData,InputData,logFile);
}


/*!****************************************************************************/
/*!****************************************************************************/
void MicroRVE::postProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	defaultPostProcessing(myDatabase,problemData,InputData,logFile);

}




