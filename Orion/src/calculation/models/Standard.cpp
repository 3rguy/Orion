/*
 * Standard.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <models/Standard.h>

/*!****************************************************************************/
/*!****************************************************************************/
Standard::Standard(Database& myDB, DataContainer* problemData,
			 InputFileData* InputData, ofstream& logFile){

	cout << "Model Calculation: Standard" << endl;
	logFile << "Model Calculation: Standard" << endl;

	myDatabase = myDB;

	problemData->setValue("surfaceLoadList",dbMatrix());

}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::preProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	initPreProcessing(myDatabase,problemData,InputData,logFile);

	defaultPreProcessing(myDatabase,problemData,InputData,logFile);
}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::loadSelectedData(DataContainer* problemData, Database& myDatabase,
			InputFileData* InputData, ofstream& logFile){


	intVector& supportDataID = problemData->getIntVector("supportDataID");

	dbMatrix& dataParametersList = problemData->getDbMatrix("dataParametersList");
	dataParametersList.resize(supportDataID.size());

	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");
	stepHistoryList.resize(supportDataID.size());

	dbMatrix& surfaceForceList = problemData->getDbMatrix("surfaceLoadList");
	surfaceForceList.resize(supportDataID.size(),dbVector());

	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "-----------------------------" << endl;
		logFile << "[" << i <<"]Loading Data file: " << mainData.getFolderName() << endl;
		logFile << "-----------------------------" << endl;

		cout << "-----------------------------" << endl;
		cout << "[" << i <<"]Loading Data file: " << mainData.getFolderName() << endl;
		cout << "-----------------------------" << endl;

		dataParametersList[i] = mainData.getParamValuesVec();

		//! The Data read Displacement and Mesh file
		logFile << "Reading mesh file" << endl;
		mainData.readMeshDataFile(InputData, logFile);

		// Read result from file
		mainData.readResultFile(InputData, logFile);

		mainData.insertZeroResultFields(InputData, logFile);

		//! Free memory
//		mainData.delMeshData();

		string fileName = mainData.getFolderName() + "load_deformation-node42-DOF2.grf";
		dbMatrix defVolLoadMat;
		mainData.readGraphFile_grfFormat(fileName, defVolLoadMat, logFile);
		surfaceForceList[i] = defVolLoadMat[1];

		// Extract step-values for step standardisation
		stepHistoryList[i] = mainData.getStepValueVec();

		// Standardise the degrees of Freedom
		if (InputData->getValue("standardiseDOF") == 1){
			standardiseDOF(mainData, problemData, InputData, logFile);
		}
	}

	cout << "-----------------------------" << endl;
	cout << endl;


}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::reducedCalculation(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	calculation(myDatabase,problemData,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::postProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	defaultPostProcessing(myDatabase,problemData,InputData,logFile);

}

/*!************************************************************************/
/*!************************************************************************/
void Standard::preROMCalculationFunctions(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

		vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));
		intVector& supportDataID = problemData->getIntVector("supportDataID");
		dbVector surfaceLoad;
		for (int i = 0; i < supportDataID.size(); i++) {
			surfaceLoad =problemData->getDbMatrix("surfaceLoadList")[i];
			dbMatrix slMat(1, surfaceLoad);
			problemData->getDbMatrixVec("resultList").push_back(slMat);
			printVector(surfaceLoad,"surfaceLoadList",logFile);
			logFile<<"slMat: " << slMat.size() << " x " << slMat[0].size() << endl;
			for(int j=0; j<slMat.size();j++)
				printVector(slMat[0],"",logFile);
		}

		cout << "Standard::preROMCalculationFunctions" << endl;
		logFile << "Standard::preROMCalculationFunctions" << endl;

		StepStandard->standardise(problemData, InputData, logFile);

		dbMatrix surfaceLoadList(supportDataID.size(), dbVector());
		for (int i = 0; i < supportDataID.size(); i++) {
			problemData->getDbMatrix("surfaceLoadList")[i] =
					problemData->getDbMatrixVec("resultList")[i][0];

			surfaceLoadList[i] =
					problemData->getDbMatrixVec("resultList")[i][0];

			printVector(surfaceLoadList[i],"surfaceLoadList[i]", logFile);
		}
		problemData->getDbMatrix("surfaceLoadList") = surfaceLoadList;

		dbVector interplSurfaceLoad;
		for (int i = 0; i < surfaceLoadList[0].size(); i++) {
			double sum = 0;
			for (int j = 0; j < surfaceLoadList.size(); j++) {
				sum += surfaceLoadList[j][i] * myData.getInterpolants()[j];
			}
			interplSurfaceLoad.push_back(sum);
		}

		vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));

		myData.setGraph("surfaceLoad",interplSurfaceLoad);


}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::additionalPostProcessingFunctions(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	int graphX = 0; // Displacement
	int graphNode = 42;
	int graphDOF = 2;

	string resultNameX, resultNameY;
	dbVector graphX_data = myData.getGraphData(graphX, graphNode, graphDOF,
			resultNameX, InputData, logFile);
	dbVector graphY_data = myData.getGraph("surfaceLoad");


	ostringstream nodeConvert;   // stream used for the conversion
	nodeConvert << graphNode;

	ostringstream DOFConvert;   // stream used for the conversion
	DOFConvert << graphDOF;

	string graphFileName = "load_deformation-node"
			+ nodeConvert.str() + "-DOF" + DOFConvert.str() + ".grf";
	myData.saveGraphResultsToFile_grf_format(graphFileName, graphX_data,
			graphY_data, logFile);

}

