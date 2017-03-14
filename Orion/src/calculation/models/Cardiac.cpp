/*
 * Cardiac.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <models/Cardiac.h>

Cardiac::Cardiac(Database& myDB, DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	cout << "Model Calculation: Cardiac" << endl;
	logFile << "Model Calculation: Cardiac" << endl;

	myDatabase = myDB;

	problemData->setValue("leftCavityVolumesList",dbMatrix());
	problemData->setValue("leftCavityPressuresList",dbMatrix());
	problemData->setValue("rightCavityVolumesList",dbMatrix());
	problemData->setValue("rightCavityPressuresList",dbMatrix());
	problemData->setValue("leftCavityVolumes",dbMatrix());
	problemData->setValue("rightCavityVolumes",dbMatrix());

}


/*!****************************************************************************/
/*!****************************************************************************/
void Cardiac::preProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	initPreProcessing(myDatabase,problemData,InputData,logFile);

	defaultPreProcessing(myDatabase,problemData,InputData,logFile);

}


/*!****************************************************************************/
/*!****************************************************************************/
void Cardiac::reducedCalculation(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	calculation(myDatabase,problemData,InputData,logFile);
}


/*!****************************************************************************/
/*!****************************************************************************/
void Cardiac::postProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	defaultPostProcessing(myDatabase,problemData,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Cardiac::loadSelectedData(DataContainer* problemData,
		Database& myDatabase, InputFileData* InputData, ofstream& logFile){

	cout << "Cardiac::loadSelectedData" << endl;

	int isRV = InputData->getValue("isRightVentriclePresent");

	intVector& supportDataID = problemData->getIntVector("supportDataID");

	dbMatrix& dataParametersList = problemData->getDbMatrix("dataParametersList");
	dataParametersList.resize(supportDataID.size());

	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");
	stepHistoryList.resize(supportDataID.size(),dbVector());

	dbMatrix& leftCavityVolumesList = problemData->getDbMatrix("leftCavityVolumesList");
	leftCavityVolumesList.resize(supportDataID.size(),dbVector());

	dbMatrix& leftCavityPressuresList = problemData->getDbMatrix("leftCavityPressuresList");
	leftCavityPressuresList.resize(supportDataID.size(),dbVector());

	dbMatrix& rightCavityVolumesList = problemData->getDbMatrix("rightCavityVolumesList");
	rightCavityVolumesList.resize(supportDataID.size(),dbVector());

	dbMatrix& rightCavityPressuresList = problemData->getDbMatrix("rightCavityPressuresList");
	rightCavityPressuresList.resize(supportDataID.size(),dbVector());

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "-----------------------------" << endl;
		logFile << "[" << i <<"]Loading Data folder: " << mainData.getFolderName() << endl;
		logFile << "-----------------------------" << endl;

		cout << "-----------------------------" << endl;
		cout << "[" << i <<"]Loading Data folder: " << mainData.getFolderName() << endl;
		cout << "-----------------------------" << endl;

		dataParametersList[i] = mainData.getParamValuesVec();

		//! The Data read Displacement and Mesh file
		logFile << "Reading mesh file" << endl;
		mainData.readMeshDataFile(InputData, logFile);

		// Read result from file
		mainData.readResultFile(InputData, logFile);

		mainData.insertZeroResultFields(InputData, logFile);
		InputData->setValue("insertZeroResultFields",1);

		mainData.readLeftVentriclePVGraphResultFile(InputData,logFile);
		mainData.calcLeftCavityVolumes(InputData,logFile);
		leftCavityVolumesList[i] = mainData.getLeftCavityVolumes();
		leftCavityPressuresList[i] = mainData.getLeftCavityPressures();

		if(isRV == 1){
			mainData.readRightVentriclePVGraphResultFile(InputData,logFile);
			mainData.calcRightCavityVolumes(InputData,logFile);
			rightCavityVolumesList[i] = mainData.getRightCavityVolumes();
			rightCavityPressuresList[i] = mainData.getRightCavityPressures();
		}

		if(InputData->getValue("standardiseStep") == 1)
			mainData.syncCardiacTimeStepsAndResults(InputData,logFile);

		//! Free memory
//		mainData.delMeshData();

		string fileName = mainData.getFolderName() + "load_time-loadID_0.grf";
		dbMatrix defVolLoadMat;
		mainData.readGraphFile_grfFormat(fileName, defVolLoadMat, logFile);
		mainData.getStepValueVec() = defVolLoadMat[0];

		// Extract step-values for step standardisation
		stepHistoryList[i] = mainData.getStepValueVec();

		// Standardise the degrees of Freedom
//		if (InputData->getValue("standardiseDOF") == 1){
//			standardiseDOF(mainData, problemData, InputData, logFile);
//		}
	}

	cout << "-----------------------------" << endl;

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Reading selected datasets completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "Reading selected datasets completed in: "
			<< elapsed_seconds.count() << " sec" << endl;


	if(InputData->getValue("standardiseDOF") == 1){

//		cout << endl;
//		logFile << endl;

		cout << "DOF Standardisation activated" << endl;
		logFile << "DOF Standardisation activated" << endl;

		logFile << "******* DOF Standardisation  *******" << endl;

		if(InputData->getValue("gridNodesType") == 1 ){

			DOFStandard->setData(myData);

			DOFStandard->initDOFStandardisation(problemData, myDatabase,InputData,
					logFile);

			myGrid = DOFStandard->getGrid();
		}

		for (int i = 0; i < supportDataID.size(); i++) {
			Data& mainData = myDatabase.getDataId(supportDataID[i]);

			cout << "Standardising Data: " << i << endl;
			logFile << "Standardising Data: " << i << endl;

			standardiseDOF(mainData, problemData, InputData, logFile);
		}
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
//void Cardiac::defaultPostProcessing(Database& myDatabase,
//		DataContainer* problemData, InputFileData* InputData,
//		ofstream& logFile) {
//
//	logFile << "Saving to File" << endl;
//	cout << "Saving to File" << endl;
//
//#ifdef _ROMCalcDebugMode_
//	logFile << "******* New Calculated Result *******" << endl;
//	vector<string>& commonResultsNameList = myData.getResultNameList();
//	for (int i = 0; i < commonResultsNameList.size(); i++) {
//
//		logFile << "Result Type: " << commonResultsNameList[i] << endl;
//		logFile << "Num of DOFs = "
//				<< myData.getResultDOF(commonResultsNameList[i].c_str())
//				<< endl;
//		printMatrix(myData.getResult(commonResultsNameList[i].c_str()),
//				commonResultsNameList[i].c_str(), logFile);
//	}
//#endif
//
//	if (InputData->getValue("standardiseDOF") == 1) {
//
//		FEMGeometryExt* gridFEMGeo = myGrid->getGridGeometry(logFile);
//
//		string oldFolderName = myData.getFolderName();
//		myData.getFolderName() = string("postProc");
//		vector<Data*> dataList(1);
//		dataList.at(0) = &myData;
//
//		if(InputData->getValue("gridNodesResultPlot") == 1)
//			saveResultInGridNodesFormat(problemData, dataList, InputData, logFile);
//
//		myData.getFolderName() = oldFolderName;
//
//		// -------------------------------------------------------------------------
//		// Loading input parameters for the interpolation calculation
//		InputData->setValue("interpolantionType",
//				InputData->getValue("gInterpolantionType"));
//		InputData->setValue("influenceRangeFactor",
//				InputData->getValue("gInfluenceRangeFactor"));
//		InputData->setValue("MLSCalculationType",
//				InputData->getValue("gMLSCalculationType"));
//		InputData->setValue("MLSPolynomialDegree",
//				InputData->getValue("gMLSPolynomialDegree"));
//		InputData->setValue("MLSWeightFunc",
//				InputData->getValue("gMLSWeightFunc"));
//		InputData->setValue("parameterPolynomialDegree",
//				InputData->getValue("gparameterPolynomialDegree"));
//		myGrid->initCalcResultOnParticles(myData, InputData, logFile);
//
//
//		myData.setMeshData(gridFEMGeo);
//
//		vector<string>& resultNameList = myData.getResultNameList();
//		for (int i = 0; i < resultNameList.size(); i++) {
//
//			myData.assignResultToParticles(resultNameList[i].c_str(), InputData,
//					logFile);
//
//			dbMatrix resultMat = myGrid->interpolateResultOnGridPoint(myData,
//					InputData, logFile);
//
//			myGrid->resetNodesStepDOFMat();
//
//			myData.deleteResult(resultNameList[i].c_str());
//			myData.setResult(resultNameList[i].c_str(), resultMat);
//
//		}
//
//		myData.setMeshData(myGrid->getGridGeometry(logFile));
//		myGrid->setGridGeometry(gridFEMGeo);
//
//		myData.saveResultsToFile(logFile);
//		system("cp -f fem_orion.res fem_orion_map.res");
//
//		myData.delMeshData();
//		myData.readMeshDataFile(InputData, logFile);
//		myData.getMeshData()->writeMeshFile("fem_orion.msh",InputData,logFile);
//	}
//
//	myData.saveResultsToFile(logFile);
//	myData.calcLeftAndRightCavityVolumes(InputData,logFile);
//}

/*!****************************************************************************/
/*!****************************************************************************/
void Cardiac::additionalPostProcessingFunctions(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile){

	myData.calcLeftCavityVolumes(InputData,logFile);

	if(InputData->getValue("isRightVentriclePresent") == 1){
		myData.calcRightCavityVolumes(InputData,logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void Cardiac::preROMCalculationFunctions(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {


	// Standardise the left and right cavity pressures

	// =========================================================================
	// Left cavity pressure
	// =========================================================================
	logFile << "Standardising the left cavity pressure." << endl;
	vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));
	intVector& supportDataID = problemData->getIntVector("supportDataID");
	dbVector leftPressure;
	for (int i = 0; i < supportDataID.size(); i++) {
		leftPressure =
				myDatabase.getDataId(supportDataID[i]).getLeftCavityPressures();
		dbMatrix lpMat(1, leftPressure);
		problemData->getDbMatrixVec("resultList").push_back(lpMat);
	}

	if(InputData->getValue("standardiseStep") == 1)
		StepStandard->standardise(problemData, InputData, logFile);

	dbMatrix leftCavityPressuresList(supportDataID.size(), dbVector());
	for (int i = 0; i < supportDataID.size(); i++) {
		myDatabase.getDataId(supportDataID[i]).getLeftCavityPressures() =
				problemData->getDbMatrixVec("resultList")[i][0];

		leftCavityPressuresList[i] =
				problemData->getDbMatrixVec("resultList")[i][0];

		printVector(
				myDatabase.getDataId(supportDataID[i]).getLeftCavityPressures(),
				"getLeftCavityPressures", logFile);
	}

	dbVector interplLeftCavityPressures;
	for (int i = 0; i < leftCavityPressuresList[0].size(); i++) {
		double sum = 0;
		for (int j = 0; j < leftCavityPressuresList.size(); j++) {
			sum += leftCavityPressuresList[j][i] * myData.getInterpolants()[j];
		}
		interplLeftCavityPressures.push_back(sum);
	}
	myData.getLeftCavityPressures() = interplLeftCavityPressures;

	printVector(myData.getLeftCavityPressures(),
			"myData.getLeftCavityPressures()", logFile);

	// =========================================================================
	// Right cavity pressure
	// =========================================================================
	int isRV = InputData->getValue("isRightVentriclePresent");
	if(isRV == 1){
		logFile << "Standardising the right cavity pressure." << endl;
		vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));
		dbVector rightPressure;
		for (int i = 0; i < supportDataID.size(); i++) {
			rightPressure =
					myDatabase.getDataId(supportDataID[i]).getRightCavityPressures();
			dbMatrix rpMat(1, rightPressure);
			problemData->getDbMatrixVec("resultList").push_back(rpMat);
		}

		if(InputData->getValue("standardiseStep") == 1)
			StepStandard->standardise(problemData, InputData, logFile);

		dbMatrix rightCavityPressuresList(supportDataID.size(), dbVector());
		for (int i = 0; i < supportDataID.size(); i++) {
			myDatabase.getDataId(supportDataID[i]).getRightCavityPressures() =
					problemData->getDbMatrixVec("resultList")[i][0];

			rightCavityPressuresList[i] =
					problemData->getDbMatrixVec("resultList")[i][0];
			printVector(
					myDatabase.getDataId(supportDataID[i]).getRightCavityPressures(),
					"getRightCavityPressures", logFile);
		}

		dbVector interplrightCavityPressures;
		for (int i = 0; i < rightCavityPressuresList[0].size(); i++) {
			double sum = 0;
			for (int j = 0; j < rightCavityPressuresList.size(); j++) {
				sum += rightCavityPressuresList[j][i] * myData.getInterpolants()[j];
			}
			interplrightCavityPressures.push_back(sum);
		}
		myData.getRightCavityPressures() = interplrightCavityPressures;

		printVector(myData.getRightCavityPressures(),
				"myData.getRightCavityPressures()", logFile);
	}

	vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));

}
