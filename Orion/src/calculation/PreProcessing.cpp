/*
 * PreProcessing.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <PreProcessing.h>

//PreProcessing::PreProcessing(Data& yourData,GridNodes* yourGrid)
//										:myData(yourData),myGrid(yourGrid),
//										 DOFStandard(NULL){};

/*!****************************************************************************/
/*!****************************************************************************/
void PreProcessing::initPreProcessing(Database& myDatabase,
		DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile) {

	using namespace std;

	problemData->setValue("supportDataID", intVector());
	problemData->setValue("myParameters", dbVector());
	problemData->setValue("parameterRadii", dbVector());
	problemData->setValue("resultList", vector<dbMatrix>());
	problemData->setValue("dataParametersList", dbMatrix());
	problemData->setValue("calcResultList", vector<dbMatrix>());
	problemData->setValue("standardStepsList", dbMatrix());
	problemData->setValue("stepHistoryList", dbMatrix());

	if (InputData->getValue("standardiseDOF") == 1) {

		DOFStandard = new DOFStandardisation(myData);

		// Generated grid points
		if (InputData->getValue("gridNodesType") == 1) {
			myData.setAnchorPoint(InputData->getValue("anchorPoint"));
			myData.readMeshDataFile(InputData, logFile);
			myData.getMeshData()->writeMeshFile("fem_orion.msh", InputData,
					logFile);
		// geometry specific grid points (e.g heart)
		} else if (InputData->getValue("gridNodesType") == 2) {
			myData.readTransformedMeshDataFile(InputData, logFile);
			myData.getMeshData()->writeMeshFile("fem_orion_map.msh", InputData,
					logFile);

			DOFStandard->initDOFStandardisation(problemData, myDatabase,
					InputData, logFile);
			myGrid = DOFStandard->getGrid();

		// Arbitrary grid points (e.g cube)
		} else if (InputData->getValue("gridNodesType") == 3) {
			myData.readMeshDataFile(InputData, logFile);
			myData.getMeshData()->writeMeshFile("fem_orion.msh", InputData,
					logFile);
			DOFStandard->initDOFStandardisation(problemData, myDatabase,
					InputData, logFile);
			myGrid = DOFStandard->getGrid();
		}

	} else {
		myData.readMeshDataFile(InputData, logFile);
		myData.getMeshData()->writeMeshFile("fem_orion.msh", InputData,
				logFile);
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
void PreProcessing::defaultPreProcessing(Database& myDatabase,
	   DataContainer* problemData, InputFileData* InputData, ofstream& logFile){
	// Search database and extract the neighbouring nodes

	searchDatabase(problemData, myDatabase, InputData, logFile);

	// Carry out a preliminary check of the quality of MLS interpolants
	// based on the selected datasets
	preliminaryMLSInterpolantsCalc(myDatabase, problemData, InputData, logFile);

	dbVector& myParameters = problemData->getDbVector("myParameters");
	myData.setParamValuesVec(myParameters);

	loadSelectedData(problemData, myDatabase, InputData, logFile);

	settingCommonResultNames(myDatabase, problemData, InputData, logFile);

	if(InputData->getValue("standardiseStep") == 1)
		standardiseStep(problemData, InputData, logFile);
	else
		assembleResultsOnly(problemData, InputData, logFile);

	if(InputData->getValue("savePreProcessedData") == 1){
		savePreProcessedData(problemData,InputData,logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Extract the parameters for which the calculation will be carried out from
//! the input file
void PreProcessing::searchDatabase(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	// Extracting the main parameters from the inputfile
	extractParameters(problemData, myDatabase, InputData, logFile);

	// Searching supports for main parameters from database
	DatabaseQuery* DBQuery = new DatabaseQuery(problemData, myDatabase,
															InputData, logFile);
	delete DBQuery;

	// Print out selected datasets details
	intVector& supportDataID = problemData->getIntVector("supportDataID");
	logFile << "Number of datasets selected: " << supportDataID.size() << endl;
	cout << "Number of datasets selected: " << supportDataID.size() << endl;

	logFile << "List of Data selected" << endl;
	logFile << "*********************" << endl;
	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "(" << i << "): " << "Data[" << mainData.getId() << "]-> "
				<< mainData.getFolderName() << endl;

	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Extract the parameters for which the calculation will be carried out from
//! the input file
void PreProcessing::extractParameters(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	//! Check that only one type of material is present
	vector<map<string, double> >& materialsSet = InputData->getMaterials();
	if (materialsSet.size() > 1) {
		logFile << "Geometry with multiple material types is not yet supported"
				<< endl;
		cout << "Geometry with multiple material types is not yet supported"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//! Extract the material's parameters
	map<string, double>& material = materialsSet[0];

	//! Extract the names of the parameters from Database
	vector<string> paramNamesVec = myDatabase.getParamNamesVec();

	logFile << " Parameters in Database:" << endl;
	for (int i = 0; i < paramNamesVec.size(); i++) {
		logFile << paramNamesVec[i] << endl;
	}

	//! Extract the parameters' value
	vector<string> missingParameters;
	dbVector& myParameters = problemData->getDbVector("myParameters");
	myParameters.resize(paramNamesVec.size());
	for (int i = 0; i < paramNamesVec.size(); i++) {

		//! Search for specific parameter in the material's parameters list
		map<string, double>::iterator it_material = material.find(
				paramNamesVec[i]);

		if (it_material == material.end()) {
			// If specific parameter is not found, record it in
			// missing parameters vector
			missingParameters.resize(missingParameters.size() + 1,
					paramNamesVec[i]);
		} else {
			// If specific parameter is found, record it in parameters vector
			myParameters[i] = it_material->second;
		}
	}

	// If parameters are missing from the input file, print them out and
	// abort program
	if (missingParameters.size() > 0) {
		logFile << "ERROR: The following parameters are missing from the "
				"input file: " << endl;
		cout << "ERROR: The following parameters are missing from the "
				"input files: " << endl;

		for (int i = 0; i < missingParameters.size(); i++) {
			logFile << missingParameters[i] << " ";
			cout << missingParameters[i] << " ";
		}
		logFile << endl;
		cout << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

#ifdef _ROMCalcDebugMode_
	logFile << "******* Parameters Extracted *******" << endl;
	for (int j = 0; j < myParameters.size(); j++) {
		logFile << "[" << j << "] " << myParameters[j] << endl;
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void PreProcessing::loadSelectedData(DataContainer* problemData,
		Database& myDatabase, InputFileData* InputData, ofstream& logFile){


	intVector& supportDataID = problemData->getIntVector("supportDataID");

	dbMatrix& dataParametersList = problemData->getDbMatrix("dataParametersList");
	dataParametersList.resize(supportDataID.size());

	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");
	stepHistoryList.resize(supportDataID.size());


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

		//! Free memory
//		mainData.delMeshData();

		// Extract step-values for step standardisation
		stepHistoryList[i] = mainData.getStepValueVec();

		// Standardise the degrees of Freedom
//		if (InputData->getValue("standardiseDOF") == 1){
//			standardiseDOF(mainData, problemData, InputData, logFile);
//		}

	}

	if(InputData->getValue("standardiseDOF") == 1 &&
				InputData->getValue("gridNodesType") == 1 ){

			DOFStandard->initDOFStandardisation(problemData, myDatabase,InputData,
					logFile);
			myGrid = DOFStandard->getGrid();
	}

	for (int i = 0; i < supportDataID.size(); i++) {
		Data& mainData = myDatabase.getDataId(supportDataID[i]);
		standardiseDOF(mainData, problemData, InputData, logFile);
	}


	cout << "-----------------------------" << endl;
	cout << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void PreProcessing::standardiseDOF(Data& mainData, DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	logFile << "Starting DOF standardisation of Data " << mainData.getId()
			<< endl;

	if(InputData->getValue("gridNodesType") == 2){

		// Read the mapped geometry(used by point-in-polygon algorithm)
		mainData.readTransformedMeshDataFile(InputData, logFile);

		// Output results in map nodes format
		string mapResFileName = mainData.getFolderName() + "fem_map.res";
		mainData.saveAllResultsToFile_res_format(mapResFileName.c_str(), InputData, logFile);

		string mapMshFileName = mainData.getFolderName() + "fem_map.msh";
		mainData.getMeshData()->writeMeshFile(mapMshFileName.c_str(), InputData,
				logFile);
	}

	// -----------------------------------------------------------------
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	// Standardise all results
	DOFStandard->standardiseResultDOF(problemData, mainData, InputData,
			logFile);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "DOF standardisation of Data completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "DOF standardisation of Data completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	// -----------------------------------------------------------------

	if(InputData->getValue("gridNodesType") == 2){
	// Output results in grid nodes format
	string gridResFileName = mainData.getFolderName() + "fem_grid.res";
	mainData.saveAllResultsToFile_res_format(gridResFileName.c_str(), InputData, logFile);

	string gridMshFileName = mainData.getFolderName() + "fem_grid.msh";
	myGrid->getGridGeometry(logFile)->writeMeshFile(gridMshFileName.c_str(),
			InputData, logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// TODO: Debug
//! Test MLS interpolants
void PreProcessing::preliminaryMLSInterpolantsCalc(Database& myDatabase,
		DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile){

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	dbMatrix dataParametersList(supportDataID.size(),dbVector());

	logFile << "List of Data selected" << endl;
	logFile << "*********************" << endl;
	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		//! The displacement, parameters and step vector are extracted and
		//! stored
		dataParametersList[i] = mainData.getParamValuesVec();
	}

	dbVector& myParameters = problemData->getDbVector("myParameters");
	dbVector& parameterRadii = problemData->getDbVector("parameterRadii");

	printVector(parameterRadii,"parameterRadii",logFile);

	dbVector interpolants;

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

	Interpolation* MLSInterpolate = new Interpolation(myParameters,
					dataParametersList, parameterRadii, interpolants, InputData,
					logFile);

	myData.setInterpolants(interpolants);

	printVector(interpolants,"****** Interpolants ******",logFile);
		oPType intpolantSum = 0;
		for(int i=0; i<interpolants.size();i++) intpolantSum += interpolants[i];
		logFile << "Sum of interpolants: " << intpolantSum << endl;
		cout << "Interpolants quality(sum): " << intpolantSum << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
void PreProcessing::settingCommonResultNames(Database& myDatabase,
		DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile){

	vector < string > commonResultsNameList
						= findCommonResultName(problemData, myDatabase,
														InputData, logFile);

	myData.setResultNameList(commonResultsNameList);

}

/*!****************************************************************************/
/*!****************************************************************************/
vector<string> PreProcessing::findCommonResultName(DataContainer* problemData,
		Database& myDatabase, InputFileData* InputData, ofstream& logFile) {

	vector<vector<string> > commonResultNameListMatrix;
	vector<string> commonResultNameList;

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	for (int i = 0; i < supportDataID.size(); i++) {
		commonResultNameListMatrix.push_back(
				myDatabase.getDataId(supportDataID[i]).getResultNameList());
	}

	vector<vector<string> > comparisonMatrix = commonResultNameListMatrix;
	for (int i = 0; i < commonResultNameListMatrix.size(); i++) {
		for (int j = 0; j < commonResultNameListMatrix[i].size(); j++) {

			int counter = 0;

			for (int k = 0; k < comparisonMatrix.size(); k++) {
				for (int l = 0; l < comparisonMatrix[k].size(); l++) {

					if (commonResultNameListMatrix[i][j]
							== comparisonMatrix[k][l]) {
						comparisonMatrix[k].erase(
								comparisonMatrix[k].begin() + l);
						counter++;
						break;
					}
				}

			}

			if (counter == commonResultNameListMatrix.size())
				commonResultNameList.push_back(
						commonResultNameListMatrix[i][j]);
		}
	}

	return commonResultNameList;
}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void PreProcessing::standardiseStep(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	cout << "Step value standardisation" << endl;
	logFile << "Step value standardisation" << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	vector<string> commonResultsNameList = myData.getResultNameList();

	StepStandard = new StepStandardisation(problemData, InputData, logFile);

	// Find the common steps for calculation
	StepStandard->initStepStandardisation(problemData, InputData, logFile);

	recordStandardisedSteps(problemData, InputData, logFile);

	// Some additional calculations -> function is overloaded by each model
	preROMCalculationFunctions(problemData, InputData, logFile);

	// Standardise the step of all commonResult name
	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");
	intVector& supportDataID = problemData->getIntVector("supportDataID");
	for (int i = 0; i < commonResultsNameList.size(); i++) {

		logFile << "***** Standardising result: "
				<< commonResultsNameList[i] << " *****" << endl;

		myData.setResultDOF(commonResultsNameList[i].c_str(),
				myDatabase.getDataId(supportDataID[0]).getResultDOF(
						commonResultsNameList[i].c_str()));

		// Compile list of specific result result
		for (int j = 0; j < supportDataID.size(); j++) {
			resultList.push_back(
					myDatabase.getDataId(supportDataID[j]).getResult(
							commonResultsNameList[i].c_str()));
			myDatabase.getDataId(supportDataID[j]).deleteResult(
					commonResultsNameList[i].c_str());
		}

//		printMatrix(resultList[0], "Non-standardise result fields:",logFile);

		StepStandard->standardise(problemData, InputData, logFile);

		// ---------------------------------------------------------------------
		// Record results

		// save the step standardised results
		for (int j = 0; j < supportDataID.size(); j++) {
			myDatabase.getDataId(supportDataID[j]).setResult(
					commonResultsNameList[i].c_str(), resultList[j]);
		}

		vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));

	}

	delete StepStandard;


	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Step value standardisation completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "Step value standardisation completed in: "
			<< elapsed_seconds.count() << " sec" << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void PreProcessing::assembleResultsOnly(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	cout << "results assembly" << endl;
	logFile << "results assembly" << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	vector<string> commonResultsNameList = myData.getResultNameList();

	// Some additional calculations -> function is overloaded by each model
	preROMCalculationFunctions(problemData, InputData, logFile);

	interpolateStandardSteps(problemData, InputData, logFile);

	// Standardise the step of all commonResult name
	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");
	intVector& supportDataID = problemData->getIntVector("supportDataID");
	for (int i = 0; i < commonResultsNameList.size(); i++) {

		logFile << "***** Assembling result: " << commonResultsNameList[i] << " *****" << endl;

		myData.setResultDOF(commonResultsNameList[i].c_str(),
				myDatabase.getDataId(supportDataID[0]).getResultDOF(
						commonResultsNameList[i].c_str()));

		// Compile list of specific result result
		for (int j = 0; j < supportDataID.size(); j++) {
			resultList.push_back(
					myDatabase.getDataId(supportDataID[j]).getResult(
							commonResultsNameList[i].c_str()));
			myDatabase.getDataId(supportDataID[j]).deleteResult(
					commonResultsNameList[i].c_str());
		}

//		printMatrix(resultList[0], "Non-standardise result fields:",logFile);

//		StepStandard->standardise(problemData, InputData, logFile);

		// ---------------------------------------------------------------------
		// Record results

		// save the step standardised results
		for (int j = 0; j < supportDataID.size(); j++) {
			myDatabase.getDataId(supportDataID[j]).setResult(
					commonResultsNameList[i].c_str(), resultList[j]);
		}

		vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));

	}


	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "results assembly completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "results assembly completed in: "
			<< elapsed_seconds.count() << " sec" << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void PreProcessing::interpolateStandardSteps(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile){

	intVector& supportDataID = problemData->getIntVector("supportDataID");

	dbMatrix stepList(supportDataID.size(),dbVector());

	for (int i=0; i< supportDataID.size(); i++){
		stepList[i] = myDatabase.getDataId(supportDataID[i]).getStepValueVec();
	}

	dbVector interpolants = myData.getInterpolants();

	if(interpolants.size() != stepList.size()){
		logFile << "ERROR: In PreProcessing::interpolateStandardSteps, "
				"interpolants and steplist size does not match" << endl;

		cout << "ERROR: In PreProcessing::interpolateStandardSteps, "
				"interpolants and steplist size does not match" << endl;
	}

	dbVector interpolatedSteps(stepList[0].size(),0);

	for(int i = 0 ; i < interpolatedSteps.size(); i++){
		for (int j = 0; j < supportDataID.size(); j++){
			interpolatedSteps[i] += stepList[j][i]*interpolants[j];
		}
	}

	myData.setStepValueVec(interpolatedSteps);
	problemData->setValue("PODIStepValueVec", myData.getStepValueVec());
}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void PreProcessing::recordStandardisedSteps(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	// Record the standardised step values of the problem at hand
	dbVector standardSteps = StepStandard->getStandardSteps(
				myData.getInterpolants(), problemData, InputData, logFile);
	myData.setStepValueVec(standardSteps);
	problemData->setValue("PODIStepValueVec", myData.getStepValueVec());


	// Record the new standardised step values of each dataset
	vector<dbMatrix> phaseList(StepStandard->getNumPhase(),dbMatrix());

	for(int i=0; i<phaseList.size(); i++){
		string ss = "standardStepPhaseList" + std::to_string(i);
		phaseList[i] = problemData->getDbMatrix(ss.c_str());
	}

	dbMatrix combinedPhaseList =
		StepStandard->combinePhaseList(phaseList,problemData,InputData,logFile);

	printMatrix(combinedPhaseList,"combinedPhaseList",logFile);

//	logFile << "PreProcessing::recordStandardisedSteps" << endl;
//		cout << "PreProcessing::recordStandardisedSteps" << endl;
//		MPI_Abort(MPI_COMM_WORLD, 1);

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	for(int i=0; i<supportDataID.size();i++){
		myDatabase.getDataId(supportDataID[i]).setStepValueVec(combinedPhaseList[i]);
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void PreProcessing::savePreProcessedData(DataContainer* problemData,
								InputFileData* InputData, ofstream& logFile){

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	for(int i=0; i<supportDataID.size();i++){
		Data& sData = myDatabase.getDataId(supportDataID[i]);

		sData.saveResultsToFile(InputData,logFile);
		sData.calcLeftAndRightCavityVolumes(InputData,logFile);
	}

	logFile << "Saving standardised result data completed" << endl;
	cout << "Saving standardised result data completed" << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);

}
