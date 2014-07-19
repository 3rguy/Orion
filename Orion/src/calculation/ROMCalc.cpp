#include "ROMCalc.h"

ROMCalc::ROMCalc(InputFileData* InputData, ofstream& logFile) :
		myData(Data()) {

	//! Setting-up database
	string databaseName = "myDatabase.csv";
	Database myDatabase(databaseName, logFile);

	DataContainer* problemData = new DataContainer();

	problemData->setValue("supportDataID",intVector());
	problemData->setValue("myParameters",dbVector());
	problemData->setValue("parameterRadii",dbVector());
	problemData->setValue("resultList",vector<dbMatrix>());
	problemData->setValue("dataParametersList",dbMatrix());
	problemData->setValue("calcResultList",vector<dbMatrix>());
	problemData->setValue("standardSteps",dbVector());

	myData.readMeshDataFile(InputData, logFile);
	system("cp -f fem.msh fem_orion.msh");

	//! Carry out preprocessing
	preProcessing(problemData, myDatabase,InputData, logFile);

	//! PODI Calculation
	ROMCalculationSet(problemData, myDatabase,InputData, logFile);

	//! Saving results to file
	postProcessing(problemData, myDatabase, InputData, logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Preprocessing:
void ROMCalc::preProcessing(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	dbMatrix& dataParametersList = problemData->getDbMatrix("dataParametersList");
	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");

	// Extracting the main parameters from the inputfile
	extractParameters(problemData, myDatabase, InputData, logFile);

	// Searching supports for main parameters from database
	searchDatabase(problemData, myDatabase, InputData, logFile);

	// Load-up all geometry and displacement data files
	dataParametersList.resize(supportDataID.size());

	dbMatrix stepHistoryList(supportDataID.size(), dbVector());

	logFile << "List of Data selected" << endl;
	logFile << "*********************" << endl;
	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "(" << i << "): " << "Data[" << mainData.getId() << "]-> "
				<< mainData.getFolderName() << endl;

		//! The Data read Displacement and Mesh file
		logFile << "Reading mesh file" << endl;
		mainData.readMeshDataFile(InputData, logFile);

		//! The displacement, parameters and step vector are extracted and
		//! stored
		dataParametersList[i] = mainData.getParamValuesVec();

		stepHistoryList[i] = mainData.getStepValueVec();
		printVector(stepHistoryList[i],"stepHistoryList",logFile);

	}

	problemData->setValue("stepHistoryList",stepHistoryList);

	dbVector& myParameters = problemData->getDbVector("myParameters");
	myData.setParamValuesVec(myParameters);
	myData.readMeshDataFile(InputData, logFile);
	myData.setAnchorPoint(InputData->getValue("anchorPoint"));

	standardisationProcess(problemData, myDatabase, InputData, logFile);

#ifdef _ROMCalcDebugMode_
	logFile << "******* stepHistoryList *******" << endl;
	for (int i = 0; i < stepHistoryList.size(); i++) {
		logFile << "[" << i << "]: ";
		for (int j = 0; j < stepHistoryList[i].size(); j++) {
			logFile << stepHistoryList[i][j] << ", ";
		}
		logFile << endl;
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Extract the parameters for which the calculation will be carried out from
//! the input file
void ROMCalc::extractParameters(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	//! Check that only one type of material is present
	vector<map<string, double> > materialsSet = InputData->getMaterials();
	if (materialsSet.size() > 1) {
		logFile << "Geometry with multiple material types is not yet supported"
				<< endl;
		cout << "Geometry with multiple material types is not yet supported"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//! Extract the material's parameters
	map<string, double> material = materialsSet[0];

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
//! Search for the influencing particle in the database for the
//! interpolation process
void ROMCalc::searchDatabase(DataContainer* problemData, Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	cout << "Searching database" << endl;
	logFile << "Searching database" << endl;

	dbVector& myParameters = problemData->getDbVector("myParameters");
	logFile << "myParameters: " <<myParameters.size() << endl;
	intVector& supportDataID = problemData->getIntVector("supportDataID");
	logFile << "supportDataID: " <<supportDataID.size() << endl;
	dbVector& parameterRadii = problemData->getDbVector("parameterRadii");
	logFile << "parameterRadii: " <<parameterRadii.size() << endl;

	// Extract the radius factor which will determine the influence range
	double radiusFactor = InputData->getValue("influenceRangeFactor");

	// Determine the max and min value of each parameter(i.e Range) and store
	// them in the matrix "myParamaterInfluenceRange"
	int counter = 0;
	dbMatrix myParamaterInfluenceRange(myParameters.size(), dbVector(2));
	parameterRadii.resize(myParameters.size());
	for (int i = 0; i < myParameters.size(); i++) {
		myParamaterInfluenceRange[counter][0] = myParameters[i]
				- (myParameters[i] * radiusFactor); 	// Min value
		myParamaterInfluenceRange[counter][1] = myParameters[i]
				+ (myParameters[i] * radiusFactor);	// Max value
		parameterRadii[counter] = myParameters[i] * radiusFactor;

		logFile << "Range: " << myParamaterInfluenceRange[counter][0]
				<< " < Parameter[" << i << "] < "
				<< myParamaterInfluenceRange[counter][1] << endl;

		counter++;
	}

	// Extract the list of data
	vector<Data> dataList = myDatabase.getDataList();

	// Determine if the data parameters are within the influence range
	dbVector dataParametersVec;
	int parameterInsideCounter = 0;
	for (int i = 0; i < myDatabase.size(); i++) {

		// Extract parameters of of a particular data
		dataParametersVec = dataList[i].getParamValuesVec();

		// Loop over each parameter and check if they are within the range
		// specified. If one of them is not, the rest of the parameters are
		// skipped.
		parameterInsideCounter = 0;
		for (int j = 0; j < myParamaterInfluenceRange.size(); j++) {

			logFile << endl << myParamaterInfluenceRange[j][0] << " <= "
					<< dataParametersVec[j] << " <= "
					<< myParamaterInfluenceRange[j][1] << endl;

			if (myParamaterInfluenceRange[j][0] <= dataParametersVec[j]
					&& dataParametersVec[j] <= myParamaterInfluenceRange[j][1])
				parameterInsideCounter++;
			else
				// No need to continue the comparison process if one of the
				// parameters is out of range
				break;
		}

		// If the all parameters are inside the influence range, record ID of
		// Data
		if (parameterInsideCounter == myParamaterInfluenceRange.size()) {
			supportDataID.resize(supportDataID.size() + 1);
			supportDataID[supportDataID.size() - 1] = dataList[i].getId();
			logFile << "----> Taken" << endl;
		}
	}

#ifdef _ROMCalcDebugMode_
	logFile << "******* Supporting Data *******" << endl;
	for (int j = 0; j < supportDataID.size(); j++) {
		logFile << j << "). ID = " << supportDataID[j] << endl;
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void ROMCalc::ROMCalculationSet(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	intVector& supportDataID = problemData->getIntVector("supportDataID");

//	if (InputData->getValue("standardiseDOF") == 1){
//
//		for (int i = 0; i < supportDataID.size(); i++) {
//			Data& mainData = myDatabase.getDataId(supportDataID[i]);
//
//			string inputFileName = mainData.getFolderName()
//										+ "fem_standardDOF.res";
//			mainData.readResultFile_resFormat(inputFileName,logFile);
//		}
//	}


	// Find the list of results available for calculation
	vector<string> commonResultsNameList =
			findCommonResultName(problemData,myDatabase,InputData,logFile);
	myData.setResultNameList(commonResultsNameList);

	// Find the common steps for calculation
	dbVector& standardSteps = problemData->getDbVector("standardSteps");
	standardSteps = findCommonStepsList(problemData,myDatabase,InputData,logFile);
	myData.setStepValueVec(standardSteps);


	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");
	problemData->setValue("calcResultList",dbMatrix());

	cout << "############# Reduced Order Calculation ############# " << endl;
	logFile << "############# Reduced Order Calculation ############# " << endl;

	for (int i = 0; i < commonResultsNameList.size(); i++) {

		cout << "*****************************************************" << endl;
		cout << "Processing result: " << commonResultsNameList[i] << endl;

		logFile << "*****************************************************"
				<< endl;

		logFile << "Processing result: " << commonResultsNameList[i] << endl;

		for (int j = 0; j < supportDataID.size(); j++) {
			resultList.push_back(
					myDatabase.getDataId(supportDataID[j]).getResult(
							commonResultsNameList[i].c_str()));
			myDatabase.getDataId(supportDataID[j]).deleteResult(
					commonResultsNameList[i].c_str());
		}

		stepStandardisation(problemData, InputData, logFile);

		reducedOrderMethodCalc(problemData, InputData, logFile);

		//myData.deleteResult(commonResultsNameList[i].c_str());
		myData.setResult(commonResultsNameList[i].c_str(),
				problemData->getDbMatrix("calcResultList"));
		dbMatrix().swap(problemData->getDbMatrix("calcResultList"));
		myData.setResultDOF(commonResultsNameList[i].c_str(),
				myDatabase.getDataId(supportDataID[0]).getResultDOF(
						commonResultsNameList[i].c_str()));

		vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));

	}

#ifdef _ROMCalcDebugMode_
	logFile << "******* Calculated Result *******" << endl;
	for (int i = 0; i < commonResultsNameList.size(); i++) {

		logFile << "Result Type: " << commonResultsNameList[i] << endl;
		logFile << "Num of DOFs = "
				<< myData.getResultDOF(commonResultsNameList[i].c_str())
				<< endl;
		printMatrix(myData.getResult(commonResultsNameList[i].c_str()),
				commonResultsNameList[i].c_str(), logFile);
	}
#endif

	cout << "##################################################### " << endl;
	logFile << "##################################################### " << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
vector<string> ROMCalc::findCommonResultName(DataContainer* problemData,
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
dbVector ROMCalc::findCommonStepsList(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile){

	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");
	dbVector standardSteps;

	int position = 0;
	for (int i = 0; i < stepHistoryList.size(); i++) {
		for (int j = 0; j < stepHistoryList[i].size(); j++) {
			position = findDoubleVecPos(stepHistoryList[i][j], 0,
					standardSteps.size(), standardSteps);
			if (position == -1) {
				standardSteps.push_back(stepHistoryList[i][j]);
			}
		}
	}
#ifdef _ROMCalcDebugMode_
	logFile << "******* Unsorted Standard Steps *******" << endl;

	logFile << "Steps are: ";
	for (int j = 0; j < standardSteps.size(); j++) {
		logFile << standardSteps[j] << ", ";
	}
	logFile << endl;
#endif

	sortDoubleVector(standardSteps, 0, standardSteps.size() - 1);

#ifdef _ROMCalcDebugMode_
	logFile << "******* Sorted Standard Steps *******" << endl;

	logFile << "Steps are: ";
	for (int j = 0; j < standardSteps.size(); j++) {
		logFile << standardSteps[j] << ", ";
	}
	logFile << endl;
#endif

	return standardSteps;

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void ROMCalc::reducedOrderMethodCalc(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	cout << "PODI Calculation started " << endl;
	logFile << "PODI Calculation started " << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	//! ROM Calculation
	PODICalc* Podi = new PODICalc(problemData->getDbVector("myParameters"),
			problemData->getDbVector("parameterRadii"),
			problemData->getIntVector("supportDataID"),
			problemData->getDbMatrix("dataParametersList"),
			problemData->getDbMatrixVec("resultList"),
			problemData->getDbMatrix("calcResultList"),
			InputData, logFile);

	delete Podi;

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "PODI Calculation completed in: " << elapsed_seconds.count()
			<< " sec" << endl;
	logFile << "PODI Calculation completed in: " << elapsed_seconds.count()
			<< " sec" << endl;

#ifdef _ROMCalcDebugMode_
//	printMatrix(problemData->getDbMatrix("resultingDisplacementMatrix"),
//			"***** Resulting Displacement Matrix *****", logFile);
	printMatrix(problemData->getDbMatrix("calcResultList"),
				"***** Calculated Result Displacement List *****", logFile);
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Post Processing: Write the calculated displacement to file
void ROMCalc::postProcessing(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	logFile << "Saving to File" << endl;
	cout << "Saving to File" << endl;

#ifdef _ROMCalcDebugMode_
	logFile << "******* New Calculated Result *******" << endl;
	vector<string>& commonResultsNameList = myData.getResultNameList();
	for (int i = 0; i < commonResultsNameList.size(); i++) {

		logFile << "Result Type: " << commonResultsNameList[i] << endl;
		logFile << "Num of DOFs = "
				<< myData.getResultDOF(commonResultsNameList[i].c_str())
				<< endl;
		printMatrix(myData.getResult(commonResultsNameList[i].c_str()),
				commonResultsNameList[i].c_str(), logFile);
	}
#endif

	if(InputData->getValue("standardiseDOF") == 1){

		FEMGeometryExt* gridFEMGeo = myGrid->getGridGeometry(logFile);

		myGrid->initCalcResultOnParticles(myData,InputData,logFile);
		myData.setMeshData(gridFEMGeo);

		vector<string>& resultNameList = myData.getResultNameList();
		for(int i = 0 ; i < resultNameList.size() ; i++){

			myData.assignResultToParticles(resultNameList[i].c_str(),InputData,logFile);

			dbMatrix resultMat = myGrid->interpolateResultOnGridPoint(myData,
					InputData, logFile);

			myGrid->resetNodesStepDOFMat();

			myData.deleteResult(resultNameList[i].c_str());
			myData.setResult(resultNameList[i].c_str(),resultMat);

		}

		myData.setMeshData(myGrid->getGridGeometry(logFile));
		myGrid->setGridGeometry(gridFEMGeo);
	}

//	dbVector sumStepValueVec(myDatabase.getDataId(supportDataID[0]).getStepValueVec().size(),0);
//	dbVector avgStepValueVec = sumStepValueVec;
//	for(int i = 0; i < supportDataID.size(); i++){
//		dbVector tempVec = myDatabase.getDataId(supportDataID[i]).getStepValueVec();
//		logFile << "tempVec size(" << tempVec.size() <<") | sumStepValueVec(" << sumStepValueVec.size() << ")" << endl;
//		for(int j = 0; j < tempVec.size(); j++){
//			sumStepValueVec[j] += tempVec[j];
//			logFile << "tempVec[j]: " << tempVec[j] << "|j:"<<j<<"/"<<tempVec.size()<<"|i:"<<i<<"/"<< supportDataID.size() <<endl;
//		}
//	}
//
//	for(int i = 0; i < sumStepValueVec.size(); i++){postProcessing
//		avgStepValueVec[i] = sumStepValueVec[i]/supportDataID.size();
//	}

	// Data myData(myParameters);
	// myData.setStepValueVec(avgStepValueVec);

//	myGrid->
//		setDispOnNodes(myData,resultingDisplacementMatrix,InputData,logFile);
//
//	resultingDisplacementMatrix.clear();
//	dbMatrix().swap(resultingDisplacementMatrix);
//
//	resultingDisplacementMatrix =
//			myGrid->calcDispOnParticles(myData, InputData, logFile);

//#ifdef _ROMCalcDebugMode_
//	printMatrix(problemData->getDbMatrix("resultingDisplacementMatrix"),
//			"***** Resulting Displacement Matrix prior to saving *****",logFile);
//#endif

//	cout << endl << "--------------------------------" << endl;
//	MPI_Abort(MPI_COMM_WORLD, 1);

//	myData.setStepValueVec(problemData->getDbVector("standardSteps"));
//	myData.setDisplacement(problemData->getDbMatrix("resultingDisplacementMatrix"));
//	myData.saveDispFile(logFile);

	myData.saveResultsToFile(logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! findSupports: Finding the supports of a particular point
void ROMCalc::findSupports(dbVector& iPoint, dbMatrix& pointList,
		intVector& supportsList, dbVector& pointRadii, InputFileData* InputData,
		ofstream& logFile) {

	double radiusFactor = InputData->getValue("influenceRangeFactor");
	radiusFactor = 1.5;		// Temporary solution -> need to find a better one

	// Determine the max and min value of each parameter(i.e Range) and store
	// them in the matrix "myParamaterInfluenceRange"
	int counter = 0;
	dbMatrix pointInfluenceRange(iPoint.size(), dbVector(2));
	pointRadii.resize(iPoint.size());
	for (int i = 0; i < iPoint.size(); i++) {
		pointInfluenceRange[counter][0] = iPoint[i]
				- (iPoint[i] * radiusFactor); 	// Min value
		pointInfluenceRange[counter][1] = iPoint[i]
				+ (iPoint[i] * radiusFactor);	// Max value
		pointRadii[counter] = iPoint[i] * radiusFactor;
		counter++;
	}

	// Determine if the data parameters are within the influence range
	dbVector dataParametersVec;
	int pointInsideCounter = 0;
	for (int i = 0; i < pointList.size(); i++) {

		// Extract parameters of of a particular data
		dbVector pointVec = pointList[i];

		// Loop over each parameter and check if they are within the range
		// specified. If one of them is not, the rest of the parameters are
		// skipped.
		pointInsideCounter = 0;
		for (int j = 0; j < pointInfluenceRange.size(); j++) {

			logFile << pointInfluenceRange[j][0] << " < " << pointVec[j]
					<< " < " << pointInfluenceRange[j][1] << endl;

			if (pointInfluenceRange[j][0] <= pointVec[j]
					&& pointVec[j] <= pointInfluenceRange[j][1]) {
				pointInsideCounter++;
				logFile << "Taken" << endl << endl;

			} else
				// No need to continue the comparison process if one of the
				// parameters is out of range
				break;
		}

		// If the all parameters are inside the influence range, record ID of
		// Data
		if (pointInsideCounter == pointInfluenceRange.size()) {
			supportsList.resize(supportsList.size() + 1);
			supportsList[supportsList.size() - 1] = i;
		}
	}

	if (supportsList.size() < 2) {
		logFile << "ERROR: Point does not have enough of supports" << endl;
		cout << "ERROR: Point does not have enough of supports" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

#ifdef _ROMCalcDebugMode_
	logFile << "******* Supporting Points *******" << endl;

	logFile << "Main point: ";
	for (int i = 0; i < iPoint.size(); i++) {
		logFile << iPoint[i] << ",";
	}
	logFile << endl;

	logFile << "Supports are: " << endl;
	for (int j = 0; j < supportsList.size(); j++) {
		logFile << "[" << j << "] : ";
		for (int k = 0; k < pointList[supportsList[j]].size(); k++) {
			logFile << pointList[supportsList[j]][k] << ",";
		}
		logFile << endl;
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void ROMCalc::standardisationProcess(DataContainer* problemData,
		Database& myDatabase, InputFileData* InputData, ofstream& logFile){

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	if (InputData->getValue("standardiseDOF") == 1) {
		initDOFStandardisation(problemData, myDatabase, InputData, logFile);
	}


	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");
	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "-----------------------------" << endl;
		logFile << "Loading Data file: " << mainData.getFolderName() << endl;
		logFile << "-----------------------------" << endl;

		cout << "-----------------------------" << endl;
		cout << "Loading Data file: " << mainData.getFolderName() << endl;
		cout << "-----------------------------" << endl;

		// Read result from file
		mainData.readResultFile(InputData, logFile);

		// Extract step-values for step standardisation
		stepHistoryList[i] = mainData.getStepValueVec();

		// Standardise the degrees of Freedom
		if (InputData->getValue("standardiseDOF") == 1){

			logFile << "Starting DOF standardisation of Data "
					<< supportDataID[i] << endl;

			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();


			standardiseResultDOF(problemData,mainData,InputData,logFile);

			string inputFileName = mainData.getFolderName()
							+ "fem_standardDOF.res";
			mainData.saveAllResultsToFile_res_format(inputFileName.c_str(),logFile);


			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;

			cout << "DOF standardisation of Data completed in: "
				<< elapsed_seconds.count() << " sec" << endl;
			logFile << "DOF standardisation of Data completed in: "
					<< elapsed_seconds.count() << " sec" << endl;
		}

	}

	cout << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! findSupports: Finding the supports of a particular point
void ROMCalc::stepStandardisation(DataContainer* problemData,
		InputFileData* InputData,ofstream& logFile) {

	// Define a list of standard steps

	cout << "Step value standardisation" << endl;
	logFile << "Step value standardisation" << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	dbVector& standardSteps = problemData->getDbVector("standardSteps");
	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");

	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");


	InputData->setValue("dispMatrixRearrange", 0);
//	vector<dbMatrix> newDisplacementList(supportDataID.size(), dbMatrix());

	//! TODO: (1)Implement method to skip the following code in case the steps are constants
	//! TODO: (2)Check that interpolated displacement are not used for other interpolation calculation (??)
	int totalNumDofs = resultList[0].size();

	for (int i = 0; i < supportDataID.size(); i++) {
		dbMatrix tempDispMatrix(totalNumDofs, dbVector(standardSteps.size()));
		for (int j = 0; j < standardSteps.size(); j++) {
			int position = findDoubleVecPos(standardSteps[j], 0,
					stepHistoryList[i].size(), stepHistoryList[i]);
			if (position != -1) {
				for (int k = 0; k < totalNumDofs; k++)
					tempDispMatrix[k][j] = resultList[i][k][position];
			} else {
				dbVector iPoint;
				iPoint.push_back(standardSteps[j]);

				dbMatrix pointList(stepHistoryList[i].size(), dbVector(1));
				for (int k = 0; k < stepHistoryList[i].size(); k++) {
					pointList[k][0] = stepHistoryList[i][k];
				}

				dbVector pointRadii;
				intVector supportsList;

				findSupports(iPoint, pointList, supportsList, pointRadii,
						InputData, logFile);

				dbMatrix supportPointList(supportsList.size(), dbVector(1));
				for (int i = 0; i < supportsList.size(); i++) {
					supportPointList[i] = pointList[supportsList[i]];
				}

				vector<dbMatrix> tempDispList(1,
						dbMatrix(resultList[i].size(),
								dbVector(supportsList.size(), 0)));

				for (int k = 0; k < supportsList.size(); k++) {
					for (int l = 0; l < resultList[i].size(); l++) {
						tempDispList[0][l][k] =
								resultList[i][l][supportsList[k]];
					}
				}

				dbMatrix tempResultingDispMatrix;

#ifdef _ROMCalcDebugMode_
				printVector(iPoint, "iPoint:", logFile);
				printVector(pointRadii, "pointRadii:", logFile);
				printVector(supportsList, "supportsList:", logFile);
				printMatrix(supportPointList, "supportPointList:", logFile);
				printMatrix(pointList, "pointList:", logFile);
#endif

				PODICalc* Podi = new PODICalc(iPoint, pointRadii, supportsList,
						supportPointList, tempDispList, tempResultingDispMatrix,
						InputData, logFile);

				for (int l = 0; l < totalNumDofs; l++)
					tempDispMatrix[l][j] = tempResultingDispMatrix[l][0];

				delete Podi;

			}
		}

		resultList[i] = tempDispMatrix;
	}

	InputData->setValue("dispMatrixRearrange", 1);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Step value standardisation completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "Step value standardisation completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
}

/*!****************************************************************************/
/*!****************************************************************************/
void ROMCalc::initDOFStandardisation(DataContainer* problemData,
		Database& myDatabase, InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	cout << "DOF Standardisation activated" << endl;
	logFile << "******* DOF Standardisation activated *******" << endl;

	intVector& supportDataID = problemData->getIntVector("supportDataID");

	dbMatrix anchorPointList(supportDataID.size(), dbVector(0));
	dbMatrix maxCoordRange;

	vector<Data> dataList(supportDataID.size());

	// Preliminary stuff
	// -----------------
	//! For each data set: Assign displacement field to every particles and reset
	//! their coordinate with respect to the defined anchor point
	if (InputData->getValue("automaticGridNodeSetup") == 1) {
		for (int i = 0; i < supportDataID.size(); i++) {

			// reset particle's coordinates
			logFile << "Coordinate setup of Data: " << supportDataID[i] << endl;
			coordinateSetup(myDatabase.getDataId(supportDataID[i]),
					maxCoordRange, InputData, logFile);
			myDatabase.getDataId(supportDataID[i]).getMeshData()->writeMeshFile(
					InputData, logFile);

		}

#ifdef _ROMCalcDebugMode_
		logFile << "******* New coordinates *******"
				<< endl;
		for (int i = 0; i < supportDataID.size(); i++) {
			logFile << "Supporting Data: " << supportDataID[i] + 1 << endl;

			vector<ParticleExt>& ptcls =
					myDatabase.getDataId(supportDataID[i]).getMeshData()->getNodesVec();
			for (int j = 0; j < ptcls.size(); j++) {
				logFile << "Particle[" << j + 1 << "] New Coordinate: ";
				for (int k = 0; k < ptcls[j].getCoords().size(); k++) {
					logFile << ptcls[j].getCoords()[k] << " ";
				}
				logFile << endl;
			}
		}
		logFile << "******************************" << endl;
#endif

		//! Reset particle's coordinates for the unknown geometry
		logFile << "Coordinate setup for myData" << endl;
		coordinateSetup(myData, maxCoordRange, InputData, logFile);
//	myData.meshDataAccess()->writeMeshFile(InputData, logFile);

		// Set maxCoordRange
		if (maxCoordRange.size() == 3) {

			maxCoordRange.resize(3, dbVector(2));

			InputData->setValue("GridNodeMinX", maxCoordRange[0][0]);
			InputData->setValue("GridNodeMaxX", maxCoordRange[0][1]);

			InputData->setValue("GridNodeMinY", maxCoordRange[1][0]);
			InputData->setValue("GridNodeMaxY", maxCoordRange[1][1]);

			InputData->setValue("GridNodeMinZ", maxCoordRange[2][0]);
			InputData->setValue("GridNodeMaxZ", maxCoordRange[2][1]);

		} else {
			cout << " In ROMCalc::DOFStandardisation, dimension greater "
					"than 3 is not supported yet " << endl;
			logFile << " In ROMCalc::DOFStandardisation, dimension greater "
					"than 3 is not supported yet " << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	//! Setup the reference grid and define the nodes
	cout << "Setting up benchmark grid" << endl;
	logFile << "Setting up benchmark grid" << endl;
	myGrid = new GridNodes(InputData, logFile);


	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;


	cout << "DOFStandardisation initialisation completed in "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "DOFStandardisation initialisation completed in "
			<< elapsed_seconds.count() << " sec" << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
void ROMCalc::standardiseResultDOF(DataContainer* problemData, Data& mainData,
		InputFileData* InputData, ofstream& logFile) {

	intVector& supportDataID = problemData->getIntVector("supportDataID");

	//! Interpolate each displacement field on the reference grid
	dbMatrix standardResult;
	string choice;

	vector<string>& resultNameList = mainData.getResultNameList();

	myGrid->interpolantSetup(mainData, InputData, logFile);

	for (int i = 0; i < resultNameList.size(); i++) {

		// assign displacement field to particles
		mainData.assignResultToParticles(resultNameList[i].c_str(), InputData,
				logFile);

		standardResult = myGrid->interpolateResultOnGridPoint(mainData,
				InputData, logFile);

		myGrid->resetNodesStepDOFMat();

		mainData.getMeshData()->writeMeshFile(InputData, logFile);

		mainData.deleteResult(resultNameList[i].c_str());

		mainData.setResult(resultNameList[i].c_str(), standardResult);
	}

	myGrid->resetNodes();

}

/*!****************************************************************************/
/*!****************************************************************************/
void ROMCalc::coordinateSetup(Data& sData, dbMatrix& maxCoordRange,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	vector<ParticleExt>& particles = sData.getMeshData()->getNodesVec();
	int anchorPoint = sData.getAnchorPoint();

	//!Extract coordinates of anchor point
	dbVector anchorCoords;
	if(anchorPoint > -1 && anchorPoint < particles.size()+1)
		anchorCoords = particles[anchorPoint-1].getCoords();
	else{
		cout << "Anchor Point does not exist or has not been specified in"
				" the input file" << endl;
		logFile << "Anchor Point does not exist or has not been specified in"
				" the input file" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

#ifdef _ROMCalcDebugMode_
	logFile << "Coordinates of AnchorPoint: ";
	for(int i=0; i<anchorCoords.size();i++){
		logFile << anchorCoords[i] << ", ";
	}
	logFile << endl;
#endif

	for (int j = 0; j < particles.size(); j++) {
		for (int k = 0; k < particles[j].getCoords().size(); k++) {

			// Reset coordinates with respect to anchor point
			logFile << "[" << j << "] " << particles[j].getCoords()[k]
			        << " - " << anchorCoords[k];
			particles[j].getCoords()[k] = particles[j].getCoords()[k]
					- anchorCoords[k];
			logFile << " = " << particles[j].getCoords()[k] << endl;

			if (maxCoordRange.size() > k) {
				if (maxCoordRange[k][0] > particles[j].getCoords()[k]) // Find smallest coordinate
					maxCoordRange[k][0] = particles[j].getCoords()[k];
				else if (maxCoordRange[k][1] < particles[j].getCoords()[k]) // Find largest coordinate
					maxCoordRange[k][1] = particles[j].getCoords()[k];
			} else
				maxCoordRange.resize(maxCoordRange.size() + 1, dbVector(2, 0));
		}
	}
}

