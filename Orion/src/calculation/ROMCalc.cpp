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
	problemData->setValue("graphResultList",vector<dbMatrix>());
	problemData->setValue("dataParametersList",dbMatrix());
	problemData->setValue("calcResultList",vector<dbMatrix>());
	problemData->setValue("calcGraphResultList",dbMatrix());
	problemData->setValue("standardSteps",dbVector());

	if (InputData->getValue("standardiseDOF") == 1){
		myData.readTransformedMeshDataFile(InputData,logFile);
		myData.getMeshData()->writeMeshFile("fem_orion_map.msh",InputData,logFile);
	}
	else{
		myData.readMeshDataFile(InputData, logFile);
		myData.getMeshData()->writeMeshFile("fem_orion.msh",InputData,logFile);
	}

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
	searchDatabase_three(problemData, myDatabase, InputData, logFile);

	// Load-up all geometry and displacement data files
	logFile << "Number of datasets selected: " << supportDataID.size() << endl;
	cout << "Number of datasets selected: " << supportDataID.size() << endl;
	dataParametersList.resize(supportDataID.size());

//	test_MLSInterpolantsCalc(problemData, myDatabase, InputData,logFile);

	dbMatrix stepHistoryList(supportDataID.size(), dbVector());

	logFile << "List of Data selected" << endl;
	logFile << "*********************" << endl;
	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "(" << i << "): " << "Data[" << mainData.getId() << "]-> "
				<< mainData.getFolderName() << endl;

		//! The Data read Displacement and Mesh file
		logFile << "Reading mesh file" << endl;
		//mainData.readMeshDataFile(InputData, logFile);

		//! The displacement, parameters and step vector are extracted and
		//! stored
		dataParametersList[i] = mainData.getParamValuesVec();

		stepHistoryList[i] = mainData.getStepValueVec();
		printVector(stepHistoryList[i],"stepHistoryList",logFile);

	}

	problemData->setValue("stepHistoryList",stepHistoryList);

	dbVector& myParameters = problemData->getDbVector("myParameters");
	myData.setParamValuesVec(myParameters);
	//myData.readMeshDataFile(InputData, logFile);
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
//! Search for the influencing particle in the database for the interpolation
//! process
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
	double radiusFactor = InputData->getValue("dbInfluenceRangeFactor");

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

#ifdef _ROMCalcDebugMode_ROMCalculationSet
	logFile << "******* Supporting Data *******" << endl;
	for (int j = 0; j < supportDataID.size(); j++) {
		logFile << j << "). ID = " << supportDataID[j] << endl;
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Search for the influencing particle in the database for the interpolation
//! process
void ROMCalc::searchDatabase_two(DataContainer* problemData, Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	cout << "Searching database" << endl;
	logFile << "Searching database" << endl;

	dbVector& myParameters = problemData->getDbVector("myParameters");
	logFile << "myParameters: " <<myParameters.size() << endl;
	intVector& supportDataID = problemData->getIntVector("supportDataID");
	logFile << "supportDataID: " <<supportDataID.size() << endl;
	dbVector& parameterRadii = problemData->getDbVector("parameterRadii");
	logFile << "parameterRadii: " <<parameterRadii.size() << endl;

	// Read influence factors from input file
	string inputParam_start = "dbInfluenceRangeFactor";
	dbVector radiusFactorVec(myParameters.size(),0);
	for(int i=0; i < myParameters.size(); i ++){

		ostringstream convert;   // stream used for the conversion
		convert << i;
		string inputParam = inputParam_start + convert.str();

		radiusFactorVec[i] = InputData->getValue(inputParam.c_str());

		logFile << "In searchDatabase_two, inputParam = " << inputParam
				<< " -> " << radiusFactorVec[i] << endl;
	}


	// Determine the max and min value of each parameter(i.e Range) and store
	// them in the matrix "myParamaterInfluenceRange"
	int counter = 0;
	dbMatrix myParamaterInfluenceRange(myParameters.size(), dbVector(2));
	parameterRadii.resize(myParameters.size());
	for (int i = 0; i < myParameters.size(); i++) {

		// record the limits of the selected parameters
		myParamaterInfluenceRange[counter][0] = myParameters[i]
				- (myParameters[i] * radiusFactorVec[i]); 	// Min value
		myParamaterInfluenceRange[counter][1] = myParameters[i]
				+ (myParameters[i] * radiusFactorVec[i]);	// Max value

		// record the influence radii
		parameterRadii[counter] = myParameters[i] * radiusFactorVec[i];

		logFile << "Range: " << myParamaterInfluenceRange[counter][0]
				<< " < Parameter[" << i << "]: " << myParameters[i] << " < "
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
			logFile << "Considering: " << endl;
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

		// If all parameters are inside the influence range, record ID of
		// Data
		if (parameterInsideCounter == myParamaterInfluenceRange.size()) {
			supportDataID.resize(supportDataID.size() + 1);
			supportDataID[supportDataID.size() - 1] = dataList[i].getId();
			logFile << "----> Taken" << endl;
		}
	}

	if(supportDataID.size() < 2){
		logFile << "Number of selected data is less than 2.\n"
				"Influence radius is probably too small." << endl;
		cout << "Number of selected data is less than 2.\n"
				"Influence radius is probably too small." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


	double extFactor = InputData->getValue("dbInfluenceRadExtFactor");
	parameterRadii = dbVector(parameterRadii.size(),0);

	for (int i = 0; i < supportDataID.size(); i++) {
		for (int j = 0; j < dataList.size(); j++) {
			if (supportDataID[i] == dataList[j].getId()) {
				dbVector& params = dataList[j].getParamValuesVec();

				oPType dist = 0;
				for (int l = 0; l < params.size(); l++) {
					dist = abs(params[l] - myParameters[l]) * extFactor;

					if(dist > parameterRadii[l]){
						parameterRadii[l] = dist;
					}
				}

			}
		}
	}

	printVector(parameterRadii,"parameterRadii",logFile);

	// List all the supporting data that has been selected
	logFile << "******* Supporting Data *******" << endl;
	for (int j = 0; j < supportDataID.size(); j++) {
		logFile << j << "). ID = " << supportDataID[j] << ", " << myDatabase.getDataId(supportDataID[j]).getFolderName() << endl;
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Search for the influencing particle in the database for the interpolation
//! process
void ROMCalc::searchDatabase_three(DataContainer* problemData, Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	cout << "Searching database" << endl;
	logFile << "Searching database" << endl;

	dbVector& myParameters = problemData->getDbVector("myParameters");
	logFile << "myParameters: " <<myParameters.size() << endl;
	intVector& supportDataID = problemData->getIntVector("supportDataID");
	logFile << "supportDataID: " <<supportDataID.size() << endl;
	dbVector& parameterRadii = problemData->getDbVector("parameterRadii");
	logFile << "parameterRadii: " <<parameterRadii.size() << endl;

	// Read influence factors from input file
	string inputParam_start = "dbInfluenceRangeFactor";
	dbVector radiusFactorVec(myParameters.size(),0);
	for(int i=0; i < myParameters.size(); i ++){

		ostringstream convert;   // stream used for the conversion
		convert << i;
		string inputParam = inputParam_start + convert.str();

		radiusFactorVec[i] = InputData->getValue(inputParam.c_str());

		logFile << "In searchDatabase_two, inputParam = " << inputParam
				<< " -> " << radiusFactorVec[i] << endl;
	}


	// Determine the max and min value of each parameter(i.e Range) and store
	// them in the matrix "myParamaterInfluenceRange"
	int counter = 0;
	dbMatrix myParamaterInfluenceRange(myParameters.size(), dbVector(2));
	parameterRadii.resize(myParameters.size());
	for (int i = 0; i < myParameters.size(); i++) {

		// record the limits of the selected parameters
		myParamaterInfluenceRange[counter][0] = myParameters[i]
				- (myParameters[i] * radiusFactorVec[i]); 	// Min value
		myParamaterInfluenceRange[counter][1] = myParameters[i]
				+ (myParameters[i] * radiusFactorVec[i]);	// Max value

		// record the influence radii
		parameterRadii[counter] = myParameters[i] * radiusFactorVec[i];

		logFile << "Range: " << myParamaterInfluenceRange[counter][0]
				<< " < Parameter[" << i << "]: " << myParameters[i] << " < "
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
			logFile << "Considering: " << endl;
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

		// If all parameters are inside the influence range, record ID of
		// Data
		if (parameterInsideCounter == myParamaterInfluenceRange.size()) {
			supportDataID.resize(supportDataID.size() + 1);
			supportDataID[supportDataID.size() - 1] = dataList[i].getId();
			logFile << "----> Taken" << endl;
		}
	}

	if(supportDataID.size() < 2){
		logFile << "Number of selected data is less than 2.\n"
				"Influence radius is probably too small." << endl;
		cout << "Number of selected data is less than 2.\n"
				"Influence radius is probably too small." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


	//
	oPType maxRad = 0;
	int maxRadId;
	for (int i = 0; i < supportDataID.size(); i++) {
		for (int j = 0; j < dataList.size(); j++) {
			if (supportDataID[i] == dataList[j].getId()) {
				dbVector& params = dataList[j].getParamValuesVec();

				oPType dist = 0;
				for (int l = 0; l < params.size(); l++) {
					dist += pow(params[l] - myParameters[l],2);
				}

				if (maxRad < dist) {
					maxRad = dist;
					maxRadId = i;
				}

			}
		}
	}

	double extFactor = InputData->getValue("dbInfluenceRadExtFactor");

	dbVector& params = dataList[supportDataID[maxRadId]].getParamValuesVec();
	for (int m = 0; m < params.size(); m++) {
		parameterRadii[m] = abs(params[m] - myParameters[m]) * extFactor;
	}

	printVector(parameterRadii,"parameterRadii",logFile);

	// List all the supporting data that has been selected
	logFile << "******* Supporting Data *******" << endl;
	for (int j = 0; j < supportDataID.size(); j++) {
		logFile << j << "). ID = " << supportDataID[j] << ", " << myDatabase.getDataId(supportDataID[j]).getFolderName() << endl;
	}

}


/*!****************************************************************************/
/*!****************************************************************************/
//! Test MLS interpolants
void ROMCalc::test_MLSInterpolantsCalc(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile){

	intVector& supportDataID = problemData->getIntVector("supportDataID");
	dbMatrix& dataParametersList = problemData->getDbMatrix("dataParametersList");

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

	printVector(interpolants,"****** Interpolants ******",logFile);
		oPType intpolantSum = 0;
		for(int i=0; i<interpolants.size();i++) intpolantSum += interpolants[i];
		logFile << "Sum of interpolants: " << intpolantSum << endl;
		cout << "Sum of interpolants: " << intpolantSum << endl;

	cout << "The end !" << endl;

	MPI_Abort(MPI_COMM_WORLD, 1);
}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void ROMCalc::ROMCalculationSet(DataContainer* problemData,Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	intVector& supportDataID = problemData->getIntVector("supportDataID");

	// Find the list of results available for calculation
	vector<string> commonResultsNameList =
			findCommonResultName(problemData,myDatabase,InputData,logFile);
	myData.setResultNameList(commonResultsNameList);

	// Find the common steps for calculation
	dbVector& standardSteps = problemData->getDbVector("standardSteps");
	standardSteps = findCommonStepsList(problemData,myDatabase,InputData,logFile);
	myData.setStepValueVec(standardSteps);
	problemData->setValue("PODIStepValueVec",standardSteps);


	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");
	problemData->setValue("calcResultList",dbMatrix());

	cout << "############# Reduced Order Calculation ############# " << endl;
	logFile << "############# Reduced Order Calculation ############# " << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	for (int i = 0; i < commonResultsNameList.size(); i++) {


		myData.setResultDOF(commonResultsNameList[i].c_str(),
				myDatabase.getDataId(supportDataID[0]).getResultDOF(
						commonResultsNameList[i].c_str()));

		cout << "*****************************************************" << endl;
		cout << "Processing result: " << commonResultsNameList[i] << endl;

		logFile << "*****************************************************" << endl;
		logFile << "Processing result: " << commonResultsNameList[i] << endl;

		// Compile list of specific result result
		for (int j = 0; j < supportDataID.size(); j++) {
			resultList.push_back(
					myDatabase.getDataId(supportDataID[j]).getResult(
							commonResultsNameList[i].c_str()));
			myDatabase.getDataId(supportDataID[j]).deleteResult(
					commonResultsNameList[i].c_str());
		}

		stepStandardisation(problemData, InputData, logFile);

		// ---------------------------------------------------------------------
		// Loading input parameters
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

		problemData->setValue("PODIResultName",commonResultsNameList[i]);
		problemData->setValue("PODIDofPerNode",
				myData.getResultDOF(commonResultsNameList[i].c_str()));

		// ---------------------------------------------------------------------
		// Reduced Order Calculation
		reducedOrderMethodCalc(problemData, InputData, logFile);

		// ---------------------------------------------------------------------
		// Record results
		myData.setResult(commonResultsNameList[i].c_str(),
				problemData->getDbMatrix("calcResultList"));
		dbMatrix().swap(problemData->getDbMatrix("calcResultList"));

		vector<dbMatrix>().swap(problemData->getDbMatrixVec("resultList"));

		problemData->deleteString("PODIResultName");
		problemData->deleteInt("PODIDofPerNode");

	}

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Reduced Order Calculation completed in: " << elapsed_seconds.count()
			<< " sec" << endl;
	logFile << "Reduced Order Calculation completed in: " << elapsed_seconds.count()
			<< " sec" << endl;


	// -------------------------------------------------------------------------
	// Graph result processing

//	vector<dbMatrix>& graphResultList=
//			problemData->getDbMatrixVec("graphResultList");
//
//	int numGraphResults = graphResultList[0].size();
//	int numSteps = graphResultList[0][0].size();
//
//	//Transpose each matrix
//	vector<dbMatrix> transGraphResultList(graphResultList.size(),dbMatrix());
//	for(int i = 0 ; i <graphResultList.size(); i++ ){
//		dbMatrix tempMat(numSteps,dbVector(numGraphResults,0));
//		for(int j = 0 ; j < numGraphResults; j++){
//			for(int k = 0 ; k < numSteps; k++){
//				tempMat[k][j] = graphResultList[i][j][k];
//			}
//		}
//		transGraphResultList[i] = tempMat;
//	}
//
//	string PODIResultName = "graph";
//	problemData->setValue("PODIResultName",PODIResultName);
//
//	problemData->deleteDbVector("PODIStepValueVec");
//	dbVector graphStepVec(numSteps,0);
//	for(int i=0; i<numSteps ; i++){
//		graphStepVec[i] = i;
//	}
//	problemData->setValue("PODIStepValueVec",graphStepVec);
//
//	int graphDof = 1;
//	problemData->setValue("PODIDofPerNode",graphDof);
//
//	dbMatrix calcGraphResultListMatrix;
//	PODICalc* PodiGraph = new PODICalc(problemData->getDbVector("myParameters"),
//				problemData->getDbVector("parameterRadii"),
//				problemData->getIntVector("supportDataID"),
//				problemData->getDbMatrix("dataParametersList"),
//				transGraphResultList,calcGraphResultListMatrix, problemData,
//				InputData, logFile);
//
//	problemData->getDbMatrix("calcGraphResultList").
//			resize(calcGraphResultListMatrix[0].size(),
//					dbVector(calcGraphResultListMatrix.size()));
//	for(int i = 0; i < calcGraphResultListMatrix.size(); i++){
//		for(int j = 0; j < calcGraphResultListMatrix[0].size(); j++){
//			problemData->getDbMatrix("calcGraphResultList")[j][i]
//			                                 = calcGraphResultListMatrix[i][j];
//		}
//	}
//
//
//	delete PodiGraph;

//	problemData->deleteString("PODIResultName");
//	problemData->deleteInt("PODIDofPerNode");
	problemData->deleteDbVector("PODIStepValueVec");


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
			problemData, InputData, logFile);

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
void ROMCalc::postProcessing(DataContainer* problemData, Database& myDatabase,
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

	if (InputData->getValue("standardiseDOF") == 1) {

		FEMGeometryExt* gridFEMGeo = myGrid->getGridGeometry(logFile);

		string oldFolderName = myData.getFolderName();
		myData.getFolderName() = string("postProc");
		vector<Data*> dataList(1);
		dataList.at(0) = &myData;

		if(InputData->getValue("gridNodesResultPlot") == 1)
			saveResultInGridNodesFormat(problemData, dataList, InputData, logFile);

		myData.getFolderName() = oldFolderName;

		// -------------------------------------------------------------------------
		// Loading input parameters for the interpolation calculation
		InputData->setValue("interpolantionType",
				InputData->getValue("gInterpolantionType"));
		InputData->setValue("influenceRangeFactor",
				InputData->getValue("gInfluenceRangeFactor"));
		InputData->setValue("MLSCalculationType",
				InputData->getValue("gMLSCalculationType"));
		InputData->setValue("MLSPolynomialDegree",
				InputData->getValue("gMLSPolynomialDegree"));
		InputData->setValue("MLSWeightFunc",
				InputData->getValue("gMLSWeightFunc"));
		InputData->setValue("parameterPolynomialDegree",
				InputData->getValue("gparameterPolynomialDegree"));
		myGrid->initCalcResultOnParticles(myData, InputData, logFile);


		myData.setMeshData(gridFEMGeo);

		vector<string>& resultNameList = myData.getResultNameList();
		for (int i = 0; i < resultNameList.size(); i++) {

			myData.assignResultToParticles(resultNameList[i].c_str(), InputData,
					logFile);

			dbMatrix resultMat = myGrid->interpolateResultOnGridPoint(myData,
					InputData, logFile);

			myGrid->resetNodesStepDOFMat();

			myData.deleteResult(resultNameList[i].c_str());
			myData.setResult(resultNameList[i].c_str(), resultMat);

		}

		myData.setMeshData(myGrid->getGridGeometry(logFile));
		myGrid->setGridGeometry(gridFEMGeo);

		myData.saveResultsToFile(logFile);
		system("cp -f fem_orion.res fem_orion_map.res");
//
		myData.delMeshData();
		myData.readMeshDataFile(InputData, logFile);
		myData.getMeshData()->writeMeshFile("fem_orion.msh",InputData,logFile);
	}

//	myData.setGraphResultList(problemData->getDbMatrix("calcGraphResultList"));
//	myData.saveGraphResultsToFile(logFile);
	myData.saveResultsToFile(logFile);
	myData.calcLeftAndRightCavityVolumes(InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! findSupports: Finding the supports of a particular point
void ROMCalc::findSupports(dbVector& iPoint, dbMatrix& pointList,
		intVector& supportsList, dbVector& pointRadii, InputFileData* InputData,
		ofstream& logFile) {

//	double radiusFactor = InputData->getValue("influenceRangeFactor");
	double radiusFactor = 1.5;		// Temporary solution -> need to find a better one

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
	vector<dbMatrix>& graphResultList = problemData->getDbMatrixVec("graphResultList");
	for (int i = 0; i < supportDataID.size(); i++) {

		Data& mainData = myDatabase.getDataId(supportDataID[i]);

		logFile << "-----------------------------" << endl;
		logFile << "[" << i <<"]Loading Data file: " << mainData.getFolderName() << endl;
		logFile << "-----------------------------" << endl;

		cout << "-----------------------------" << endl;
		cout << "[" << i <<"]Loading Data file: " << mainData.getFolderName() << endl;
		cout << "-----------------------------" << endl;

		//! The Data read Displacement and Mesh file
		logFile << "Reading mesh file" << endl;
		mainData.readMeshDataFile(InputData, logFile);

		// Read result from file
		mainData.readResultFile(InputData, logFile);

		// Read result from graph file
		mainData.readGraphResultFile(InputData,logFile);
		graphResultList.push_back(mainData.getGraphResultList());

		// Re-calculate the cavity volumes and write to graph file
		mainData.calcLeftAndRightCavityVolumes(InputData, logFile);

		//! Free memory
		mainData.delMeshData();

		// Extract step-values for step standardisation
		stepHistoryList[i] = mainData.getStepValueVec();

		// Standardise the degrees of Freedom
		if (InputData->getValue("standardiseDOF") == 1){

			logFile << "Starting DOF standardisation of Data "
					<< supportDataID[i] << endl;

			// Read the mapped goemetry(used by point-in-polygon algorithm)
			mainData.readTransformedMeshDataFile(InputData, logFile);

			// Output results in map nodes format
			string mapResFileName = mainData.getFolderName() + "fem_map.res";
			mainData.saveAllResultsToFile_res_format(mapResFileName.c_str(),logFile);

			string mapMshFileName = mainData.getFolderName() + "fem_map.msh";
			mainData.getMeshData()->writeMeshFile(mapMshFileName.c_str(),InputData,logFile);

			// -----------------------------------------------------------------
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();

			// Standardise all results
			standardiseResultDOF(problemData,mainData,InputData,logFile);

			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;

			cout << "DOF standardisation of Data completed in: "
				 << elapsed_seconds.count() << " sec" << endl;
			logFile << "DOF standardisation of Data completed in: "
					<< elapsed_seconds.count() << " sec" << endl;
			// -----------------------------------------------------------------

			// Output results in grid nodes format
			string gridResFileName = mainData.getFolderName() + "fem_grid.res";
			mainData.saveAllResultsToFile_res_format(gridResFileName.c_str(),logFile);

			string gridMshFileName = mainData.getFolderName() + "fem_grid.msh";
			myGrid->getGridGeometry(logFile)->
					writeMeshFile(gridMshFileName.c_str(),InputData,logFile);

		}

//		cout << "CONTROLLED ENDING" << endl;
//		logFile << "CONTROLLED ENDING" << endl;
//		MPI_Abort(MPI_COMM_WORLD, 1);

	}

	cout << endl;

	// Save results in GridNodes format (superseeded by fem_grid outputs)
//	if (InputData->getValue("standardiseDOF") == 1
//			&& InputData->getValue("gridNodesResultPlot") == 1) {
//
//		vector<Data*> dataList(supportDataID.size());
//		for (int i = 0; i < supportDataID.size(); i++) {
//			dataList.at(i) = &myDatabase.getDataId(supportDataID[i]);
//		}
//
//		saveResultInGridNodesFormat(problemData, dataList, InputData, logFile);
//	}

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

	vector<dbMatrix>& graphResultList = problemData->getDbMatrixVec("graphResultList");

	//! TODO: (1)Implement method to skip the following code in case the steps are constants
	//! TODO: (2)Check that interpolated displacement are not used for other interpolation calculation (??)
	int totalNumDofs = resultList[0].size();
	int numGraphResult = graphResultList[0].size();

	for (int i = 0; i < supportDataID.size(); i++) {
		dbMatrix tempDispMatrix(totalNumDofs, dbVector(standardSteps.size()));
		dbMatrix tempGraphMatrix(numGraphResult,dbVector(standardSteps.size(),0));

		for (int j = 0; j < standardSteps.size(); j++) {


//			int position = findDoubleVecPos(standardSteps[j], 0,
//					stepHistoryList[i].size(), stepHistoryList[i]);

			int position = -1;
			for(int k = 0 ; k < stepHistoryList[i].size(); k++){
				if(fabs(stepHistoryList[i][k]-standardSteps[j]) < 1e-10)
					position = k;
			}


			if (position != -1) {
				for (int k = 0; k < totalNumDofs; k++)
					tempDispMatrix[k][j] = resultList[i][k][position];

				for(int k = 0; k < numGraphResult; k++)
					tempGraphMatrix[k][j] = graphResultList[i][k][position];

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

				// Loading input parameters
				InputData->setValue("interpolantionType",
						InputData->getValue("dbInterpolantionType"));
				InputData->setValue("influenceRangeFactor",1);
				InputData->setValue("MLSCalculationType",
						InputData->getValue("dbMLSCalculationType"));
				InputData->setValue("MLSPolynomialDegree",
						InputData->getValue("dbMLSPolynomialDegree"));
				InputData->setValue("MLSWeightFunc",
						InputData->getValue("dbMLSWeightFunc"));
				InputData->setValue("parameterPolynomialDegree",
						InputData->getValue("dbparameterPolynomialDegree"));

				PODICalc* Podi = new PODICalc(iPoint, pointRadii, supportsList,
						supportPointList, tempDispList, tempResultingDispMatrix,
						problemData, InputData, logFile);

				// -------------------------------------------------------------
				dbVector interpolants = Podi->getPODInterpolants();

				for(int k = 0; k < numGraphResult; k++){
					for(int m = 0; m < supportsList.size(); m++){
						tempGraphMatrix[k][j] +=
								interpolants[m] * graphResultList[i][k][supportsList[m]];
					}
				}

				// -------------------------------------------------------------

				for (int l = 0; l < totalNumDofs; l++)
					tempDispMatrix[l][j] = tempResultingDispMatrix[l][0];

				delete Podi;

			}
		}

		resultList[i] = tempDispMatrix;
		graphResultList[i] = tempGraphMatrix;
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


	// Preliminary stuff
	// -----------------
	if (InputData->getValue("automaticGridNodeSetup") == 1) {

		//! --------------------------------------------------------------------
		//! For each data set: Assign displacement field to every particles and
		//! reset their coordinate with respect to the defined anchor point

		cout << "Setting up automatic benchmark grid" << endl;
		logFile << "Setting up automatic benchmark grid" << endl;

		intVector& supportDataID = problemData->getIntVector("supportDataID");
		dbMatrix maxCoordRange;

		for (int i = 0; i < supportDataID.size(); i++) {

			// reset particle's coordinates
			logFile << "Coordinate setup of Data: " << supportDataID[i] << endl;
			coordinateSetup(myDatabase.getDataId(supportDataID[i]),
					maxCoordRange, InputData, logFile);
			myDatabase.getDataId(supportDataID[i]).getMeshData()->writeMeshFile(
					"fem2.res",InputData, logFile);

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
	else{

		//! --------------------------------------------------------------------
		//! Setup the reference grid and define the nodes
		cout << "Loading up benchmark grid" << endl;
		logFile << "Loading up benchmark grid" << endl;

		myGrid = new GridNodes(InputData, logFile);
	}


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

	// -------------------------------------------------------------------------
	// Loading input parameters for the interpolation calculation
	InputData->setValue("interpolantionType",
			InputData->getValue("gInterpolantionType"));
	InputData->setValue("influenceRangeFactor",
			InputData->getValue("gInfluenceRangeFactor"));
	InputData->setValue("MLSCalculationType",
			InputData->getValue("gMLSCalculationType"));
	InputData->setValue("MLSPolynomialDegree",
			InputData->getValue("gMLSPolynomialDegree"));
	InputData->setValue("MLSWeightFunc",
			InputData->getValue("gMLSWeightFunc"));
	InputData->setValue("parameterPolynomialDegree",
			InputData->getValue("gparameterPolynomialDegree"));

	// Setting up the interpolants
	myGrid->interpolantSetup(mainData, InputData, logFile);
	// -------------------------------------------------------------------------

	// Interpolating all results + deleting previous ones + saving the new ones
	vector<string>& resultNameList = mainData.getResultNameList();
	for (int i = 0; i < resultNameList.size(); i++) {

		// assign displacement field to particles
		mainData.assignResultToParticles(resultNameList[i].c_str(), InputData,
				logFile);

		// Interpolate result on gridNodes
		dbMatrix standardResult = myGrid->interpolateResultOnGridPoint(mainData,
				InputData, logFile);

		// Reset all dofs in myGrid
		myGrid->resetNodesStepDOFMat();

		// Delete previous results and save the new one
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

			particles[j].getCoords()[k] = particles[j].getCoords()[k]
								- anchorCoords[k];

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

/*!****************************************************************************/
/*!****************************************************************************/
void ROMCalc::saveResultInGridNodesFormat(DataContainer* problemData,
		vector<Data*>& dataList, InputFileData* InputData, ofstream& logFile){

	using namespace std;

	Data* gridData = new Data();
	system("cd gridNodes");
	//gridData->setFileName(string("gridNodes/"));
	int old_choice = (int) InputData->getValue("FEMGeometrySetupType");
	InputData->setValue("FEMGeometrySetupType",1);
	gridData->setFolderName(string("gridNodes/"));
	gridData->readMeshDataFile(InputData,logFile);
	InputData->setValue("FEMGeometrySetupType",old_choice);

	string filename;

	for(int i = 0; i < dataList.size(); i++){

		gridData->getResultNameList() = dataList.at(i)->getResultNameList();
		gridData->getResultDOFList() = dataList.at(i)->getResultDOFList();
		gridData->getResultList() = dataList.at(i)->getResultList();
		gridData->getStepValueVec() = dataList.at(i)->getStepValueVec();

		filename = "gridNode_" + dataList.at(i)->getFolderName();

		replace(filename.begin(),filename.end(),'/','_');

		string resfile = filename + ".res";
		gridData->saveAllResultsToFile_res_format(resfile.c_str(),logFile);

		string mshfile = filename + ".msh";
		string command = "cp fem.msh " + mshfile;
		system(command.c_str());

		gridData->getResultNameList().clear();
		gridData->getResultDOFList().clear();
		gridData->getResultList().clear();
		gridData->getStepValueVec().clear();
	}

	system("cd ..");
}
