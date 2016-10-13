/*
 * DatabaseQuery.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <DatabaseQuery.h>

DatabaseQuery::DatabaseQuery(DataContainer* problemData, Database& myDatabase,
		InputFileData* InputData, ofstream& logFile) {

	int choice = InputData->getValue("dbInfluenceRangeType");

	switch(choice){
	case 1:
		databaseQueryAlgorithm_one(problemData,myDatabase,InputData,logFile);
		break;

	case 2:
		databaseQueryAlgorithm_two(problemData,myDatabase,InputData,logFile);
		break;

	case 3:
		databaseQueryAlgorithm_three(problemData,myDatabase,InputData,logFile);
		break;

	case 4:
		databaseQueryAlgorithm_four(problemData,myDatabase,InputData,logFile);
		break;

	default:
		logFile << "ERROR: Incorrect DatabaseQuery algorithm selected. Valid options: 1, 2, 3." << endl;
		cout << "ERROR: Incorrect DatabaseQuery algorithm selected. Valid options: 1, 2, 3." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Search for the influencing particle in the database for the interpolation
//! process
void DatabaseQuery::databaseQueryAlgorithm_one(DataContainer* problemData, Database& myDatabase,
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
void DatabaseQuery::databaseQueryAlgorithm_two(DataContainer* problemData, Database& myDatabase,
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
void DatabaseQuery::databaseQueryAlgorithm_three(DataContainer* problemData, Database& myDatabase,
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
//! Search for the influencing particle in the database for the interpolation
//! process
void DatabaseQuery::databaseQueryAlgorithm_four(DataContainer* problemData, Database& myDatabase,
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
				- radiusFactorVec[i]; 	// Min value
		myParamaterInfluenceRange[counter][1] = myParameters[i]
				+ radiusFactorVec[i];	// Max value

		// record the influence radii
		parameterRadii[counter] = radiusFactorVec[i];

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
