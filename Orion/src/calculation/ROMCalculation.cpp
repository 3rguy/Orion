/*
 * ROMCalculation.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: rama
 */

#include <ROMCalculation.h>

//ROMCalculation::ROMCalculation(Data& yourData):myData(yourData){};

/*!****************************************************************************/
/*!****************************************************************************/
void ROMCalculation::calculation(Database& myDatabase, DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	// Find the list of results available for calculation
	vector < string > commonResultsNameList = myData.getResultNameList();

	// Find the common steps for calculation
	problemData->setValue("PODIStepValueVec", myData.getStepValueVec());

	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");
	problemData->setValue("calcResultList", dbMatrix());
	intVector& supportDataID = problemData->getIntVector("supportDataID");

	cout << "############# Reduced Order Calculation ############# " << endl;
	logFile << "############# Reduced Order Calculation ############# " << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	if(commonResultsNameList.size() == 0){
		logFile << "In ROMCalculation::calculation, commonResultsNameList.size() = 0"
				<< endl;
		cout << "In ROMCalculation::calculation, commonResultsNameList.size() = 0"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	for (int i = 0; i < commonResultsNameList.size(); i++) {

		cout << "*****************************************************" << endl;
		cout << "Processing result: " << commonResultsNameList[i] << endl;

		logFile << "*****************************************************"
				<< endl;
		logFile << "Processing result: " << commonResultsNameList[i] << endl;

		// Compile list of specific result result
		for (int j = 0; j < supportDataID.size(); j++) {
			resultList.push_back(
					myDatabase.getDataId(supportDataID[j]).getResult(
							commonResultsNameList[i].c_str()));
			myDatabase.getDataId(supportDataID[j]).deleteResult(
					commonResultsNameList[i].c_str());
		}

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

		InputData->setValue("PODMeanCalculation",
				InputData->getValue("dbPODMeanCalculation"));
		InputData->setValue("PODCalculationType",
				InputData->getValue("dbPODCalculationType"));
		InputData->setValue("PODEnergyLevel",
				InputData->getValue("dbPODEnergyLevel"));
		InputData->setValue("PODICalculationType",
				InputData->getValue("dbPODICalculationType"));
		InputData->setValue("PODIPolynomialDegree",
				InputData->getValue("dbPODIPolynomialDegree"));

		problemData->setValue("PODIResultName", commonResultsNameList[i]);
		problemData->setValue("PODIDofPerNode",
				myData.getResultDOF(commonResultsNameList[i].c_str()));

		// ---------------------------------------------------------------------
		// Reduced Order Calculation
		reducedOrderCalculation(problemData, InputData, logFile);

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

	cout << "Reduced Order Calculation completed in: "
			<< elapsed_seconds.count() << " sec" << endl;
	logFile << "Reduced Order Calculation completed in: "
			<< elapsed_seconds.count() << " sec" << endl;

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
void ROMCalculation::reducedOrderCalculation(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	cout << "PODI Calculation started " << endl;
	logFile << "PODI Calculation started " << endl;

//	vector<dbMatrix>& rList = problemData->getDbMatrixVec("resultList");
//	logFile << "Printing the whole of resultList:" << endl;
//	for(int i=0; i<rList.size(); i++){
//		logFile << "rList[" << i << "]: " << endl;
//		printMatrix(rList[i],"rList",logFile);
//	}

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
//!
vector<string> ROMCalculation::findCommonResultName(DataContainer* problemData,
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

