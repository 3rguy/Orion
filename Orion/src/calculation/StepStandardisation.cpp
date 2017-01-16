/*
 * StepStandardisation.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <StepStandardisation.h>

/*!****************************************************************************/
/*!****************************************************************************/
//!
StepStandardisation::StepStandardisation(DataContainer* problemData,
		InputFileData* InputData,ofstream& logFile){

	nPhase = 0;
	defineNormaliseSteps(problemData,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void StepStandardisation::defineNormaliseSteps(DataContainer* problemData,
		InputFileData* InputData,ofstream& logFile){

	double standStepSpacing = InputData->getValue("standardStepSpacing");
	int nStandStep = floor(1/standStepSpacing) + 1;
	normStandStep = dbVector (nStandStep);
	for(int i=0; i<nStandStep-1; i++){
		normStandStep[i] = 0 + (i*standStepSpacing);
	}
	normStandStep[nStandStep-1] = 1;

}

/*!****************************************************************************/
/*!****************************************************************************/
void StepStandardisation::standardise(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	vector<dbMatrix>& resultList = problemData->getDbMatrixVec("resultList");
	int nData = resultList.size();
	int nSteps = problemData->getDbVector("PODIStepValueVec").size();

	for (int i = 0; i < nData; i++) {

		logFile << "***** nData[" << i << "] *****" << endl;

		dbMatrix& result = resultList[i];
		dbMatrix newResult(result.size(),dbVector(nSteps,0));

		for (int j = 0; j < nPhase; j++) {

			logFile << "***** nPhase[" << j << "] *****" << endl;

			// *****************************************************************
			// Get All info related to the timesteps
			// *****************************************************************
			// Extract the standard steps
			string ss = "standardStepPhaseList" + std::to_string(j);
			dbVector& standardStepPhase = problemData->getDbMatrix(ss.c_str())[i];

			problemData->setValue("standardStepPhase", standardStepPhase);
			printVector(standardStepPhase, "standardStepPhase", logFile);

			// Extract the non-standard steps
			string ns = "nonStandardStepPhaseList" + std::to_string(j);
			dbVector& nonStandardStepPhase = problemData->getDbMatrix(
					ns.c_str())[i];
			problemData->setValue("nonStandardStepPhase", nonStandardStepPhase);


			// *****************************************************************
			// Get all info related to the result fields
			// *****************************************************************
			// Extract the non-standard results fields
			dbMatrix nonStandardResultPhase(result.size(), dbVector());
			for (int k = 0; k < result.size(); k++) {
				for (int l = phaseIndex[i][j]; l < phaseIndex[i][j + 1] + 1;
						l++) {
					nonStandardResultPhase[k].push_back(result[k][l]);
				}
			}
			problemData->setValue("nonStandardResultPhase",
					nonStandardResultPhase);

			// Standardise the results
			stepStandardisationCalc(problemData, InputData, logFile);
			dbMatrix standardresultPhase = problemData->getDbMatrix("nonStandardResultPhase");

//			if(i==2 && j==2){
//				cout << "Controlled ending" << endl;
//				logFile << "Controlled ending" << endl;
//				MPI_Abort(MPI_COMM_WORLD, 1);
//			}

			if(j==0){
				newResult = standardresultPhase;
			}
			else{
//				cout << "building newResultMatrix" << endl;
//				logFile << "building newResultMatrix" << endl;
				for (int k = 0; k < standardresultPhase.size(); k++) {
					for (int l = 1; l < standardresultPhase[k].size(); l++) { // Intentionally skip first entry since it has already been recorded from the previous phase
					newResult[k].push_back(standardresultPhase[k][l]);
					}
				}

			}

		}

		resultList[i].swap(newResult);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void StepStandardisation::stepStandardisationCalc(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {
	// Define a list of standard steps

	dbVector& standardStepPhase = problemData->getDbVector("standardStepPhase");
	dbVector& nonStandardStepPhase = problemData->getDbVector("nonStandardStepPhase");
	dbMatrix& nonStandardResultPhase = problemData->getDbMatrix("nonStandardResultPhase");

	InputData->setValue("dispMatrixRearrange", 0);

	int totalNumDofs = nonStandardResultPhase.size();

#ifdef _StepStandardisationDebugMode_
	printMatrix(nonStandardResultPhase,"nonStandardResultPhase",logFile);
#endif

	dbMatrix standardResultPhase(totalNumDofs, dbVector(standardStepPhase.size()));

	printVector(standardStepPhase,"standardStepPhase",logFile);
	InputData->setValue("printPodiInfo",0);

	for (int j = 0; j < standardStepPhase.size(); j++) {

		int position = -1;
		for (int k = 0; k < nonStandardStepPhase.size(); k++) {
			if (fabs(nonStandardStepPhase[k] - standardStepPhase[j]) < 1e-10)
				position = k;
		}

		// If the actual step is found among the standardsteps,
		if (position != -1) {
			for (int k = 0; k < totalNumDofs; k++)
				standardResultPhase[k][j] = nonStandardResultPhase[k][position];

		} else {
			dbVector iPoint;
			iPoint.push_back(standardStepPhase[j]);

			dbMatrix pointList(nonStandardStepPhase.size(), dbVector(1));
			for (int k = 0; k < nonStandardStepPhase.size(); k++) {
				pointList[k][0] = nonStandardStepPhase[k];
			}

			dbVector pointRadii;
			intVector supportsList;

			findSupports_two(iPoint, pointList, supportsList, pointRadii, InputData,
					logFile);

			dbMatrix supportPointList(supportsList.size(), dbVector(1));
			for (int i = 0; i < supportsList.size(); i++) {
				supportPointList[i] = pointList[supportsList[i]];
			}

			vector<dbMatrix> tempDispList(1,
					dbMatrix(nonStandardResultPhase.size(),
							dbVector(supportsList.size(), 0)));

			for (int k = 0; k < supportsList.size(); k++) {
				for (int l = 0; l < nonStandardResultPhase.size(); l++) {
					tempDispList[0][l][k] = nonStandardResultPhase[l][supportsList[k]];
				}
			}

#ifdef _StepStandardisationDebugMode_
			printVector(iPoint, "iPoint:", logFile);
			printMatrix(supportPointList, "supportPointList:", logFile);

			for(int p=0; p<tempDispList.size();p++){
				logFile << "----------- p[" << p << "] -----------" << endl;
				printMatrix(tempDispList[p],"tempDispList",logFile);
			}

			printVector(iPoint, "iPoint:", logFile);
			printVector(pointRadii, "pointRadii:", logFile);
			printVector(supportsList, "supportsList:", logFile);
			printMatrix(supportPointList, "supportPointList:", logFile);
			printMatrix(pointList, "pointList:", logFile);
#endif

			// Loading input parameters
			InputData->setValue("interpolantionType",
					InputData->getValue("dbInterpolantionType"));
			InputData->setValue("influenceRangeFactor", 1);
			InputData->setValue("MLSCalculationType",
					InputData->getValue("dbMLSCalculationType"));
			InputData->setValue("MLSPolynomialDegree",
					InputData->getValue("dbMLSPolynomialDegree"));
			InputData->setValue("MLSWeightFunc",
					InputData->getValue("dbMLSWeightFunc"));
			InputData->setValue("parameterPolynomialDegree",
					InputData->getValue("standStepstepPODIPolynomialDegree"));

			InputData->setValue("PODMeanCalculation",
					InputData->getValue("standStepPODMeanCalculation"));
			InputData->setValue("PODCalculationType",
					InputData->getValue("standStepstepPODCalculationType"));
			InputData->setValue("PODEnergyLevel",
					InputData->getValue("standStepstepPODEnergyLevel"));
			InputData->setValue("PODICalculationType",
					InputData->getValue("standStepstepPODICalculationType"));
			InputData->setValue("PODIPolynomialDegree",
					InputData->getValue("standStepstepPODIPolynomialDegree"));

			double plotPOVandPOM = InputData->getValue("PlotPOMsAndPOVs");
			InputData->setValue("PlotPOMsAndPOVs",0);

			dbMatrix tempResultingDispMatrix;

			PODICalc* Podi = new PODICalc(iPoint, pointRadii, supportsList,
					supportPointList, tempDispList, tempResultingDispMatrix,
					problemData, InputData, logFile);

			InputData->setValue("PlotPOMsAndPOVs",plotPOVandPOM);

			for (int l = 0; l < totalNumDofs; l++)
				standardResultPhase[l][j] = tempResultingDispMatrix[l][0];

			delete Podi;

#ifdef _StepStandardisationDebugMode_
			printMatrix(tempResultingDispMatrix,"tempResultingDispMatrix",logFile);
			printMatrix(standardResultPhase,"standardResultPhase",logFile);
#endif

		}
	}

#ifdef _StepStandardisationDebugMode_
	printMatrix(standardResultPhase,"standardResultPhase(final)",logFile);
#endif


	nonStandardResultPhase = standardResultPhase;

	InputData->setValue("printPodiInfo",1);
	InputData->setValue("dispMatrixRearrange", 1);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! findSupports: Finding the supports of a particular point
void StepStandardisation::initStepStandardisation(DataContainer* problemData,
		InputFileData* InputData,ofstream& logFile){

	dbMatrix standardStepsList;

	if(InputData->getValue("modelType") != 2)
		findCommonStepsList(problemData,InputData,logFile);
	else
		findCommonCardiacStepsList(problemData,InputData,logFile);
}

///*!****************************************************************************/
///*!****************************************************************************/
////!
//dbMatrix StepStandardisation::findCommonStepsList(DataContainer* problemData,
//		InputFileData* InputData, ofstream& logFile){
//
//	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");
//	dbVector standardSteps;
//
//	int position = 0;
//	for (int i = 0; i < stepHistoryList.size(); i++) {
//		for (int j = 0; j < stepHistoryList[i].size(); j++) {
//			position = findDoubleVecPos(stepHistoryList[i][j], 0,
//					standardSteps.size(), standardSteps);
//			if (position == -1) {
//				standardSteps.push_back(stepHistoryList[i][j]);
//			}
//		}
//	}
//#ifdef _StepStandardisationDebugMode_
//	logFile << "******* Unsorted Standard Steps *******" << endl;
//
//	logFile << "Steps are: ";
//	for (int j = 0; j < standardSteps.size(); j++) {
//		logFile << standardSteps[j] << ", ";
//	}
//	logFile << endl;
//#endif
//
//	sortDoubleVector(standardSteps, 0, standardSteps.size() - 1);
//
//#ifdef _StepStandardisationDebugMode_
//	logFile << "******* Sorted Standard Steps *******" << endl;
//
//	logFile << "Steps are: ";
//	for (int j = 0; j < standardSteps.size(); j++) {
//		logFile << standardSteps[j] << ", ";
//	}
//	logFile << endl;
//#endif
//
//	int nStepHistoryList = stepHistoryList.size();
//	dbMatrix standardStepsList(nStepHistoryList,standardSteps);
//
//	return standardStepsList;
//}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void StepStandardisation::findCommonStepsList(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){



	int nData = problemData->getIntVector("supportDataID").size();

	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");


	phaseIndex = intMatrix(nData,intVector());

	for(int i = 0; i < nData; i++){

		phaseIndex[i].push_back(0);
		phaseIndex[i].push_back(stepHistoryList[i].size()-1);

		printMatrix(phaseIndex,"phaseIndex",logFile);
	}


	// Retain only the common cardiac phases (simple implementation)
	int minPhaseIndex = INT_MAX;
	for(int i = 0 ; i < phaseIndex.size(); i++){
		if(minPhaseIndex > phaseIndex[i].size())
			minPhaseIndex = phaseIndex[i].size();
	}
	if(minPhaseIndex <= 0){
		cout << "ERROR: In StepStandardisation::findCommonCardiacStepsList, the"
				" minimum of cardiac phases recorded is zero." << endl;
		logFile << "ERROR: In StepStandardisation::findCommonCardiacStepsList, the"
				" minimum of cardiac phases recorded is zero." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else{
		for(int i = 0 ; i < phaseIndex.size(); i++){
			phaseIndex[i].resize(minPhaseIndex);
		}
	}

	nPhase = minPhaseIndex - 1;

	vector<dbMatrix> standardStepPhaseListMat, nonStandardStepPhaseListMat;
	createPhaseList(stepHistoryList,standardStepPhaseListMat,
			nonStandardStepPhaseListMat, problemData, InputData, logFile);

	for(int i = 0; i < nPhase; i++){
		string ss = "standardStepPhaseList" + std::to_string(i);
		problemData->setValue(ss.c_str(),standardStepPhaseListMat[i]);

		string ns = "nonStandardStepPhaseList" + std::to_string(i);
		problemData->setValue(ns.c_str(),nonStandardStepPhaseListMat[i]);

		printMatrix(standardStepPhaseListMat[i],"standardStepPhaseList",logFile);
		printMatrix(nonStandardStepPhaseListMat[i],"nonStandardStepPhaseList",logFile);
	}

}


/*!****************************************************************************/
/*!****************************************************************************/
//! findSupports: Finding the supports of a particular point
void StepStandardisation::findSupports(dbVector& iPoint, dbMatrix& pointList,
		intVector& supportsList, dbVector& pointRadii, InputFileData* InputData,
		ofstream& logFile) {

	int findSupport_choice = 0;

	double dist_sum = 0, sqr_dist = 0;
	dbVector distVec(pointList.size());
	for(int i=0; i<pointList.size(); i++){
		sqr_dist = 0;
		for(int j=0; j<pointList[i].size(); j++){
			sqr_dist += pow(pointList[i][j]-iPoint[j],2);
		}
		distVec[i] = sqrt(sqr_dist);
	}

	printVector(distVec,"distVec (before sorting)",logFile);
	std::sort(distVec.begin(),distVec.end());
	printVector(distVec,"distVec (after sorting)",logFile);

//	double radiusFactor = 1.1;
	double radiusFactor = InputData->getValue("stepStandRadiusFactor");

	int counter = 0;
	dbMatrix pointInfluenceRange(iPoint.size(), dbVector(2));
	pointRadii.resize(iPoint.size());
	for (int i = 0; i < iPoint.size(); i++) {
		pointInfluenceRange[counter][0] = iPoint[i]
				- (distVec[1] * radiusFactor); 	// Min value
		pointInfluenceRange[counter][1] = iPoint[i]
				+ (distVec[1] * radiusFactor);	// Max value
		pointRadii[counter] = distVec[1] * radiusFactor;
		counter++;
	}

	printMatrix(pointInfluenceRange,"pointInfluenceRange: min(left column) max(right column)",logFile);
	printVector(pointRadii,"pointRadii",logFile);

	// =========================================================================

	// Temporary solution -> need to find a better one
//	double radiusFactor = InputData->getValue("stepStandardRadiusFactor");
//
//	// Determine the max and min value of each parameter(i.e Range) and store
//	// them in the matrix "myParamaterInfluenceRange"
//	int counter = 0;
//	dbMatrix pointInfluenceRange(iPoint.size(), dbVector(2));
//	pointRadii.resize(iPoint.size());
//	for (int i = 0; i < iPoint.size(); i++) {
//		pointInfluenceRange[counter][0] = iPoint[i]
//				- (iPoint[i] * radiusFactor); 	// Min value
//		pointInfluenceRange[counter][1] = iPoint[i]
//				+ (iPoint[i] * radiusFactor);	// Max value
//		pointRadii[counter] = iPoint[i] * radiusFactor;
//		counter++;
//	}
	// =========================================================================

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

//			logFile << pointInfluenceRange[j][0] << " < " << pointVec[j]
//					<< " < " << pointInfluenceRange[j][1] << endl;

			if (pointInfluenceRange[j][0] <= pointVec[j]
					&& pointVec[j] <= pointInfluenceRange[j][1]) {
				pointInsideCounter++;
//				logFile << "Taken" << endl << endl;

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

#ifdef _StepStandardisationDebugMode_
	logFile << "******* Supporting Points *******" << endl;

	logFile << "Main point: ";
	for (int i = 0; i < iPoint.size(); i++) {
		logFile << iPoint[i] << ",";
	}
	logFile << endl;

	printMatrix(pointList,"pointList",logFile);

	printVector(supportsList,"supportsList",logFile);

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


void StepStandardisation::findSupports_two(dbVector& iPoint, dbMatrix& pointList,
		intVector& supportsList, dbVector& pointRadii, InputFileData* InputData,
		ofstream& logFile) {

	int findSupport_choice = 0;

	double dist_sum = 0, sqr_dist = 0;
	dbVector distVec(pointList.size());
	for(int i=0; i<pointList.size(); i++){
		sqr_dist = 0;
		for(int j=0; j<pointList[i].size(); j++){
			sqr_dist += pow(pointList[i][j]-iPoint[j],2);
		}
		distVec[i] = sqrt(sqr_dist);
	}

//	printVector(distVec,"distVec (before sorting)",logFile);
	std::sort(distVec.begin(),distVec.end());
//	printVector(distVec,"distVec (after sorting)",logFile);

//	double radiusFactor = 1.1;
	double radiusFactor = InputData->getValue("stepStandRadiusFactor");

	int counter = 0;
	dbMatrix pointInfluenceRange(iPoint.size(), dbVector(2));
	pointRadii.resize(iPoint.size());
	for (int i = 0; i < iPoint.size(); i++) {
		pointInfluenceRange[counter][0] = iPoint[i]
				- (distVec[1] * radiusFactor); 	// Min value
		pointInfluenceRange[counter][1] = iPoint[i]
				+ (distVec[1] * radiusFactor);	// Max value
		pointRadii[counter] = distVec[1] * radiusFactor;
		counter++;
	}

//	printMatrix(pointInfluenceRange,"pointInfluenceRange (1st col [min], 2nd col [max])",logFile);
//	printVector(pointRadii,"pointRadii",logFile);


	// =========================================================================

	// Temporary solution -> need to find a better one
//	double radiusFactor = InputData->getValue("stepStandardRadiusFactor");
//
//	// Determine the max and min value of each parameter(i.e Range) and store
//	// them in the matrix "myParamaterInfluenceRange"
//	int counter = 0;
//	dbMatrix pointInfluenceRange(iPoint.size(), dbVector(2));
//	pointRadii.resize(iPoint.size());
//	for (int i = 0; i < iPoint.size(); i++) {
//		pointInfluenceRange[counter][0] = iPoint[i]
//				- (iPoint[i] * radiusFactor); 	// Min value
//		pointInfluenceRange[counter][1] = iPoint[i]
//				+ (iPoint[i] * radiusFactor);	// Max value
//		pointRadii[counter] = iPoint[i] * radiusFactor;
//		counter++;
//	}
	// =========================================================================

	// Determine if the data parameters are within the influence range
//	dbVector dataParametersVec;
//	int pointInsideCounter = 0;
//	for (int i = 0; i < pointList.size(); i++) {
//
//		// Extract parameters of of a particular data
//		dbVector pointVec = pointList[i];
//
//		// Loop over each parameter and check if they are within the range
//		// specified. If one of them is not, the rest of the parameters are
//		// skipped.
//		pointInsideCounter = 0;
//		for (int j = 0; j < pointInfluenceRange.size(); j++) {
//
////			logFile << pointInfluenceRange[j][0] << " < " << pointVec[j]
////					<< " < " << pointInfluenceRange[j][1] << endl;
//
//			if (pointInfluenceRange[j][0] <= pointVec[j]
//					&& pointVec[j] <= pointInfluenceRange[j][1]) {
//				pointInsideCounter++;
////				logFile << "Taken" << endl << endl;
//
//			} else
//				// No need to continue the comparison process if one of the
//				// parameters is out of range
//				break;
//		}
//
//		// If the all parameters are inside the influence range, record ID of
//		// Data
//		if (pointInsideCounter == pointInfluenceRange.size()) {
//			supportsList.resize(supportsList.size() + 1);
//			supportsList[supportsList.size() - 1] = i;
//		}
//	}

	vector<int>().swap(supportsList);
	for (int i = 0; i < pointList.size(); i++) {
		supportsList.resize(supportsList.size() + 1);
		supportsList[supportsList.size() - 1] = i;
	}

	if (supportsList.size() < 2) {
		logFile << "ERROR: Point does not have enough of supports" << endl;
		cout << "ERROR: Point does not have enough of supports" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (supportsList.size() > 2 && InputData->getValue("StepStandardisationSupportFilter")==1) {

		if (pointList[0].size() == 1) {

//			logFile
//					<< "In StepStandardisation::findSupports, retaining only the two closest points."
//					<< endl;

			double dist_sum = 0, sqr_dist = 0;
			int postiviteID = -1, negativeID = -1;
			double dist = 0, maxDist = 0;
			double dist_positive = DBL_MAX;
			double dist_negative = -DBL_MAX;
			for (int i = 0; i < supportsList.size(); i++) {
				dist = 0;
				for (int j = 0; j < pointList[supportsList[i]].size(); j++) {

					dist += pointList[supportsList[i]][j] - iPoint[j];

//					logFile << "Comparing points: " << pointList[supportsList[i]][j] <<" and " << iPoint[j] << " with dist: " << dist << endl;
				}

				if(dist > 0){
//					logFile << "dist[" << dist << "] < [" << dist_positive << "]dist_positive: ";
					if(dist < dist_positive){
						dist_positive = dist;
						postiviteID = supportsList[i];
//						logFile << "Recorded with postiviteID: " << postiviteID;
					}
//					logFile << endl;
				}
				else{
//					logFile << "dist[" << dist << "] > [" << dist_negative << "]dist_negative: ";
					if(dist > dist_negative){
						dist_negative = dist;
						negativeID = supportsList[i];
//						logFile << "Recorded with negativeID: " << negativeID;
					}
//					logFile << endl;
				}

			}

			vector<int>().swap(supportsList);

			supportsList.push_back(postiviteID);
			supportsList.push_back(negativeID);

			double radiusFactor = InputData->getValue("stepStandRadiusFactor");

//			logFile << "dist_positive: " << dist_positive << endl;
//			logFile << "dist_negative: " << dist_negative << endl;
			if(fabs(dist_negative) > dist_positive)
				maxDist = fabs(dist_negative);
			else
				maxDist = dist_positive;

//			logFile << "maxDist: " << maxDist << endl;

			int counter = 0;
			pointRadii.clear();
			pointRadii.resize(iPoint.size());
			for (int i = 0; i < iPoint.size(); i++) {
				pointInfluenceRange[counter][0] = iPoint[i]
						- (maxDist * radiusFactor); 	// Min value
				pointInfluenceRange[counter][1] = iPoint[i]
						+ (maxDist * radiusFactor);	// Max value
				pointRadii[counter] = maxDist * radiusFactor;
				counter++;
			}

		}

	}

#ifdef _StepStandardisationDebugMode_
	logFile << "******* Supporting Points *******" << endl;

	logFile << "Main point: ";
	for (int i = 0; i < iPoint.size(); i++) {
		logFile << iPoint[i] << ",";
	}
	logFile << endl;

	printMatrix(pointList,"pointList",logFile);

	printVector(supportsList,"supportsList",logFile);

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

void StepStandardisation::findSupports_three(dbVector& iPoint,
		dbMatrix& pointList, intVector& supportsList, dbVector& pointRadii,
		InputFileData* InputData, ofstream& logFile) {

	vector<int>().swap(supportsList);
	for (int i = 0; i < pointList.size(); i++) {
		supportsList.resize(supportsList.size() + 1);
		supportsList[supportsList.size() - 1] = i;
	}

	if (supportsList.size() < 2) {
		logFile << "ERROR: Point does not have enough of supports" << endl;
		cout << "ERROR: Point does not have enough of supports" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (pointList[0].size() == 1) {

		logFile
				<< "In StepStandardisation::findSupports, retaining only the two closest points."
				<< endl;

		double dist_sum = 0, sqr_dist = 0;
		int postiviteID = -1, negativeID = -1;
		double dist = 0, maxDist = 0;
		double dist_positive = DBL_MAX;
		double dist_negative = -DBL_MAX;
		for (int i = 0; i < supportsList.size(); i++) {
			dist = 0;
			for (int j = 0; j < pointList[supportsList[i]].size(); j++) {

				dist += pointList[supportsList[i]][j] - iPoint[j];

				logFile << "Comparing points: " << pointList[supportsList[i]][j]
						<< " and " << iPoint[j] << " with dist: " << dist
						<< endl;
			}

			if (dist > 0) {
				logFile << "dist[" << dist << "] < [" << dist_positive
						<< "]dist_positive: ";
				if (dist < dist_positive) {
					dist_positive = dist;
					postiviteID = supportsList[i];
					logFile << "Recorded with postiviteID: " << postiviteID;
				}
				logFile << endl;
			} else {
				logFile << "dist[" << dist << "] > [" << dist_negative
						<< "]dist_negative: ";
				if (dist > dist_negative) {
					dist_negative = dist;
					negativeID = supportsList[i];
					logFile << "Recorded with negativeID: " << negativeID;
				}
				logFile << endl;
			}

		}

		vector<int>().swap(supportsList);

		supportsList.push_back(postiviteID);
		supportsList.push_back(negativeID);

		double radiusFactor = InputData->getValue("stepStandRadiusFactor");

		logFile << "dist_positive: " << dist_positive << endl;
		logFile << "dist_negative: " << dist_negative << endl;
		if (fabs(dist_negative) > dist_positive)
			maxDist = fabs(dist_negative);
		else
			maxDist = dist_positive;

		logFile << "maxDist: " << maxDist << endl;

		int counter = 0;
		pointRadii.clear();
		pointRadii.resize(iPoint.size());
		dbMatrix pointInfluenceRange(iPoint.size(), dbVector(2));
		for (int i = 0; i < iPoint.size(); i++) {
			pointInfluenceRange[counter][0] = iPoint[i]
					- (maxDist * radiusFactor); 	// Min value
			pointInfluenceRange[counter][1] = iPoint[i]
					+ (maxDist * radiusFactor);	// Max value
			pointRadii[counter] = maxDist * radiusFactor;
			counter++;
		}

	} else {

		logFile
				<< "Error: StepStandardisation::findSupports_three only support 1D points"
				<< endl;
		cout
				<< "Error: StepStandardisation::findSupports_three only support 1D points"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

#ifdef _StepStandardisationDebugMode_
	logFile << "******* Supporting Points *******" << endl;

	logFile << "Main point: ";
	for (int i = 0; i < iPoint.size(); i++) {
		logFile << iPoint[i] << ",";
	}
	logFile << endl;

	printMatrix(pointList,"pointList",logFile);

	printVector(supportsList,"supportsList",logFile);

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
//!
void StepStandardisation::findCommonCardiacStepsList(
	   DataContainer* problemData, InputFileData* InputData, ofstream& logFile){

	int nData = problemData->getIntVector("supportDataID").size();

	dbMatrix& leftCavityVolumesList =
			problemData->getDbMatrix("leftCavityVolumesList");
	dbMatrix leftCavityPressuresList =
				problemData->getDbMatrix("leftCavityPressuresList");

	dbMatrix& rightCavityVolumesList =
			problemData->getDbMatrix("rightCavityVolumesList");
	dbMatrix rightCavityPressuresList =
			problemData->getDbMatrix("rightCavityPressuresList");

	logFile << "leftCavityVolumesList size: " << leftCavityVolumesList.size() << endl;
	logFile << "leftCavityPressuresList size: " << leftCavityVolumesList.size() << endl;
	logFile << "rightCavityVolumesList size: " << leftCavityVolumesList.size() << endl;
	logFile << "rightCavityPressuresList size: " << leftCavityVolumesList.size() << endl;

	dbMatrix& stepHistoryList = problemData->getDbMatrix("stepHistoryList");


	phaseIndex = intMatrix(nData,intVector());

	for(int i = 0; i < nData; i++){

		intVector cPi;

		getCardiacLoopDetails(leftCavityVolumesList[i], leftCavityPressuresList[i],
				stepHistoryList[i],cPi,problemData,InputData,logFile);

		phaseIndex[i] = cPi;

		printMatrix(phaseIndex,"phaseIndex",logFile);
	}


	// Retain only the common cardiac phases (simple implementation)
	int minPhaseIndex = INT_MAX;
	for(int i = 0 ; i < phaseIndex.size(); i++){
		if(minPhaseIndex > phaseIndex[i].size())
			minPhaseIndex = phaseIndex[i].size();
	}
	if(minPhaseIndex <= 0){
		cout << "ERROR: In StepStandardisation::findCommonCardiacStepsList, the"
				" minimum of cardiac phases recorded is zero." << endl;
		logFile << "ERROR: In StepStandardisation::findCommonCardiacStepsList, the"
				" minimum of cardiac phases recorded is zero." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else{
		for(int i = 0 ; i < phaseIndex.size(); i++){
			phaseIndex[i].resize(minPhaseIndex);
		}
	}

	nPhase = minPhaseIndex - 1;


//	double startTime, endTime, incrementTime;
//	for(int i = 0; i < nPhase; i++){
//
//		logFile << "Phase: " << i << endl;
//		logFile << "-------" << endl;
//		dbMatrix standardStepPhaseList(nData);
//		dbMatrix nonStandardStepPhaseList(nData);
//
//		for(int j = 0 ; j < nData; j++){
//
//			dbVector& dataStepHistory = stepHistoryList[j];
//
//			logFile << "Data[" << j << "]:" << endl;
//			printVector(dataStepHistory,"dataStepHistory",logFile);
//
//			dbVector::const_iterator phaseFirst = dataStepHistory.begin() + phaseIndex[j][i];
//			dbVector::const_iterator phaselast = dataStepHistory.begin() + phaseIndex[j][i+1]+1;
//			dbVector nonStandardStepPhase(phaseFirst,phaselast);
//
//			startTime = dataStepHistory[phaseIndex[j][i]];
//			endTime = dataStepHistory[phaseIndex[j][i+1]];
//			incrementTime = endTime - startTime;
//
//			dbVector standardStepPhase = normStandStep;
//			for(int k=0; k<normStandStep.size(); k++){
//				standardStepPhase[k] = (standardStepPhase[k]*incrementTime) + startTime;
//			}
//
//			standardStepPhaseList[j] = standardStepPhase;
//			nonStandardStepPhaseList[j] = nonStandardStepPhase;
//
//		}
//
//		string ss = "standardStepPhaseList" + std::to_string(i);
//		problemData->setValue(ss.c_str(),standardStepPhaseList);
//
//		string ns = "nonStandardStepPhaseList" + std::to_string(i);
//		problemData->setValue(ns.c_str(),nonStandardStepPhaseList);
//
//		printMatrix(standardStepPhaseList,"standardStepPhaseList",logFile);
//		printMatrix(nonStandardStepPhaseList,"nonStandardStepPhaseList",logFile);
//	}

	vector<dbMatrix> standardStepPhaseListMat, nonStandardStepPhaseListMat;
	createPhaseList(stepHistoryList,standardStepPhaseListMat,
			nonStandardStepPhaseListMat, problemData, InputData, logFile);

	for(int i = 0; i < nPhase; i++){
		string ss = "standardStepPhaseList" + std::to_string(i);
		problemData->setValue(ss.c_str(),standardStepPhaseListMat[i]);

		string ns = "nonStandardStepPhaseList" + std::to_string(i);
		problemData->setValue(ns.c_str(),nonStandardStepPhaseListMat[i]);

		printMatrix(standardStepPhaseListMat[i],"standardStepPhaseList",logFile);
		printMatrix(nonStandardStepPhaseListMat[i],"nonStandardStepPhaseList",logFile);
	}

//	vector<dbMatrix> standardLPPhaseListMat, nonStandardLPPhaseListMat;
//	createPhaseList(leftCavityPressuresList, standardLPPhaseListMat,
//			nonStandardLPPhaseListMat, problemData, InputData, logFile);
//
//	for (int i = 0; i < nPhase; i++) {
//		string ss = "standardLPPhaseList" + std::to_string(i);
//		problemData->setValue(ss.c_str(), standardLPPhaseListMat[i]);
//
//		string ns = "nonStandardLPPhaseList" + std::to_string(i);
//		problemData->setValue(ns.c_str(), standardLPPhaseListMat[i]);
//
//		printMatrix(standardLPPhaseListMat[i], "standardLPPhaseList", logFile);
//		printMatrix(nonStandardLPPhaseListMat[i], "nonStandardLPPhaseList",
//				logFile);
//	}
//
//	vector<dbMatrix> standardRPPhaseListMat, nonStandardRPPhaseListMat;
//	createPhaseList(rightCavityPressuresList, standardRPPhaseListMat,
//			nonStandardRPPhaseListMat, problemData, InputData, logFile);
//
//	for (int i = 0; i < nPhase; i++) {
//		string ss = "standardRPPhaseList" + std::to_string(i);
//		problemData->setValue(ss.c_str(), standardRPPhaseListMat[i]);
//
//		string ns = "nonStandardRPPhaseList" + std::to_string(i);
//		problemData->setValue(ns.c_str(), standardRPPhaseListMat[i]);
//
//		printMatrix(standardRPPhaseListMat[i], "standardRPPhaseList", logFile);
//		printMatrix(nonStandardRPPhaseListMat[i], "nonStandardRPPhaseList",
//				logFile);
//	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void StepStandardisation::createPhaseList(dbMatrix& stepList,
		vector<dbMatrix>& standardStepPhaseListMat,
		vector<dbMatrix>& nonStandardStepPhaseListMat,
		DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile) {

	int nData = stepList.size();

	double startTime, endTime, incrementTime;
	for (int i = 0; i < nPhase; i++) {

		logFile << "Phase: " << i << endl;
		logFile << "-------" << endl;
		dbMatrix standardStepPhaseList(nData);
		dbMatrix nonStandardStepPhaseList(nData);

		for (int j = 0; j < nData; j++) {

			dbVector& dataStepHistory = stepList[j];

			logFile << "Data[" << j << "]:" << endl;
			printVector(dataStepHistory, "dataStepHistory", logFile);

			dbVector::const_iterator phaseFirst = dataStepHistory.begin()
					+ phaseIndex[j][i];
			dbVector::const_iterator phaselast = dataStepHistory.begin()
					+ phaseIndex[j][i + 1] + 1;
			dbVector nonStandardStepPhase(phaseFirst, phaselast);

			startTime = dataStepHistory[phaseIndex[j][i]];
			endTime = dataStepHistory[phaseIndex[j][i + 1]];
			incrementTime = endTime - startTime;

			dbVector standardStepPhase = normStandStep;
			for (int k = 0; k < normStandStep.size(); k++) {
				standardStepPhase[k] = (standardStepPhase[k] * incrementTime)
						+ startTime;
			}

			standardStepPhaseList[j] = standardStepPhase;
			nonStandardStepPhaseList[j] = nonStandardStepPhase;

		}

		standardStepPhaseListMat.push_back(standardStepPhaseList);
		nonStandardStepPhaseListMat.push_back(nonStandardStepPhaseList);

//		printMatrix(standardStepPhaseList, "standardStepPhaseList", logFile);
//		printMatrix(nonStandardStepPhaseList, "nonStandardStepPhaseList",
//				logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
void StepStandardisation::getCardiacLoopDetails(
		dbVector& volumesList,dbVector& pressuresList, dbVector& stepList,
		intVector& cPi,
		DataContainer* problemData, InputFileData* InputData, ofstream& logFile){


	intVector eDi; // End-diastole index
	intVector eIi; // End-isovolumetric contraction index
	intVector eEi; // End-ejection index
	intVector eRi; // End-relaxation index

	//Check
	if(volumesList.size() != pressuresList.size()){
		cout << "ERROR: In StepStandardisation::getCardiacLoopDetails(), volumesList and pressureList are not of the same size." << endl;
		logFile << "ERROR: In StepStandardisation::getCardiacLoopDetails(), volumesList and pressureList are not of the same size." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	cPi.push_back(0);

	// normalise volume and pressure vector
	double vNorm = computeNorm(volumesList,2,logFile);
	double pNorm = computeNorm(pressuresList,2,logFile);

#ifdef _StepStandardisationDebugMode_
	logFile << "volumesList.size(): " << volumesList.size() << endl;
	logFile << "pressuresList.size(): " << pressuresList.size() << endl;

	logFile <<"Volume and pressure list (before normalisation)" << endl;
	printVector(volumesList,"volumesList",logFile);
	printVector(pressuresList,"pressuresList",logFile);

	logFile << "vNorm = " << vNorm << endl;
	logFile << "pNorm = " << pNorm << endl;
#endif

	for(int i=0; i<volumesList.size();i++){
		volumesList[i] = volumesList[i]/vNorm;
		pressuresList[i] = pressuresList[i]/pNorm;
	}

	double maxVol=0, maxPress=0;
	for(int i=0; i<volumesList.size();i++){

		if(maxVol<volumesList[i]){
			maxVol = volumesList[i];
		}

		if(maxPress<pressuresList[i]){
			maxPress = pressuresList[i];
		}

	}

	for(int i=0; i<volumesList.size();i++){
			volumesList[i] = volumesList[i]/maxVol;
			pressuresList[i] = pressuresList[i]/maxPress;
	}

#ifdef _StepStandardisationDebugMode_
	logFile <<"Volume and pressure list (after normalisation)" << endl;
	printVector(volumesList,"volumesList",logFile);
	printVector(pressuresList,"pressuresList",logFile);
#endif

	bool eDFound = false, eIFound = true, eEFound = true, eRFound = true, phaseJustFound = false;
	double gradientPrevious=0, gradientNext=0, angle=0;
	for(int i=1; i < volumesList.size(); i++){

		gradientPrevious = (volumesList[i] - volumesList[i-1])/(pressuresList[i] - pressuresList[i-1]);

		if(i != (volumesList.size()-1))
			gradientNext = (volumesList[i+1] - volumesList[i])/(pressuresList[i+1] - pressuresList[i]);

		angle = atan((gradientPrevious-gradientNext)/(1+(gradientPrevious*gradientNext)));

#ifdef _StepStandardisationDebugMode_
		logFile << "[" << i << "]"
		<< "\t PV Gradient(previous): " << gradientPrevious
		<< "\t PV Gradient(next): " << gradientNext
		<< "\t angle: " << angle
		<< "\t Volume(previous): " << volumesList[i-1]
        << "\t Volume(current): " << volumesList[i]
		<< "\t Pressure(previous): " << pressuresList[i-1]
		<< "\t Pressure(current): " << pressuresList[i]
		<< endl;
		logFile << "---------------------------------------" << endl;
#endif

		// Find the end-diastole point
		if(eDFound == false){
			if((gradientNext < 1 && phaseJustFound == false) || i==(volumesList.size()-1)){

#ifdef _StepStandardisationDebugMode_
				logFile << "!! FOUND: End-diastole point" << endl;
#endif

				eDi.push_back(i);
				cPi.push_back(i);
				eDFound = true;
				eIFound = false;
				phaseJustFound = true;
				continue;
			}
		}

		// Find the end-isovolumetric contraction point
		if(eIFound == false){
			if((angle > 1.2e-02 && phaseJustFound == false) || i==(volumesList.size()-1)){

#ifdef _StepStandardisationDebugMode_
				logFile << "!! FOUND: End-isovolumetric contraction point" << endl;
#endif

				eIi.push_back(i);
				cPi.push_back(i);
				eIFound = true;
				eEFound = false;
				phaseJustFound = true;
				continue;
			}
		}

		// Find the end-ejection point
		if(eEFound == false){
			if((gradientPrevious < 0 && gradientNext < 0 && angle < 1.0e-5 && phaseJustFound == false) || i==(volumesList.size()-1)){

#ifdef _StepStandardisationDebugMode_
				logFile << "!! FOUND: End-ejection contraction point" << endl;
#endif

				eEi.push_back(i);
				cPi.push_back(i);
				eEFound = true;
				eRFound = false;
				phaseJustFound = true;
				continue;
			}
		}

		// Find the end-relaxation point
		if (eRFound == false) {
			if (gradientPrevious > 99999 && phaseJustFound == false) {

#ifdef _StepStandardisationDebugMode_
				logFile << "!! FOUND: End-relaxation point" << endl;
#endif

				eRi.push_back(i);
				cPi.push_back(i);
				eRFound = true;
				eDFound = false;
				phaseJustFound = true;
				continue;
			}
		}

		phaseJustFound = false;

	}

	if(eDFound == true && eIFound == true && eEFound == true && eRFound == false && volumesList.back() != eEi.back()){

#ifdef _StepStandardisationDebugMode_
		logFile << "!! FOUND: End-relaxation point" << endl;
#endif

		eRi.push_back(volumesList.size() - 1);
		cPi.push_back(volumesList.size() - 1);
		eRFound = true;
	}


	// Print the heart cycle details
	logFile << "========================" << endl;
	logFile << "Heart cycle information " << endl;
	logFile << "========================" << endl;

	logFile << "End-diastole: " << endl;
	logFile << "--------------" << endl;
	for(int i=0; i<eDi.size(); i++){
		logFile << "[" << i << "]" << endl;
		logFile << "eDi: " << eDi[i] << endl;
		logFile << "Volume = "<< volumesList[eDi[i]]  << endl;
		logFile << "Pressure = " << pressuresList[eDi[i]] << endl;
		logFile << "Timestep = " << stepList[eDi[i]] << endl;
	}

	logFile << "End-isovolumetric contraction:" << endl;
	logFile << "------------------------------" << endl;
	for(int i=0; i<eIi.size(); i++){
		logFile << "[" << i << "]" << endl;
		logFile << "eIi: " << eIi[i] << endl;
		logFile << "Volume = "<< volumesList[eIi[i]] << endl;
		logFile << "Pressure = " << pressuresList[eIi[i]] << endl;
		logFile << "Timestep = " << stepList[eIi[i]] << endl;
	}

	logFile << "End-ejection: " << endl;
	logFile << "--------------" << endl;
	for(int i=0; i<eEi.size(); i++){
		logFile << "[" << i << "]" << endl;
		logFile << "eEi: " << eEi[i] << endl;
		logFile << "Volume = "<< volumesList[eEi[i]]  << endl;
		logFile << "Pressure = " << pressuresList[eEi[i]] << endl;
		logFile << "Timestep = " << stepList[eEi[i]] << endl;
	}

	logFile << "End-relaxation: " << endl;
	logFile << "----------------" << endl;
	for(int i=0; i<eRi.size(); i++){
		logFile << "[" << i << "]" << endl;
		logFile << "eRi: " << eRi[i] << endl;
		logFile << volumesList[eRi[i]] << endl;
		logFile << "Pressure = " << pressuresList[eRi[i]] << endl;
		logFile << "Timestep = " << stepList[eRi[i]] << endl;
	}

//	cPi.clear();
//	cPi.push_back(0);
//	cPi.push_back(volumesList.size()-1);


}

/*!****************************************************************************/
/*!****************************************************************************/
//!
dbVector StepStandardisation::getStandardSteps(dbVector& interpolants,
	   DataContainer* problemData, InputFileData* InputData, ofstream& logFile){

#ifdef _StepStandardisationDebugMode_
	logFile << "In StepStandardisation::getStandardSteps" << endl;
#endif

	dbVector interpolatedSteps;

	for(int i=0; i<nPhase; i++){
		string ss = "standardStepPhaseList" + std::to_string(i);
		dbMatrix& standardStepPhaseList = problemData->getDbMatrix(ss.c_str());

#ifdef _StepStandardisationDebugMode_
		printMatrix(standardStepPhaseList,"standardStepPhaseList",logFile);
#endif

		int nSteps = standardStepPhaseList[0].size();
		dbVector interpolatedStepPhase(nSteps,0);
		int nData = interpolants.size();

		if(nData != standardStepPhaseList.size()){
			logFile << "ERROR: In StepStandardisation::getStandardSteps, "
					"interpolants size and standardStepPhaseList size are not"
					" the same." << endl;
			cout << "ERROR: In StepStandardisation::getStandardSteps, "
					"interpolants size and standardStepPhaseList size are not"
					" the same." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for(int j=0; j<nSteps; j++){
			for(int k=0; k<nData; k++){
				interpolatedStepPhase[j] += interpolants[k] * standardStepPhaseList[k][j];
			}
		}

#ifdef _StepStandardisationDebugMode_
		printVector(interpolatedStepPhase,"interpolatedStepPhase",logFile);
#endif

		if(interpolatedSteps.size() == 0)
			interpolatedSteps = interpolatedStepPhase;
		else{
			interpolatedSteps.pop_back();
			interpolatedSteps.insert(interpolatedSteps.end(),
					interpolatedStepPhase.begin(),interpolatedStepPhase.end());
		}

#ifdef _StepStandardisationDebugMode_
		printVector(interpolatedSteps,"interpolatedSteps",logFile);
#endif

	}

	return interpolatedSteps;

}

/*!****************************************************************************/
/*!****************************************************************************/
//!
dbVector StepStandardisation::interpolateStandardSteps(const char* namePhaseList,
		dbVector& interpolants, DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	logFile << "In StepStandardisation::getStandardSteps" << endl;

	dbVector interpolatedSteps;

	for(int i=0; i<nPhase; i++){
		string ss = namePhaseList + std::to_string(i);
		dbMatrix& standardStepPhaseList = problemData->getDbMatrix(ss.c_str());

		printMatrix(standardStepPhaseList,"standardStepPhaseList",logFile);

		int nSteps = standardStepPhaseList[0].size();
		dbVector interpolatedStepPhase(nSteps,0);
		int nData = interpolants.size();

		if(nData != standardStepPhaseList.size()){
			logFile << "ERROR: In StepStandardisation::getStandardSteps, "
					"interpolants size and standardStepPhaseList size are not"
					" the same." << endl;
			cout << "ERROR: In StepStandardisation::getStandardSteps, "
					"interpolants size and standardStepPhaseList size are not"
					" the same." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for(int j=0; j<nSteps; j++){
			for(int k=0; k<nData; k++){
				interpolatedStepPhase[j] += interpolants[k] * standardStepPhaseList[k][j];
			}
		}

		printVector(interpolatedStepPhase,"interpolatedStepPhase",logFile);

		if(interpolatedSteps.size() == 0)
			interpolatedSteps = interpolatedStepPhase;
		else{
			interpolatedSteps.pop_back();
			interpolatedSteps.insert(interpolatedSteps.end(),
					interpolatedStepPhase.begin(),interpolatedStepPhase.end());

		}
		printVector(interpolatedSteps,"interpolatedSteps",logFile);

	}

	return interpolatedSteps;
}

/*!****************************************************************************/
/*!****************************************************************************/
//!
dbMatrix StepStandardisation::combinePhaseList(vector<dbMatrix> phaseList,
		DataContainer* problemData,InputFileData* InputData, ofstream& logFile){

	int numPhase = phaseList.size();

	int nData;
	if(numPhase > 0)
		nData = phaseList[0].size();
	else{
		nData = 0 ;
		return dbMatrix();
	}

	dbMatrix combinedPhase(nData);

	for(int i=0; i < nData; i++){
		dbVector combinedDataPhase;
		for(int j=0; j<numPhase; j++){
			dbVector selectedDataPhase = phaseList[j][i];
			for(int k=0; k<selectedDataPhase.size();k++){
				combinedDataPhase.push_back(selectedDataPhase[k]);
			}
			if(j != numPhase-1){
				combinedDataPhase.pop_back();
			}
		}
		combinedPhase[i] = combinedDataPhase;
	}


	return combinedPhase;

}
