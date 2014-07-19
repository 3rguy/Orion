#include "PODICalc.h"

using namespace std;

/*!****************************************************************************/
/*!****************************************************************************/
PODICalc::PODICalc(InputFileData* InputData,ofstream& logFile){

	logFile << "PODICalc::PODICalc(InputFileData* InputData,ofstream& logFile)"
			" has not yet been implemented" << endl;
	MPI_Abort(MPI_COMM_WORLD,1);

}

/*!****************************************************************************/
/*!****************************************************************************/
PODICalc::PODICalc(dbVector& myParameters, dbVector& parameterRadii,
		intVector& supportDataID, dbMatrix& dataParametersList,
		vector<dbMatrix>& displacementList,
		dbMatrix& resultingDisplacementMatrix,InputFileData* InputData,
		ofstream& logFile){

	// Calculate the interpolants used during the interpolation process
	dbVector interpolants;
	interpolantsCalc(myParameters,dataParametersList,parameterRadii,
			interpolants,InputData,logFile);

	vector<dbMatrix> rearrangedDisplacementList;
	if(InputData->getValue("dispMatrixRearrange") == 1 )
		rearrangeDisplacementMatrix(displacementList,rearrangedDisplacementList,
				InputData,logFile);
	else
		rearrangedDisplacementList = displacementList;

	int choice = InputData->getValue("PODICalculationType");
	switch(choice){
	case 1:
		PODInterpolation(rearrangedDisplacementList,interpolants,
				resultingDisplacementMatrix,InputData,logFile);
		break;
	case 2:
		PODInterpolationEnhanced(rearrangedDisplacementList,interpolants,
				resultingDisplacementMatrix,InputData,logFile);
		break;
	default:
		logFile << "PODICalculationType: " << choice << " does not exist.\n"
				" The valid options are <1,2>." << endl;
		cout << "PODICalculationType: " << choice << " does not exist.\n"
						"Valid options are <1,2>." << endl;
		MPI_Abort(MPI_COMM_WORLD,1);

	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void PODICalc::interpolantsCalc(dbVector& myParameters,
		dbMatrix& dataParametersList,dbVector& parametersRadii,
		dbVector& interpolants,	InputFileData* InputData,ofstream& logFile){

	int choice = InputData->getValue("interpolantionType");
	//Interpolation* interpolationSpace;

	switch(choice){

	case 1:

		// Constant Interpolation function
		cout << "ERROR: Constant Interpolation not yet implemented" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
		break;

	case 2:

	{
		Interpolation* MLSInterpolate = new Interpolation(myParameters,
				dataParametersList, parametersRadii, interpolants, InputData,
				logFile);

		delete MLSInterpolate;

	}
		break;

	default:
		cout << "ERROR: Interpolation Type not implemented" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
		break;

	}
}

/*!****************************************************************************/
/*!****************************************************************************/
void PODICalc::rearrangeDisplacementMatrix(vector<dbMatrix>& displacementList,
		vector<dbMatrix>& rearrangeDisplacementList, InputFileData* InputData,
		ofstream& logFile){


#ifdef _PODICalcDebugMode_
	logFile << "********** PODICalc::rearrangeResultMatrix **********"
			<< endl;
	logFile << "------- Result Matrix: Before rearrangement  -------"
			<< endl;
	for(int i=0; i < displacementList.size(); i++){
		printMatrix(displacementList[i],"",logFile);
		logFile << endl;
	}
#endif

	rearrangeDisplacementList.resize(displacementList[0][0].size(),
			dbMatrix(displacementList[0].size(),
					dbVector(displacementList.size())));

	for(int i = 0; i < displacementList[0][0].size(); i++){
		for(int j = 0; j < displacementList.size(); j++){
			for(int k = 0; k < displacementList[0].size(); k++){
				rearrangeDisplacementList[i][k][j] = displacementList[j][k][i];
			}
		}
	}

#ifdef _PODICalcDebugMode_
	logFile << "------- Result Matrix: After rearrangement  -------"
			<< endl;
	for(int i=0; i < rearrangeDisplacementList.size(); i++){
			printMatrix(rearrangeDisplacementList[i],"",logFile);
			logFile << endl;
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void PODICalc::PODInterpolation(vector<dbMatrix>& rearrangeDisplacementList,
		dbVector& interpolants, dbMatrix& resultingDisplacementMatrix,
		InputFileData* InputData, ofstream& logFile) {


#ifdef _PODICalcCheckMode_
	for(int i = 0; i < rearrangeDisplacementList.size(); i++) {
		if(rearrangeDisplacementList[0][0].size() != interpolants.size()) {
			cout << "ERROR: In PODICalc::PODInterpolation, the number of "
			"columns in the displacement matrix does not corresponds "
			"to the number of interpolants" << endl;
			logFile << "ERROR: In PODICalc::PODInterpolation, the number of "
			"columns in the displacement matrix does not corresponds"
			" to the number of interpolants" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
#endif

	int interpolation_choice = InputData->getValue("directInterpolation");
	if (interpolation_choice == 1) {

		logFile << "WARNING: POD CALCULATION WAS DISABLED AND DIRECT"
				" INTERPOLATION HAS BEEN ACTIVATED" << endl;
		cout << "WARNING: POD CALCULATION WAS DISABLED AND DIRECT"
				" INTERPOLATION HAS BEEN ACTIVATED" << endl;
	}


	resultingDisplacementMatrix.resize(rearrangeDisplacementList[0].size(),
			dbVector());


	double energyLevel = InputData->getValue("PODEnergyLevel");
	double totalPOVsConserved = 0;
	double totalEnergyConverserved = 0;
	for (int i = 0; i < rearrangeDisplacementList.size(); i++) {

		dbMatrix reduceDisplacementMatrix;
		dbVector interpolatedReducedDispVec, fullInterpolatedDispVec;

		if (interpolation_choice == 0) {

			//! ----------------------------------------------------------------
			//! POD Reduction + direction interpolation

			//! Setup the POD space and compress displacement matrix
			PODCalc* PODSpace = new PODCalc(rearrangeDisplacementList[i],
					reduceDisplacementMatrix, energyLevel, InputData, logFile);

			logFile << " ---------------------------------------------------"
					<< endl;
			logFile << "Carrying out the interpolation on the reduced matrix"
					<< endl;

			interpolatedReducedDispVec.resize(reduceDisplacementMatrix.size(),
					0);
			for (int j = 0; j < reduceDisplacementMatrix.size(); j++) {
				for (int k = 0; k < reduceDisplacementMatrix[j].size(); k++) {
					interpolatedReducedDispVec[j] +=
							reduceDisplacementMatrix[j][k] * interpolants[k];
				}
			}

			logFile << "Number of POVs conserved: "
					<< PODSpace->getNumPOVConserved() << endl;

			logFile << "Number of POVs conserved: "
					<< PODSpace->getEnergyConserved() << endl;

			totalPOVsConserved += PODSpace->getNumPOVConserved();
			totalEnergyConverserved += PODSpace->getEnergyConserved();

			//! Convert vector to high dimensional space
			PODSpace->expandVector(interpolatedReducedDispVec,
					fullInterpolatedDispVec, InputData, logFile);

			delete PODSpace;

		}
		else if (interpolation_choice == 1) {

			//! ----------------------------------------------------------------
			//! Direction interpolation

			fullInterpolatedDispVec.resize(rearrangeDisplacementList[i].size(),0);
			for (int j = 0; j < rearrangeDisplacementList[i].size(); j++) {
				for (int k = 0; k < rearrangeDisplacementList[i][j].size();
						k++) {
					fullInterpolatedDispVec[j] +=
							rearrangeDisplacementList[i][j][k]
									* interpolants[k];
				}
			}
		}
		else{
			logFile << "In PODICalc::PODInterpolation, interpolation method "
					"does not exist.\n directInterpolation(0,1)" << endl;
			cout << "In PODICalc::PODInterpolation, interpolation method "
					"does not exist.\n directInterpolation(0,1)" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);

		}

		//! Store the full interpolated displacement vector
		for (int j = 0; j < resultingDisplacementMatrix.size(); j++) {
			resultingDisplacementMatrix[j].resize(
					resultingDisplacementMatrix[j].size() + 1);
			resultingDisplacementMatrix[j][resultingDisplacementMatrix[j].size()
					- 1] = fullInterpolatedDispVec[j];
		}
	}

	// Info display
	if (interpolation_choice == 0) {

		logFile << "Average POVs conserved = "
				<< totalPOVsConserved/rearrangeDisplacementList.size() << endl;
		cout << "Average POVs conserved = "
			 << totalPOVsConserved/rearrangeDisplacementList.size() << endl;

		logFile << "Targeted energy level: " << energyLevel
				<< "\t Average Energy conserved: "
				<< totalEnergyConverserved/rearrangeDisplacementList.size()
				<< endl;
		cout << "Targeted energy level: " << energyLevel
			 << "\t Average Energy conserved: "
			 << totalEnergyConverserved/rearrangeDisplacementList.size()
			 << endl;
	}



#ifdef _PODICalcDebugMode_

	logFile << "********* PODICalc::PODInterpolation *********" << endl;

	logFile << "---------- Result Matrix List ----------" << endl;
	for(int i = 0; i < rearrangeDisplacementList.size(); i++) {
		logFile << "Matrix No: " << i << endl;
		for(int j = 0; j < rearrangeDisplacementList[i].size(); j++) {
			for(int k = 0; k < rearrangeDisplacementList[i][j].size(); k++) {
				logFile << rearrangeDisplacementList[i][j][k] << " ";
			}
			logFile << endl;
		}
	}

	logFile << "------- PODI Calculated Matrix -------" << endl;
	for(int i = 0; i<resultingDisplacementMatrix.size(); i++ ) {
		for(int j = 0; j<resultingDisplacementMatrix[i].size(); j++ ) {
			logFile << resultingDisplacementMatrix[i][j] << " ";
		}
		logFile << endl;
	}

#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
void PODICalc::PODInterpolationEnhanced(vector<dbMatrix>& rearrangeDisplacementList,
		dbVector& interpolants, dbMatrix& resultingDisplacementMatrix,
		InputFileData* InputData, ofstream& logFile) {

	logFile << "PODInterpolationEnhanced Algorithm used" << endl;


#ifdef _PODICalcCheckMode_
	for(int i = 0; i < rearrangeDisplacementList.size(); i++) {
		if(rearrangeDisplacementList[0][0].size() != interpolants.size()) {
			cout << "ERROR: In PODICalc::PODInterpolation, the number of "
			"columns in the displacement matrix does not corresponds "
			"to the number of interpolants" << endl;
			logFile << "ERROR: In PODICalc::PODInterpolation, the number of "
			"columns in the displacement matrix does not corresponds"
			" to the number of interpolants" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
#endif

	int interpolation_choice = InputData->getValue("directInterpolation");
	if (interpolation_choice == 1) {

		logFile << "WARNING: POD CALCULATION WAS DISABLED AND DIRECT"
				" INTEPOLATION HAS BEEN ACTIVATED" << endl;
		cout << "WARNING: POD CALCULATION WAS DISABLED AND DIRECT"
				" INTEPOLATION HAS BEEN ACTIVATED" << endl;
	}

#ifdef _PODICalcDebugMode_
	logFile << "********* PODICalc::PODInterpolation *********" << endl;
	logFile << "---------- Result Matrix List ----------" << endl;
	for(int i = 0; i < rearrangeDisplacementList.size(); i++) {
		logFile << "Matrix No: " << i << endl;
		for(int j = 0; j < rearrangeDisplacementList[i].size(); j++) {
			for(int k = 0; k < rearrangeDisplacementList[i][j].size(); k++) {
				logFile << rearrangeDisplacementList[i][j][k] << " ";
			}
			logFile << endl;
		}
	}
#endif


	resultingDisplacementMatrix.resize(rearrangeDisplacementList[0].size(),
			dbVector());

	double energyLevel = InputData->getValue("PODEnergyLevel");
	double totalPOVsConserved = 0;
	double totalEnergyConverserved = 0;
	for (int i = 0; i < rearrangeDisplacementList.size(); i++) {

		// Find the zeroes entries
		intVector zeroEntriesVec;
		dbMatrix& fullDispMat = rearrangeDisplacementList[i];
		for(int j=0; j<fullDispMat.size(); j++){

			int counter = 0;
			for(int k=0; k<fullDispMat[j].size(); k++){
				if(fullDispMat[j][k] != 0)
					break;
				else
					counter++;
			}

			if(counter == fullDispMat[j].size()){
				resizeArray(zeroEntriesVec,zeroEntriesVec.size()+1);
				zeroEntriesVec[zeroEntriesVec.size()-1] = j;
			}
		}

		//Erase the zeroes entries
		logFile << "Number of zero rows found = " << zeroEntriesVec.size() << endl;

		printMatrix(fullDispMat,"fullDispMat -> Full",logFile);

		for(int j=zeroEntriesVec.size()-1; j>0 ; j--){
//			logFile << "j: " << j << endl;
//			logFile << zeroEntriesVec[j] << ", ";
			fullDispMat.erase(fullDispMat.begin()+zeroEntriesVec[j]);
		}
		logFile << endl;

		printMatrix(fullDispMat,"fullDispMat -> Reduced",logFile);

		dbMatrix reduceDisplacementMatrix;
		dbVector interpolatedReducedDispVec, fullInterpolatedDispVec;

		if (interpolation_choice == 0) {

			//! ----------------------------------------------------------------
			//! POD Reduction + direction interpolation

			//! Setup the POD space and compress displacement matrix
			PODCalc* PODSpace = new PODCalc(rearrangeDisplacementList[i],
					reduceDisplacementMatrix, energyLevel, InputData, logFile);

			logFile << " ---------------------------------------------------"
					<< endl;
			logFile << "Carrying out the interpolation on the reduced matrix"
					<< endl;

			interpolatedReducedDispVec.resize(reduceDisplacementMatrix.size(),
					0);
			for (int j = 0; j < reduceDisplacementMatrix.size(); j++) {
				for (int k = 0; k < reduceDisplacementMatrix[j].size(); k++) {
					interpolatedReducedDispVec[j] +=
							reduceDisplacementMatrix[j][k] * interpolants[k];
				}
			}

			logFile << "Number of POVs conserved: "
					<< PODSpace->getNumPOVConserved() << endl;

			logFile << "Number of POVs conserved: "
					<< PODSpace->getEnergyConserved() << endl;

			totalPOVsConserved += PODSpace->getNumPOVConserved();
			totalEnergyConverserved += PODSpace->getEnergyConserved();

			//! Convert vector to high dimensional space
			PODSpace->expandVector(interpolatedReducedDispVec,
					fullInterpolatedDispVec, InputData, logFile);

			delete PODSpace;



		}
		else if (interpolation_choice == 1) {

			//! ----------------------------------------------------------------
			//! Direction interpolation
			fullInterpolatedDispVec.resize(rearrangeDisplacementList[i].size(),0);
			for (int j = 0; j < rearrangeDisplacementList[i].size(); j++) {
				for (int k = 0; k < rearrangeDisplacementList[i][j].size();
						k++) {
					fullInterpolatedDispVec[j] +=
							rearrangeDisplacementList[i][j][k]
									* interpolants[k];
				}
			}
		}
		else{
			logFile << "In PODICalc::PODInterpolation, interpolation method "
					"does not exist.\n directInterpolation(0,1)" << endl;
			cout << "In PODICalc::PODInterpolation, interpolation method "
					"does not exist.\n directInterpolation(0,1)" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);

		}

		// Restore the zeroes vector
//		for (int j = 0; j < zeroEntriesVec.size(); j++) {
//			fullDispMat.insert(fullDispMat.begin()+zeroEntriesVec[j],dbVector(fullDispMat[0].size(),0));
//		}
		for (int j = 0; j<zeroEntriesVec.size(); j++){
			fullInterpolatedDispVec.insert(fullInterpolatedDispVec.begin()+zeroEntriesVec[j],0);
		}

		//! Store the full interpolated displacement vector
		for (int j = 0; j < resultingDisplacementMatrix.size(); j++) {
			resultingDisplacementMatrix[j].resize(
					resultingDisplacementMatrix[j].size() + 1);
			resultingDisplacementMatrix[j][resultingDisplacementMatrix[j].size()
					- 1] = fullInterpolatedDispVec[j];
		}
	}

	// Info display
	if (interpolation_choice == 0) {

		logFile << "Average POVs conserved = "
				<< totalPOVsConserved/rearrangeDisplacementList.size() << endl;
		cout << "Average POVs conserved = "
			 << totalPOVsConserved/rearrangeDisplacementList.size() << endl;

		logFile << "Targeted energy level: " << energyLevel
				<< "\t Average Energy conserved: "
				<< totalEnergyConverserved/rearrangeDisplacementList.size()
				<< endl;
		cout << "Targeted energy level: " << energyLevel
			 << "\t Average Energy conserved: "
			 << totalEnergyConverserved/rearrangeDisplacementList.size()
			 << endl;
	}


#ifdef _PODICalcDebugMode_
	logFile << "------- PODI Calculated Matrix -------" << endl;
	for(int i = 0; i<resultingDisplacementMatrix.size(); i++ ) {
		for(int j = 0; j<resultingDisplacementMatrix[i].size(); j++ ) {
			logFile << resultingDisplacementMatrix[i][j] << " ";
		}
		logFile << endl;
	}
#endif

}


