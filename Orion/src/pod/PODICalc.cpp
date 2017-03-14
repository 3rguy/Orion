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
		dbMatrix& resultingDisplacementMatrix,DataContainer* problemData,
		InputFileData* InputData,
		ofstream& logFile){

	// Calculate the interpolants used during the interpolation process
	interpolantsCalc(myParameters,dataParametersList,parameterRadii,
			PODInterpolants,InputData,logFile);


	// Rearrange displacement matrix
	vector<dbMatrix> rearrangedDisplacementList;
	if(InputData->getValue("dispMatrixRearrange") == 1 )
		rearrangeDisplacementMatrix(displacementList,rearrangedDisplacementList,
				InputData,logFile);
	else
		rearrangedDisplacementList = displacementList;


	// Carry out the PODI calculation
	int choice = InputData->getValue("PODICalculationType");
	switch(choice){
	case 1:
		PODInterpolation(rearrangedDisplacementList,PODInterpolants,
				resultingDisplacementMatrix,problemData,InputData,logFile);
		break;
	case 2:
		PODInterpolationEnhanced(rearrangedDisplacementList,PODInterpolants,
				resultingDisplacementMatrix,problemData,InputData,logFile);
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
		dbMatrix& dataParametersList, dbVector& parametersRadii,
		dbVector& interpolants, InputFileData* InputData, ofstream& logFile) {

	int choice = InputData->getValue("interpolantionType");
	//Interpolation* interpolationSpace;

	switch (choice) {

	case 1:

		// Constant Interpolation function
		cout << "ERROR: Constant Interpolation has not yet been implemented"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
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
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;

	}

	printVector(interpolants,"****** Interpolants ******",logFile);
	oPType intpolantSum = 0;
	for(int i=0; i<interpolants.size();i++) intpolantSum += interpolants[i];
	logFile << "Sum of interpolants: " << intpolantSum << endl;

	if(InputData->getValue("printPodiInfo") == 1) {
		cout <<  "Sum of interpolants: " << intpolantSum << endl;
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
		DataContainer* problemData,InputFileData* InputData, ofstream& logFile)
{

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

	int PlotPOMs_choice = InputData->getValue("PlotPOMsAndPOVs");

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

			// Save POMs to file
			if (PlotPOMs_choice == 1) {
				cout << "Saving POMs to file" << endl;
				logFile << "Saving POMs to file" << endl;
				savePOMsToFile_ResFormat(PODSpace->getPOMs(), PODSpace->getMeanVec(),
						i, problemData, InputData, logFile);
				savePOVsToFile_GrfFormat(PODSpace->getPOVs(), i, problemData,
						InputData, logFile);
				saveConservedPOVsToFile_GrfFormat(PODSpace->getNumPOVConserved(),
						PODSpace->getEnergyConserved(), i, problemData,
						InputData, logFile);
			}

			delete PODSpace;

		}
		else if (interpolation_choice == 1) {

			//! ----------------------------------------------------------------
			//! Direction interpolation

			logFile << "Direction interpolation used" << endl;

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
	if (interpolation_choice == 0 && InputData->getValue("printPodiInfo") == 1) {

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
		DataContainer* problemData, InputFileData* InputData, ofstream& logFile)
{

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
		logFile << "Size: " << rearrangeDisplacementList[0].size() << " x " << rearrangeDisplacementList[0][0].size() << endl;
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
	int count_interpolChange = 0;
	for (int i = 0; i < rearrangeDisplacementList.size(); i++) {

		PODCalc* PODSpace = new PODCalc(rearrangeDisplacementList[i],
													InputData, logFile);

		// Find the zeroes entries
		intVector zeroEntriesVec;
		dbMatrix& fullDispMat = PODSpace->getDataMatrix();

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

		int PlotPOMs_choice = InputData->getValue("PlotPOMsAndPOVs");

		//Erase the zeroes entries
		if(zeroEntriesVec.size() > 0){

			PlotPOMs_choice = 0;

			logFile << "Number of zero rows found = " << zeroEntriesVec.size() << endl;
			for(int j=zeroEntriesVec.size()-1; j>-1 ; j--){
				fullDispMat.erase(fullDispMat.begin()+zeroEntriesVec[j]);
			}
			logFile << endl;
		}

		// If all entries has been deleted from the fullDispMat, the
		// interpolation_choice is switched to 1
		int oldInterpolation_choice = interpolation_choice;
		bool zeroDispMat = false;
		if(fullDispMat.size() == 0){
			zeroDispMat = true;
			interpolation_choice = 1;
			count_interpolChange++;
		}

		dbMatrix reduceDisplacementMatrix;
		dbVector interpolatedReducedDispVec, fullInterpolatedDispVec;

		if (interpolation_choice == 0) {

			//! ----------------------------------------------------------------
			//! POD Reduction + direction interpolation

			// Define the energy level to be conserved
			PODSpace->setEnergyLevel(energyLevel);

			// Calculate the POMs and POVs
			PODSpace->setPOVsandPOMs(InputData,logFile);

			// Determine the number of POMs needed to meet the minimum energy level
			// requirement
			PODSpace->energyConservCalc(logFile);

			// Reduce Matrix's size
			PODSpace->compressMatrix(reduceDisplacementMatrix,InputData,logFile);


			logFile << " ---------------------------------------------------"
					<< endl;
			logFile << "Interpolating the reduced matrix"
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

			logFile << "Energy conserved: "
					<< PODSpace->getEnergyConserved() << endl;

			totalPOVsConserved += PODSpace->getNumPOVConserved();
			totalEnergyConverserved += PODSpace->getEnergyConserved();

			//! Convert vector to high dimensional space
			PODSpace->expandVector(interpolatedReducedDispVec,
					fullInterpolatedDispVec, InputData, logFile);

			// Save POMs to file
			if (PlotPOMs_choice == 1) {
				savePOMsToFile_ResFormat(PODSpace->getPOMs(),
						PODSpace->getMeanVec(), i, problemData, InputData,
						logFile);

				printVector(PODSpace->getPOVs(),"Print vector, POVs: ",logFile);
				savePOVsToFile_GrfFormat(PODSpace->getPOVs(), i, problemData,
						InputData, logFile);
				saveConservedPOVsToFile_GrfFormat(PODSpace->getNumPOVConserved(),
						PODSpace->getEnergyConserved(), i, problemData,
						InputData, logFile);
			}
		}
		else if (interpolation_choice == 1) {

			//! ----------------------------------------------------------------
			//! Direction interpolation

			logFile << "Direction interpolation used" << endl;

			fullInterpolatedDispVec.resize(fullDispMat.size(),0);
			for (int j = 0; j < fullDispMat.size(); j++) {
				for (int k = 0; k < fullDispMat[j].size();
						k++) {
					fullInterpolatedDispVec[j] +=
							fullDispMat[j][k]
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

//		logFile << "---------- Before adding zero row ----------" << endl;
//		printVector(fullInterpolatedDispVec,"fullInterpolatedDispVec",logFile);

		//! Restore the zero rows
		for (int j = 0; j<zeroEntriesVec.size(); j++){
			fullInterpolatedDispVec.insert(fullInterpolatedDispVec.begin()+zeroEntriesVec[j],0);
		}

		// If a direct interpolation is carried out, the mean of the ensemble matrix is
		// not automatically added.
		PODSpace->addMean(fullInterpolatedDispVec,InputData,logFile);

		delete PODSpace;

		if(zeroDispMat == true)
					interpolation_choice = oldInterpolation_choice;

//		logFile << "---------- After adding zero row ----------" << endl;
//		printVector(fullInterpolatedDispVec,"fullInterpolatedDispVec",logFile);

		//! Store the full interpolated displacement vector
		for (int j = 0; j < resultingDisplacementMatrix.size(); j++) {
			resultingDisplacementMatrix[j].resize(
					resultingDisplacementMatrix[j].size() + 1);
			resultingDisplacementMatrix[j][resultingDisplacementMatrix[j].size()
					- 1] = fullInterpolatedDispVec[j];
		}
	}

	// PODI info display
	if (interpolation_choice == 0 && InputData->getValue("printPodiInfo") == 1) {

		logFile << "Average POVs conserved = "
				<< totalPOVsConserved/(rearrangeDisplacementList.size()-count_interpolChange) << endl;
		cout << "Average POVs conserved = "
			 << totalPOVsConserved/(rearrangeDisplacementList.size()-count_interpolChange) << endl;

		logFile << "Targeted energy level: " << energyLevel
				<< "\t Average Energy conserved: "
				<< totalEnergyConverserved/(rearrangeDisplacementList.size()-count_interpolChange)
				<< endl;
		cout << "Targeted energy level: " << energyLevel
			 << "\t Average Energy conserved: "
			 << totalEnergyConverserved/(rearrangeDisplacementList.size()-count_interpolChange)
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

/*!****************************************************************************/
/*!****************************************************************************/
void PODICalc::savePOMsToFile_ResFormat(dbMatrix& POMs, dbVector& meanVec,
		int stepValueID, DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile) {

	string& resultName = problemData->getString("PODIResultName");
	int& nDOFPerNode = problemData->getInt("PODIDofPerNode");

	string fileName = "POMs_" + resultName + ".res";

	bool fileExist = true;

	struct stat buffer;
	if (stat(fileName.c_str(), &buffer) != 0) { // File exists
		fileExist = false;
	}

	ofstream writeToFile(fileName, ios::out | ios::app);

	if (fileExist == false) {
		writeToFile << "GiD Post Results File 1.0" << endl;
	}

	string analysis_name = " ";
	string my_location = "onNodes";

	string DofType;
	if (nDOFPerNode == 1)
		DofType = "scalar";
	else if (nDOFPerNode > 1)
		DofType = "vector";
	else {
		cout << "In PODICalc::savePOMsToFile_ResFormat, nDOFPerNode is invalid"
				<< endl;
		logFile
				<< "In PODICalc::savePOMsToFile_ResFormat, nDOFPerNode is invalid"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	string componentNames = "ComponentNames ";
	for (int i = 0; i < nDOFPerNode; i++) {

		ostringstream convert;   // stream used for the conversion
		convert << i + 1;

		componentNames = componentNames + "\"DOF_" + convert.str() + "\", ";
	}

	int numDOFs = POMs.size();
	int numPOMs = POMs[0].size();

	logFile << "------------------------------------------------------" << endl;
//	printMatrix(POMs,"POMs",logFile);

	for (int i = 0; i < numPOMs; i++) {
		writeToFile << "Result " << "\"" << "POM_" << i << "\" \""
				<< analysis_name << "\" " << stepValueID << " " << DofType
				<< " " << my_location << endl;

		writeToFile << componentNames << endl;
		writeToFile << "Values" << endl;

		int DOFCounter = 0;
		int nodeCounter = 1;
		while (DOFCounter < numDOFs) {

//			logFile << "writeToFile[nodeCounter] " << nodeCounter << endl;
			writeToFile << nodeCounter << " ";

			for (int k = 0; k < nDOFPerNode; k++) {
//				logFile << "Accessing: " << DOFCounter << ", " << i << endl;
				writeToFile << POMs[DOFCounter][i] << " ";
				DOFCounter++;
			}
			writeToFile << endl;
			nodeCounter++;
		}

		writeToFile << "End Values" << endl;
	}

	// -------------------------------------------------------------------------
	// Write mean vector to file

	if (meanVec.size() > 0) {
		writeToFile << "Result " << "\"" << "Mean\" \"" << analysis_name
				<< "\" " << stepValueID << " " << DofType << " " << my_location
				<< endl;

		writeToFile << componentNames << endl;
		writeToFile << "Values" << endl;

		int DOFCounter = 0;
		int nodeCounter = 1;
		while (DOFCounter < numDOFs) {

			writeToFile << nodeCounter << " ";

			for (int k = 0; k < nDOFPerNode; k++) {
				writeToFile << meanVec[DOFCounter] << " ";
				DOFCounter++;
			}
			writeToFile << endl;
			nodeCounter++;
		}

		writeToFile << "End Values" << endl;
	}

	writeToFile.close();
}

/*!****************************************************************************/
/*!****************************************************************************/
void PODICalc::savePOVsToFile_GrfFormat(dbVector& POVs, int stepValueID,
		DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile){

	string& resultName = problemData->getString("PODIResultName");
	int& nDOFPerNode = problemData->getInt("PODIDofPerNode");

	string fileName = "POVs_" + resultName + ".grf";
	string fileNamePer = "POVs_Perctage_" + resultName + ".grf";

	bool fileExist = true;

	struct stat buffer;
	if (stat(fileName.c_str(), &buffer) != 0) { // File exists
		fileExist = false;
	}

	ofstream writeToFile(fileName,ios::out | ios::app);
	ofstream writeToFilePer(fileNamePer, ios::out | ios::app);

	if (fileExist == false) {
		writeToFile << "#POV_id POV_value" << endl;
	}

	double sumPOVs = 0 ;
	for(int i=0; i < POVs.size(); i++){
		sumPOVs += POVs[i];
	}


	for(int i=0; i < POVs.size(); i++){
		writeToFile << POVs[i] << "\t";
		writeToFilePer << (POVs[i]/sumPOVs)*100 << "\t";
	}
	writeToFile  << endl;
	writeToFilePer  << endl;

	writeToFile.close();
	writeToFilePer.close();
}

/*!****************************************************************************/
/*!****************************************************************************/
void PODICalc::saveConservedPOVsToFile_GrfFormat(int& numConsrvPOV,
		double& consrvPOVPercentage, int stepValueID, DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	string& resultName = problemData->getString("PODIResultName");
	int& nDOFPerNode = problemData->getInt("PODIDofPerNode");

	string fileName = "NumConservedPOVs_" + resultName + ".grf";
	string fileNamePer = "conservedPOVs_Perctage_" + resultName + ".grf";

	// Check if file exist
	bool fileExist = true;
	bool fileExistPer = true;

	struct stat buffer;
	if (stat(fileName.c_str(), &buffer) != 0) { // File exists
		fileExist = false;
	}
	if (stat(fileNamePer.c_str(), &buffer) != 0) { // File exists
		fileExistPer = false;
	}

	// Create/Open file
	ofstream writeToFile(fileName, ios::out | ios::app);
	ofstream writeToFilePer(fileNamePer, ios::out | ios::app);

	// Write header line if file has been created
	if (fileExist == false) {
		writeToFile << "#Iteration_ID POV_value" << endl;
	}
	if (fileExist == false) {
		writeToFilePer << "#Iteration_ID POV_energy" << endl;
	}

	// Write content to files
	writeToFile << stepValueID << "\t" << numConsrvPOV << endl;
	writeToFilePer << stepValueID << "\t" << consrvPOVPercentage << endl;

	// Close files access
	writeToFile.close();
	writeToFilePer.close();
}


