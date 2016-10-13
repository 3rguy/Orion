/*
 * DOFStandardisation.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <DOFStandardisation.h>

DOFStandardisation::DOFStandardisation(Data& yourData) {

	myData = yourData;
	myGrid = NULL;
}

/*!****************************************************************************/
/*!****************************************************************************/
void DOFStandardisation::initDOFStandardisation(DataContainer* problemData,
		Database& myDatabase, InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();


	// Preliminary stuff
	// -----------------
	if (InputData->getValue("gridNodesType") == 1) {

		//! --------------------------------------------------------------------
		//! For each data set: Assign displacement field to every particles and
		//! reset their coordinates with respect to the defined anchor point
		//!
		cout << "Loading up cube grid" << endl;
		logFile << "Loading up cube grid" << endl;

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
		logFile << "Coordinate setup of myData" << endl;
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

		//! --------------------------------------------------------------------
		//! Setup the reference grid and define the nodes

		myGrid = new GridNodes(InputData, logFile);
	}
	else{

		//! --------------------------------------------------------------------
		//! Setup the reference grid and define the nodes
		cout << "Loading up template grid" << endl;
		logFile << "Loading up template grid" << endl;

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
void DOFStandardisation::standardiseResultDOF(DataContainer* problemData, Data& mainData,
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
void DOFStandardisation::coordinateSetup(Data& sData, dbMatrix& maxCoordRange,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	vector<ParticleExt>& particles = sData.getMeshData()->getNodesVec();
	int anchorPoint = sData.getAnchorPoint();

	//!Extract coordinates of anchor point
	dbVector anchorCoords;
	if(anchorPoint == 0)
		anchorCoords = dbVector(particles[0].getCoords().size(),0);
	else if(anchorPoint > 0 && anchorPoint < particles.size()+1)
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
