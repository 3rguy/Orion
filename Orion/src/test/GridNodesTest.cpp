/*
 * GridNodesTest.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: rama
 */

#include <GridNodesTest.h>

GridNodesTest::GridNodesTest(InputFileData* InputData, ofstream& logFile){

	int GridNodesTest_choice = InputData->getValue("GridNodesTestMode");

	switch(GridNodesTest_choice){
	case 1:
	{
		DOFStandardisation(InputData,logFile);
		break;

	}
	default:
		cout << "GridNodesTestMode is invalid. Check input file." << endl;
		logFile << "GridNodesTestMode is invalid. Check input file." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodesTest::DOFStandardisation(InputFileData* InputData, ofstream& logFile){

	string templateFolder = "gridNodes/";
	string sourceFolder = "source/";

	Data sData, tData;
	sData.setFolderName(sourceFolder);
	tData.setFolderName(templateFolder);

	DataContainer* problemData = new DataContainer();

	// Setup the template grid nodes
	initDOFStandardisation(problemData,InputData,logFile);


	logFile << "-----------------------------" << endl;
	logFile << "Loading Data file: " << sData.getFolderName() << endl;
	logFile << "-----------------------------" << endl;

	cout << "-----------------------------" << endl;
	cout << "Loading Data file: " << sData.getFolderName() << endl;
	cout << "-----------------------------" << endl;

	cout << "Reading mesh.dat file ... ";
	sData.readMeshDataFile(InputData,logFile);
	cout << "OK" << endl;

	// Read result from file
	cout << "Reading fem.res file ... ";
	sData.readResultFile(InputData, logFile);
	cout << "OK" << endl;

	// Standardise the results
	ResultStandardisation(sData,problemData,InputData,logFile);


	logFile << "-----------------------------" << endl;
	logFile << "Loading Data file: " << sData.getFolderName() << endl;
	logFile << "-----------------------------" << endl;

	cout << "-----------------------------" << endl;
	cout << "Loading Data file: " << sData.getFolderName() << endl;
	cout << "-----------------------------" << endl;

	cout << "Reading mesh.dat file ... ";
	tData.readMeshDataFile(InputData,logFile);
	cout << "OK" << endl;

	tData.getResultList() = sData.getResultList();
	tData.getResultDOFList() = sData.getResultDOFList();
	tData.getStepValueVec() = sData.getStepValueVec();
	tData.getResultNameList() = sData.getResultNameList();

	cout << "Writting to fem.msh file ... ";
	string fileName = tData.getFolderName() + "fem.msh";
	tData.writeMeshToMSHFile(fileName,InputData,logFile);
	cout << "OK" << endl;

	cout << "Writting to fem.res file ... ";
	tData.saveResultsToFile(InputData,logFile);
	cout << "OK" << endl;
}

/*!****************************************************************************/
/*!****************************************************************************/
void GridNodesTest::initDOFStandardisation(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	cout << "DOF Standardisation" << endl;
	logFile << "******* DOF Standardisation *******" << endl;

	//! --------------------------------------------------------------------
	//! Setup the reference grid and define the nodes
	cout << "Loading up benchmark grid" << endl;
	logFile << "Loading up benchmark grid" << endl;

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
void GridNodesTest::ResultStandardisation(Data& sData,
		DataContainer* problemData,InputFileData* InputData, ofstream& logFile){

	//! Interpolate each displacement field on the reference grid
		dbMatrix standardResult;
		string choice;

		// Loading input parameters
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

		vector<string>& resultNameList = sData.getResultNameList();

		logFile << "Iterpolants setup" << endl;
		cout << "Iterpolants setup" << endl;
		myGrid->interpolantSetup(sData, InputData, logFile);


		logFile << "Iterpolanting results ... ";
		cout << "Iterpolanting results ... ";
		for (int i = 0; i < resultNameList.size(); i++) {

			// assign displacement field to particles
			sData.assignResultToParticles(resultNameList[i].c_str(), InputData,
					logFile);

			standardResult = myGrid->interpolateResultOnGridPoint(sData,
					InputData, logFile);

			myGrid->resetNodesStepDOFMat();

			sData.getMeshData()->writeMeshFile("fem2.msh",InputData, logFile);

			sData.deleteResult(resultNameList[i].c_str());

			sData.setResult(resultNameList[i].c_str(), standardResult);
		}

		myGrid->resetNodes();

		logFile << "OK" << endl;
		cout << "OK" << endl;

}
