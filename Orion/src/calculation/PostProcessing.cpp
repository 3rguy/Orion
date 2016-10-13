/*
 * PostProcessing.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <PostProcessing.h>

void PostProcessing::defaultPostProcessing(Database& myDatabase,
		DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile) {

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

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

//	InputData->setValue("gridNodesType", 2);

	if (InputData->getValue("standardiseDOF") == 1) {

		std::chrono::time_point<std::chrono::system_clock> start_stdDOF, end_stdDOF;
		start_stdDOF = std::chrono::system_clock::now();

		if (InputData->getValue("gridNodesType") == 1) {

			FEMGeometryExt* gridFEMGeo = myGrid->getGridGeometry(logFile);

			string oldFolderName = myData.getFolderName();
			myData.getFolderName() = string("postProc");
			vector<Data*> dataList(1);
			dataList.at(0) = &myData;

//			if (InputData->getValue("gridNodesResultPlot") == 1)
//				saveResultInGridNodesFormat(problemData, dataList, InputData,
//						logFile);

			myData.getFolderName() = oldFolderName;

			// -----------------------------------------------------------------
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

				myData.assignResultToParticles(resultNameList[i].c_str(),
						InputData, logFile);

				dbMatrix resultMat = myGrid->interpolateResultOnGridPoint(
						myData, InputData, logFile);

				myGrid->resetNodesStepDOFMat();

				logFile << "resultNameList: " << resultNameList[i] << endl;
				printMatrix(resultMat, "resultMat", logFile);

				myData.deleteResult(resultNameList[i].c_str());
				myData.setResult(resultNameList[i].c_str(), resultMat);

			}

			myData.setMeshData(myGrid->getGridGeometry(logFile));
			myGrid->setGridGeometry(gridFEMGeo);

//			myData.saveResultsToFile(InputData, logFile);
//			system("cp -f fem_orion.res fem_orion_map.res");

			myData.getMeshData()->writeMeshFile("fem_orion.msh", InputData,
					logFile);
		}

		else if (InputData->getValue("gridNodesType") == 2) {

			FEMGeometryExt* gridFEMGeo = myGrid->getGridGeometry(logFile);

			if (InputData->getValue("gridNodesResultPlot") == 1){
				string oldFolderName = myData.getFolderName();
				myData.getFolderName() = string("postProc");
				vector<Data*> dataList(1);
				dataList.at(0) = &myData;

				saveResultInGridNodesFormat(problemData, dataList, InputData,
						logFile);

				myData.getFolderName() = oldFolderName;

			}

			// ---------------------------------------------------------------------
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

				myData.assignResultToParticles(resultNameList[i].c_str(),
						InputData, logFile);

				dbMatrix resultMat = myGrid->interpolateResultOnGridPoint(
						myData, InputData, logFile);

				myGrid->resetNodesStepDOFMat();

				myData.deleteResult(resultNameList[i].c_str());
				myData.setResult(resultNameList[i].c_str(), resultMat);

			}

			myData.setMeshData(myGrid->getGridGeometry(logFile));
			myGrid->setGridGeometry(gridFEMGeo);

			myData.saveResultsToFile(InputData, logFile);
			system("cp -f fem_orion.res fem_orion_map.res");

			myData.delMeshData();
			myData.readMeshDataFile(InputData, logFile);
			myData.getMeshData()->writeMeshFile("fem_orion.msh", InputData,
					logFile);
		} else if (InputData->getValue("gridNodesType") == 3) {

			FEMGeometryExt* gridFEMGeo = myGrid->getGridGeometry(logFile);

//			if (InputData->getValue("gridNodesResultPlot") == 1){
//				string oldFolderName = myData.getFolderName();
//				myData.getFolderName() = string("postProc");
//				vector<Data*> dataList(1);
//				dataList.at(0) = &myData;
//
//				saveResultInGridNodesFormat(problemData, dataList, InputData,
//						logFile);
//
//				myData.getFolderName() = oldFolderName;
//
//			}

			// ---------------------------------------------------------------------
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

				myData.assignResultToParticles(resultNameList[i].c_str(),
						InputData, logFile);

				dbMatrix resultMat = myGrid->interpolateResultOnGridPoint(
						myData, InputData, logFile);

				myGrid->resetNodesStepDOFMat();

				myData.deleteResult(resultNameList[i].c_str());
				myData.setResult(resultNameList[i].c_str(), resultMat);

			}

			myData.setMeshData(myGrid->getGridGeometry(logFile));
			myGrid->setGridGeometry(gridFEMGeo);

			myData.saveResultsToFile(InputData, logFile);
			system("cp -f fem_orion.res fem_orion_map.res");

			myData.delMeshData();
			myData.readMeshDataFile(InputData, logFile);
			myData.getMeshData()->writeMeshFile("fem_orion.msh", InputData,
					logFile);
		} else if (InputData->getValue("gridNodesType") == 4) {

			// ---------------------------------------------------------------------
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
			myGrid->interpolantSetup(myData, InputData, logFile);
			// ---------------------------------------------------------------------

			// Interpolating all results + deleting previous ones + saving the new ones
			vector<string>& resultNameList = myData.getResultNameList();
			for (int i = 0; i < resultNameList.size(); i++) {

				// assign displacement field to particles
				myData.assignResultToParticles(resultNameList[i].c_str(),
						InputData, logFile);

				// Interpolate result on gridNodes
				dbMatrix standardResult = myGrid->interpolateResultOnGridPoint(
						myData, InputData, logFile);

				// Reset all dofs in myGrid
				myGrid->resetNodesStepDOFMat();

				// Delete previous results and save the new one
				myData.deleteResult(resultNameList[i].c_str());
				myData.setResult(resultNameList[i].c_str(), standardResult);
			}

			myGrid->resetNodes();

			myData.saveResultsToFile(InputData, logFile);
			system("cp -f fem_orion.res fem_orion_map.res");

			myData.getMeshData()->writeMeshFile("fem_orion.msh", InputData,
					logFile);
		}

		end_stdDOF = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds_stdDOF = end_stdDOF - start_stdDOF;

		cout << "DOF standardisation time: " << elapsed_seconds_stdDOF.count() << " sec"
		     << endl;
		logFile << "DOF standardisation time: " << elapsed_seconds_stdDOF.count() << " sec"
				<< endl;

	}

	myData.saveResultsToFile(InputData, logFile);
	myData.plotPostProcessGraph(InputData, logFile);

	additionalPostProcessingFunctions(problemData, InputData, logFile);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Total postprocessing time: " << elapsed_seconds.count() << " sec"
	     << endl;
	logFile << "Total postprocessing time: " << elapsed_seconds.count() << " sec"
			<< endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
void PostProcessing::saveResultInGridNodesFormat(DataContainer* problemData,
		vector<Data*>& dataList, InputFileData* InputData, ofstream& logFile) {

	using namespace std;

	Data* gridData = new Data();
	system("cd gridNodes");
//gridData->setFileName(string("gridNodes/"));
	int old_choice = (int) InputData->getValue("FEMGeometrySetupType");
	InputData->setValue("FEMGeometrySetupType", 1);
	gridData->setFolderName(string("gridNodes/"));
	gridData->readMeshDataFile(InputData, logFile);
	InputData->setValue("FEMGeometrySetupType", old_choice);

	string filename;

	for (int i = 0; i < dataList.size(); i++) {

		gridData->getResultNameList() = dataList.at(i)->getResultNameList();
		gridData->getResultDOFList() = dataList.at(i)->getResultDOFList();
		gridData->getResultList() = dataList.at(i)->getResultList();
		gridData->getStepValueVec() = dataList.at(i)->getStepValueVec();

		filename = "gridNode_" + dataList.at(i)->getFolderName();

		replace(filename.begin(), filename.end(), '/', '_');

		string resfile = filename + ".res";
		gridData->saveAllResultsToFile_res_format(resfile.c_str(), InputData,
				logFile);

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

