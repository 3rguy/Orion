/*
 * PostProcessing.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <PostProcessing.h>

void PostProcessing::defaultPostProcessing(Database& myDatabase, DataContainer* problemData,
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

		myData.delMeshData();
		myData.readMeshDataFile(InputData, logFile);
		myData.getMeshData()->writeMeshFile("fem_orion.msh",InputData,logFile);
	}

	myData.saveResultsToFile(logFile);
	myData.plotPostProcessGraph(InputData,logFile);

	additionalPostProcessingFunctions(problemData,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void PostProcessing::saveResultInGridNodesFormat(DataContainer* problemData,
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

