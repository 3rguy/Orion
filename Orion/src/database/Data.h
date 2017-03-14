/*
 * Data.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef DATA_H_
#define DATA_H_


#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "defs.h"
#include "FEMGeometryExt.h"
#include "ParticleExt.h"
#include "InputFileData.h"
#include "DataContainer.h"

using namespace std;

class Data
{
public:
	Data();
    Data(dbVector paramVec);
//    ~Data(){delMeshData();};
    ~Data(){};

    void setNdofs(int val){nDofs=val;};
    int& getNdofs(void){return nDofs;};

    void setFolderName(string name){folderName=name;};
    string& getFolderName(void){return folderName;};

//    void setDisplacement(dbMatrix disp){dispMatrixList=disp;};
//    dbMatrix& getDisplacement(void){return dispMatrixList;};

    void setCoords(dbMatrix coords){coordsList=coords;};
    dbMatrix& getCoords(void){return coordsList;};

    void setId(int i){id=i;};
    int getId(void){return id;};

    void setParamValuesVec(dbVector paramVec){paramValues=paramVec;};
    dbVector& getParamValuesVec(){return paramValues;};

    void setStepValueVec(dbVector valueVec){step_value_vec=valueVec;};
    dbVector& getStepValueVec(){return step_value_vec;};

    void delMeshData(void);
    FEMGeometryExt* getMeshData(){return meshData;};
    void setMeshData(FEMGeometryExt* p){meshData = p;};


    void readMeshDataFile(InputFileData* InputData, ofstream& logFile);
    void readMeshDataFile(string& meshFileName,string& inputFileName,
    		InputFileData* InputData, ofstream& logFile);
    void readTransformedMeshDataFile(InputFileData* InputData, ofstream& logFile);
    void writeMeshToMSHFile(string& fileName,InputFileData* InputData,
    		ofstream& logFile);

    // Reading displacement matrix from files
    void readResultFile(InputFileData* InputData,ofstream& logFile);
    void readResultFile_txtFormat(ofstream& logFile);
//    void readResultFile_resFormat(ofstream& logFile);
    void readResultFile_resFormat_resultwise(string& inputFileName,ofstream& logFile);
	bool readResultFile_resFormat_specificResult(string& inputFileName,
			string& result_name, dbMatrix& resultMatrix, int& numDofs,
			dbVector& resultStepVector, ofstream& logFile);
	vector<vector<string> > readResultFile_resFormat_HeadersOnly(
			string& inputFileName, ofstream& logFile);

	void readResultFile_resFormat(string& inputFileName,ofstream& logFile);
	void readResultFile_resFormat_allResult(string& inputFileName,
			vector<string>& resultNameList,
			map<double, map<string,dbMatrix> >& resultsData,
			ofstream& logFile);

	void readGraphFile_grfFormat(std::string& fileName, dbMatrix& grfMatrix,
			ofstream& logFile);

	void readVentriclesPVGraphResultFile(InputFileData* InputData,
															ofstream& logFile);
	void readLeftVentriclePVGraphResultFile(InputFileData* InputData,
															ofstream& logFile);
	void readRightVentriclePVGraphResultFile(InputFileData* InputData,
															ofstream& logFile);

	void saveGraphResultsToFile_grf_format(std::string outputFileName,
			dbVector& graphVecOne,dbVector& graphVecTwo,ofstream& logFile);

//    Writing generated displacement matrix to files
//    void saveDispFile(ofstream& logFile);
//    void saveDispFile_res_format(const char* outputFileName,ofstream& logFile);

    void saveResultsToFile(InputFileData* InputData,ofstream& logFile);
    void saveResultsToFile_step(InputFileData* InputData,ofstream& logFile);
    void saveAllResultsToFile_res_format(const char* outputFileName,
    		InputFileData* InputData,ofstream& logFile);

    void setAnchorPoint(int aPoint){anchorPoint = aPoint;};
    int& getAnchorPoint(){return anchorPoint;};

//	void assignDispToParticles(InputFileData* InputData, ofstream& logFile);
	void saveResultsToFile_res_format(const char* outputFileName,
				vector<string>& saveResultNameList,InputFileData* InputData,
				ofstream& logFile);
	void assignResultToParticles(const char* resultName,
			InputFileData* InputData, ofstream& logFile);

    void setSupportingNodesList(intMatrix matList){supportingNodesList = matList;};
    intMatrix& getSupportingNodesList(){return supportingNodesList;};

    void setResult(const char* resultName,dbMatrix result);
    dbMatrix& getResult(const char* resultName);
    dbVector& getResultAcrossSteps(const char* resultName,int dofID);
    dbVector& getResultAcrossDOFs(const char* resultName,int step);
    map<string,dbMatrix>& getResultList(){return resultList;};
    void deleteResult(const char* resultName);

    void setResultDOF(const char* resultName,int& dof);
    int& getResultDOF(const char* resultName);
    map<string,int>& getResultDOFList(){return resultDofList;};
    void deleteResultDOF(const char* resultName);

    vector<string>& getResultNameList(){return resultNameList;};
    void setResultNameList(vector<string> nameList){resultNameList = nameList;};

    void calcCavityVolumes(InputFileData* InputData,ofstream& logFile);

    void calcLeftAndRightCavityVolumes(InputFileData* InputData, ofstream& logFile);
    void calcLeftCavityVolumes(InputFileData* InputData, ofstream& logFile);
    void calcLeftCavityVolumes_step(InputFileData* InputData, ofstream& logFile);
    void calcRightCavityVolumes(InputFileData* InputData, ofstream& logFile);
    void calcRightCavityVolumes_step(InputFileData* InputData, ofstream& logFile);

    dbVector& getLeftCavityVolumes(){return leftCavityVolumes;};
    dbVector& getRightCavityVolumes(){return rightCavityVolumes;};

    dbVector& getLeftCavityPressures(){return leftCavityPressures;};
    dbVector& getRightCavityPressures(){return rightCavityPressures;};

    void setInterpolants(dbVector& vals){interpolants = vals;};
    dbVector& getInterpolants(){return interpolants;};

    void setGraph(const char* graphName,dbVector graphResult);
    dbVector& getGraph(const char* graphName);
    void deleteGraph(const char* graphName);

    void syncCardiacTimeStepsAndResults(InputFileData* InputData, ofstream& logFile);

    void insertZeroResultFields(InputFileData* InputData, ofstream& logFile);
    void removeZeroResultFields(InputFileData* InputData, ofstream& logFile);

    void plotPostProcessGraph(InputFileData* InputData, ofstream& logFile);
    dbVector& getGraphData(int graphType, int node, int DOF, string& resultName,
    		InputFileData* InputData, ofstream& logFile);

private:

    // Unique identity of Data
    int id;

    // ?? -> no longer useful
    int nDofs;
    int nNodes;
    int nDims;

    // Data specific info
    dbVector paramValues;
    string folderName;

    // Results extracted from result file
    dbVector step_value_vec;

    vector<vector<string> > allResultsNameList;
    vector<string> resultNameList;

    map<string,dbMatrix> resultList;
    map<string,int> resultDofList;

    map<string,dbVector> graphList;

    // DOFStandardisation stuff
    FEMGeometryExt* meshData;
    dbMatrix coordsList;
    int anchorPoint;
    intMatrix supportingNodesList;

    dbVector leftCavityVolumes;
    dbVector leftCavityPressures;

    dbVector rightCavityVolumes;
    dbVector rightCavityPressures;

    dbVector interpolants;

};

#endif /* DATA_H_ */
