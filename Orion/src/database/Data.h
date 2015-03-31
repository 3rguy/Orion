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

    void setFileName(string name){folderName=name;};
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

    void setGraphResultList(dbMatrix graphMatList){graphResultList = graphMatList;};
    dbMatrix& getGraphResultList(){return graphResultList;};

//    void printDispMatrix(ofstream& logFile);
//    void plotNode(int dof,ofstream& logFile);


    void readMeshDataFile(InputFileData* InputData, ofstream& logFile);
    //void deleteMeshDataFile(InputFileData* InputData, ofstream& logFile);

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

	void readGraphResultFile(InputFileData* InputData, ofstream& logFile);
	void readGraphFile_grfFormat(std::string& fileName, dbMatrix& grfMatrix,
			ofstream& logFile);
	void saveGraphResultsToFile(ofstream& logFile);
	void saveGraphResultsToFile_grf_format(std::string outputFileName,
			dbVector& graphVecOne,dbVector& graphVecTwo,ofstream& logFile);

    // Writing generated displacement matrix to files
//    void saveDispFile(ofstream& logFile);
//    void saveDispFile_res_format(const char* outputFileName,ofstream& logFile);

    void saveResultsToFile(ofstream& logFile);
    void saveAllResultsToFile_res_format(const char* outputFileName,
    		ofstream& logFile);

    void setAnchorPoint(int aPoint){anchorPoint = aPoint;};
    int& getAnchorPoint(){return anchorPoint;};

//	void assignDispToParticles(InputFileData* InputData, ofstream& logFile);
	void saveResultsToFile_res_format(const char* outputFileName,
				vector<string>& saveResultNameList, ofstream& logFile);
	void assignResultToParticles(const char* resultName,
			InputFileData* InputData, ofstream& logFile);

    void setSupportingNodesList(intMatrix matList){supportingNodesList = matList;};
    intMatrix& getSupportingNodesList(){return supportingNodesList;};

    void setResult(const char* resultName,dbMatrix result);
    dbMatrix& getResult(const char* resultName);
    map<string,dbMatrix>& getResultList(){return resultList;};
    void deleteResult(const char* resultName);

    void setResultDOF(const char* resultName,int& dof);
    int& getResultDOF(const char* resultName);
    map<string,int>& getResultDOFList(){return resultDofList;};
    void deleteResultDOF(const char* resultName);

    vector<string>& getResultNameList(){return resultNameList;};
    void setResultNameList(vector<string> nameList){resultNameList = nameList;};

    void calcCavityVolumes(InputFileData* InputData,ofstream& logFile);


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
    vector<vector<string> > allResultsNameList;
    vector<string> resultNameList;
    dbVector step_value_vec;
    map<string,dbMatrix> resultList;
    map<string,int> resultDofList;

    dbMatrix graphResultList;

    // DOFStandardisation stuff
    FEMGeometryExt* meshData;
    dbMatrix coordsList;
    int anchorPoint;
    intMatrix supportingNodesList;

};

#endif /* DATA_H_ */
