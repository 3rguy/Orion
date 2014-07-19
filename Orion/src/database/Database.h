/*
 * Database.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef DATABASE_H_
#define DATABASE_H_


#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <cmath>

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "defs.h"
#include "Data.h"

using namespace std;

class Database{
public:
    Database(string& databaseName,ofstream& logFile);
    ~Database(){}
    void addEntry(Data entry);
    void deleteEntry(int id);
    void readDatabase(string& databaseName, vector<string>& paramNameVec,
    		dbMatrix& paramMatrix, intVector& anchorPointList,
    		vector<string> &fileNameVec, ofstream& logFile);
    int size(){return dataList.size();}
    int getNumOfParam(){return paramNamesVec.size();}
    vector<string> getParamNamesVec(){return paramNamesVec;}
    void printDatabase(ofstream& logFile);
    Data& getDataEntry(int id){return dataList[id];};
    Data& getDataId(int id){return dataList[findData(id)];}
    int findData(int val);

    vector<Data>& getDataList(){return dataList;}

    void importDisp(InputFileData* InputData,int val,ofstream& logFile);
    void printDataInfo(int val,ofstream& logFile);

private:
    vector<Data> dataList;
    vector<string> paramNamesVec;
};


#endif /* DATABASE_H_ */
