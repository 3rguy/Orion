/*
 * DataContainer.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef DATACONTAINER_H_
#define DATACONTAINER_H_


#include "commonFunctions.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "commonTypedefs.h"
#include "defs.h"
#include "mpi.h"

using namespace std;

struct DataContainer {

public:

	DataContainer();
	DataContainer(std::ofstream& logFile);

	void setValue(const char* name,int val);
	bool findInt(const char* name);
	int& getInt(const char* name);
	void deleteInt(const char* name);

	void setValue(const char* name,double val);
	bool findDouble(const char* name);
	double& getDouble(const char* name);
	void deleteDouble(const char* name);

	void setValue(const char* name,intVector vec);
	bool findIntVector(const char* name);
	intVector& getIntVector(const char* name);
	void deleteIntVector(const char* name);

	void setValue(const char* name,dbVector vec);
	bool findDbVector(const char* name);
	dbVector& getDbVector(const char* name);
	void deleteDbVector(const char* name);

	void setValue(const char* name,intMatrix mat);
	bool findIntMatrix(const char* name);
	intMatrix& getIntMatrix(const char* name);
	void deleteIntMatrix(const char* name);

	void setValue(const char* name,dbMatrix mat);
	bool findDbMatrix(const char* name);
	dbMatrix& getDbMatrix(const char* name);
	void deleteDbMatrix(const char* name);

	void setValue(const char* name,vector<dbMatrix> mat);
	bool findDbMatrixVec(const char* name);
	vector<dbMatrix>& getDbMatrixVec(const char* name);
	void deleteDbMatrixVec(const char* name);

	void setValue(const char* name,string val);
	bool findString(const char* name);
	string& getString(const char* name);
	void deleteString(const char* name);

	void setValue(const char* name,vector<string> val);
	bool findStringVec(const char* name);
	vector<string>& getStringVec(const char* name);
	void deleteStringVec(const char* name);

	map<string,int> intList;
	map<string,double> doubleList;

	map<string,intVector> intVectorList;
	map<string,dbVector> dbVectorList;

	map<string,intMatrix> intMatrixList;
	map<string,dbMatrix> dbMatrixList;

	map<string,vector<dbMatrix> > dbMatrixVecList;

	map<string,string> stringList;
	map<string,vector<string>> stringVectorList;

};

#endif /* DATACONTAINER_H_ */
