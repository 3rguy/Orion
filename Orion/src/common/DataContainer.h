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

	void setValue(const char* name,intVector vec);
	bool findIntVectorList(const char* name);
	intVector& getIntVector(const char* name);

	void setValue(const char* name,dbVector vec);
	bool findDbVectorList(const char* name);
	dbVector& getDbVector(const char* name);

	void setValue(const char* name,dbMatrix mat);
	bool findDbMatrixList(const char* name);
	dbMatrix& getDbMatrix(const char* name);

	void setValue(const char* name,vector<dbMatrix> mat);
	bool findDbMatrixVecList(const char* name);
	vector<dbMatrix>& getDbMatrixVec(const char* name);

	void setValue(const char* name,intMatrix mat);
	bool findIntMatrixList(const char* name);
	intMatrix& getIntMatrix(const char* name);


	map<string,intVector> intVectorList;
	map<string,dbVector> dbVectorList;

	map<string,intMatrix> intMatrixList;
	map<string,dbMatrix> dbMatrixList;

	map<string,vector<dbMatrix> > dbMatrixVecList;

};

#endif /* DATACONTAINER_H_ */
