/*
 * DatabaseQuery.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef DATABASEQUERY_H_
#define DATABASEQUERY_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>
#include <map>

#include <chrono>
#include <ctime>

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "DataContainer.h"
#include "defs.h"
#include "Database.h"
#include "Data.h"
#include "ParticleExt.h"
#include "ParticleDistributionX.h"
#include "InputFileData.h"
#include "GridNodes.h"
#include "PODICalc.h"

class DatabaseQuery {
public:
	DatabaseQuery(DataContainer* problemData, Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	~DatabaseQuery(){};

	/*!************************************************************************/
	/*!************************************************************************/
	void databaseQueryAlgorithm_one(DataContainer* problemData,
			Database& myDatabase, InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void databaseQueryAlgorithm_two(DataContainer* problemData,
			Database& myDatabase, InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void databaseQueryAlgorithm_three(DataContainer* problemData,
			Database& myDatabase, InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void databaseQueryAlgorithm_four(DataContainer* problemData,
			Database& myDatabase, InputFileData* InputData, ofstream& logFile);
};

#endif /* DATABASEQUERY_H_ */
