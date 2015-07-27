/*
 * GridNodesTest.h
 *
 *  Created on: Mar 31, 2015
 *      Author: rama
 */

#ifndef GRIDNODESTEST_H_
#define GRIDNODESTEST_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <cmath>

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
#include "Data.h"
#include "GridNodes.h"


class GridNodesTest {
public:

	GridNodesTest(InputFileData* InputData, ofstream& logFile);
	~GridNodesTest(){};

	/*!************************************************************************/
	/*!************************************************************************/
	void DOFStandardisation(InputFileData* InputData, ofstream& logFile);


	/*!************************************************************************/
	/*!************************************************************************/
	void initDOFStandardisation(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);


	/*!************************************************************************/
	/*!************************************************************************/
	void ResultStandardisation(Data& sData,	DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);




	GridNodes* myGrid;


};

#endif /* GRIDNODESTEST_H_ */
