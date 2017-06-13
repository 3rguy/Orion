/*
 * MicroRVE.h
 *
 *  Created on: 14 Mar 2017
 *      Author: rama
 */

#ifndef MICRORVE_H_
#define MICRORVE_H_

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

#include "defs.h"
#include "DataContainer.h"
#include "InputFileData.h"
#include "commonTypedefs.h"
#include "commonFunctions.h"

#include "Database.h"
#include "Data.h"
#include "PreProcessing.h"
#include "ROMCalculation.h"
#include "PostProcessing.h"

class MicroRVE: public PreProcessing,public ROMCalculation,public PostProcessing {
public:
	MicroRVE(Database& myDB, DataContainer* problemData,
			 InputFileData* InputData, ofstream& logFile);

	~MicroRVE(){};

	/*!************************************************************************/
	/*!************************************************************************/
	void preProcessing(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void defaultPreProcessing(Database& myDatabase,
		   DataContainer* problemData, InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void loadSelectedData(DataContainer* problemData, Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void preROMCalculationFunctions(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile){};

	/*!************************************************************************/
	/*!************************************************************************/
	void reducedCalculation(DataContainer* problemData,	InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void postProcessing(DataContainer* problemData,	InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
//	void defaultPostProcessing(Database& myDatabase, DataContainer* problemData,
//			InputFileData* InputData, ofstream& logFile);
	void additionalPostProcessingFunctions(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile){};
};

#endif /* MICRORVE_H_ */
