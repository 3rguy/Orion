/*
 * Plasticity.h
 *
 *  Created on: Nov 13, 2015
 *      Author: rama
 */

#ifndef PLASTICITY_H_
#define PLASTICITY_H_

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

class Plasticity :public PreProcessing,public ROMCalculation,public PostProcessing{
public:
	Plasticity(Database& myDB, DataContainer* problemData,
								 InputFileData* InputData, ofstream& logFile);
	~Plasticity(){};

	/*!************************************************************************/
	/*!************************************************************************/
	void preProcessing(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void loadSelectedData(DataContainer* problemData, Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

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
	void preROMCalculationFunctions(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!****************************************************************************/
	/*!****************************************************************************/
	void additionalPostProcessingFunctions(DataContainer* problemData,
				InputFileData* InputData, ofstream& logFile);

};

#endif /* PLASTICITY_H_ */
