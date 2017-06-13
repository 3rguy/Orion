/*
 * ROMCalculation.h
 *
 *  Created on: Jul 29, 2015
 *      Author: rama
 */

#ifndef ROMCALCULATION_H_
#define ROMCALCULATION_H_

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
#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "InputFileData.h"
#include "DataContainer.h"

#include "Database.h"
#include "Data.h"

#include "PODICalc.h"
#include "Model.h"

class ROMCalculation: virtual public Model{
public:

//	ROMCalculation(Data& yourData);
	~ROMCalculation(){};

	void setData(Data& yourData){myData = yourData;};
	Data& getData(){return myData;};

	/*!************************************************************************/
	/*!************************************************************************/
	void calculation(Database& myDatabase, DataContainer* problemData,
				InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void reducedOrderCalculation(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	vector<string> findCommonResultName(DataContainer* problemData,
			Database& myDatabase, InputFileData* InputData, ofstream& logFile);

};

#endif /* ROMCALCULATION_H_ */
