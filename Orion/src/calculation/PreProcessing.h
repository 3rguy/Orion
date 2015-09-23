/*
 * PreProcessing.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef PREPROCESSING_H_
#define PREPROCESSING_H_

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
#include "GridNodes.h"

#include "Model.h"
#include "DatabaseQuery.h"
#include "Interpolation.h"
#include "DOFStandardisation.h"
#include "StepStandardisation.h"

class PreProcessing: virtual public Model{
public:

//	PreProcessing(Data& yourData,GridNodes* yourGrid);
	~PreProcessing(){};

	void setData(Data& yourData){myData = yourData;};
	Data& getData(){return myData;};

	void setGrid(GridNodes* yourGrid){yourGrid = myGrid;};
	GridNodes* getGrid(){return myGrid;};

	/*!************************************************************************/
	/*!************************************************************************/
	void initPreProcessing(Database& myDatabase, DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void defaultPreProcessing(Database& myDatabase, DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void extractParameters(DataContainer* problemData, Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	virtual void loadSelectedData(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void searchDatabase(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void standardiseDOF(Data& mainData, DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void preliminaryMLSInterpolantsCalc(Database& myDatabase,
			DataContainer* problemData,	InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void settingCommonResultNames(Database& myDatabase,
			DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	vector<string> findCommonResultName(DataContainer* problemData,
			Database& myDatabase, InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void stepStandardisation(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	virtual void preROMCalculationFunctions(DataContainer* problemData,
								InputFileData* InputData, ofstream& logFile){};

	/*!************************************************************************/
	/*!************************************************************************/
	void recordStandardisedSteps(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void savePreProcessedData(DataContainer* problemData,
								 InputFileData* InputData, ofstream& logFile);


protected:
	DOFStandardisation* DOFStandard;
	StepStandardisation* StepStandard;
};

#endif /* PREPROCESSING_H_ */
