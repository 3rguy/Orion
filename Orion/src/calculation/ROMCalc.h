/*
 * ROMCalc.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef ROMCALC_H_
#define ROMCALC_H_

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

class ROMCalc
{
public:
	ROMCalc(InputFileData* InputData,ofstream& logFile);
	~ROMCalc(){};

	/*!************************************************************************/
	/*!************************************************************************/
	//! Extract the parameters for which the calculation will be carried out
	//! from the input file
	void extractParameters(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData,ofstream& logFile);

	void test_MLSInterpolantsCalc(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	//! Search for the influencing particle in the database for the
	//! interpolation process
	void searchDatabase(DataContainer* problemData, Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	void searchDatabase_two(DataContainer* problemData, Database& myDatabase,
				InputFileData* InputData, ofstream& logFile);

	void searchDatabase_three(DataContainer* problemData, Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	// Preprocessing:
	void preProcessing(DataContainer* problemData,Database& myDatabase,
				InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void ROMCalculationSet(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	vector<string> findCommonResultName(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	dbVector findCommonStepsList(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	void reducedOrderMethodCalc(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void findSupports(dbVector& iPoint, dbMatrix& pointList,
			intVector& supportsList, dbVector& parameterRadii,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void standardisationProcess(DataContainer* problemData,
			Database& myDatabase, InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void initDOFStandardisation(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	void coordinateSetup(Data& sData, dbMatrix& maxCoordRange,
			InputFileData* InputData, ofstream& logFile);

	void standardiseResultDOF(DataContainer* problemData,
			Data& myData, InputFileData* InputData, ofstream& logFile);


	/*!************************************************************************/
	/*!************************************************************************/
	void stepStandardisation(DataContainer* problemData,
			InputFileData* InputData,ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void postProcessing(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void saveResultInGridNodesFormat(DataContainer* problemData,
			std::vector<Data*>& dataList, InputFileData* InputData,
			ofstream& logFile);

private:
	Data myData;
	GridNodes* myGrid;

};

#endif /* ROMCALC_H_ */
