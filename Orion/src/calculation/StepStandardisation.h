/*
 * StepStandardisation.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef CALCULATION_STEPSTANDARDISATION_H_
#define CALCULATION_STEPSTANDARDISATION_H_

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
#include <climits>

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

class StepStandardisation {
public:

	StepStandardisation(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);
	~StepStandardisation(){};


	void setNumPhase(int& val){nPhase = val;};
	int& getNumPhase(){return nPhase;};

	/*!************************************************************************/
	/*!************************************************************************/
	void defineNormaliseSteps(DataContainer* problemData,
			InputFileData* InputData,ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void standardise(DataContainer* problemData, InputFileData* InputData,
															ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void stepStandardisationCalc(DataContainer* problemData,
			InputFileData* InputData,ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void initStepStandardisation( DataContainer* problemData,
			InputFileData* InputData,ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void findCommonStepsList(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void findSupports(dbVector& iPoint, dbMatrix& pointList,
			intVector& supportsList, dbVector& pointRadii,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void findCommonCardiacStepsList(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void createPhaseList(dbMatrix& stepList,
			vector<dbMatrix>& standardStepPhaseListMat,
			vector<dbMatrix>& nonStandardStepPhaseListMat,
			DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void getCardiacLoopDetails(dbVector& volumesList,dbVector& pressureList,
			dbVector& stepList, intVector& cPi,
			DataContainer* problemData,InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	dbVector getStandardSteps(dbVector& interpolants,DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	dbVector interpolateStandardSteps( const char* namePhaseList,
			dbVector& interpolants, DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	dbMatrix combinePhaseList(vector<dbMatrix> phaseList,
			DataContainer* problemData,InputFileData* InputData, ofstream& logFile);

	dbVector normStandStep;
	intMatrix phaseIndex;
	int nPhase;
};

#endif /* CALCULATION_STEPSTANDARDISATION_H_ */
