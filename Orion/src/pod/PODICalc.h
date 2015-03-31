/*
 * PODICalc.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef PODICALC_H_
#define PODICALC_H_


#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>
#include <numeric>

#include <sys/stat.h>

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "PODCalc.h"
#include "defs.h"
#include "Database.h"
#include "Data.h"
#include "ParticleX.h"
//#include "ParticleDistributionX.h"
#include "InputFileData.h"
#include "GridNodes.h"
#include "Interpolation.h"

class PODICalc
{
public:

	PODICalc(InputFileData* InputData, ofstream& logFile);

	dbVector& getPODInterpolants(){return PODInterpolants;};
	void setPODInterpolants(dbVector interpolants){PODInterpolants = interpolants;};

	/*!************************************************************************/
	/*!************************************************************************/
	PODICalc(dbVector& myParameters, dbVector& parameterRadii,
			intVector& supportDataID, dbMatrix& dataParametersList,
			vector<dbMatrix>& displacementList,
			dbMatrix& resultingDisplacementMatrix, DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void interpolantsCalc(dbVector& myParameters, dbMatrix& dataParametersList,
			dbVector& parameterRadii, dbVector& interpolants,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void rearrangeDisplacementMatrix(vector<dbMatrix>& displacementList,
			vector<dbMatrix>& rearrangeDisplacementList,
			InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void PODInterpolation(vector<dbMatrix>& rearrangedisplacementList,
			dbVector& interpolants, dbMatrix& resultingDisplacementMatrix,
			DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void PODInterpolationEnhanced(vector<dbMatrix>& rearrangedisplacementList,
			dbVector& interpolants, dbMatrix& resultingDisplacementMatrix,
			DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void savePOMsToFile_ResFormat(dbMatrix& POMs, int stepValueID,
			DataContainer* problemData,	InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void savePOVsToFile_ResFormat(dbVector& POVs, int stepValueID,
			DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);



private:

	dbVector PODInterpolants;

};

#endif /* PODICALC_H_ */
