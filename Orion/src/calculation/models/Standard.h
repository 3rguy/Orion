/*
 * Standard.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef STANDARD_H_
#define STANDARD_H_

#include <fstream>

#include "defs.h"
#include "DataContainer.h"
#include "InputFileData.h"

#include "PreProcessing.h"
#include "ROMCalculation.h"
#include "PostProcessing.h"


class Standard:public PreProcessing,public ROMCalculation,public PostProcessing{
public:

	Standard(Database& myDB, DataContainer* problemData,
								 InputFileData* InputData, ofstream& logFile);
	~Standard(){};

	/*!************************************************************************/
	/*!************************************************************************/
	void preProcessing(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void reducedCalculation(DataContainer* problemData,	InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void postProcessing(DataContainer* problemData,	InputFileData* InputData,
			ofstream& logFile);

};

#endif /* STANDARD_H_ */
