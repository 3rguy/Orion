/*
 * ErrorCalc.h
 *
 *  Created on: Oct 21, 2014
 *      Author: rama
 */

#ifndef ERRORCALC_H_
#define ERRORCALC_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "defs.h"
#include "InputFileData.h"
#include "Data.h"
#include "DataContainer.h"

class ErrorCalc {

public:
	ErrorCalc(InputFileData* InputData, ofstream& logFile);
	~ErrorCalc(){};

	void readingErrorFile(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	void resErrorCalc(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	void grfErrorCalc(DataContainer* problemData, InputFileData* InputData,
			ofstream& logFile);

	void compareTimeSteps(Data* exactData, Data* approxData,
			InputFileData* InputData, ofstream& logFile);

	vector<string> findCommonDataResults(Data* exactData, Data* approxData,
			InputFileData* InputData, ofstream& logFile);

	void calculateErrors_res_Format(Data* exactData, Data* approxData,
			vector<string>& commonNameList, InputFileData* InputData,
			ofstream& logFile);

	void calculateErrors_grf_Format(dbMatrix& exactGrfMatrix,
			dbMatrix& approxGrfMatrix,InputFileData* InputData,
			ofstream& logFile);

	static double relativeL2norm(dbMatrix& exactMat, dbMatrix& approxMat,
			InputFileData* InputData, ofstream& logFile);
};

#endif /* ERRORCALC_H_ */
