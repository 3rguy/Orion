/*
 * DOFStandardisation.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef DOFSTANDARDISATION_H_
#define DOFSTANDARDISATION_H_

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

class DOFStandardisation {
public:
	DOFStandardisation(Data& yourData);
	~DOFStandardisation(){};

	void setData(Data& yourData){myData = yourData;};
	Data& getData(){return myData;};

	void setGrid(GridNodes* yourGrid){yourGrid = myGrid;};
	GridNodes* getGrid(){return myGrid;};

	void initDOFStandardisation(DataContainer* problemData,Database& myDatabase,
			InputFileData* InputData, ofstream& logFile);

	void standardiseResultDOF(DataContainer* problemData, Data& mainData,
			InputFileData* InputData, ofstream& logFile);

	void coordinateSetup(Data& sData, dbMatrix& maxCoordRange,
			InputFileData* InputData, ofstream& logFile);

private:
	Data myData;
	GridNodes* myGrid;
};

#endif /* DOFSTANDARDISATION_H_ */
