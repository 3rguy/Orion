/*
 * PostProcessing.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef POSTPROCESSING_H_
#define POSTPROCESSING_H_

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

#include <time.h>
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
#include "FEMGeometryExt.h"
#include "Model.h"

class PostProcessing: virtual public Model{
public:

//	PostProcessing(Data& yourData,GridNodes* yourGrid);
	virtual ~PostProcessing(){};

	void setData(Data& yourData){myData = yourData;};
	Data& getData(){return myData;};

	void setGrid(GridNodes* yourGrid){yourGrid = myGrid;};
	GridNodes* getGrid(){return myGrid;};

	/*!************************************************************************/
	/*!************************************************************************/
	virtual void defaultPostProcessing(Database& myDatabase, DataContainer* problemData,
				InputFileData* InputData, ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	void saveResultInGridNodesFormat(DataContainer* problemData,
			std::vector<Data*>& dataList, InputFileData* InputData,
			ofstream& logFile);

	/*!************************************************************************/
	/*!************************************************************************/
	virtual void additionalPostProcessingFunctions(DataContainer* problemData,
			InputFileData* InputData, ofstream& logFile){};
};

#endif /* POSTPROCESSING_H_ */
