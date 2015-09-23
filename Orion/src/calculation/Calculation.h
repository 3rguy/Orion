/*
 * Calculation.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef CALCULATION_H_
#define CALCULATION_H_

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

#include "ModelSet.h"

class Calculation {
public:
	Calculation(InputFileData* InputData, ofstream& logFile);
	~Calculation(){};

};

#endif /* CALCULATION_H_ */
