/*
 * ErrorCalculation.h
 *
 *  Created on: Aug 19, 2015
 *      Author: rama
 */

#ifndef CALCULATION_ERRORCALCULATION_H_
#define CALCULATION_ERRORCALCULATION_H_


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

class ErrorCalculation: virtual public Model{
public:
	ErrorCalculation();
	virtual ~ErrorCalculation();
};

#endif /* CALCULATION_ERRORCALCULATION_H_ */
