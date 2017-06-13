/*
 * Model.h
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#ifndef MODELSET_H_
#define MODELSET_H_

#include <fstream>

#include "defs.h"
#include "DataContainer.h"
#include "InputFileData.h"

#include "Database.h"
#include "Standard.h"
#include "Cardiac.h"
#include "Plasticity.h"
#include "MicroRVE.h"

class ModelSet{
public:
	ModelSet(DataContainer* problemData, InputFileData* InputData,
															ofstream& logFile);
	~ModelSet(){};

};

#endif /* MODELSET_H_ */
