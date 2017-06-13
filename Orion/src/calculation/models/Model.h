/*
 * Model.h
 *
 *  Created on: 09 Aug 2015
 *      Author: ritesh
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Database.h"
#include "Data.h"
#include "GridNodes.h"

class Model {
public:
	Model();
	~Model(){};

	Database myDatabase;
	Data myData;
	GridNodes* myGrid;

};

#endif /* MODEL_H_ */
