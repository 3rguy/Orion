/*
 * Node.h
 *
 *  Created on: 03 Jul 2014
 *      Author: ritesh
 */

#ifndef NODE_H_
#define NODE_H_


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

class Node {
public:
	Node(dbVector coordinates);
	~Node(){}

	dbVector& getCoords(){return coords;}
	void setCoords(dbVector coordinates){coords=coordinates;}
	int getId(){return id;}
	void setId(int val){id=val;}
	dbVector& getInterpolants(){return interpolants;}
	void setInterpolants(dbVector interpolantList)
								{interpolants = interpolantList;}

	intVector& getSPtlcs(){return sPtcls;}
	void setSPtcls(intVector& ptcls){sPtcls = ptcls;}

	void setStepDOFMat(dbMatrix mat){mat = stepDOFMat;}
	dbMatrix& getStepDOFMat(){return stepDOFMat;}

	void reset();
	void resetStepDOFMat();

	int& getVolumeElement(){return volumeElement;}
	void setVolumeElement(int volumeElem){this->volumeElement = volumeElem;}

private:

	int id;
	dbVector coords;

	int volumeElement;

	intVector sPtcls;
	dbVector interpolants;

	dbMatrix stepDOFMat;

};

#endif /* NODE_H_ */
