/*
 * Node.cpp
 *
 *  Created on: 03 Jul 2014
 *      Author: ritesh
 */

#include "Node.h"

Node::Node(dbVector coordinates):sPtcls(intVector(0)),interpolants(dbVector(0)),
stepDOFMat(dbMatrix(0,dbVector(0))){

	volumeElement = -1;
	id = -1;
	coords=coordinates;

}

void Node::reset(){
	volumeElement = -1;
	sPtcls = intVector(0);
	interpolants = dbVector(0);
	this->resetStepDOFMat();
}


void Node::resetStepDOFMat(){
	stepDOFMat = dbMatrix(0,dbVector(0));
}

