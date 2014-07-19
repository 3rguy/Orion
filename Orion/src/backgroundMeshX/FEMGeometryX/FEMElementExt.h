/*
 * FEMElementExt.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef FEMELEMENTEXT_H_
#define FEMELEMENTEXT_H_

#include <vector>

#include "commonTypedefs.h"
#include "FEMElement.h"

class FEMElementExt : public FEMElement {

  public:
	FEMElementExt(int usedDOF);
	FEMElementExt(FEMElement& element);
    ~FEMElementExt() {};

    void setMaterialID(int matID){this->getMaterialID() =matID;};

    void setNodes(intVector nodesVec){this->getNodes() = nodesVec;};

    // FEM element type.
    void setElemType(int value) {this->getElemType() = value;};
    void setElemOrder(int value) {this->getElemOrder() = value;};

    // Manipulating surface elements
    intVector& getSurfaceElems(){return surfaceElems;};
    void setSurfaceElems(intVector surfElemVec){surfaceElems = surfElemVec;};

  private:
    intVector surfaceElems;

};

#endif /* FEMELEMENTEXT_H_ */
