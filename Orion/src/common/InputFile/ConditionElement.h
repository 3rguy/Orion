// Stores all properties of an element where loading or Dirichlet condition are applied
// Note: The current magnitude of the condition is computed as
// p(t) = lambda(t)*p0 where p0 denotes the reference value and lambda
// the current condition factor. This factor is incremented in time.

// tractionLoad, lineForceLoad, pointForceLoad, surfacePressureLoad,
// surfaceMomentLoad, surfaceElectricChargeLoad
// displacementConstraint, rotationConstraint, electricConstraint,
// depolarisationConstraint, stressConstraint, microConstraint

#ifndef ConditionElement_h_
#define ConditionElement_h_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"

class ConditionElement {

  public:

    ConditionElement() : ID(0),motherElemID(0) {};
    ~ConditionElement() {};

    int& getID() { return ID; };

    // Return the global ID of the corresponding mother element
    int& getMotherElementID() { return motherElemID; };

    // Return element nodes.
    intVector& getNodes() { return nodes; };

    // Return element Gauss points.
    intVector& getGaussPts() { return gaussPts; };
    
  private:

    int ID;

    intVector nodes;
    intVector gaussPts;
    int motherElemID; // surface or line elements belong to a volume element

};

#endif
