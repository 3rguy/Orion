// Calculate for a given point a shapefunction set. - one for each
// particle.

#ifndef ShapefunctionSet_h_
#define ShapefunctionSet_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "petscsys.h"

class ShapefunctionSet {

  public:
    ShapefunctionSet() {};
    virtual ~ShapefunctionSet() {};

    dbVector& getShapefunctions() { return shapefunctions; };
    dbMatrix& getFirstDerivShapes() { return firstDerivShapes; };
    dbMatrix& getSecondDerivShapes() { return secondDerivShapes; };

  protected:
    dbVector shapefunctions;
    dbMatrix firstDerivShapes;
    dbMatrix secondDerivShapes;

 private:

    // Calculate at a gauss point for all its supporting particles their 
    // shape functions.
    //virtual void calcShapes() {};
 
};

#endif
