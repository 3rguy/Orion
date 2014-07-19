// Stores all properties of a single point.

#ifndef Point_h_
#define Point_h_

#include <iostream>
#include <vector>

#include "commonTypedefs.h"


class Point {

  public:

    Point(int usedDims) { 
      coords = dbVector(usedDims); 
    };
    ~Point() {};
   
    dbVector& getCoords() { return coords; };
    intVector& getSupportPtcls() { return supportingPtcls; };
    dbVector& getShapeFuncs() { return shapeFuncs; };

    dbVector& getMLSvec() { return MLSvec; };
    dbMatrix& getMLSmat() { return MLSmat; };

  private:

    dbVector coords;
    intVector supportingPtcls;
    dbVector shapeFuncs;

    dbVector MLSvec;
    dbMatrix MLSmat;

};

#endif
