// Stores all properties of a customized symmetric/asymmetric spline.

#ifndef CustomSpline_h_
#define CustomSpline_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"  
#include "commonTypedefs.h"


class CustomSpline {


  public:

    CustomSpline() : splineIsSet(false),splineIntegral(1.0) {};
    ~CustomSpline() {};

    bool splineIsSet;

    double& getSplineIntegral() { return splineIntegral; };
    dbMatrix3& getSplineCoefficients() { return splineCoefficients; };
    dbMatrix& getSplineKnots() { return splineKnots; };

  private:

    double splineIntegral;
    dbMatrix splineKnots;
    dbMatrix3 splineCoefficients;


};

#endif
