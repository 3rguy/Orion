// Calculate the window function ordinate of each particle supporting at a
// given point.

#ifndef WindowFunctionSet_h_
#define WindowFunctionSet_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"

class WindowFunctionSet {

  public:

    virtual ~WindowFunctionSet() {};

    virtual dbVector& getWindowFuncs() = 0;

    virtual dbVector& getXDerivWinFuncs() = 0;
    virtual dbVector& getYDerivWinFuncs() = 0;
    virtual dbVector& getZDerivWinFuncs() = 0;

    virtual dbVector& getXXDerivWinFuncs() = 0;
    virtual dbVector& getYYDerivWinFuncs() = 0;
    virtual dbVector& getZZDerivWinFuncs() = 0;
    virtual dbVector& getXYDerivWinFuncs() = 0;
    virtual dbVector& getYZDerivWinFuncs() = 0;
    virtual dbVector& getZXDerivWinFuncs() = 0;

};

#endif
