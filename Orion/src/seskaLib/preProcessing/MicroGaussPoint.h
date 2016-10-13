// Stores all properties of a single micro Gauss point.

#ifndef MicroGaussPoint_h_
#define MicroGaussPoint_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "InputFileData.h"

#include "IntegrationPoint.h"

class MicroGaussPoint : public virtual IntegrationPoint {

  public:

    MicroGaussPoint() {};
    ~MicroGaussPoint() {};


  private:
  
    int microSpaceID;

};

#endif
