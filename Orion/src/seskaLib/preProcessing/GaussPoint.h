// Stores all properties of a single gauss point.

#ifndef GaussPoint_h_
#define GaussPoint_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"  
#include "commonTypedefs.h"

#include "IntegrationPoint.h"
#include "MicroSpace.h"

class GaussPoint : public virtual IntegrationPoint {


 public:

  GaussPoint() {};
  ~GaussPoint() {};

  std::vector<MicroSpace>& getMicroSpaces() { return microSpaces; };

  // micro boundary conditions
  blVector& getMicroBoundDOF(int ID);
  dbVector& getMicroBoundConds(int ID);
  blVector& getMicroBoundDOF() { return getMicroBoundDOF(0); };
  dbVector& getMicroBoundConds() { return getMicroBoundConds(0); };
  dbVector& getDeltaMicroBoundConds(int ID);
  dbVector& getInitialMicroBoundConds(int ID);    

 private:  

  // micro boundary conditions
  blMatrix microBoundDOF;
  dbMatrix microBoundConds;
  dbMatrix deltaMicroBoundConds;
  dbMatrix initialMicroBoundConds;

  std::vector<MicroSpace> microSpaces;


};

#endif
