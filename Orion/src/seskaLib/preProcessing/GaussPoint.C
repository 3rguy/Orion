// Stores all properties of a single Gauss point.

#include "GaussPoint.h"


// micro boundary conditions

blVector& GaussPoint::getMicroBoundDOF(int ID) { 

  if(ID < microBoundDOF.size())
    return microBoundDOF[ID]; 

  else {
    microBoundDOF.resize(ID+1);
    return microBoundDOF[ID]; 
  }

}

dbVector& GaussPoint::getMicroBoundConds(int ID) { 

  if(ID < microBoundConds.size())
    return microBoundConds[ID]; 

  else {
    microBoundConds.resize(ID+1);
   return microBoundConds[ID]; 
  }

}

dbVector& GaussPoint::getDeltaMicroBoundConds(int ID) { 

  if(ID < deltaMicroBoundConds.size())
    return deltaMicroBoundConds[ID]; 

  else {
    deltaMicroBoundConds.resize(ID+1);
    return deltaMicroBoundConds[ID]; 
  }

}

dbVector& GaussPoint::getInitialMicroBoundConds(int ID) { 

  if(ID < initialMicroBoundConds.size())
    return initialMicroBoundConds[ID]; 

  else {
    initialMicroBoundConds.resize(ID+1);
    return initialMicroBoundConds[ID]; 
  }

}

// void& GaussPoint::getMicroSpaces() { 

//   return mSpaces;

// }

