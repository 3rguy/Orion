/*
 * GaussPointX.cpp
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#include "GaussPointX.h"

// micro boundary conditions

blVector& GaussPointX::getMicroBoundDOF(int ID) {

  if(ID < microBoundDOF.size())
    return microBoundDOF[ID];

  else {
    microBoundDOF.resize(ID+1);
    return microBoundDOF[ID];
  }

}

dbVector& GaussPointX::getMicroBoundConds(int ID) {

  if(ID < microBoundConds.size())
    return microBoundConds[ID];

  else {
    microBoundConds.resize(ID+1);
   return microBoundConds[ID];
  }

}

dbVector& GaussPointX::getDeltaMicroBoundConds(int ID) {

  if(ID < deltaMicroBoundConds.size())
    return deltaMicroBoundConds[ID];

  else {
    deltaMicroBoundConds.resize(ID+1);
    return deltaMicroBoundConds[ID];
  }

}

dbVector& GaussPointX::getInitialMicroBoundConds(int ID) {

  if(ID < initialMicroBoundConds.size())
    return initialMicroBoundConds[ID];

  else {
    initialMicroBoundConds.resize(ID+1);
    return initialMicroBoundConds[ID];
  }

}

// void& GaussPointX::getMicroSpaces() {

//   return mSpaces;

// }

