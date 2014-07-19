/*
 * GaussPointX.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef GAUSSPOINTX_H_
#define GAUSSPOINTX_H_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"

#include "IntegrationPointX.h"
//#include "MicroSpace.h"

class GaussPointX : public virtual IntegrationPointX {


 public:

  GaussPointX() {};
  ~GaussPointX() {};

  //std::vector<MicroSpace>& getMicroSpaces() { return microSpaces; };

  // Global index datum manipulating
  //  void setGlobalID(int idx);
  //  int& getGlobalID() { return globalID; };
  //  int& getMaterialID() { return materialID; };

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

 // std::vector<MicroSpace> microSpaces;


};

#endif /* GAUSSPOINTX_H_ */
