/*
 * GaussPointSetsX.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef GAUSSPOINTSETSX_H_
#define GAUSSPOINTSETSX_H_

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"

struct GaussPointSetX {

  dbMatrix coord;
  dbVector weight;

  virtual ~GaussPointSetX() {};

};

struct GaussSetLine1X : virtual public GaussPointSetX {

  GaussSetLine1X();
  ~GaussSetLine1X() {};

};

struct GaussSetLine2X : virtual public GaussPointSetX {

  GaussSetLine2X();
  ~GaussSetLine2X() {};

};

struct GaussSetLine3X : virtual public GaussPointSetX {

  GaussSetLine3X();
  ~GaussSetLine3X() {};

};

struct GaussSetLine4X : virtual public GaussPointSetX {

  GaussSetLine4X();
  ~GaussSetLine4X() {};

};

struct GaussSetLine5X : virtual public GaussPointSetX {

  GaussSetLine5X();
  ~GaussSetLine5X() {};

};

/***********************************************************************/
struct GaussSetTria1X : virtual public GaussPointSetX {

  GaussSetTria1X();
  ~GaussSetTria1X() {};

};

struct GaussSetTetra1X : virtual public GaussPointSetX {

  GaussSetTetra1X();
  ~GaussSetTetra1X() {};

};

/***********************************************************************/

struct GaussSetTria3X : virtual public GaussPointSetX {

  GaussSetTria3X();
  ~GaussSetTria3X() {};

};

struct GaussSetTetra4X : virtual public GaussPointSetX {

  GaussSetTetra4X();
  ~GaussSetTetra4X() {};

};

/***********************************************************************/

struct GaussSetTria4X : virtual public GaussPointSetX {

  GaussSetTria4X();
  ~GaussSetTria4X() {};

};

struct GaussSetTetra5X : virtual public GaussPointSetX {

  GaussSetTetra5X();
  ~GaussSetTetra5X() {};

};

/***********************************************************************/
struct GaussSetRect1X : virtual public GaussPointSetX {

  GaussSetRect1X();
  ~GaussSetRect1X() {};

};

struct GaussSetCube1X : virtual public GaussPointSetX {

  GaussSetCube1X();
  ~GaussSetCube1X() {};

};

/***********************************************************************/
struct GaussSetRect4X : virtual public GaussPointSetX {

  GaussSetRect4X();
  ~GaussSetRect4X() {};

};

struct GaussSetCube8X : virtual public GaussPointSetX {

  GaussSetCube8X();
  ~GaussSetCube8X() {};

};

/***********************************************************************/
struct GaussSetRect9X : virtual public GaussPointSetX {

  GaussSetRect9X();
  ~GaussSetRect9X() {};

};

struct GaussSetCube27X : virtual public GaussPointSetX {

  GaussSetCube27X();
  ~GaussSetCube27X() {};

};

/***********************************************************************/
struct GaussSetRect16X : virtual public GaussPointSetX {

  GaussSetRect16X();
  ~GaussSetRect16X() {};

};

struct GaussSetCube64X : virtual public GaussPointSetX {

  GaussSetCube64X();
  ~GaussSetCube64X() {};

};

/***********************************************************************/
struct GaussSetRect25X : virtual public GaussPointSetX {

  GaussSetRect25X();
  ~GaussSetRect25X() {};

};

struct GaussSetCube125X : virtual public GaussPointSetX {

  GaussSetCube125X();
  ~GaussSetCube125X() {};

};

#endif /* GAUSSPOINTSETSX_H_ */
