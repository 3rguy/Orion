#ifndef GaussPointSets_h_
#define GaussPointSets_h_

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"

struct GaussPointSet {

  dbMatrix coord;
  dbVector weight;

  static const int Gauss_Point_Positioning = 1; // 0: Zienkiewcz, 1: GiD

  virtual ~GaussPointSet() {};

};

struct GaussSetLine1 : virtual public GaussPointSet {

  GaussSetLine1();
  ~GaussSetLine1() {};

};

struct GaussSetLine2 : virtual public GaussPointSet {

  GaussSetLine2();
  ~GaussSetLine2() {};

};

struct GaussSetLine3 : virtual public GaussPointSet {

  GaussSetLine3();
  ~GaussSetLine3() {};

};

struct GaussSetLine4 : virtual public GaussPointSet {

  GaussSetLine4();
  ~GaussSetLine4() {};

};

struct GaussSetLine5 : virtual public GaussPointSet {

  GaussSetLine5();
  ~GaussSetLine5() {};

};

/***********************************************************************/
struct GaussSetTria1 : virtual public GaussPointSet {

  GaussSetTria1();
  ~GaussSetTria1() {};

};

struct GaussSetTetra1 : virtual public GaussPointSet {

  GaussSetTetra1();
  ~GaussSetTetra1() {};

};

/***********************************************************************/

struct GaussSetTria3 : virtual public GaussPointSet {

  GaussSetTria3();
  ~GaussSetTria3() {};

};

struct GaussSetTetra4 : virtual public GaussPointSet {

  GaussSetTetra4();
  ~GaussSetTetra4() {};

};

/***********************************************************************/

struct GaussSetTria4 : virtual public GaussPointSet {

  GaussSetTria4();
  ~GaussSetTria4() {};

};

struct GaussSetTetra5 : virtual public GaussPointSet {

  GaussSetTetra5();
  ~GaussSetTetra5() {};

};

/***********************************************************************/
struct GaussSetRect1 : virtual public GaussPointSet {

  GaussSetRect1();
  ~GaussSetRect1() {};

};

struct GaussSetCube1 : virtual public GaussPointSet {

  GaussSetCube1();
  ~GaussSetCube1() {};

};

/***********************************************************************/
struct GaussSetRect4 : virtual public GaussPointSet {

  GaussSetRect4();
  ~GaussSetRect4() {};

};

struct GaussSetCube8 : virtual public GaussPointSet {

  GaussSetCube8();
  ~GaussSetCube8() {};

};

/***********************************************************************/
struct GaussSetRect9 : virtual public GaussPointSet {

  GaussSetRect9();
  ~GaussSetRect9() {};

};

struct GaussSetCube27 : virtual public GaussPointSet {

  GaussSetCube27();
  ~GaussSetCube27() {};

};

/***********************************************************************/
struct GaussSetRect16 : virtual public GaussPointSet {

  GaussSetRect16();
  ~GaussSetRect16() {};

};

struct GaussSetCube64 : virtual public GaussPointSet {

  GaussSetCube64();
  ~GaussSetCube64() {};

};

/***********************************************************************/
struct GaussSetRect25 : virtual public GaussPointSet {

  GaussSetRect25();
  ~GaussSetRect25() {};

};

struct GaussSetCube125 : virtual public GaussPointSet {

  GaussSetCube125();
  ~GaussSetCube125() {};

};

#endif
