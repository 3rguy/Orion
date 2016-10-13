/*
 * SurfaceElementTemplatesX.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef SURFACEELEMENTTEMPLATESX_H_
#define SURFACEELEMENTTEMPLATESX_H_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "ElementTemplateX.h"
#include "ParticleX.h"

#include "mpi.h"
#include "petscsys.h"


struct Tria3ElementTemplateX : public ElementTemplateX {

  ~Tria3ElementTemplateX() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);

  Tria3ElementTemplateX();

  // Calculate the transformation and metric factor for the integration
  // of a surface element's area.
  double metricFactor(std::vector<ParticleX>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile);

};

/***********************************************************************/

struct Rect4ElementTemplateX : public ElementTemplateX {

  ~Rect4ElementTemplateX() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);
  double N4(dbVector& coords);

  Rect4ElementTemplateX();

  // Calculate the transformation and metric factor for the integration
  // of a surface element's area.
  double metricFactor(std::vector<ParticleX>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile);

};

/***********************************************************************/

struct Rect9ElementTemplateX : public ElementTemplateX {

  // ordering: (0,0), (-1,1),  (-1,-1), (1,-1), (0,1), (-1,0), (0,-1),
  //           (1,0), (0,0)
  ~Rect9ElementTemplateX() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);
  double N4(dbVector& coords);
  double N5(dbVector& coords);
  double N6(dbVector& coords);
  double N7(dbVector& coords);
  double N8(dbVector& coords);
  double N9(dbVector& coords);

  Rect9ElementTemplateX();

  // Calculate the transformation and metric factor for the integration
  // of a surface element's area.
  double metricFactor(std::vector<ParticleX>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile);

};

#endif /* SURFACEELEMENTTEMPLATESX_H_ */
