#ifndef SurfaceElementTemplates_h_
#define SurfaceElementTemplates_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "ElementTemplate.h"
#include "Particle.h"

#include "mpi.h"
#include "petscsys.h"


struct Tria3ElementTemplate : public ElementTemplate {

  ~Tria3ElementTemplate() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);

  Tria3ElementTemplate();

  // Calculate the transformation and metric factor for the integration
  // of a surface element's area.
  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile);

};

/***********************************************************************/

struct Rect4ElementTemplate : public ElementTemplate {

  ~Rect4ElementTemplate() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);
  double N4(dbVector& coords);

  Rect4ElementTemplate();

  // Calculate the transformation and metric factor for the integration
  // of a surface element's area.
  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile);

};

/***********************************************************************/

struct Rect9ElementTemplate : public ElementTemplate {

  // ordering: (0,0), (-1,1),  (-1,-1), (1,-1), (0,1), (-1,0), (0,-1), 
  //           (1,0), (0,0)
  ~Rect9ElementTemplate() {};

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

  Rect9ElementTemplate();

  // Calculate the transformation and metric factor for the integration
  // of a surface element's area.
  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile);

};

#endif
