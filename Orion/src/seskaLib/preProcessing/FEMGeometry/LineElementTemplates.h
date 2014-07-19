#ifndef LineElementTemplates_h_
#define LineElementTemplates_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "ElementTemplate.h"
#include "Particle.h"

#include "mpi.h"
#include "petscsys.h"


struct Line2ElementTemplate : public ElementTemplate {

  Line2ElementTemplate();
  ~Line2ElementTemplate() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);

  // Calculate the transformation and metric factor for the integration
  // of a line element's length.
  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,std::ofstream& logFile);

};

/***********************************************************************/

struct Line3ElementTemplate : public ElementTemplate {

  // ordering: (1), (-1), (0)
  Line3ElementTemplate();
  ~Line3ElementTemplate() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);

  // Calculate the transformation and metric factor for the integration
  // of a line element's length.
  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,std::ofstream& logFile);

};

#endif
