/*
 * LineElementTemplatesX.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef LINEELEMENTTEMPLATESX_H_
#define LINEELEMENTTEMPLATESX_H_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "ElementTemplateX.h"
#include "ParticleX.h"

#include "mpi.h"
#include "petscsys.h"


struct Line2ElementTemplateX : public ElementTemplateX {

  Line2ElementTemplateX();
  ~Line2ElementTemplateX() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);

  // Calculate the transformation and metric factor for the integration
  // of a line element's length.
  double metricFactor(std::vector<ParticleX>& particles,
		      intVector& nodes,
		      dbVector& coords,std::ofstream& logFile);

};

/***********************************************************************/

struct Line3ElementTemplateX : public ElementTemplateX {

  // ordering: (1), (-1), (0)
  Line3ElementTemplateX();
  ~Line3ElementTemplateX() {};

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);

  // Calculate the transformation and metric factor for the integration
  // of a line element's length.
  double metricFactor(std::vector<ParticleX>& particles,
		      intVector& nodes,
		      dbVector& coords,std::ofstream& logFile);

};

#endif /* LINEELEMENTTEMPLATESX_H_ */
