#ifndef ElementTemplate_h_
#define ElementTemplate_h_

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonTypedefs.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"


struct ElementTemplate  {

  ~ElementTemplate() {};

  static const int Gauss_Point_Positioning = 1; // 0: Zienkiewcz, 1: GiD

  dbMatrix nodalCoords;
  dbVector nodalWeights;

  virtual double N(int func,dbVector& coords) { return 0.0; };
  virtual double dN(int func,int deriv,dbVector& coords) { return 0.0; };

  // Compute Jacobian dx_i/dL_j
  virtual void jacobian(std::vector<Particle>& particles,
			intVector& nodes,
			dbVector& coords,
			dbMatrix& jacobian,
			std::ofstream& logFile) {};


  // Calculate the transformation and metric factor for the integration
  // of an element's length, surface or volume.
  virtual double metricFactor(std::vector<Particle>& particles,
			      intVector& nodes,
			      dbVector& coords,
			      std::ofstream& logFile) 
  { return 0.0; };

  // Calculate the length of a edge of volume element to integrate
  // over this edge.
  virtual void lineMetric(std::vector<Particle>& particles,
			  intVector& nodes,
			  dbVector& coords,double& length,
			  std::ofstream& logFile) {};

  // Calculate surface normal vector on volume element's side and its 
  // length, which is the absolute value of the surface, in order to 
  // integrate over this volume element's side.
  virtual void surfaceMetric(std::vector<Particle>& particles,
			     intVector& nodes,
			     dbVector& coords,dbVector& surfaceNormal,
			     double& surface,std::ofstream& logFile)
  {};

  // Calculate the transformation and metric factor for the integration
  // of an element's volume.
  virtual double metricFactorVolume(std::vector<Particle>& particles,
				    intVector& nodes,
				    dbVector& coords,
				    std::ofstream& logFile) 
  { return 0.0; };



  dbMatrix& getNodalCoords() { return nodalCoords; };

};

#endif
