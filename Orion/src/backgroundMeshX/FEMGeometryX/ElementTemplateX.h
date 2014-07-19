#ifndef ElementTemplateX_h_
#define ElementTemplateX_h_

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "commonTypedefs.h"
#include "mpi.h"
#include "ParticleX.h"
#include "petscsys.h"


struct ElementTemplateX  {

  ~ElementTemplateX() {};

  dbMatrix nodalCoords;
  dbVector nodalWeights;

  virtual double N(int func,dbVector& coords) { return 0.0; };
  virtual double dN(int func,int deriv,dbVector& coords) { return 0.0; };

  // Compute Jacobi matrix dx_i/dL_j
  virtual void jacobiMatrix(std::vector<ParticleX>& particles,
			    intVector& nodes,
			    dbVector& coords,
			    dbMatrix& jacobiMatrix,
			    std::ofstream& logFile) {};


  // Calculate the transformation and metric factor for the integration
  // of an element's length, surface or volume.
  virtual double metricFactor(std::vector<ParticleX>& particles,
			      intVector& nodes,
			      dbVector& coords,
			      std::ofstream& logFile) 
  { return 0.0; };

  // Calculate the length of a edge of volume element to integrate
  // over this edge.
  virtual void lineMetric(std::vector<ParticleX>& particles,
			  intVector& nodes,
			  dbVector& coords,double& length,
			  std::ofstream& logFile) {};

  // Calculate surface normal vector on volume element's side and its 
  // length, which is the absolute value of the surface, in order to 
  // integrate over this volume element's side.
  virtual void surfaceMetric(std::vector<ParticleX>& particles,
			     intVector& nodes,
			     dbVector& coords,dbVector& surfaceNormal,
			     double& surface,std::ofstream& logFile)
  {};

  // Calculate the transformation and metric factor for the integration
  // of an element's volume.
  virtual double metricFactorVolume(std::vector<ParticleX>& particles,
				    intVector& nodes,
				    dbVector& coords,
				    std::ofstream& logFile) 
  { return 0.0; };



  dbMatrix& getNodalCoords() { return nodalCoords; };

};

#endif
