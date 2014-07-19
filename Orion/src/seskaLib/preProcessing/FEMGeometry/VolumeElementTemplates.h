#ifndef VolumeElementTemplates_h_
#define VolumeElementTemplates_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "ElementTemplate.h"
#include "Particle.h"

#include "mpi.h"
#include "petscsys.h"


struct Tetra4ElementTemplate : public ElementTemplate {


  ~Tetra4ElementTemplate() {};

  // NOTE: Node order of coordinates' vector must be as follows:


  Tetra4ElementTemplate();

  // The set of shape functions.
  double N(int func,dbVector& coords);
  double N1(dbVector& coords);
  double N2(dbVector& coords);
  double N3(dbVector& coords);
  double N4(dbVector& coords);

  double dN(int func,int deriv,dbVector& coords);
  double dN1(int deriv,dbVector& coords);
  double dN2(int deriv,dbVector& coords);
  double dN3(int deriv,dbVector& coords);
  double dN4(int deriv,dbVector& coords);

  // Compute Jacobian dx_i/dL_j
  void jacobian(std::vector<Particle>& particles,
		intVector& nodes,
		dbVector& coords,
		dbMatrix& Jmat,
		std::ofstream& logFile) {};

  // Calculate the length of a edge of volume element to integrate
  // over this edge.
  void lineMetric(std::vector<Particle>& particles,
		  intVector& nodes,
		  dbVector& coords,double& length,
		  std::ofstream& logFile);

  // Calculate surface normal vector on volume element's side and its 
  // length, which is the absolute value of the surface, in order to 
  // integrate over this volume element's side.
  void surfaceMetric(std::vector<Particle>& particles,
		     intVector& nodes,
		     dbVector& coords,dbVector& surfaceNormal,
		     double& surface,std::ofstream& logFile);

  // Calculate the transformation and metric factor for the integration
  // of a tetrahedra element's volume.
  double metricFactorVolume(std::vector<Particle>& particles,
			    intVector& nodes,
			    dbVector& coords,
			    std::ofstream& logFile);

  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile) {

    return metricFactorVolume(particles,nodes,coords,logFile); };

};

/***********************************************************************/

struct Cube8ElementTemplate : public ElementTemplate {

  ~Cube8ElementTemplate() {};

  // NOTE: Node order of coordinates' vector must be as follows: 
  // (1,1,1),(-1,1,1),(-1,-1,1),(1,-1,1),
  // (1,1,-1),(-1,1,-1),(-1,-1,-1),(1,-1,-1)

  Cube8ElementTemplate();

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

  double dN(int func,int deriv,dbVector& coords);
  double dN1(int deriv,dbVector& coords);
  double dN2(int deriv,dbVector& coords);
  double dN3(int deriv,dbVector& coords);
  double dN4(int deriv,dbVector& coords);
  double dN5(int deriv,dbVector& coords);
  double dN6(int deriv,dbVector& coords);
  double dN7(int deriv,dbVector& coords);
  double dN8(int deriv,dbVector& coords);

  // Compute Jacobian dx_i/dL_j
  void jacobian(std::vector<Particle>& particles,
		intVector& nodes,
		dbVector& coords,
		dbMatrix& Jmat,
		std::ofstream& logFile);

  // Calculate the length of a edge of volume element to integrate
  // over this edge.
  void lineMetric(std::vector<Particle>& particles,
		  intVector& nodes,
		  dbVector& coords,double& length,
		  std::ofstream& logFile);

  // Calculate surface normal vector on volume element's side and its 
  // length, which is the absolute value of the surface, in order to 
  // integrate over this volume element's side.
  void surfaceMetric(std::vector<Particle>& particles,
		     intVector& nodes,
		     dbVector& coords,
		     dbVector& surfaceNormal, double& surface,
		     std::ofstream& logFile);

  // Calculate the transformation and metric factor for the integration
  // of a hexahedra element's volume.
  double metricFactorVolume(std::vector<Particle>& particles,
			    intVector& nodes,
			    dbVector& coords,std::ofstream& logFile);

  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile) {

    return metricFactorVolume(particles,nodes,coords,logFile); };


};

/***********************************************************************/

struct Cube27ElementTemplate : public ElementTemplate {

  ~Cube27ElementTemplate() {};

  // NOTE: Node order of coordinates' vector must be as follows: 


  Cube27ElementTemplate();

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
  double N10(dbVector& coords);
  double N11(dbVector& coords);
  double N12(dbVector& coords);
  double N13(dbVector& coords);
  double N14(dbVector& coords);
  double N15(dbVector& coords);
  double N16(dbVector& coords);
  double N17(dbVector& coords);
  double N18(dbVector& coords);
  double N19(dbVector& coords);
  double N20(dbVector& coords);
  double N21(dbVector& coords);
  double N22(dbVector& coords);
  double N23(dbVector& coords);
  double N24(dbVector& coords);
  double N25(dbVector& coords);
  double N26(dbVector& coords);
  double N27(dbVector& coords);

  double dN(int func,int deriv,dbVector& coords) ;
  double dN1(int deriv,dbVector& coords) ;
  double dN2(int deriv,dbVector& coords) ;
  double dN3(int deriv,dbVector& coords) ;
  double dN4(int deriv,dbVector& coords) ;
  double dN5(int deriv,dbVector& coords) ;
  double dN6(int deriv,dbVector& coords) ;
  double dN7(int deriv,dbVector& coords) ;
  double dN8(int deriv,dbVector& coords) ;
  double dN9(int deriv,dbVector& coords) ;
  double dN10(int deriv,dbVector& coords) ;
  double dN11(int deriv,dbVector& coords) ;
  double dN12(int deriv,dbVector& coords) ;
  double dN13(int deriv,dbVector& coords) ;
  double dN14(int deriv,dbVector& coords) ;
  double dN15(int deriv,dbVector& coords) ;
  double dN16(int deriv,dbVector& coords) ;
  double dN17(int deriv,dbVector& coords) ;
  double dN18(int deriv,dbVector& coords) ;
  double dN19(int deriv,dbVector& coords) ;
  double dN20(int deriv,dbVector& coords) ;
  double dN21(int deriv,dbVector& coords) ;
  double dN22(int deriv,dbVector& coords) ;
  double dN23(int deriv,dbVector& coords) ;
  double dN24(int deriv,dbVector& coords) ;
  double dN25(int deriv,dbVector& coords) ;
  double dN26(int deriv,dbVector& coords) ;
  double dN27(int deriv,dbVector& coords) ;

  // Compute Jacobian dx_i/dL_j
  void jacobian(std::vector<Particle>& particles,
		intVector& nodes,
		dbVector& coords,
		dbMatrix& Jmat,
		std::ofstream& logFile);


  // Calculate the length of a edge of volume element to integrate
  // over this edge.
  void lineMetric(std::vector<Particle>& particles,
		  intVector& nodes,
		  dbVector& coords,double& length,
		  std::ofstream& logFile);

  // Calculate surface normal vector on volume element's side and its 
  // length, which is the absolute value of the surface, in order to 
  // integrate over this volume element's side.
  void surfaceMetric(std::vector<Particle>& particles,
		     intVector& nodes,
		     dbVector& coords,
		     dbVector& surfaceNormal, double& surface,
		     std::ofstream& logFile);

  // Calculate the transformation and metric factor for the integration
  // of a hexahedra element's volume.
  double metricFactorVolume(std::vector<Particle>& particles,
			    intVector& nodes,
			    dbVector& coords,
			    std::ofstream& logFile);

  double metricFactor(std::vector<Particle>& particles,
		      intVector& nodes,
		      dbVector& coords,
		      std::ofstream& logFile) {

    return metricFactorVolume(particles,nodes,coords,logFile); };
};

#endif
