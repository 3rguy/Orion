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


struct ElementTemplate {

    virtual ~ElementTemplate() {};

    static const int Gauss_Point_Positioning = 1; // 0: Zienkiewicz, 1: GiD

    dbMatrix nodalCoords;
    dbVector nodalWeights;

    virtual double N(int func,dbVector& coords) {
      return 0.0;
    };
    virtual double dN(int func,int deriv,dbVector& coords) {
      return 0.0;
    };

    /// Compute Jacobian dx_i/dL_j
    virtual void getJacobian(std::vector<Particle>& particles,intVector& nodes,
                             dbVector& coords,dbMatrix& jacobian,
                             std::ofstream& logFile) { };

    /// Calculate the transformation and metric factor for the integration
    /// of an element's length, surface or volume.
    virtual double getMetricFactor(std::vector<Particle>& particles,
                                   intVector& nodes,dbVector& coords,
                                   std::ofstream& logFile) {
      return 0.0;
    };

    /// Calculate the metrics
    virtual void getMetrics(std::vector<Particle>& particles,intVector& nodes,
                            dbVector& coords,dbMatrix& metrics,double& metric,
                            std::ofstream& logFile) {};

    dbMatrix& getNodalCoords() {
      return nodalCoords;
    };

};

#endif
