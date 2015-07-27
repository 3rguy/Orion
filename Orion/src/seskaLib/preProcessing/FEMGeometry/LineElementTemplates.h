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

    /// The set of shape functions.
    double N(int func,dbVector& coords);
    double N1(dbVector& coords);
    double N2(dbVector& coords);

    double dN(int func,int deriv,dbVector& coords);
    double dN1(int deriv,dbVector& coords);
    double dN2(int deriv,dbVector& coords);
    double dN3(int deriv,dbVector& coords);

    /// Calculate the metrics: tangent vector
    void getMetrics(std::vector<Particle>& particles,intVector& nodes,
                    dbVector& coords,dbMatrix& metrics,double& metric,
                    std::ofstream& logFile);

    /// Calculate the transformation and metric factor for the integration
    /// of a line element's length.
    double getMetricFactor(std::vector<Particle>& particles,intVector& nodes,
                           dbVector& coords,std::ofstream& logFile);

};

/***********************************************************************/

struct Line3ElementTemplate : public ElementTemplate {

    /// ordering: (1), (-1), (0)
    Line3ElementTemplate();
    ~Line3ElementTemplate() {};

    /// The set of shape functions.
    double N(int func,dbVector& coords);
    double N1(dbVector& coords);
    double N2(dbVector& coords);
    double N3(dbVector& coords);

    double dN(int func,int deriv,dbVector& coords) {
      return 0;
    };
    double dN1(int deriv,dbVector& coords) {
      return 0;
    };
    double dN2(int deriv,dbVector& coords) {
      return 0;
    };
    double dN3(int deriv,dbVector& coords) {
      return 0;
    };

    /// Calculate the metrics: tangent vector
    void getMetrics(std::vector<Particle>& particles,intVector& nodes,
                    dbVector& coords,dbMatrix& metrics,double& metric,
                    std::ofstream& logFile);

    /// Calculate the transformation and metric factor for the integration
    /// of a line element's length.
    double getMetricFactor(std::vector<Particle>& particles,intVector& nodes,
                           dbVector& coords,std::ofstream& logFile);

};

#endif
