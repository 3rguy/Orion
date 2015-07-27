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

    Tria3ElementTemplate();
    ~Tria3ElementTemplate() {};

    /// The set of shape functions.
    double N(int func,dbVector& coords);
    double N1(dbVector& coords);
    double N2(dbVector& coords);
    double N3(dbVector& coords);

    double dN(int func,int deriv,dbVector& coords) { return 0; };
    double dN1(int deriv,dbVector& coords) { return 0; };
    double dN2(int deriv,dbVector& coords) { return 0; };
    double dN3(int deriv,dbVector& coords) { return 0;};

    /// Calculate the metrics: normal vector
    void getMetrics(std::vector<Particle>& particles,intVector& nodes,
                       dbVector& coords,dbMatrix& metrics,double& metric,
                       std::ofstream& logFile);

    /// Calculate the transformation and metric factor for the integration
    /// of a surface element's area.
    double getMetricFactor(std::vector<Particle>& particles,intVector& nodes,
                           dbVector& coords,std::ofstream& logFile);

};

/***********************************************************************/

struct Rect4ElementTemplate : public ElementTemplate {

    Rect4ElementTemplate();
   ~Rect4ElementTemplate() {};

    /// The set of shape functions.
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

    /// Calculate the metrics: normal vector
    void getMetrics(std::vector<Particle>& particles,intVector& nodes,
                       dbVector& coords,dbMatrix& metrics,double& metric,
                       std::ofstream& logFile);

    /// Calculate the transformation and metric factor for the integration
    /// of a surface element's area.
    double getMetricFactor(std::vector<Particle>& particles,intVector& nodes,
                           dbVector& coords,std::ofstream& logFile);

};

/***********************************************************************/

struct Rect9ElementTemplate : public ElementTemplate {

    /// ordering: (0,0), (-1,1),  (-1,-1), (1,-1), (0,1), (-1,0), (0,-1),
    ///           (1,0), (0,0)
    Rect9ElementTemplate();
    ~Rect9ElementTemplate() {};

    /// The set of shape functions.
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

    double dN(int func,int deriv,dbVector& coords) { return 0; };

    /// Calculate the metrics: normal vector
    void getMetrics(std::vector<Particle>& particles,intVector& nodes,
                    dbVector& coords,dbMatrix& metrics,double& metric,
                    std::ofstream& logFile);

    /// Calculate the transformation and metric factor for the integration
    /// of a surface element's area.
    double getMetricFactor(std::vector<Particle>& particles,intVector& nodes,
                           dbVector& coords,std::ofstream& logFile);

};

#endif
