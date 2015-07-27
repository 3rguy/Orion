#include "LineElementTemplates.h"

Line2ElementTemplate::Line2ElementTemplate() {

  allocateArray(nodalCoords,2,1);

  nodalCoords[0][0] = 1;
  nodalCoords[1][0] = -1;
}

/************************************************************************/
double Line2ElementTemplate::N(int func,dbVector& coords) {

  using namespace std;

  double result;

  switch(func) {
  case 0:
    result = N1(coords);
    break;
  case 1:
    result = N2(coords);
    break;
  default:
    cerr << "Line2ElementTemplate has no function N" << func + 1 << "!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Line2ElementTemplate::N1(dbVector& coords) {
  return (1.0 / 2.0 * (1 + coords[0]));
}
double Line2ElementTemplate::N2(dbVector& coords) {
  return (1.0 / 2.0 * (1 - coords[0]));
}

/************************************************************************/
double Line2ElementTemplate::dN(int func,int deriv,dbVector& coords) {

  using namespace std;

  double result;

  switch(func) {
  case 0:
    result = dN1(deriv,coords);
    break;
  case 1:
    result = dN2(deriv,coords);
    break;
  default:
    cerr << "Line2ElementTemplate has no function dN" << func + 1 << "!"
        << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Line2ElementTemplate::dN1(int deriv,dbVector& coords) {

  using namespace std;

  return (1.0 / 2.0);
}
double Line2ElementTemplate::dN2(int deriv,dbVector& coords) {

  using namespace std;

  return (-1.0 / 2.0);
}

/************************************************************************/
// Calculate the metrics: tangent vector
void Line2ElementTemplate::getMetrics(std::vector<Particle>& particles,
                                       intVector& nodes,dbVector& coords,
                                       dbMatrix& metrics,
                                       double& metric,
                                       std::ofstream& logFile) {

  using namespace std;

  if(metrics.size() == 0) allocateArray(metrics,1,3);
  clearArray(metrics);

  // t_i = dx_i/dL
  //
  // t_i = dN^I/dL*x^I_i

  // Loop over all 3 dimensions
  for(int i = 0;i < metrics[0].size();i++) {

    for(int I = 0;I < nodes.size();I++) {
      metrics[0][i] += dN(I,0,coords) * particles[nodes[I] - 1].getCoord(i);
    }

  }

  // Calculate the surface area and normalize normal vector
  metric = computeNorm(metrics[0],2,logFile);

  for(int i=0;i<metrics[0].size();i++)
    metrics[0][i] /= metric;

}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a line element's length.
double Line2ElementTemplate::getMetricFactor(std::vector<Particle>& particles,
                                             intVector& nodes,dbVector& coords,
                                             std::ofstream& logFile) {


  using namespace std;

  double metric;
  dbMatrix metrics;

  getMetrics(particles,nodes,coords,metrics,metric,
             logFile);
  return (metric);
}

/************************************************************************/
/************************************************************************/
Line3ElementTemplate::Line3ElementTemplate() {

  allocateArray(nodalCoords,3,1);

  nodalCoords[0][0] = 1;
  nodalCoords[1][0] = -1;
  nodalCoords[2][0] = 0;
}

/************************************************************************/
double Line3ElementTemplate::N(int func,dbVector& coords) {

  using namespace std;

  double result;

  switch(func) {
  case 0:
    result = N1(coords);
    break;
  case 1:
    result = N2(coords);
    break;
  case 2:
    result = N3(coords);
    break;
  default:
    cerr << "Line3ElementTemplate has no function N" << func + 1 << "!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Line3ElementTemplate::N1(dbVector& coords) {
  return (1.0 / 2.0 * (pow(coords[0],2) + coords[0]));
}
double Line3ElementTemplate::N2(dbVector& coords) {
  return (1.0 / 2.0 * (pow(coords[0],2) - coords[0]));
}
double Line3ElementTemplate::N3(dbVector& coords) {
  return ( -pow(coords[0],2) + 1.0);
}

/************************************************************************/
// Calculate the metrics: tangent vector
void Line3ElementTemplate::getMetrics(std::vector<Particle>& particles,
                                       intVector& nodes,dbVector& coords,
                                       dbMatrix& metrics,
                                       double& metric,
                                       std::ofstream& logFile) {

  using namespace std;

  if(metrics.size() == 0) allocateArray(metrics,1,3);
  clearArray(metrics);

  // t_i = dx_i/dL
  //
  // t_i = dN^I/dL*x^I_i

  // Loop over all 3 dimensions
  for(int i = 0;i < metrics[0].size();i++) {

    for(int I = 0;I < nodes.size();I++) {
      metrics[0][i] += dN(I,0,coords) * particles[nodes[I] - 1].getCoord(i);
    }

  }

  // Calculate the surface area and normalize normal vector
  metric = computeNorm(metrics[0],2,logFile);

  for(int i=0;i<metrics[0].size();i++)
    metrics[0][i] /= metric;

}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a line element's length.
double Line3ElementTemplate::getMetricFactor(std::vector<Particle>& particles,
                                             intVector& nodes,dbVector& coords,
                                             std::ofstream& logFile) {

  using namespace std;

  double metric;
  dbMatrix metrics;

  getMetrics(particles,nodes,coords,metrics,metric,
             logFile);
  return (metric);
}
