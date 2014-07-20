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
    cerr <<"Line2ElementTemplate has no function N"<<func+1<<"!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Line2ElementTemplate::N1(dbVector& coords) {
  return (1.0/2.0*(1+coords[0]));
}
double Line2ElementTemplate::N2(dbVector& coords) {
  return (1.0/2.0*(1-coords[0]));
}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a line element's length.
double Line2ElementTemplate::metricFactor(std::vector<Particle>& particles,
					  intVector& nodes,
					  dbVector& coords,
					  std::ofstream& logFile) {

  return(fabs(0.5*particles[nodes[0]-1].getCoord(0)-
	      0.5*particles[nodes[1]-1].getCoord(0)));

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
    cerr <<"Line3ElementTemplate has no function N"<<func+1<<"!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Line3ElementTemplate::N1(dbVector& coords) {
  return (1.0/2.0*(pow(coords[0],2)+coords[0]));
}
double Line3ElementTemplate::N2(dbVector& coords) {
  return (1.0/2.0*(pow(coords[0],2)-coords[0]));
}
double Line3ElementTemplate::N3(dbVector& coords) {
  return (-pow(coords[0],2)+1.0);
}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a line element's length.
double Line3ElementTemplate::metricFactor(std::vector<Particle>& particles,
					  intVector& nodes,
					  dbVector& coords,
					  std::ofstream& logFile) {

  return((coords[0]+1.0)*particles[nodes[0]-1].getCoord(0)-
	 2.0*coords[0]*particles[nodes[1]-1].getCoord(0)+
	 (coords[0]-1.0)*particles[nodes[2]-1].getCoord(0));

}