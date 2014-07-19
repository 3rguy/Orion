#include "SurfaceElementTemplates.h"



Tria3ElementTemplate::Tria3ElementTemplate() {

  allocateArray(nodalCoords,3,3);

  nodalCoords[0][0] = 1;
  nodalCoords[0][1] = 0;
  nodalCoords[0][2] = 0;

  nodalCoords[1][0] = 0;
  nodalCoords[1][1] = 1;
  nodalCoords[1][2] = 0;

  nodalCoords[2][0] = 0;
  nodalCoords[2][1] = 0;
  nodalCoords[2][2] = 1;

}

/***********************************************************************/
double Tria3ElementTemplate::N(int func,dbVector& coords) {

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
    cerr <<"Tria3ElementTemplate has no function N"<<func+1<<"!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Tria3ElementTemplate::N1(dbVector& coords) {
  return (coords[0]);
}
double Tria3ElementTemplate::N2(dbVector& coords) {
  return (coords[1]);
}
double Tria3ElementTemplate::N3(dbVector& coords) {
  return (coords[2]);
}

/***********************************************************************/
// Calculate the transformation and metric factor for the integration
// of a surface element's area.
double Tria3ElementTemplate::metricFactor(std::vector<Particle>& particles,
					  intVector& nodes,
					  dbVector& coords,
					  std::ofstream& logFile) {
  
  double det;

  det = 0.5*(particles[nodes[1]-1].getCoord(0)
	     *particles[nodes[2]-1].getCoord(1)
	    +particles[nodes[0]-1].getCoord(0)
	     *particles[nodes[1]-1].getCoord(1)
	    +particles[nodes[0]-1].getCoord(1)
	     *particles[nodes[2]-1].getCoord(0)
	    -particles[nodes[0]-1].getCoord(0)
	     *particles[nodes[2]-1].getCoord(1)
	    -particles[nodes[1]-1].getCoord(1)
	     *particles[nodes[2]-1].getCoord(0)
	    -particles[nodes[0]-1].getCoord(1)
	     *particles[nodes[1]-1].getCoord(0));

  return(fabs(det));
}

/************************************************************************/
/************************************************************************/
Rect4ElementTemplate::Rect4ElementTemplate() {

  allocateArray(nodalCoords,4,2);

  nodalCoords[0][0] = 1;
  nodalCoords[0][1] = 1;

  nodalCoords[1][0] = -1;
  nodalCoords[1][1] = 1;

  nodalCoords[2][0] = -1;
  nodalCoords[2][1] = -1;

  nodalCoords[3][0] = 1;
  nodalCoords[3][1] = -1;

}

/************************************************************************/
double Rect4ElementTemplate::N(int func,dbVector& coords) {

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
  case 3:
    result = N4(coords);
    break;
  default:
    cerr <<"Rect4ElementTemplate has no function N"<<func+1<<"!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Rect4ElementTemplate::N1(dbVector& coords) {
  return (1.0/4.0*(1+coords[0])*(1+coords[1]));
}
double Rect4ElementTemplate::N2(dbVector& coords) {
  return (1.0/4.0*(1-coords[0])*(1+coords[1]));
}
double Rect4ElementTemplate::N3(dbVector& coords) {
  return (1.0/4.0*(1-coords[0])*(1-coords[1]));
}
double Rect4ElementTemplate::N4(dbVector& coords) {
  return (1.0/4.0*(1+coords[0])*(1-coords[1]));
}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a surface element's area.
double Rect4ElementTemplate::metricFactor(std::vector<Particle>& particles,
					  intVector& nodes,
					  dbVector& coords,
					  std::ofstream& logFile) {
  
  double det;

  det =(1.0/4.0*(1+coords[1])*particles[nodes[0]-1].getCoord(0)
	-1.0/4.0*(1+coords[1])*particles[nodes[1]-1].getCoord(0)
	-1.0/4.0*(1-coords[1])*particles[nodes[2]-1].getCoord(0)
	+1.0/4.0*(1-coords[1])*particles[nodes[3]-1].getCoord(0))
      *(1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(1)
	+1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(1)
	-1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(1)
	-1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(1))
      -(1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(0)
	+1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(0)
	-1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(0)
	-1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(0))
      *(1.0/4.0*(1+coords[1])*particles[nodes[0]-1].getCoord(1)
	-1.0/4.0*(1+coords[1])*particles[nodes[1]-1].getCoord(1)
	-1.0/4.0*(1-coords[1])*particles[nodes[2]-1].getCoord(1)
	+1.0/4.0*(1-coords[1])*particles[nodes[3]-1].getCoord(1));

  return(fabs(det));
}

/************************************************************************/
/************************************************************************/
Rect9ElementTemplate::Rect9ElementTemplate() {

  allocateArray(nodalCoords,9,2);

  nodalCoords[0][0] = 1;
  nodalCoords[0][1] = 1;

  nodalCoords[1][0] = -1;
  nodalCoords[1][1] = 1;

  nodalCoords[2][0] = -1;
  nodalCoords[2][1] = -1;

  nodalCoords[3][0] = 1;
  nodalCoords[3][1] = -1;

  nodalCoords[4][0] = 0;
  nodalCoords[4][1] = 1;

  nodalCoords[5][0] = -1;
  nodalCoords[5][1] = 0;

  nodalCoords[6][0] = 0;
  nodalCoords[6][1] = -1;

  nodalCoords[7][0] = 1;
  nodalCoords[7][1] = 0;

  nodalCoords[8][0] = 0;
  nodalCoords[8][1] = 0;

}

/************************************************************************/
double Rect9ElementTemplate::N(int func,dbVector& coords) {

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
  case 3:
    result = N4(coords);
    break;
  case 4:
    result = N5(coords);
    break;
  case 5:
    result = N6(coords);
    break;
  case 6:
    result = N7(coords);
    break;
  case 7:
    result = N8(coords);
    break;
  case 8:
    result = N9(coords);
    break;
  default:
    cerr <<"Rect9ElementTemplate has no function N"<<func+1<<"!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Rect9ElementTemplate::N1(dbVector& coords) {
  return(1.0/4.0*(pow(coords[0],2)+coords[0])*(pow(coords[1],2)+coords[1]));
}
double Rect9ElementTemplate::N2(dbVector& coords) {
  return(1.0/4.0*(pow(coords[0],2)-coords[0])*(pow(coords[1],2)+coords[1]));
}
double Rect9ElementTemplate::N3(dbVector& coords) {
  return(1.0/4.0*(pow(coords[0],2)-coords[0])*(pow(coords[1],2)-coords[1]));
}
double Rect9ElementTemplate::N4(dbVector& coords) {
  return(1.0/4.0*(pow(coords[0],2)+coords[0])*(pow(coords[1],2)-coords[1]));
}
double Rect9ElementTemplate::N5(dbVector& coords) {
  return(-1.0/2.0*(pow(coords[0],2)-1.0)*(pow(coords[1],2)+coords[1]));
}
double Rect9ElementTemplate::N6(dbVector& coords) {
  return(-1.0/2.0*(pow(coords[0],2)-coords[0])*(pow(coords[1],2)-1.0));
}
double Rect9ElementTemplate::N7(dbVector& coords) {
  return(-1.0/2.0*(pow(coords[0],2)-1.0)*(pow(coords[1],2)-coords[1]));
}
double Rect9ElementTemplate::N8(dbVector& coords) {
  return(-1.0/2.0*(pow(coords[0],2)+coords[0])*(pow(coords[1],2)-1.0));
}
double Rect9ElementTemplate::N9(dbVector& coords) {
  return((pow(coords[0],2)-1.0)*(pow(coords[1],2)-1.0));
}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a surface element's area.
double Rect9ElementTemplate::metricFactor(std::vector<Particle>& particles,
					  intVector& nodes,
					  dbVector& coords,
					  std::ofstream& logFile) {
  
  double det;

  det = 

    (1.0/4.0*(2.0*coords[0]+1.0)*(pow(coords[1],2)+coords[1])*particles[nodes[0]-1].getCoord(0)
     -coords[0]*(pow(coords[1],2)+coords[1])*particles[nodes[1]-1].getCoord(0)
     +1.0/4.0*(2.0*coords[0]-1.0)*(pow(coords[1],2)+coords[1])*particles[nodes[2]-1].getCoord(0)
     -1.0/2.0*(2.0*coords[0]-1.0)*(pow(coords[1],2)-1.0)*particles[nodes[3]-1].getCoord(0)
     +1.0/4.0*(2.0*coords[0]-1.0)*(pow(coords[1],2)-coords[1])*particles[nodes[4]-1].getCoord(0)
     -coords[0]*(pow(coords[1],2)-coords[1])*particles[nodes[5]-1].getCoord(0)
     +1.0/4.0*(2.0*coords[0]+1.0)*(pow(coords[1],2)-coords[1])*particles[nodes[6]-1].getCoord(0)
     -1.0/2.0*(2.0*coords[0]+1.0)*(pow(coords[1],2)-1.0)*particles[nodes[7]-1].getCoord(0)
     +2.0*coords[0]*(pow(coords[1],2)-1.0)*particles[nodes[8]-1].getCoord(0)) *

    (1.0/4.0*(pow(coords[0],2)+coords[0])*(2.0*coords[1]+1.0)*particles[nodes[0]-1].getCoord(1)
     -1.0/2.0*(pow(coords[0],2)-1.0)*(2.0*coords[1]+1.0)*particles[nodes[1]-1].getCoord(1)
     +1.0/4.0*(pow(coords[0],2)-coords[0])*(2.0*coords[1]+1.0)*particles[nodes[2]-1].getCoord(1)
     -(pow(coords[0],2)-coords[0])*coords[1]*particles[nodes[3]-1].getCoord(1)
     +1.0/4.0*(pow(coords[0],2)-coords[0])*(2.0*coords[1]-1.0)*particles[nodes[4]-1].getCoord(1)
     -1.0/2.0*(pow(coords[0],2)-1.0)*(2.0*coords[1]-1.0)*particles[nodes[5]-1].getCoord(1)
     +1.0/4.0*(pow(coords[0],2)+coords[0])*(2.0*coords[1]-1.0)*particles[nodes[6]-1].getCoord(1)
     -(pow(coords[0],2)+coords[0])*coords[1]*particles[nodes[7]-1].getCoord(1)
     +2.0*(pow(coords[0],2)-1.0)*coords[1]*particles[nodes[8]-1].getCoord(1)) -

    (1.0/4.0*(pow(coords[0],2)+coords[0])*(2.0*coords[1]+1.0)*particles[nodes[0]-1].getCoord(0)
     -1.0/2.0*(pow(coords[0],2)-1.0)*(2.0*coords[1]+1.0)*particles[nodes[1]-1].getCoord(0)
     +1.0/4.0*(pow(coords[0],2)-coords[0])*(2.0*coords[1]+1.0)*particles[nodes[2]-1].getCoord(0)
     -(pow(coords[0],2)-coords[0])*coords[1]*particles[nodes[3]-1].getCoord(0)
     +1.0/4.0*(pow(coords[0],2)-coords[0])*(2.0*coords[1]-1.0)*particles[nodes[4]-1].getCoord(0)
     -1.0/2.0*(pow(coords[0],2)-1.0)*(2.0*coords[1]-1.0)*particles[nodes[5]-1].getCoord(0)
     +1.0/4.0*(pow(coords[0],2)+coords[0])*(2.0*coords[1]-1.0)*particles[nodes[6]-1].getCoord(0)
     -(pow(coords[0],2)+coords[0])*coords[1]*particles[nodes[7]-1].getCoord(0)
     +2.0*(pow(coords[0],2)-1.0)*coords[1]*particles[nodes[8]-1].getCoord(0)) *

    (1.0/4.0*(2.0*coords[0]+1.0)*(pow(coords[1],2)+coords[1])*particles[nodes[0]-1].getCoord(1)
     -coords[0]*(pow(coords[1],2)+coords[1])*particles[nodes[1]-1].getCoord(1)
     +1.0/4.0*(2.0*coords[0]-1.0)*(pow(coords[1],2)+coords[1])*particles[nodes[2]-1].getCoord(1)
     -1.0/2.0*(2.0*coords[0]-1.0)*(pow(coords[1],2)-1.0)*particles[nodes[3]-1].getCoord(1)
     +1.0/4.0*(2.0*coords[0]-1.0)*(pow(coords[1],2)-coords[1])*particles[nodes[4]-1].getCoord(1)
     -coords[0]*(pow(coords[1],2)-coords[1])*particles[nodes[5]-1].getCoord(1)
     +1.0/4.0*(2.0*coords[0]+1.0)*(pow(coords[1],2)-coords[1])*particles[nodes[6]-1].getCoord(1)
     -1.0/2.0*(2.0*coords[0]+1.0)*(pow(coords[1],2)-1.0)*particles[nodes[7]-1].getCoord(1)
     +2.0*coords[0]*(pow(coords[1],2)-1.0)*particles[nodes[8]-1].getCoord(1));

  return(fabs(det));
}
