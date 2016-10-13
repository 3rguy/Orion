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
    cerr << "Tria3ElementTemplate has no function N" << func + 1 << "!" << endl;
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

/************************************************************************/
// Calculate the metrics: normal vector
void Tria3ElementTemplate::getMetrics(std::vector<Particle>& particles,
                                         intVector& nodes,dbVector& coords,
                                         dbMatrix& metrics,double& metric,
                                         std::ofstream& logFile) {

  using namespace std;

  if(metrics.size() == 0) allocateArray(metrics,1,3);
  clearArray(metrics);

  // Calculate the surface normal vector
  metrics[0][0] = (particles[nodes[0] - 1].getCoord(1)
    * particles[nodes[1] - 1].getCoord(2)
    - particles[nodes[0] - 1].getCoord(2) * particles[nodes[1] - 1].getCoord(1))
    +

    (particles[nodes[1] - 1].getCoord(1) * particles[nodes[2] - 1].getCoord(2)
      - particles[nodes[1] - 1].getCoord(2)
        * particles[nodes[2] - 1].getCoord(1))
    +

    (particles[nodes[2] - 1].getCoord(1) * particles[nodes[0] - 1].getCoord(2)
      - particles[nodes[2] - 1].getCoord(2)
        * particles[nodes[0] - 1].getCoord(1));

  metrics[0][1] = (particles[nodes[0] - 1].getCoord(2)
    * particles[nodes[1] - 1].getCoord(0)
    - particles[nodes[0] - 1].getCoord(0) * particles[nodes[1] - 1].getCoord(2))
    +

    (particles[nodes[1] - 1].getCoord(2) * particles[nodes[2] - 1].getCoord(0)
      - particles[nodes[1] - 1].getCoord(0)
        * particles[nodes[2] - 1].getCoord(2))
    +

    (particles[nodes[2] - 1].getCoord(2) * particles[nodes[0] - 1].getCoord(0)
      - particles[nodes[2] - 1].getCoord(0)
        * particles[nodes[0] - 1].getCoord(2));

  metrics[0][2] = (particles[nodes[0] - 1].getCoord(0)
    * particles[nodes[1] - 1].getCoord(1)
    - particles[nodes[0] - 1].getCoord(1) * particles[nodes[1] - 1].getCoord(0))
    +

    (particles[nodes[1] - 1].getCoord(0) * particles[nodes[2] - 1].getCoord(1)
      - particles[nodes[1] - 1].getCoord(1)
        * particles[nodes[2] - 1].getCoord(0))
    +

    (particles[nodes[2] - 1].getCoord(0) * particles[nodes[0] - 1].getCoord(1)
      - particles[nodes[2] - 1].getCoord(1)
        * particles[nodes[0] - 1].getCoord(0));

  // Calculate the surface area and normalize normal vector
  metric = 0.5*computeNorm(metrics[0],2,logFile);

  for(int i=0;i<metrics[0].size();i++)
    metrics[0][i] /= (2 * metric);

}

/***********************************************************************/
// Calculate the transformation and metric factor for the integration
// of a surface element's area.
double Tria3ElementTemplate::getMetricFactor(std::vector<Particle>& particles,
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

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    intVector idx(4);

    idx[0] = 2;
    idx[1] = 3;
    idx[2] = 0;
    idx[3] = 1;

    reorderVector(nodalCoords,idx);
  }

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
    cerr << "Rect4ElementTemplate has no function N" << func + 1 << "!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Rect4ElementTemplate::N1(dbVector& coords) {
  return (1.0 / 4.0 * (1 + coords[0]) * (1 + coords[1]));
}
double Rect4ElementTemplate::N2(dbVector& coords) {
  return (1.0 / 4.0 * (1 - coords[0]) * (1 + coords[1]));
}
double Rect4ElementTemplate::N3(dbVector& coords) {
  return (1.0 / 4.0 * (1 - coords[0]) * (1 - coords[1]));
}
double Rect4ElementTemplate::N4(dbVector& coords) {
  return (1.0 / 4.0 * (1 + coords[0]) * (1 - coords[1]));
}

/************************************************************************/
double Rect4ElementTemplate::dN(int func,int deriv,dbVector& coords) {

  using namespace std;

  double result;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    switch(func) {
    case 2:
      result = dN1(deriv,coords);
      break;
    case 3:
      result = dN2(deriv,coords);
      break;
    case 0:
      result = dN3(deriv,coords);
      break;
    case 1:
      result = dN4(deriv,coords);
      break;
    default:
      cerr << "Rect4ElementTemplate has no function dN" << func + 1 << "!"
          << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }
  // default Zienkiewicz nodal numbering
  else {

    switch(func) {
    case 0:
      result = dN1(deriv,coords);
      break;
    case 1:
      result = dN2(deriv,coords);
      break;
    case 2:
      result = dN3(deriv,coords);
      break;
    case 3:
      result = dN4(deriv,coords);
      break;
    default:
      cerr << "Rect4ElementTemplate has no function dN" << func + 1 << "!"
          << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  return result;
}

double Rect4ElementTemplate::dN1(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (1.0 / 4.0 * (1 + coords[1]));
    break;

  case 1:
    return (1.0 / 4.0 * (1 + coords[0]));
    break;

  default:
    cerr << "In Rect4ElementTemplate::dN1/dL" << deriv << " not available!"
        << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Rect4ElementTemplate::dN2(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return ( -1.0 / 4.0 * (1 + coords[1]));
    break;

  case 1:
    return (1.0 / 4.0 * (1 - coords[0]));
    break;

  default:
    cerr << "In Rect4ElementTemplate::dN1/dL" << deriv << " not available!"
        << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Rect4ElementTemplate::dN3(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return ( -1.0 / 4.0 * (1 - coords[1]));
    break;

  case 1:
    return ( -1.0 / 4.0 * (1 - coords[0]));
    break;

  default:
    cerr << "In Rect4ElementTemplate::dN3/dL" << deriv << " not available!"
        << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Rect4ElementTemplate::dN4(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (1.0 / 4.0 * (1 - coords[1]));
    break;

  case 1:
    return ( -1.0 / 4.0 * (1 + coords[0]));
    break;

  default:
    cerr << "In Rect4ElementTemplate::dN4/dL" << deriv << " not available!"
        << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/************************************************************************/
// Calculate the metrics: normal vector
void Rect4ElementTemplate::getMetrics(std::vector<Particle>& particles,
                                         intVector& nodes,dbVector& coords,
                                         dbMatrix& metrics,
                                         double& metric,
                                         std::ofstream& logFile) {

  using namespace std;

  if(metrics.size() == 0) allocateArray(metrics,1,3);
  clearArray(metrics);

  intMatrix e = getPermutations(3);

  // Calculate the surface normal vector
  //
  // n_k = eps_{ijk} dx_i/dL_1*dx_j/dL_2
  //
  // n_k = eps_{ijk} dN^I/dL_1*x^I_i dN^J/dL_2*x^J_j

  // Loop 2 times over all three rotational degrees of freedom
  // -> 3 odd and 3 even permutations!
  for(int p = 0;p < e.size();p++) {

    for(int I = 0;I < nodes.size();I++) {

      for(int J = 0;J < nodes.size();J++) {
        metrics[0][e[p][2]] += e[p][3] * dN(I,0,coords)
          * particles[nodes[I] - 1].getCoord(e[p][0]) * dN(J,1,coords)
          * particles[nodes[J] - 1].getCoord(e[p][1]);
      }

    }

  }

  // Calculate the surface area and normalize normal vector
  metric = computeNorm(metrics[0],2,logFile);

  for(int i=0;i<metrics[0].size();i++)
    metrics[0][i] /= metric;

}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a surface element's area.
double Rect4ElementTemplate::getMetricFactor(std::vector<Particle>& particles,
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
    cerr << "Rect9ElementTemplate has no function N" << func + 1 << "!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Rect9ElementTemplate::N1(dbVector& coords) {
  return (1.0 / 4.0 * (pow(coords[0],2) + coords[0])
    * (pow(coords[1],2) + coords[1]));
}
double Rect9ElementTemplate::N2(dbVector& coords) {
  return (1.0 / 4.0 * (pow(coords[0],2) - coords[0])
    * (pow(coords[1],2) + coords[1]));
}
double Rect9ElementTemplate::N3(dbVector& coords) {
  return (1.0 / 4.0 * (pow(coords[0],2) - coords[0])
    * (pow(coords[1],2) - coords[1]));
}
double Rect9ElementTemplate::N4(dbVector& coords) {
  return (1.0 / 4.0 * (pow(coords[0],2) + coords[0])
    * (pow(coords[1],2) - coords[1]));
}
double Rect9ElementTemplate::N5(dbVector& coords) {
  return ( -1.0 / 2.0 * (pow(coords[0],2) - 1.0)
    * (pow(coords[1],2) + coords[1]));
}
double Rect9ElementTemplate::N6(dbVector& coords) {
  return ( -1.0 / 2.0 * (pow(coords[0],2) - coords[0])
    * (pow(coords[1],2) - 1.0));
}
double Rect9ElementTemplate::N7(dbVector& coords) {
  return ( -1.0 / 2.0 * (pow(coords[0],2) - 1.0)
    * (pow(coords[1],2) - coords[1]));
}
double Rect9ElementTemplate::N8(dbVector& coords) {
  return ( -1.0 / 2.0 * (pow(coords[0],2) + coords[0])
    * (pow(coords[1],2) - 1.0));
}
double Rect9ElementTemplate::N9(dbVector& coords) {
  return ((pow(coords[0],2) - 1.0) * (pow(coords[1],2) - 1.0));
}

/************************************************************************/
/// Calculate the metrics: normal vector
void Rect9ElementTemplate::getMetrics(std::vector<Particle>& particles,
                                      intVector& nodes,dbVector& coords,
                                      dbMatrix& metrics,
                                      double& metric,
                                      std::ofstream& logFile) {

  using namespace std;

  if(metrics.size() == 0) allocateArray(metrics,1,3);
  clearArray(metrics);

  // Calculate the surface normal vector

  double xi = coords[0];
  double eta = coords[1];

  double x1 = particles[nodes[0] - 1].getCoord(0);
  double x2 = particles[nodes[1] - 1].getCoord(0);
  double x3 = particles[nodes[2] - 1].getCoord(0);
  double x4 = particles[nodes[3] - 1].getCoord(0);
  double x9 = particles[nodes[4] - 1].getCoord(0);
  double x10 = particles[nodes[5] - 1].getCoord(0);
  double x11 = particles[nodes[6] - 1].getCoord(0);
  double x12 = particles[nodes[7] - 1].getCoord(0);
  double x21 = particles[nodes[8] - 1].getCoord(0);

  double y1 = particles[nodes[0] - 1].getCoord(1);
  double y2 = particles[nodes[1] - 1].getCoord(1);
  double y3 = particles[nodes[2] - 1].getCoord(1);
  double y4 = particles[nodes[3] - 1].getCoord(1);
  double y9 = particles[nodes[4] - 1].getCoord(1);
  double y10 = particles[nodes[5] - 1].getCoord(1);
  double y11 = particles[nodes[6] - 1].getCoord(1);
  double y12 = particles[nodes[7] - 1].getCoord(1);
  double y21 = particles[nodes[8] - 1].getCoord(1);

  double z1 = particles[nodes[0] - 1].getCoord(2);
  double z2 = particles[nodes[1] - 1].getCoord(2);
  double z3 = particles[nodes[2] - 1].getCoord(2);
  double z4 = particles[nodes[3] - 1].getCoord(2);
  double z9 = particles[nodes[4] - 1].getCoord(2);
  double z10 = particles[nodes[5] - 1].getCoord(2);
  double z11 = particles[nodes[6] - 1].getCoord(2);
  double z12 = particles[nodes[7] - 1].getCoord(2);
  double z21 = particles[nodes[8] - 1].getCoord(2);

  metrics[0][0] = (2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta + eta) * y1
    + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta + eta) * y2
    + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta - eta) * y3
    + 2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta - eta) * y4
    - xi * (eta * eta + eta) * y9
    + 2.0 * ( -xi / 2.0 + 1.0 / 4.0) * (eta * eta - 1.0) * y10
    - xi * (eta * eta - eta) * y11
    + 2.0 * ( -xi / 2.0 - 1.0 / 4.0) * (eta * eta - 1.0) * y12
    + 2.0 * xi * (eta * eta - 1.0) * y21)
    * (2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta + 1.0) * z1
      + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta + 1.0) * z2
      + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta - 1.0) * z3
      + 2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta - 1.0) * z4
      + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta + 1.0) * z9
      + 4.0 * ( -xi * xi / 4.0 + xi / 4.0) * eta * z10
      + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta - 1.0) * z11
      + 4.0 * ( -xi * xi / 4.0 - xi / 4.0) * eta * z12
      + 4.0 * (xi * xi / 2.0 - 1.0 / 2.0) * eta * z21)
    - (2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta + eta) * z1
      + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta + eta) * z2
      + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta - eta) * z3
      + 2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta - eta) * z4
      - xi * (eta * eta + eta) * z9
      + 2.0 * ( -xi / 2.0 + 1.0 / 4.0) * (eta * eta - 1.0) * z10
      - xi * (eta * eta - eta) * z11
      + 2.0 * ( -xi / 2.0 - 1.0 / 4.0) * (eta * eta - 1.0) * z12
      + 2.0 * xi * (eta * eta - 1.0) * z21)
      * (2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta + 1.0) * y1
        + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta + 1.0) * y2
        + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta - 1.0) * y3
        + 2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta - 1.0) * y4
        + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta + 1.0) * y9
        + 4.0 * ( -xi * xi / 4.0 + xi / 4.0) * eta * y10
        + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta - 1.0) * y11
        + 4.0 * ( -xi * xi / 4.0 - xi / 4.0) * eta * y12
        + 4.0 * (xi * xi / 2.0 - 1.0 / 2.0) * eta * y21);

  metrics[0][1] = (2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta + eta) * z1
    + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta + eta) * z2
    + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta - eta) * z3
    + 2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta - eta) * z4
    - xi * (eta * eta + eta) * z9
    + 2.0 * ( -xi / 2.0 + 1.0 / 4.0) * (eta * eta - 1.0) * z10
    - xi * (eta * eta - eta) * z11
    + 2.0 * ( -xi / 2.0 - 1.0 / 4.0) * (eta * eta - 1.0) * z12
    + 2.0 * xi * (eta * eta - 1.0) * z21)
    * (2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta + 1.0) * x1
      + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta + 1.0) * x2
      + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta - 1.0) * x3
      + 2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta - 1.0) * x4
      + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta + 1.0) * x9
      + 4.0 * ( -xi * xi / 4.0 + xi / 4.0) * eta * x10
      + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta - 1.0) * x11
      + 4.0 * ( -xi * xi / 4.0 - xi / 4.0) * eta * x12
      + 4.0 * (xi * xi / 2.0 - 1.0 / 2.0) * eta * x21)
    - (2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta + eta) * x1
      + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta + eta) * x2
      + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta - eta) * x3
      + 2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta - eta) * x4
      - xi * (eta * eta + eta) * x9
      + 2.0 * ( -xi / 2.0 + 1.0 / 4.0) * (eta * eta - 1.0) * x10
      - xi * (eta * eta - eta) * x11
      + 2.0 * ( -xi / 2.0 - 1.0 / 4.0) * (eta * eta - 1.0) * x12
      + 2.0 * xi * (eta * eta - 1.0) * x21)
      * (2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta + 1.0) * z1
        + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta + 1.0) * z2
        + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta - 1.0) * z3
        + 2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta - 1.0) * z4
        + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta + 1.0) * z9
        + 4.0 * ( -xi * xi / 4.0 + xi / 4.0) * eta * z10
        + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta - 1.0) * z11
        + 4.0 * ( -xi * xi / 4.0 - xi / 4.0) * eta * z12
        + 4.0 * (xi * xi / 2.0 - 1.0 / 2.0) * eta * z21);

  metrics[0][2] = (2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta + eta) * x1
    + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta + eta) * x2
    + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta - eta) * x3
    + 2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta - eta) * x4
    - xi * (eta * eta + eta) * x9
    + 2.0 * ( -xi / 2.0 + 1.0 / 4.0) * (eta * eta - 1.0) * x10
    - xi * (eta * eta - eta) * x11
    + 2.0 * ( -xi / 2.0 - 1.0 / 4.0) * (eta * eta - 1.0) * x12
    + 2.0 * xi * (eta * eta - 1.0) * x21)
    * (2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta + 1.0) * y1
      + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta + 1.0) * y2
      + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta - 1.0) * y3
      + 2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta - 1.0) * y4
      + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta + 1.0) * y9
      + 4.0 * ( -xi * xi / 4.0 + xi / 4.0) * eta * y10
      + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta - 1.0) * y11
      + 4.0 * ( -xi * xi / 4.0 - xi / 4.0) * eta * y12
      + 4.0 * (xi * xi / 2.0 - 1.0 / 2.0) * eta * y21)
    - (2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta + eta) * y1
      + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta + eta) * y2
      + 2.0 * (xi / 4.0 - 1.0 / 8.0) * (eta * eta - eta) * y3
      + 2.0 * (xi / 4.0 + 1.0 / 8.0) * (eta * eta - eta) * y4
      - xi * (eta * eta + eta) * y9
      + 2.0 * ( -xi / 2.0 + 1.0 / 4.0) * (eta * eta - 1.0) * y10
      - xi * (eta * eta - eta) * y11
      + 2.0 * ( -xi / 2.0 - 1.0 / 4.0) * (eta * eta - 1.0) * y12
      + 2.0 * xi * (eta * eta - 1.0) * y21)
      * (2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta + 1.0) * x1
        + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta + 1.0) * x2
        + 2.0 * (xi * xi / 8.0 - xi / 8.0) * (2.0 * eta - 1.0) * x3
        + 2.0 * (xi * xi / 8.0 + xi / 8.0) * (2.0 * eta - 1.0) * x4
        + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta + 1.0) * x9
        + 4.0 * ( -xi * xi / 4.0 + xi / 4.0) * eta * x10
        + 2.0 * ( -xi * xi / 4.0 + 1.0 / 4.0) * (2.0 * eta - 1.0) * x11
        + 4.0 * ( -xi * xi / 4.0 - xi / 4.0) * eta * x12
        + 4.0 * (xi * xi / 2.0 - 1.0 / 2.0) * eta * x21);


  // Calculate the surface area and normalize normal vector
  metric = computeNorm(metrics[0],2,logFile);

  for(int i=0;i<metrics[0].size();i++)
    metrics[0][i] /= metric;

}

/************************************************************************/
// Calculate the transformation and metric factor for the integration
// of a surface element's area.
double Rect9ElementTemplate::getMetricFactor(std::vector<Particle>& particles,
                                            intVector& nodes,dbVector& coords,
                                            std::ofstream& logFile) {
  
  using namespace std;

  double metric;
  dbMatrix metrics;

  getMetrics(particles,nodes,coords,metrics,metric,
             logFile);
  return (metric);
}
