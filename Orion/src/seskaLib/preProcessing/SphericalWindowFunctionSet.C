// Calculate the window function ordinate of each particle supporting at a
// given point. The window functions have a spherical shape.

#include "SphericalWindowFunctionSet.h"

SphericalWindowFunctionSet::SphericalWindowFunctionSet(
    InputFileData* InputData,std::vector<Particle>& ptcls,intVector& sPtcls,
    dbVector& coords,int& supportSize,unsigned int derivationOrder,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  calcWinFunctionSet(InputData,ptcls,sPtcls,coords,supportSize,derivationOrder,
                     modelData,logFile);

}

/***********************************************************************/
/***********************************************************************/
// Calculation of the window functions.
void
SphericalWindowFunctionSet::calcWinFunctionSet(
    InputFileData* InputData,std::vector<Particle>& ptcls,intVector& sPtcls,
    dbVector& coords,int& supportSize,unsigned int derivationOrder,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int winFuncNorming = (int) InputData->getValue("windowfunctionNorming");
  int windowFuncType = (int) InputData->getValue("windowFunctionType");
  int radiusMode = (int) InputData->getValue("radiusDeterminationAlgorithm");

  int particlesNum = ptcls.size();
  intMatrix v = vectorToMatrix(usedDims);

  dbVector xnorm(usedDims),dOrds(usedDims),d2Ords(v.size());
  double normFactor,dnormFactor;
  int m,n;

  if(windowFuncs.size() < supportSize) windowFuncs.resize(supportSize);

  if(derivationOrder > 0 && firstDerivWinFuncs.size() < supportSize)

  allocateArray(firstDerivWinFuncs,usedDims,supportSize);

  if(derivationOrder > 1 && secondDerivWinFuncs.size() < v.size()) allocateArray(
      secondDerivWinFuncs,v.size(),supportSize);

  for(int i = 0;i < supportSize;i++) {
    m = sPtcls[i];

    for(int k = 0;k < coords.size();k++)
      xnorm[k] = (coords[k] - ptcls[m].getCoord(k));

    // Choose a window function type.
    switch(windowFuncType) {

      // constant 1:
      case 0:

        break;

        // cubic spline
      case 1:

        SphericalWindowFunctionSet::calcCubicWinFunc(ptcls[m],xnorm,
                                                     windowFuncs[i],modelData,
                                                     logFile);

        if(derivationOrder > 0) SphericalWindowFunctionSet::calcCubicWinFunc1stDerivs(
            ptcls[m],xnorm,dOrds,modelData,logFile);

        if(derivationOrder > 1) SphericalWindowFunctionSet::calcCubicWinFunc2ndDerivs(
            ptcls[m],xnorm,d2Ords,modelData,logFile);

        break;

      default:
        cerr << "Chosen window function type isn't supported!" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
        break;
    }

    // ------------------------------------------------------------------
    // norming by 3D domain of influence
    if(winFuncNorming == 1) normFactor = 1.0
      / pow(ptcls[m].getRadius(0),usedDims);

    // norming by particle weight obtained from FEM mesh
    else if(winFuncNorming == 2) normFactor = 1.0 / ptcls[m].getWeight();

    else normFactor = 1.0;

    windowFuncs[i] *= normFactor;

    // Calculation of first order derivation of the window functions
    if(derivationOrder > 0) {

      for(int k = 0;k < usedDims;k++) {
        firstDerivWinFuncs[k][i] = dOrds[k] * normFactor;
      }

    }

    // Calculation of second order derivation of the window functions
    if(derivationOrder > 1) {

      for(int k = 0;k < usedDims;k++) {
        secondDerivWinFuncs[k][i] = d2Ords[k] * normFactor;
      }
      for(int k = usedDims;k < v.size();k++) {
        secondDerivWinFuncs[k][i] = d2Ords[k] * normFactor;
      }

    }

  }

}

/***********************************************************************/
/***********************************************************************/
// Calculate the window function in 3-D.
void
SphericalWindowFunctionSet::calcCubicWinFunc(
    Particle& ptcle,dbVector& xnorm,double& wFunc,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  dbVector& radii = ptcle.getRadii();

  // s = [(x-xI)^2+(y-yI)^2+(z-zI)^2)^0.5 / r0 = r/r0

  double rad = 0;

  for(int k = 0;k < xnorm.size();k++)
    rad += pow(xnorm[k],2);

  rad = sqrt(rad);

  double s = rad / radii[0];

  if(s <= 0.5) wFunc = cubicSplineW1(s);

  else if(s > 0.5 && s <= 1.0) wFunc = cubicSplineW2(s);

  else {
    cerr << "While calculating window functions chosen point is out of "
        << "particle range!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void
SphericalWindowFunctionSet::calcCubicWinFunc1stDerivs(
    Particle& ptcle,dbVector& xnorm,dbVector& dWFuncs,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  dbVector& radii = ptcle.getRadii();

  // s = [(x-xI)^2+(y-yI)^2+(z-zI)^2)^0.5 / r0 = r/r0

  double rad = 0;

  for(int k = 0;k < xnorm.size();k++)
    rad += pow(xnorm[k],2);

  rad = sqrt(rad);

  double s = rad / radii[0];

  // dw/dx = dw/ds*ds/dr*dr/dx = dw/ds*1/r0*r^(-1)*xnorm

  if(s < DBL_EPSILON) {

    for(int k = 0;k < xnorm.size();k++)
      dWFuncs[k] = 0;

  }
  else if(s > 0 && s <= 0.5) {

    for(int k = 0;k < xnorm.size();k++)
      dWFuncs[k] = cubicSplineDW1(s) / radii[0] * 1.0 / rad * xnorm[k];

  }
  else if(s > 0.5 && s <= 1.0) {

    for(int k = 0;k < xnorm.size();k++)
      dWFuncs[k] = cubicSplineDW2(s) / radii[0] * 1.0 / rad * xnorm[k];

  }
  else {
    cerr << "While calculating window functions chosen point is out of "
        << "particle range!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void
SphericalWindowFunctionSet::calcCubicWinFunc2ndDerivs(
    Particle& ptcle,dbVector& xnorm,dbVector& d2WFuncs,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int usedDims = xnorm.size();
  intMatrix v = vectorToMatrix(xnorm.size());

  dbVector& radii = ptcle.getRadii();

  // s = [(x-xI)^2+(y-yI)^2+(z-zI)^2)^0.5 / r0 = r/r0

  double rad = 0;

  for(int k = 0;k < xnorm.size();k++)
    rad += pow(xnorm[k],2);

  rad = sqrt(rad);

  double s = rad / radii[0];

  // d^2w/dx^2 = d^2w/ds^2*ds/dr*dr/dx*ds/dr*dr/dx + dw/ds*(d^2s/dr^2*dr/dx + ds/dr*d^2r/dx^2) =
  //           = d^2w/ds^2*1/r0^2*r^(-2)*xnorm^2 + dw/ds*[0 + ds/dr*(-2r^(-2)*xnorm^2 + r^(-1))] =
  //           = d^2w/ds^2*1/r0^2*r^(-2)*xnorm^2 -2dw/ds*1/r0*r^(-2)*xnorm^2 + dw/ds*1/r0*r^(-1)
  //
  // ------------------------------------------------------------------------------------------
  //
  // d^2w/dxdy = d^2w/ds^2*ds/dy*ds/dx + dw/ds*(d^2s/dr^2*dr/dy + ds/dr*d^2r/dxdy) =
  //           = d^2w/ds^2*1/r0^2*r^(-2)*xnorm*ynorm + dw/ds*[0 + ds/dr*(-2r^(-2)*xnorm*ynorm] =
  //           = d^2w/ds^2*1/r0^2*r^(-2)*xnorm*ynorm -2dw/ds*1/r0*r^(-2)*xnorm*ynorm

  if(s <= 0.5) {

    for(int k = 0;k < xnorm.size();k++)
      d2WFuncs[k] = cubicSplineD2W1(s) * 1.0 / pow(radii[0],2.0) * 1.0
        / pow(rad,2.0) * pow(xnorm[k],2.0)
        - 2.0 * cubicSplineDW1(s) * 1.0 / radii[0] * 1.0 / pow(rad,2.0)
          * pow(xnorm[k],2.0) + cubicSplineDW1(s) * 1.0 / radii[0] * 1.0 / rad;

    for(int k = usedDims;k < v.size();k++)
      d2WFuncs[k] = cubicSplineD2W1(s) * 1.0 / pow(radii[0],2.0) * 1.0
        / pow(rad,2.0) * (xnorm[v[k][0]] * xnorm[v[k][1]])
        - 2.0 * cubicSplineDW1(s) * 1.0 / radii[0] * 1.0 / pow(rad,2.0)
          * xnorm[v[k][0]] * xnorm[v[k][1]];

  }
  else if(s > 0.5 && s <= 1.0) {

    for(int k = 0;k < xnorm.size();k++)
      d2WFuncs[k] = cubicSplineD2W2(s) * 1.0 / pow(radii[0],2.0) * 1.0
        / pow(rad,2.0) * pow(xnorm[k],2.0)
        - 2.0 * cubicSplineDW2(s) * 1.0 / radii[0] * 1.0 / pow(rad,2.0)
          * pow(xnorm[k],2.0) + cubicSplineDW2(s) * 1.0 / radii[0] * 1.0 / rad;

    for(int k = usedDims;k < v.size();k++)
      d2WFuncs[k] = cubicSplineD2W2(s) * 1.0 / pow(radii[0],2.0) * 1.0
        / pow(rad,2.0) * (xnorm[v[k][0]] * xnorm[v[k][1]])
        - 2.0 * cubicSplineDW2(s) * 1.0 / radii[0] * 1.0 / pow(rad,2.0)
          * xnorm[v[k][0]] * xnorm[v[k][1]];

  }
  else {
    cerr << "While calculating window functions chosen point is out of "
        << "particle range!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
/***********************************************************************/
// cubic spline and its derivatives
double
SphericalWindowFunctionSet::cubicSplineW1(double& s) {

  return (2.0 / 3.0 - 4.0 * pow(s,2) + 4.0 * pow(s,3));
}

double
SphericalWindowFunctionSet::cubicSplineW2(double& s) {

  return (4.0 / 3.0 - 4.0 * s + 4.0 * pow(s,2) - 4.0 / 3.0 * pow(s,3));
}

// First order derivation of the cubic spline.
double
SphericalWindowFunctionSet::cubicSplineDW1(double& s) {

  return ( -8.0 * s + 12.0 * pow(s,2));
}

double
SphericalWindowFunctionSet::cubicSplineDW2(double& s) {

  return ( -4.0 + 8.0 * s - 4.0 * pow(s,2));
}

// Second order derivation of the cubic spline.
double
SphericalWindowFunctionSet::cubicSplineD2W1(double& s) {

  return ( -8.0 + 24.0 * s);
}

double
SphericalWindowFunctionSet::cubicSplineD2W2(double& s) {

  return (8.0 - 8.0 * s);
}
