// Calculate the window function of a point in relation to all its 
// supporting parrticles considering different influence radii
// in negative and posititve coordinate direction.

#include "WindowFuncAsym.h"

WindowFuncAsym::WindowFuncAsym(InputFileData* InputData,
                               std::vector<Particle>& ptcls,intVector& sPtcls,
                               double& x,double& y,double& z,int& supportSize,
                               unsigned int derivationOrder,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile,PetscViewer& viewerSEQ) :
    numOfSegments(0), integrationOrder(0) {

  using namespace std;

  // Calculation of the window functions 
  if(derivationOrder == 0) calcWinFunctions(InputData,ptcls,sPtcls,x,y,z,
                                            supportSize,modelData,logFile,
                                            viewerSEQ);
  
  // Calculation of the window function and its first order derivatives 
  else if(derivationOrder == 1) calcWinFunction1stDerivs(InputData,ptcls,sPtcls,
                                                         x,y,z,supportSize,
                                                         modelData,logFile,
                                                         viewerSEQ);

  // Calculation of the window function, its first and second order 
  // derivatives 
  else calcWinFunction2ndDerivs(InputData,ptcls,sPtcls,x,y,z,supportSize,
                                modelData,logFile,viewerSEQ);

}

double WindowFuncAsym::cubicSpline(double s) {

  using namespace std;

  double sabs = fabs(s);

  if(sabs <= 0.5)

  return (2.0 / 3.0 - 4 * pow(sabs,2.0) + 4.0 * pow(sabs,3.0));

  else if(sabs > 0.5 && sabs <= 1.0)

  return (4.0 / 3.0 - 4.0 * sabs + 4.0 * pow(sabs,2.0)
    - 4.0 / 3.0 * pow(sabs,3.0));

  else if(sabs > 1.0)

  return (0);

  else cerr << "In WindowFuncAsym::cubicSpline failure in weight function!"
      << endl;
  
}

// First order derivation of the cubic spline.
double WindowFuncAsym::dcubicSpline(double s) {

  using namespace std;

  double sabs = fabs(s);
  double sign;

  if(sabs > 0) sign = s / sabs;

  else sign = 1.0;

  if(sabs <= 0.5)

  return ( -8.0 * s + 12.0 * sabs * s);

  else if(sabs > 0.5 && sabs <= 1.0)

  return ( -4.0 * sign + 8.0 * s - 4.0 * sabs * s);

  else if(sabs > 1.0)

  return (0);

  else cerr << "In WindowFuncAsym::dcubicSpline failure in weight function!"
      << endl;
  
}

double WindowFuncAsym::d2cubicSpline(double s) {

  using namespace std;

  double sabs = fabs(s);
  double sign;

  if(sabs > 0) sign = s / sabs;

  else sign = 1.0;

  if(sabs <= 0.5)

  return ( -8.0 + 24.0 * sabs);

  else if(sabs > 0.5 && sabs <= 1.0)

  return (8.0 - 8.0 * sabs);

  else if(sabs > 1.0)

  return (0);

  else cerr << "In WindowFuncAsym::d2cubicSpline failure in weight function!"
      << endl;

}

double WindowFuncAsym::quarticSpline(double s) {
  
  using namespace std;

  return (1.0 - 6.0 * pow(s,2.0) + 8.0 * pow(s,3.0) - 3.0 * pow(s,4.0));

}

/************************************************************************/
/************************************************************************/
// cubic basis polynomial
dbVector WindowFuncAsym::P(double x) {

  using namespace std;

  vector<double> basis(5);

  basis[0] = 1;            // 1
  basis[1] = x;            // x
  basis[2] = x * x;          // x^2
  basis[3] = x * x * x;        // x^3
  basis[4] = x * x * x * x;      // x^4

  return (basis);
}

dbVector WindowFuncAsym::dP(double x) {

  using namespace std;

  vector<double> basis(5);

  basis[0] = 0;            // 0
  basis[1] = 1;            // 1
  basis[2] = 2 * x;          // 2x
  basis[3] = 3 * x * x;        // 3x^2
  basis[4] = 4 * x * x * x;      // 4x^3

  return (basis);
}

dbVector WindowFuncAsym::d2P(double x) {

  using namespace std;

  vector<double> basis(5);

  basis[0] = 0;            // 0
  basis[1] = 0;            // 0
  basis[2] = 2;            // 2
  basis[3] = 6 * x;          // 6x
  basis[4] = 12 * x * x;       // 12x^2

  return (basis);
}

dbVector WindowFuncAsym::d3P(double x) {

  using namespace std;

  vector<double> basis(5);

  basis[0] = 0;            // 0
  basis[1] = 0;            // 0
  basis[2] = 0;            // 0
  basis[3] = 6;            // 6
  basis[4] = 24 * x;         // 24x

  return (basis);
}

dbVector WindowFuncAsym::d4P(double x) {

  using namespace std;

  vector<double> basis(5);

  basis[0] = 0;            // 0
  basis[1] = 0;            // 0
  basis[2] = 0;            // 0
  basis[3] = 0;            // 0
  basis[4] = 24;           // 24

  return (basis);
}

dbVector WindowFuncAsym::d5P(double x) {

  using namespace std;

  vector<double> basis(5);

  basis[0] = 0;            // 0
  basis[1] = 0;            // 0
  basis[2] = 0;            // 0
  basis[3] = 0;            // 0
  basis[4] = 0;            // 0
  basis[5] = 0;          // 0

  return (basis);
}

/************************************************************************/
/************************************************************************/
// cubic exponential basis polynomial
dbVector WindowFuncAsym::expP(double x) {

  using namespace std;

  vector<double> basis(4);

  double lambda = 0.5;
  
  basis[0] = 1.0;
  basis[1] = x;
  basis[2] = exp(lambda * x);
  basis[3] = exp( -lambda * x);

  return (basis);
}

dbVector WindowFuncAsym::dexpP(double x) {

  using namespace std;

  vector<double> basis(4);

  double lambda = 0.5;
  
  basis[0] = 0.0;
  basis[1] = 1.0;
  basis[2] = lambda * exp(lambda * x);
  basis[3] = -lambda * exp( -lambda * x);

  return (basis);
}

dbVector WindowFuncAsym::d2expP(double x) {

  using namespace std;

  vector<double> basis(4);

  double lambda = 0.5;
  
  basis[0] = 0.0;
  basis[1] = 0.0;
  basis[2] = lambda * lambda * exp(lambda * x);
  basis[3] = lambda * lambda * exp( -lambda * x);

  return (basis);
}

/************************************************************************/
/************************************************************************/
// compute spline coefficients (currently only cubic)
void WindowFuncAsym::setCustomPtcleSpline(
    InputFileData* InputData,Particle& ptcle,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int windowFuncType = (int) InputData->getValue("windowFunctionType");
  int winFuncNorming = (int) InputData->getValue("windowfunctionNorming");

  dbVector& radii = ptcle.getRadii();
  CustomSpline* spline = ptcle.getSpline();

  dbMatrix3& cValues = spline->getSplineCoefficients();
  dbMatrix& splineKnots = spline->getSplineKnots();

  int pSize,linEQSize;

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"*********** asymmetric spline setup ******************"<<endl;
  logFile<<"PTCLE "<<ptcle.getID()<<endl;
#endif

  /**********************************************************************/
  // loop over dimensions and set the spline leg lengths for each 
  // dimension as well as the spline segments
  dbMatrix radius;
  allocateArray(radius,usedDims,2);

  for(int dim = 0;dim < usedDims;dim++) {

    // both spline legs 
    if(radii[dim] > 0 && radii[usedDims + dim] > 0) {

      radius[dim][0] = radii[dim];
      radius[dim][1] = radii[usedDims + dim];
    }

    // positive spline leg only 
    else if(radii[dim] > 0 && radii[usedDims + dim] == 0) {

      radius[dim][0] = radii[dim];
      radius[dim][1] = radii[dim];
    }

    // negative spline leg only
    else if(radii[dim] == 0 && radii[usedDims + dim] > 0) {

      radius[dim][0] = radii[usedDims + dim];
      radius[dim][1] = radii[usedDims + dim];
    }

    else {
      logFile << "In WindowFuncAsym1::setCustomPtcleSpline influence radius\n"
          << "zero!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

  /*********************************************************************/
  // set the segments of the spline
  // standard cubic spline:
  if(windowFuncType == 1) {

    numOfSegments = 2;
    integrationOrder = 3;
    
    pSize = 0;
    linEQSize = 0;

  }
  // customized quartic spline:
  else if(windowFuncType == 2) {

    numOfSegments = 6;
    integrationOrder = 3;
    
    pSize = P(0).size();
    linEQSize = pSize * numOfSegments;

  }
  // customized exponential spline
  else if(windowFuncType == 4) {

    numOfSegments = 8;
    integrationOrder = 3;
    
    pSize = expP(0).size();
    linEQSize = pSize * numOfSegments;

  }
  else {
    logFile << "WindowFuncAsym::setCustomPtcleSpline windowfunction type "
        << windowFuncType << " is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  allocateArray(splineKnots,usedDims,numOfSegments + 1);

  for(int dim = 0;dim < usedDims;dim++) {

    for(int i = 0;i < (int) numOfSegments / 2.0;i++) {

      splineKnots[dim][(int) (numOfSegments / 2.0) - 1 - i] =

      -radius[dim][1] * (i + 1) / (numOfSegments / 2.0);

      splineKnots[dim][(int) (numOfSegments / 2.0) + 1 + i] =

      radius[dim][0] * (i + 1) / (numOfSegments / 2.0);
      
    }

  }

#ifdef _geometryDebugMode_
  logFile<<"pSize="<<pSize<<endl;
  logFile<<"linEQSize="<<linEQSize<<endl;
  for(int dim=0;dim<usedDims;dim++) {
    for(int i=0;i<splineKnots[dim].size()-1;i++) {
      logFile<<"segment "<<i<<": "
      <<splineKnots[dim][i]<<" to "
      <<splineKnots[dim][i+1]<<endl;

    }
  }
#endif

  /**********************************************************************/
  // set Gauss integration tools
  int nodesPerElem,numOfIntPoints;
  ElementTemplate* FEMSet;
  GaussPointSet* GaussSet;

  intVector data(4);
  map<string,double> params;

  data[0] = 2; // element type (line,rectangle,brick)
  data[1] = integrationOrder; // volume integration order
  data[2] = integrationOrder; // surface integration order 
  data[3] = integrationOrder; // line integration order

  getGaussQuadratureData(data,params);

  // set FEM approximation tools
  intVector nodesIdx(2);
  FEMSet = new Line2ElementTemplate();
  vector<Particle> nodes(2,Particle(0));

  numOfIntPoints = (int) params["gaussPointsPerLineElement"];

  switch(numOfIntPoints) {

  case 1:
    GaussSet = new GaussSetLine1();
    break;

  case 2:
    GaussSet = new GaussSetLine2();
    break;

  case 3:
    GaussSet = new GaussSetLine3();
    break;

  case 4:
    GaussSet = new GaussSetLine4();
    break;

  case 5:
    GaussSet = new GaussSetLine5();
    break;

  default:
    logFile << "In WindowFuncAsym1::setCustomPtcleSpline " << numOfIntPoints
        << " line Gauss points are not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /*********************************************************************/
  // assemble the equation system, i.e. set coefficient matrix and right 
  // side vector
  // continuity and boundary conditions
  vector<vector<double> > conditions(linEQSize,vector<double>(5));

  // conditions[][0]:  x-ordinate
  //
  // conditions[][1]:  != 1.0e+06: y-ordinate 
  //                   == 1.0e+06: segment(s1) = segment(s2)
  //
  // conditions[][2]:     0: phi 
  //                      1: dphi 
  //                      2: d2phi 
  //                      3: d3phi
  //
  // conditions[][3]:     segment 1 
  //
  // conditions[][4]:     segment 2 

  int idx,idx1,idx2,cnd;
  double xord,xord1,xord2;
  dbVector values(pSize);
  dbVector values1(pSize);
  dbVector values2(pSize);

  dbMatrix AValues(linEQSize,vector<double>(linEQSize));
  dbVector bValues(linEQSize);

  allocateArray(cValues,usedDims,numOfSegments,pSize);

  if(windowFuncType == 2 || windowFuncType == 4) {

    // loop over all dimensions
    for(int dim = 0;dim < usedDims;dim++) {

      clearArray(AValues);
      clearArray(bValues);

#ifdef _geometryDebugMode_
      logFile<<"****************************************************"<<endl;
      logFile<<"* dimension "<<dim<<endl;
      logFile<<"----------------------------------------------------"<<endl;
#endif

      cnd = 0;

      // quartic spline: 6 segments i.e. 5 x 6 = 30 conditions needed
      if(windowFuncType == 2) {

        // continuity conditions of phi, dphi, d2phi, d3phi:
        // (6-1)*4  = 20

        // phi_a =  phi_b

        for(int i = 0;i < splineKnots[dim].size() - 2;i++) {

          conditions[cnd][0] = splineKnots[dim][i + 1];
          conditions[cnd][1] = 1.0e+06;
          conditions[cnd][2] = 0;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i + 1;

          cnd++;
        }

        // dphi_a =  dphi_b

        for(int i = 0;i < splineKnots[dim].size() - 2;i++) {

          conditions[cnd][0] = splineKnots[dim][i + 1];
          conditions[cnd][1] = 1.0e+06;
          conditions[cnd][2] = 1;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i + 1;
          cnd++;

        }

        // d2phi_a =  d2phi_b

        for(int i = 0;i < splineKnots[dim].size() - 2;i++) {

          conditions[cnd][0] = splineKnots[dim][i + 1];
          conditions[cnd][1] = 1.0e+06;
          conditions[cnd][2] = 2;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i + 1;
          cnd++;

        }

        // d3phi_a =  d3phi_b

        for(int i = 0;i < splineKnots[dim].size() - 2;i++) {

          conditions[cnd][0] = splineKnots[dim][i + 1];
          conditions[cnd][1] = 1.0e+06;
          conditions[cnd][2] = 3;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i + 1;
          cnd++;

        }

        // -----------------------------------------------------------------

        // phi_l = standard_phi                                         (1-3)

        for(int i = 0;i < (int) splineKnots[dim].size() / 2.0 - 1;i++) {

          conditions[cnd][0] = splineKnots[dim][i];
          conditions[cnd][1] = quarticSpline(
              fabs(splineKnots[dim][i]) / radius[dim][1]);
          conditions[cnd][2] = 0;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i;

          cnd++;
        }

        // phi_r = standard_phi                                         (4-7)

        for(int i = (int) (splineKnots[dim].size() / 2.0);
            i < (int) splineKnots[dim].size();i++) {

          conditions[cnd][0] = splineKnots[dim][i];
          conditions[cnd][1] = quarticSpline(
              fabs(splineKnots[dim][i]) / radius[dim][0]);
          conditions[cnd][2] = 0;

          if(i < (int) splineKnots[dim].size() - 1) {
            conditions[cnd][3] = i;
            conditions[cnd][4] = i;
          }
          else {
            conditions[cnd][3] = i - 1;
            conditions[cnd][4] = i - 1;
          }

          cnd++;
        }

        // -----------------------------------------------------------------

        // dphi(-1) = 0                                                  (8)

        conditions[cnd][0] = splineKnots[dim][0];
        conditions[cnd][1] = 0;
        conditions[cnd][2] = 1;
        conditions[cnd][3] = 0;
        conditions[cnd][4] = 0;

        cnd++;

        // dphi_r(0) = 0                                                 (9)

        conditions[cnd][0] = 0;
        conditions[cnd][1] = 0;
        conditions[cnd][2] = 1;
        conditions[cnd][3] = splineKnots[dim].size() / 2;
        conditions[cnd][4] = splineKnots[dim].size() / 2;

        cnd++;

        // dphi(1) = 0                                                  (10)

        conditions[cnd][0] = splineKnots[dim][splineKnots[dim].size() - 1];
        conditions[cnd][1] = 0;
        conditions[cnd][2] = 1;
        conditions[cnd][3] = splineKnots[dim].size() - 2;
        conditions[cnd][4] = splineKnots[dim].size() - 2;

        cnd++;

#ifdef _geometryDebugMode_
        for(int i=0;i<conditions.size();i++) {
          logFile<<"condition "<<i<<": xord="
          <<conditions[i][0]<<" yord="
          <<conditions[i][1]<<" deriv="
          <<conditions[i][2]<<" segment1="
          <<conditions[i][3]<<" segment2="
          <<conditions[i][4]<<endl;
        }
#endif

        // =================================================================

        cnd = 0;

        // loop over all conditions except the integral one.
        for(int i = 0;i < linEQSize;i++) {

          // --------------------
          // enforce a y-ordinate

          if(conditions[i][1] != 1.0e+06) {

            xord = conditions[i][0];
            idx = (int) conditions[i][3];

            switch((int) conditions[i][2]) {
            case 0:
              values = P(xord);
              break;

            case 1:
              values = dP(xord);
              break;

            case 2:
              values = d2P(xord);
              break;

            case 3:
              values = d3P(xord);
              break;

            case 4:
              values = d4P(xord);
              break;

            case 5:
              values = d5P(xord);
              break;

            default:
              logFile
                  << "In WindowFuncAsym1::setCustomPtcleSpline polynom order"
                  << (int) conditions[i][2] << " is not supported!" << endl;
              MPI_Abort(MPI_COMM_WORLD,1);
              break;
            }

            for(int j = 0;j < values.size();j++)

              AValues[i][idx * pSize + j] = values[j];

            bValues[i] = conditions[i][1];

#ifdef _geometryDebugMode_
            logFile<<"----------------------------------------------"<<endl;
            logFile<<"xord="<<xord<<endl;
            for(int j=0;j<values.size();j++)
            logFile<<"A["<<i<<"]["<<idx*pSize+j<<"] = "
            <<AValues[i][idx*pSize+j]<<endl;
            logFile<<"----------------"<<endl;
            logFile<<"b["<<i<<"] = "<<bValues[i]<<endl;
#endif

          }

          // ---------------------------------------------
          // enforce continuity between neighbouring segments

          if(conditions[i][1] == 1.0e+06) {

            // segment 3,4
            if(conditions[i][3] > 1)

            xord1 = conditions[i][0];

            else

            xord1 = conditions[i][0];

            idx1 = (int) conditions[i][3];

            // segment 3,4
            if(conditions[i][4] > 1)

            xord2 = conditions[i][0];

            else

            xord2 = conditions[i][0];

            idx2 = (int) conditions[i][4];

            switch((int) conditions[i][2]) {
            case 0:
              values1 = P(xord1);
              values2 = P(xord2);
              break;

            case 1:
              values1 = dP(xord1);
              values2 = dP(xord2);
              break;

            case 2:
              values1 = d2P(xord1);
              values2 = d2P(xord2);
              break;

            case 3:
              values1 = d3P(xord1);
              values2 = d3P(xord2);
              break;

            case 4:
              values1 = d4P(xord1);
              values2 = d4P(xord2);
              break;

            case 5:
              values = d5P(xord1);
              values = d5P(xord2);
              break;

            default:
              logFile
                  << "In WindowFuncAsym1::setCustomPtcleSpline polynom order"
                  << (int) conditions[i][2] << " is not supported!" << endl;
              MPI_Abort(MPI_COMM_WORLD,1);
              break;
            }

            for(int j = 0;j < values1.size();j++)

              AValues[i][idx1 * pSize + j] = values1[j];

            for(int j = 0;j < values2.size();j++)

              AValues[i][idx2 * pSize + j] = ( -1.0) * values2[j];

            bValues[i] = 0;

#ifdef _geometryDebugMode_
            logFile<<"----------------------------------------------"<<endl;
            logFile<<"xord1="<<xord1<<endl;
            logFile<<"xord2="<<xord2<<endl;
            for(int j=0;j<values.size();j++)
            logFile<<"A["<<i<<"]["<<idx1*pSize+j<<"] = "
            <<AValues[i][idx1*pSize+j]<<endl;
            logFile<<"----------------"<<endl;
            for(int j=0;j<values.size();j++)
            logFile<<"A["<<i<<"]["<<idx2*pSize+j<<"] = "
            <<AValues[i][idx2*pSize+j]<<endl;
            logFile<<"----------------"<<endl;
            logFile<<"b["<<i<<"] = "<<bValues[i]<<endl;
#endif

          }

          cnd++;
        }

      }
      
      /********************************************************************/
      // exponential spline 
      // E(x) = a0 + a1*(x - xj) + a2*exp(lambdaj*(x - xj) + a3*exp(-lambdaj*(x - xj)
      //
      // lambdaj: tension parameter to oppress oscillation (-> 0 stronger)
      // xj: starting node of segment
      //
      // 8 segments i.e. 8 x 4 = 32 conditions needed
      else if(windowFuncType == 4) {

        // continuity conditions of phi, dphi, d2phi:
        // (s-1)*3

        // phi_a =  phi_b

        for(int i = 0;i < splineKnots[dim].size() - 2;i++) {

          conditions[cnd][0] = i + 1;
          conditions[cnd][1] = 1.0e+06;
          conditions[cnd][2] = 0;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i + 1;

          cnd++;
        }

        // dphi_a =  dphi_b

        for(int i = 0;i < splineKnots[dim].size() - 2;i++) {

          conditions[cnd][0] = i + 1;
          conditions[cnd][1] = 1.0e+06;
          conditions[cnd][2] = 1;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i + 1;
          cnd++;

        }

        // d2phi_a =  d2phi_b

        for(int i = 0;i < splineKnots[dim].size() - 2;i++) {

          conditions[cnd][0] = i + 1;
          conditions[cnd][1] = 1.0e+06;
          conditions[cnd][2] = 2;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i + 1;
          cnd++;

        }

        // -----------------------------------------------------------------

        // phi_l = standard_phi

        for(int i = 0;i < (int) splineKnots[dim].size() / 2.0 - 1;i++) {

          conditions[cnd][0] = splineKnots[dim][i] - splineKnots[dim][i];
          conditions[cnd][1] = cubicSpline(
              fabs(splineKnots[dim][i]) / radius[dim][1]);
          conditions[cnd][2] = 0;
          conditions[cnd][3] = i;
          conditions[cnd][4] = i;

          cnd++;
        }

        // phi_r = standard_phi

        for(int i = (int) (splineKnots[dim].size() / 2.0);
            i < (int) splineKnots[dim].size();i++) {

          if(i < (int) splineKnots[dim].size() - 1) {
            conditions[cnd][0] = splineKnots[dim][i] - splineKnots[dim][i];
            conditions[cnd][1] = cubicSpline(
                fabs(splineKnots[dim][i]) / radius[dim][0]);
            conditions[cnd][2] = 0;
            conditions[cnd][3] = i;
            conditions[cnd][4] = i;
          }
          else {
            conditions[cnd][0] = splineKnots[dim][i] - splineKnots[dim][i - 1];
            conditions[cnd][1] = cubicSpline(
                fabs(splineKnots[dim][i]) / radius[dim][0]);
            conditions[cnd][2] = 0;
            conditions[cnd][3] = i - 1;
            conditions[cnd][4] = i - 1;
          }

          cnd++;
        }

        // -----------------------------------------------------------------

        // dphi(-1)) = 0                                                 (1)

        conditions[cnd][0] = splineKnots[dim][0] - splineKnots[dim][0];
        conditions[cnd][1] = 0;
        conditions[cnd][2] = 1;
        conditions[cnd][3] = 0;
        conditions[cnd][4] = 0;

        cnd++;

        //       // dphi_r(0) = 0              x=kn(3), xj = kn(3)                (2)

        //       conditions[cnd][0] = splineKnots[dim][(int)splineKnots[dim].size()/2]
        // 	-splineKnots[dim][(int)splineKnots[dim].size()/2];
        //       conditions[cnd][1] = 0;
        //       conditions[cnd][2] = 1;
        //       conditions[cnd][3] = splineKnots[dim].size()/2;
        //       conditions[cnd][4] = splineKnots[dim].size()/2;

        //       cnd++;

        // dphi(1) = 0                                                   (2)

        conditions[cnd][0] = splineKnots[dim][splineKnots[dim].size() - 1]
          - splineKnots[dim][splineKnots[dim].size() - 2];
        conditions[cnd][1] = 0;
        conditions[cnd][2] = 1;
        conditions[cnd][3] = splineKnots[dim].size() - 2;
        conditions[cnd][4] = splineKnots[dim].size() - 2;

        cnd++;

#ifdef _geometryDebugMode_
        for(int i=0;i<conditions.size();i++) {
          logFile<<"condition "<<i<<": xord="
          <<conditions[i][0]<<" yord="
          <<conditions[i][1]<<" deriv="
          <<conditions[i][2]<<" segment1="
          <<conditions[i][3]<<" segment2="
          <<conditions[i][4]<<endl;
        }
#endif

        // =================================================================

        cnd = 0;

        // loop over all conditions except the integral one.
        for(int i = 0;i < linEQSize;i++) {

          // --------------------
          // enforce a y-ordinate

          if(conditions[i][1] != 1.0e+06) {

            xord = conditions[i][0];
            idx = (int) conditions[i][3];

            switch((int) conditions[i][2]) {
            case 0:
              values = expP(xord);
              break;

            case 1:
              values = dexpP(xord);
              break;

            case 2:
              values = d2expP(xord);
              break;

            default:
              logFile
                  << "In WindowFuncAsym1::setCustomPtcleSpline polynom order"
                  << (int) conditions[i][2] << " is not supported!" << endl;
              MPI_Abort(MPI_COMM_WORLD,1);
              break;
            }

            for(int j = 0;j < values.size();j++)

              AValues[i][idx * pSize + j] = values[j];

            bValues[i] = conditions[i][1];

#ifdef _geometryDebugMode_
            logFile<<"----------------------------------------------"<<endl;
            logFile<<"xord="<<xord<<endl;
            for(int j=0;j<values.size();j++)
            logFile<<"A["<<i<<"]["<<idx*pSize+j<<"] = "
            <<AValues[i][idx*pSize+j]<<endl;
            logFile<<"----------------"<<endl;
            logFile<<"b["<<i<<"] = "<<bValues[i]<<endl;
#endif

          }

          // ---------------------------------------------
          // enforce continuity between neighbouring segments

          if(conditions[i][1] == 1.0e+06) {

            // segment 3,4
            xord1 = splineKnots[dim][(int) conditions[i][0]]
              - splineKnots[dim][(int) conditions[i][0] - 1];
            idx1 = (int) conditions[i][3];

            // segment 3,4
            xord2 = splineKnots[dim][(int) conditions[i][0]]
              - splineKnots[dim][(int) conditions[i][0]];
            idx2 = (int) conditions[i][4];

            switch((int) conditions[i][2]) {
            case 0:
              values1 = expP(xord1);
              values2 = expP(xord2);
              break;

            case 1:
              values1 = dexpP(xord1);
              values2 = dexpP(xord2);
              break;

            case 2:
              values1 = d2expP(xord1);
              values2 = d2expP(xord2);
              break;

            default:
              logFile
                  << "In WindowFuncAsym1::setCustomPtcleSpline polynom order"
                  << (int) conditions[i][2] << " is not supported!" << endl;
              MPI_Abort(MPI_COMM_WORLD,1);
              break;
            }

            for(int j = 0;j < values1.size();j++)

              AValues[i][idx1 * pSize + j] = values1[j];

            for(int j = 0;j < values2.size();j++)

              AValues[i][idx2 * pSize + j] = ( -1.0) * values2[j];

            bValues[i] = 0;

#ifdef _geometryDebugMode_
            logFile<<"----------------------------------------------"<<endl;
            logFile<<"xord1="<<xord1<<endl;
            logFile<<"xord2="<<xord2<<endl;
            for(int j=0;j<values.size();j++)
            logFile<<"A["<<i<<"]["<<idx1*pSize+j<<"] = "
            <<AValues[i][idx1*pSize+j]<<endl;
            logFile<<"----------------"<<endl;
            for(int j=0;j<values.size();j++)
            logFile<<"A["<<i<<"]["<<idx2*pSize+j<<"] = "
            <<AValues[i][idx2*pSize+j]<<endl;
            logFile<<"----------------"<<endl;
            logFile<<"b["<<i<<"] = "<<bValues[i]<<endl;
#endif

          }

          cnd++;
        }

      }
      else {
        logFile << "WindowFuncAsym1::setCustomPtcleSpline windowfunction type "
            << windowFuncType << " is not supported!" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      // -----------------------------------------------------------------
      // solve the equation system

      vector<double> xValues(linEQSize);
      double* xCalc;

      // Initialize PETSc objects.
      KSP ksp;
      Vec b;
      Vec x;
      VecCreateSeq(PETSC_COMM_SELF,linEQSize, &b);
      VecCreateSeq(PETSC_COMM_SELF,linEQSize, &x);

      Mat A;
      MatCreateSeqDense(PETSC_COMM_SELF,linEQSize,linEQSize,PETSC_NULL, &A);
      MatSetFromOptions(A);

      vector<int> idxVec(linEQSize);

      for(int i = 0;i < linEQSize;i++)
        idxVec[i] = i;

      for(int i = 0;i < linEQSize;i++)

        MatSetValues(A,1, &i,linEQSize, &idxVec[0], &AValues[i][0],
                     INSERT_VALUES);

      MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

      VecSetValues(b,linEQSize, &idxVec[0], &bValues[0],INSERT_VALUES);

      //MatView(A,viewerSEQ);
      //VecView(b,viewerSEQ);

      // Initialize Solver.
      int its;
      PC pc;
      KSPCreate(PETSC_COMM_SELF, &ksp);
      KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
      KSPGetPC(ksp, &pc);
      PCSetType(pc,"lu"); // if one want to use lu-factorizing PCSetType(pc,"lu");
      //  PCLUSetUseInPlace(pc); // destroy the original matrix
      KSPSetType(ksp,"preonly"); // direct solving

      KSPSolve(ksp,b,x);
      VecGetArray(x, &xCalc);

      //KSPView(ksp,viewerSEQ);
      //VecView(x,viewerSEQ);

      /*********************************************************************/
      // retrieve the polynomial coefficients
      for(int i = 0;i < numOfSegments;i++)

        for(int j = 0;j < pSize;j++)

          cValues[dim][i][j] = xCalc[i * pSize + j];

      VecRestoreArray(x, &xCalc);
      
      // Destroy all petsc objects.
      destroyPETScSolver(ksp);
      destroyPETScVec(x);
      destroyPETScVec(b);
      destroyPETScMat(A);
      
#ifdef _geometryDebugMode_
      bool allZero;
      logFile<<"******************************************************"<<endl;
      for(int i=0;i<AValues.size();i++) {
        allZero = true;
        for(int j=0;j<AValues[i].size();j++) {
          if(AValues[i][j] != 0) allZero = false;
          logFile<<"A["<<i<<"]["<<j<<"] = "<<AValues[i][j]<<endl;
        }
      }
      logFile<<"----------------"<<endl;
      for(int i=0;i<bValues.size();i++) {
        logFile<<"b["<<i<<"] = "<<bValues[i]<<endl;
      }
      if(allZero) {
        logFile<<"An entire coefficient row is zero!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }
#endif

    }

    delete FEMSet,GaussSet;

#ifdef _geometryDebugMode_
    logFile<<"#######################################################"<<endl;
    logFile<<"************ final spline coefficients ****************"<<endl;
    for(int dim=0;dim<usedDims;dim++)
    for(int i=0;i<numOfSegments;i++)
    for(int j=0;j<pSize;j++)
    logFile<<"a["<<i<<"]["<<j<<"] = "<<cValues[dim][i][j]<<endl;
#endif

  }
  
  /**********************************************************************/
  // compute the window function integral
  if(winFuncNorming == 1)

  setCustomPtcleSplineIntegral(InputData,ptcle,modelData,logFile,viewerSEQ);
  
}

/************************************************************************/
/************************************************************************/
// compute the window function integral
void WindowFuncAsym::setCustomPtcleSplineIntegral(
    InputFileData* InputData,Particle& ptcle,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int radiusMethod = (int) InputData->getValue("radiusDeterminationAlgorithm");
  int shapeType = (int) InputData->getValue("shapefunctionType");
  int windowFuncType = (int) InputData->getValue("windowFunctionType");

  dbVector& radii = ptcle.getRadii();
  CustomSpline* spline = ptcle.getSpline();

  double& splineIntegral = spline->getSplineIntegral();
  dbMatrix3& cValues = spline->getSplineCoefficients();
  dbMatrix& splineKnots = spline->getSplineKnots();

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"*********** asymmetric spline integral ***************"<<endl;
#endif

  // ---------------------------------------------------------------------
  // Set the FEM approximation tools, since a FEM element is used 
  // as integration cell.

  int nodesPerElem;
  ElementTemplate* FEMSet;

  intVector data(4);
  map<string,double> params;

  data[0] = 2; // element type (line,rectangle,brick)
  data[1] = 1; // element order
  getFEMMeshData(data,params);

  switch(usedDims) {

  // one dimensional
  case 1:

    nodesPerElem = (int) params["nodesPerLineElement"];
    FEMSet = new Line2ElementTemplate();

    break;

    // two dimensional
  case 2:
    nodesPerElem = (int) params["nodesPerSurfaceElement"];
    FEMSet = new Rect4ElementTemplate();

    break;

    // three dimensional
  case 3:
    nodesPerElem = (int) params["nodesPerVolumeElement"];
    FEMSet = new Cube8ElementTemplate();

    break;

  default:
    logFile << "In WindowFuncAsym1::setCustomPtcleSplineIntegral "
        << "dimension '" << usedDims << "' is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /*********************************************************************/
  // set Gauss integration tools
  int numOfIntPoints;
  GaussPointSet* GaussSet;
  std::vector<GaussPoint> integrationPoints;

  data[0] = 2; // element type (line,rectangle,brick)
  data[1] = integrationOrder; // volume integration order
  data[2] = integrationOrder; // surface integration order 
  data[3] = integrationOrder; // line integration order

  getGaussQuadratureData(data,params);

  switch(usedDims) {

  // one dimensional
  case 1:

    numOfIntPoints = (int) params["gaussPointsPerLineElement"];
    integrationPoints.resize(numOfIntPoints);

    switch(numOfIntPoints) {

    case 1:
      GaussSet = new GaussSetLine1();
      break;
      
    case 2:
      GaussSet = new GaussSetLine2();
      break;
      
    case 3:
      GaussSet = new GaussSetLine3();
      break;
      
    case 4:
      GaussSet = new GaussSetLine4();
      break;

    case 5:
      GaussSet = new GaussSetLine5();
      break;

    default:
      logFile << "In WindowFuncAsym1::setCustomPtcleSplineIntegral "
          << numOfIntPoints << " line Gauss points are not supported!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    break;

    // two dimensional
  case 2:

    numOfIntPoints = (int) params["gaussPointsPerSurfaceElement"];
    integrationPoints.resize(numOfIntPoints);
    
    switch(numOfIntPoints) {
    
    case 1:
      GaussSet = new GaussSetRect1();
      break;

    case 4:
      GaussSet = new GaussSetRect4();
      break;

    case 9:
      GaussSet = new GaussSetRect9();
      break;
      
    case 16:
      GaussSet = new GaussSetRect16();
      break;
      
    case 25:
      GaussSet = new GaussSetRect25();
      break;
      
    default:
      logFile << "In WindowFuncAsym1::setCustomPtcleSplineIntegral "
          << numOfIntPoints << " line Gauss points are not supported!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    break;

    // three dimensional
  case 3:

    numOfIntPoints = (int) params["gaussPointsPerVolumeElement"];
    integrationPoints.resize(numOfIntPoints);

    switch(numOfIntPoints) {

    case 1:
      GaussSet = new GaussSetCube1();
      break;
      
    case 8:
      GaussSet = new GaussSetCube8();
      break;
      
    case 27:
      GaussSet = new GaussSetCube27();
      break;

    case 64:
      GaussSet = new GaussSetCube64();
      break;

    case 125:
      GaussSet = new GaussSetCube125();
      break;

    default:
      logFile << "In WindowFuncAsym1::setCustomPtcleSplineIntegral "
          << numOfIntPoints << " volume Gauss points are not supported!"
          << endl;
      MPI_Abort(MPI_COMM_WORLD,1);

      break;
    }

    break;
    
  default:
    logFile << "In WindowFuncAsym1::setCustomPtcleSplineIntegral "
        << "dimension '" << usedDims << "' is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /**********************************************************************/
  // Calculate the nodal coordinates of the FEM element
  intMatrix quadrant;
  allocateArray(quadrant,8,3);

  int m = 0;

  // loop over all 8 quadrants
  for(int q = 0;q < quadrant.size();q++) {

    switch(q) {

    // (1) x > 0, y > 0, z > 0
    case 0:
      quadrant[m][0] = 0;
      quadrant[m][1] = 0;
      quadrant[m][2] = 0;

      break;

      // (2) x < 0, y > 0, z > 0
    case 1:

      quadrant[m][0] = usedDims;
      quadrant[m][1] = 0;
      quadrant[m][2] = 0;

      break;

      // (3) x > 0, y < 0, z > 0
    case 2:
      quadrant[m][0] = 0;
      quadrant[m][1] = usedDims;
      quadrant[m][2] = 0;

      break;

      // (4) x < 0, y < 0, z > 0
    case 3:
      quadrant[m][0] = usedDims;
      quadrant[m][1] = usedDims;
      quadrant[m][2] = 0;

      break;

      // (5) x > 0, y > 0, z < 0
    case 4:
      quadrant[m][0] = 0;
      quadrant[m][1] = 0;
      quadrant[m][2] = usedDims;

      break;

      // (6) x < 0, y > 0, z < 0
    case 5:
      quadrant[m][0] = usedDims;
      quadrant[m][1] = 0;
      quadrant[m][2] = usedDims;

      break;

      // (7) x > 0, y < 0, z < 0
    case 6:
      quadrant[m][0] = 0;
      quadrant[m][1] = usedDims;
      quadrant[m][2] = usedDims;

      break;

      // (8) x < 0, y < 0, z < 0
    case 7:
      quadrant[m][0] = usedDims;
      quadrant[m][1] = usedDims;
      quadrant[m][2] = usedDims;

      break;
    }

    // check whether quadrant exists
    if(radii[quadrant[m][0] + 0] != 0 && radii[quadrant[m][1] + 1] != 0
      && radii[quadrant[m][2] + 2] != 0)

    m++;
  }

  resizeArray(quadrant,m);

  // ---------------------------------------------------------------------

  vector<Particle> nodes(quadrant.size() * nodesPerElem,Particle(0));

  // loop over all eight quadrants
  for(int q = 0;q < quadrant.size();q++) {

    for(int i = 0;i < nodesPerElem;i++) {
      dbVector& nCoords = nodes[q * nodesPerElem + i].getCoords();
      
      // loop over all dimensions
      for(int j = 0;j < usedDims;j++) {

        // (+1) - positive leg
        if(FEMSet->nodalCoords[i][j] > 0 && quadrant[q][j] == 0) nCoords[j] =
          FEMSet->nodalCoords[i][j] * radii[j];

        // (+1) - negative leg
        else if(FEMSet->nodalCoords[i][j] > 0 && quadrant[q][j] > 0) nCoords[j] =
          0;

        // (-1) - positive leg
        else if(FEMSet->nodalCoords[i][j] < 0 && quadrant[q][j] == 0) nCoords[j] =
          0;

        // (-1) - negative leg
        else if(FEMSet->nodalCoords[i][j] < 0 && quadrant[q][j] > 0) nCoords[j] =
          FEMSet->nodalCoords[i][j] * radii[quadrant[q][j] + j];

      }

    }

  }

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"***************** cube integration ******************"<<endl;
  logFile<<"integration order="<<integrationOrder<<endl;
  for(int i=0;i<radii.size();i++)
  logFile<<"radii["<<i<<"] = "<<radii[i]<<endl;
  logFile<<"************ integration cell nodes ****************"<<endl;
  for(int q=0;q<quadrant.size();q++) {
    logFile<<"Quadrant "<<q<<":"<<endl;
    for(int i=0;i<nodesPerElem;i++) {
      dbVector& coords = nodes[q*nodesPerElem+i].getCoords();
      logFile<<"node "<<i<<" coords: ";
      for(int j=0;j<coords.size();j++)
      logFile<<coords[j]<<" ";
      logFile<<endl;
    }
  }
#endif

  /**********************************************************************/
  // Calculate the shape function ordinates of all nodes at all 
  // integration points.
  dbMatrix sFuncs(numOfIntPoints,dbVector(nodesPerElem));
  
  allocateArray(sFuncs,numOfIntPoints,nodesPerElem);

  for(int i = 0;i < numOfIntPoints;i++)

    for(int j = 0;j < nodesPerElem;j++)

      sFuncs[i][j] = FEMSet->N(j,GaussSet->coord[i]);

#ifdef _geometryDebugMode_
  logFile<<"************ shape functions ****************"<<endl;
  for(int i=0;i<integrationPoints.size();i++) {
    logFile<<"gPoint "<<i<<": ";
    for(int j=0;j<nodesPerElem;j++)
    logFile<<sFuncs[i][j]<<" ";
    logFile<<endl;
  }
#endif

  // Calculate the coordinates of the integration points.

  integrationPoints.resize(numOfIntPoints * quadrant.size());

  // loop over all eight quadrants
  for(int q = 0;q < quadrant.size();q++) {

    for(int i = 0;i < numOfIntPoints;i++) {
      dbVector& gCoords = integrationPoints[q * numOfIntPoints + i].getCoords();
      gCoords.resize(usedDims);
      
      for(int j = 0;j < usedDims;j++) {
        gCoords[j] = 0;

        for(int k = 0;k < nodesPerElem;k++) {
          dbVector& nCoords = nodes[q * nodesPerElem + k].getCoords();

          gCoords[j] += sFuncs[i][k] * nCoords[j];

        }

      }

    }
    
  }
  
#ifdef _geometryDebugMode_
  logFile<<"************** integration points ******************"<<endl;
  for(int i=0;i<integrationPoints.size();i++) {
    dbVector& coords = integrationPoints[i].getCoords();
    logFile<<"integration point "<<i<<" coords: ";
    for(int j=0;j<usedDims;j++)
    logFile<<coords[j]<<" ";
    logFile<<endl;
  }
  logFile<<"----------------------------------------------------"<<endl;
#endif

  /***********************************************************************/
  // Calculate the weight of all integration points.
  intVector nodesIdx(nodesPerElem);
  
  // loop over all eight quadrants
  for(int q = 0;q < quadrant.size();q++) {

    // loop over quadrant's nodes
    for(int i = 0;i < nodesPerElem;i++)
      nodesIdx[i] = q * nodesPerElem + i + 1;

    // loop over quadrants integration points
    for(int i = 0;i < numOfIntPoints;i++) {
      
      double& weight = integrationPoints[q * numOfIntPoints + i].getWeight();
      
      weight = FEMSet->getMetricFactor(nodes,nodesIdx,GaussSet->coord[i],
                                       logFile)

      * GaussSet->weight[i];
    }

  }

#ifdef _geometryDebugMode_
  logFile<<"************** integration points ******************"<<endl;
  for(int i=0;i<integrationPoints.size();i++) {
    logFile<<"integration point "<<i<<" weight: "
    <<integrationPoints[i].getWeight()<<endl;
  }
  logFile<<"----------------------------------------------------"<<endl;
#endif

  delete FEMSet,GaussSet;

  /**********************************************************************/
  // compute the spline volume
  double ordinate,volume;
  double totalWeight = 0;

  splineIntegral = 0;

  for(int k = 0;k < integrationPoints.size();k++) {
    double& weight = integrationPoints[k].getWeight();
    dbVector& coords = integrationPoints[k].getCoords();
    
    calcSplineValue(InputData,ptcle,coords,ordinate,modelData,logFile,
                    viewerSEQ);
    splineIntegral += ordinate * weight;

#ifdef _geometryDebugMode_
    totalWeight += weight;
    logFile<<"----------"<<endl;
    logFile<<"integration point "<<k<<" coords: ";
    logFile<<coords[0]<<" ";
    logFile<<coords[1]<<" ";
    logFile<<coords[2]<<endl;
    logFile<<"weight = "<<totalWeight<<" ";
    logFile<<"ordinate = "<<ordinate<<" ";
    logFile<<"splineVolume = "<<splineIntegral<<endl;
#endif

  }

}

/************************************************************************/
/************************************************************************/
// compute the window function ordinate for local point x
void WindowFuncAsym::calcSplineValue(InputFileData* InputData,Particle& ptcle,
                                     dbVector& x,double& w,
                                     std::map<std::string,double>& modelData,
                                     std::ofstream& logFile,
                                     PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int windowFuncType = (int) InputData->getValue("windowFunctionType");

  CustomSpline* spline = ptcle.getSpline();
  dbMatrix3& cValues = spline->getSplineCoefficients();
  dbMatrix& splineKnots = spline->getSplineKnots();

  double xValue;
  dbVector yValue(3);
  dbVector values;

  // spline is not set yet
  if( !spline->splineIsSet) {

    spline->splineIsSet = true;
    setCustomPtcleSpline(InputData,ptcle,modelData,logFile,viewerSEQ);

  }

  if(windowFuncType == 1) {

    dbVector xnorm(x.size());

    if(x[0] >= 0) xnorm[0] = x[0] / ptcle.getRadius(0);

    else xnorm[0] = x[0] / ptcle.getRadius(usedDims);

    if(x[1] >= 0) xnorm[1] = x[1] / ptcle.getRadius(1);

    else xnorm[1] = x[1] / ptcle.getRadius(usedDims + 1);

    if(x[2] >= 0) xnorm[2] = x[2] / ptcle.getRadius(2);

    else xnorm[2] = x[2] / ptcle.getRadius(usedDims + 2);

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {
      yValue[dim] = cubicSpline(xnorm[dim]);
    }

  }
  else if(windowFuncType == 2) {

    int pSize = P(0).size();

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {

      yValue[dim] = 0;

      values = P(x[dim]);

      // loop over all segments and check where x lies
      for(int i = 0;i < splineKnots[dim].size() - 1;i++) { // < 9-1 = 8
        if(x[dim] >= splineKnots[dim][i] && x[dim] <= splineKnots[dim][i + 1]) { //7-8

          for(int j = 0;j < pSize;j++)
            yValue[dim] += values[j] * cValues[dim][i][j];

          break;

        }

        else if(i + 1 == splineKnots[dim].size() - 1) {
          cerr << "In WindowFuncAsym1::calcSplineValue point lies\n"
              << "outside spline!" << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }
      
    }

  }

  else if(windowFuncType == 4) {

    int pSize = expP(0).size();

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {

      yValue[dim] = 0;

      // loop over all segments and check where x lies
      for(int i = 0;i < splineKnots[dim].size() - 1;i++) { // < 9-1 = 8

        if(x[dim] >= splineKnots[dim][i] && x[dim] <= splineKnots[dim][i + 1]) { //7-8

          values = expP(x[dim] - splineKnots[dim][i]);

          for(int j = 0;j < pSize;j++)
            yValue[dim] += values[j] * cValues[dim][i][j];

          break;

        }

        else if(i + 1 == splineKnots[dim].size() - 1) {
          cerr << "In WindowFuncAsym1::calcSplineValue point lies\n"
              << "outside spline!" << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }
      
    }

  }

  else {
    logFile << "WindowFuncAsym1::calcSplineValue windowfunction type "
        << windowFuncType << " is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  w = yValue[0] * yValue[1] * yValue[2];

}

/************************************************************************/
/************************************************************************/
// compute the window function ordinate and its for first order 
// derivatives for local point x
void WindowFuncAsym::calcSplineValue(InputFileData* InputData,Particle& ptcle,
                                     dbVector& x,double& w,double& dxW,
                                     double& dyW,double& dzW,
                                     std::map<std::string,double>& modelData,
                                     std::ofstream& logFile,
                                     PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int windowFuncType = (int) InputData->getValue("windowFunctionType");

  CustomSpline* spline = ptcle.getSpline();

  //double& splineIntegral = spline->getSplineIntegral();
  dbMatrix3& cValues = spline->getSplineCoefficients();
  dbMatrix& splineKnots = spline->getSplineKnots();

  vector<double> yValue(3);
  vector<double> d1yValue(3);
  vector<double> values,d1values;

  // spline is not set yet
  if( !spline->splineIsSet) {

    spline->splineIsSet = true;
    setCustomPtcleSpline(InputData,ptcle,modelData,logFile,viewerSEQ);

  }

  if(windowFuncType == 1) {

    dbVector xnorm(x.size());

    if(x[0] >= 0) xnorm[0] = x[0] / ptcle.getRadius(0);

    else xnorm[0] = x[0] / ptcle.getRadius(usedDims);

    if(x[1] >= 0) xnorm[1] = x[1] / ptcle.getRadius(1);

    else xnorm[1] = x[1] / ptcle.getRadius(usedDims + 1);

    if(x[2] >= 0) xnorm[2] = x[2] / ptcle.getRadius(2);

    else xnorm[2] = x[2] / ptcle.getRadius(usedDims + 2);

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {
      yValue[dim] = cubicSpline(xnorm[dim]);
      d1yValue[dim] = dcubicSpline(xnorm[dim]);
    }

  }
  else if(windowFuncType == 2) {

    int pSize = P(0).size();

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {

      yValue[dim] = 0;
      d1yValue[dim] = 0;

      values = P(x[dim]);
      d1values = dP(x[dim]);

      // loop over all segments and check where x lies
      for(int i = 0;i < splineKnots[dim].size() - 1;i++) { // < 9-1 = 8

        if(x[dim] >= splineKnots[dim][i] && x[dim] <= splineKnots[dim][i + 1]) { //7-8

          for(int j = 0;j < pSize;j++) {
            yValue[dim] += values[j] * cValues[dim][i][j];
            d1yValue[dim] += d1values[j] * cValues[dim][i][j];
          }

          break;

        }

        else if(i + 1 == splineKnots[dim].size() - 1) {
          cerr << "In WindowFuncAsym1::calcSplineValue point lies\n"
              << "outside spline!" << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }
      
    }

  }

  else if(windowFuncType == 4) {

    int pSize = expP(0).size();

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {

      yValue[dim] = 0;
      d1yValue[dim] = 0;

      // loop over all segments and check where x lies
      for(int i = 0;i < splineKnots[dim].size() - 1;i++) { // < 9-1 = 8

        if(x[dim] >= splineKnots[dim][i] && x[dim] <= splineKnots[dim][i + 1]) { //7-8

          values = expP(x[dim] - splineKnots[dim][i]);
          d1values = dexpP(x[dim] - splineKnots[dim][i]);

          for(int j = 0;j < pSize;j++) {
            yValue[dim] += values[j] * cValues[dim][i][j];
            d1yValue[dim] += d1values[j] * cValues[dim][i][j];
          }

          break;

        }

        else if(i + 1 == splineKnots[dim].size() - 1) {
          cerr << "In WindowFuncAsym1::calcSplineValue point lies\n"
              << "outside spline!" << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }
      
    }

  }

  else {
    logFile << "WindowFuncAsym1::calcSplineValue windowfunction type "
        << windowFuncType << " is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  w = yValue[0] * yValue[1] * yValue[2];
  dxW = d1yValue[0] * yValue[1] * yValue[2];
  dyW = yValue[0] * d1yValue[1] * yValue[2];
  dzW = yValue[0] * yValue[1] * d1yValue[2];
}

/************************************************************************/
/************************************************************************/
// compute the window function ordinate and its for first-and second-order 
// derivatives for local point x
void WindowFuncAsym::calcSplineValue(InputFileData* InputData,Particle& ptcle,
                                     dbVector& x,double& w,double& dxW,
                                     double& dyW,double& dzW,double& dxxW,
                                     double& dyyW,double& dzzW,double& dxyW,
                                     double& dyzW,double& dzxW,
                                     std::map<std::string,double>& modelData,
                                     std::ofstream& logFile,
                                     PetscViewer& viewerSEQ) {
  
  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int windowFuncType = (int) InputData->getValue("windowFunctionType");

  CustomSpline* spline = ptcle.getSpline();

  //double& splineIntegral = spline->getSplineIntegral();
  dbMatrix3& cValues = spline->getSplineCoefficients();
  dbMatrix& splineKnots = spline->getSplineKnots();

  vector<double> yValue(3);
  vector<double> d1yValue(3);
  vector<double> d2yValue(3);
  vector<double> values,d1values,d2values;

  // spline is not set yet
  if( !spline->splineIsSet) {

    spline->splineIsSet = true;
    setCustomPtcleSpline(InputData,ptcle,modelData,logFile,viewerSEQ);

  }

  if(windowFuncType == 1) {

    dbVector xnorm(x.size());

    if(x[0] >= 0) xnorm[0] = x[0] / ptcle.getRadius(0);

    else xnorm[0] = x[0] / ptcle.getRadius(usedDims);

    if(x[1] >= 0) xnorm[1] = x[1] / ptcle.getRadius(1);

    else xnorm[1] = x[1] / ptcle.getRadius(usedDims + 1);

    if(x[2] >= 0) xnorm[2] = x[2] / ptcle.getRadius(2);

    else xnorm[2] = x[2] / ptcle.getRadius(usedDims + 2);

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {
      yValue[dim] = cubicSpline(xnorm[dim]);
      d1yValue[dim] = dcubicSpline(xnorm[dim]);
      d2yValue[dim] = d2cubicSpline(xnorm[dim]);
    }

  }
  else if(windowFuncType == 2) {

    int pSize = P(0).size();

    // loop over all 3 dimensions
    for(int dim = 0;dim < 3;dim++) {

      yValue[dim] = 0;
      d1yValue[dim] = 0;
      d2yValue[dim] = 0;

      values = P(x[dim]);
      d1values = dP(x[dim]);
      d2values = d2P(x[dim]);

      // loop over all segments and check where x lies
      for(int i = 0;i < splineKnots[dim].size() - 1;i++) { // < 9-1 = 8

        if(x[dim] >= splineKnots[dim][i] && x[dim] <= splineKnots[dim][i + 1]) { //7-8

          for(int j = 0;j < pSize;j++) {
            yValue[dim] += values[j] * cValues[dim][i][j];
            d1yValue[dim] += d1values[j] * cValues[dim][i][j];
            d2yValue[dim] += d2values[j] * cValues[dim][i][j];
          }

          break;

        }

        else if(i + 1 == splineKnots[dim].size() - 1) {
          cerr << "In WindowFuncAsym1::calcSplineValue point lies\n"
              << "outside spline!" << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }
      
    }

  }

  else if(windowFuncType == 4) {

    int pSize = expP(0).size();

    // loop over all 3 dimensions

    for(int dim = 0;dim < 3;dim++) {

      yValue[dim] = 0;
      d1yValue[dim] = 0;
      d2yValue[dim] = 0;

      // loop over all segments and check where x lies
      for(int i = 0;i < splineKnots[dim].size() - 1;i++) { // < 9-1 = 8

        if(x[dim] >= splineKnots[dim][i] && x[dim] <= splineKnots[dim][i + 1]) { //7-8

          values = expP(x[dim] - splineKnots[dim][i]);
          d1values = dexpP(x[dim] - splineKnots[dim][i]);
          d2values = d2expP(x[dim] - splineKnots[dim][i]);

          for(int j = 0;j < pSize;j++) {
            yValue[dim] += values[j] * cValues[dim][i][j];
            d1yValue[dim] += d1values[j] * cValues[dim][i][j];
            d2yValue[dim] += d2values[j] * cValues[dim][i][j];
          }

          break;

        }

        else if(i + 1 == splineKnots[dim].size() - 1) {
          cerr << "In WindowFuncAsym1::calcSplineValue point lies\n"
              << "outside spline!" << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }
      
    }

  }

  else {
    logFile << "WindowFuncAsym1::calcSplineValue windowfunction type "
        << windowFuncType << " is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  w = yValue[0] * yValue[1] * yValue[2];

  dxW = d1yValue[0] * yValue[1] * yValue[2];
  dyW = yValue[0] * d1yValue[1] * yValue[2];
  dzW = yValue[0] * yValue[1] * d1yValue[2];

  dxxW = d2yValue[0] * yValue[1] * yValue[2];
  dyyW = yValue[0] * d2yValue[1] * yValue[2];
  dzzW = yValue[0] * yValue[1] * d2yValue[2];
  dxyW = d1yValue[0] * d1yValue[1] * yValue[2];
  dyzW = yValue[0] * d1yValue[1] * d1yValue[2];
  dzxW = yValue[0] * d1yValue[1] * d1yValue[2];
}

/***********************************************************************/
/***********************************************************************/
// Calculation of the window functions 
void WindowFuncAsym::calcWinFunctions(InputFileData* InputData,
                                      std::vector<Particle>& ptcls,
                                      intVector& sPtcls,double& x,double& y,
                                      double& z,int& supportSize,
                                      std::map<std::string,double>& modelData,
                                      std::ofstream& logFile,
                                      PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int winFuncNorming = (int) InputData->getValue("windowfunctionNorming");
  int windowFuncType = (int) InputData->getValue("windowFunctionType");
  int particlesNum = ptcls.size();

  dbVector xnorm(3);
  int m;
  double normFactor;

  windowFuncs = dbVector(supportSize);

  for(int i = 0;i < supportSize;i++) {
    m = sPtcls[i];

    xnorm[0] = x - ptcls[m].getCoord(0);
    xnorm[1] = y - ptcls[m].getCoord(1);
    xnorm[2] = z - ptcls[m].getCoord(2);
    
    calcSplineValue(InputData,ptcls[m],xnorm,windowFuncs[i],modelData,logFile,
                    viewerSEQ);

    // norming by 3D domain of influence
    if(winFuncNorming == 1) {
      CustomSpline* spline = ptcls[m].getSpline();
      windowFuncs[i] /= spline->getSplineIntegral();
    }

    // norming by particle weight obtained from FEM mesh
    else if(winFuncNorming == 2) windowFuncs[i] /= ptcls[m].getWeight();
    
  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of window function and its first order derivatives
void WindowFuncAsym::calcWinFunction1stDerivs(
    InputFileData* InputData,std::vector<Particle>& ptcls,intVector& sPtcls,
    double& x,double& y,double& z,int& supportSize,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int winFuncNorming = (int) InputData->getValue("windowfunctionNorming");
  int windowFuncType = (int) InputData->getValue("windowFunctionType");
  int particlesNum = ptcls.size();

  dbVector xnorm(3),radius(3);
  int m;
  double normFactor;

  windowFuncs = dbVector(supportSize);

  xDerivWinFuncs = dbVector(supportSize);
  yDerivWinFuncs = dbVector(supportSize);
  zDerivWinFuncs = dbVector(supportSize);

  for(int i = 0;i < supportSize;i++) {
    m = sPtcls[i];

    xnorm[0] = x - ptcls[m].getCoord(0);
    xnorm[1] = y - ptcls[m].getCoord(1);
    xnorm[2] = z - ptcls[m].getCoord(2);
    
    calcSplineValue(InputData,ptcls[m],xnorm,windowFuncs[i],xDerivWinFuncs[i],
                    yDerivWinFuncs[i],zDerivWinFuncs[i],modelData,logFile,
                    viewerSEQ);

    // norming by 3D domain of influence
    if(winFuncNorming == 1 && windowFuncType != 1) {
      CustomSpline* spline = ptcls[m].getSpline();

      windowFuncs[i] /= spline->getSplineIntegral();
      xDerivWinFuncs[i] /= spline->getSplineIntegral();
      yDerivWinFuncs[i] /= spline->getSplineIntegral();
      zDerivWinFuncs[i] /= spline->getSplineIntegral();
    }
    else if(winFuncNorming == 1 && windowFuncType == 1) {
      CustomSpline* spline = ptcls[m].getSpline();
      dbVector& radii = ptcls[m].getRadii();

      if(xnorm[0] > 0) radius[0] = radii[0];
      else radius[0] = radii[usedDims];

      if(xnorm[1] > 0) radius[1] = radii[1];
      else radius[1] = radii[usedDims + 1];

      if(xnorm[2] > 0) radius[2] = radii[2];
      else radius[2] = radii[usedDims + 2];

      windowFuncs[i] /= spline->getSplineIntegral();
      xDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[0];
      yDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[1];
      zDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[2];
    }

    // norming by particle weight obtained from FEM mesh
    else if(winFuncNorming == 2 && windowFuncType != 1) {
      windowFuncs[i] /= ptcls[m].getWeight();
      xDerivWinFuncs[i] /= ptcls[m].getWeight();
      yDerivWinFuncs[i] /= ptcls[m].getWeight();
      zDerivWinFuncs[i] /= ptcls[m].getWeight();
    }

    // norming by particle weight obtained from FEM mesh
    else if(winFuncNorming == 2 && windowFuncType == 1) {
      windowFuncs[i] /= ptcls[m].getWeight();
      xDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[0];
      yDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[1];
      zDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[2];
    }

    else if(windowFuncType == 1) {
      xDerivWinFuncs[i] /= radius[0];
      yDerivWinFuncs[i] /= radius[1];
      zDerivWinFuncs[i] /= radius[2];
    }

  }

}

/***********************************************************************/
/***********************************************************************/
// Calculation of window function and its first and second order 
// derivatives
void WindowFuncAsym::calcWinFunction2ndDerivs(
    InputFileData* InputData,std::vector<Particle>& ptcls,intVector& sPtcls,
    double& x,double& y,double& z,int& supportSize,
    std::map<std::string,double>& modelData,std::ofstream& logFile,
    PetscViewer& viewerSEQ) {
  
  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int winFuncNorming = (int) InputData->getValue("windowfunctionNorming");
  int windowFuncType = (int) InputData->getValue("windowFunctionType");
  int particlesNum = ptcls.size();

  dbVector xnorm(3),radius(3);
  int m;

  windowFuncs = dbVector(supportSize);

  xDerivWinFuncs = dbVector(supportSize);
  yDerivWinFuncs = dbVector(supportSize);
  zDerivWinFuncs = dbVector(supportSize);

  xDerivWinFuncs = dbVector(supportSize);
  yDerivWinFuncs = dbVector(supportSize);
  zDerivWinFuncs = dbVector(supportSize);

  xxDerivWinFuncs = dbVector(supportSize);
  yyDerivWinFuncs = dbVector(supportSize);
  zzDerivWinFuncs = dbVector(supportSize);

  xyDerivWinFuncs = dbVector(supportSize);
  yzDerivWinFuncs = dbVector(supportSize);
  zxDerivWinFuncs = dbVector(supportSize);

  for(int i = 0;i < supportSize;i++) {
    m = sPtcls[i];

    xnorm[0] = x - ptcls[m].getCoord(0);
    xnorm[1] = y - ptcls[m].getCoord(1);
    xnorm[2] = z - ptcls[m].getCoord(2);
    
    calcSplineValue(InputData,ptcls[m],xnorm,windowFuncs[i],xDerivWinFuncs[i],
                    yDerivWinFuncs[i],zDerivWinFuncs[i],xxDerivWinFuncs[i],
                    yyDerivWinFuncs[i],zzDerivWinFuncs[i],xyDerivWinFuncs[i],
                    yzDerivWinFuncs[i],zxDerivWinFuncs[i],modelData,logFile,
                    viewerSEQ);

    CustomSpline* spline = ptcls[m].getSpline();

    // norming by 3D domain of influence
    if(winFuncNorming == 1 && windowFuncType != 1) {
      CustomSpline* spline = ptcls[m].getSpline();

      windowFuncs[i] /= spline->getSplineIntegral();
      xDerivWinFuncs[i] /= spline->getSplineIntegral();
      yDerivWinFuncs[i] /= spline->getSplineIntegral();
      zDerivWinFuncs[i] /= spline->getSplineIntegral();
      xxDerivWinFuncs[i] /= spline->getSplineIntegral();
      yyDerivWinFuncs[i] /= spline->getSplineIntegral();
      zzDerivWinFuncs[i] /= spline->getSplineIntegral();
      xyDerivWinFuncs[i] /= spline->getSplineIntegral();
      yzDerivWinFuncs[i] /= spline->getSplineIntegral();
      zxDerivWinFuncs[i] /= spline->getSplineIntegral();
    }
    else if(winFuncNorming == 1 && windowFuncType == 1) {
      CustomSpline* spline = ptcls[m].getSpline();
      dbVector& radii = ptcls[m].getRadii();

      if(xnorm[0] > 0) radius[0] = radii[0];
      else radius[0] = radii[usedDims];

      if(xnorm[1] > 0) radius[1] = radii[1];
      else radius[1] = radii[usedDims + 1];

      if(xnorm[2] > 0) radius[2] = radii[2];
      else radius[2] = radii[usedDims + 2];

      windowFuncs[i] /= spline->getSplineIntegral();
      xDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[0];
      yDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[1];
      zDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[2];
      xxDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[0] * radius[0];
      yyDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[1] * radius[1];
      zzDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[2] * radius[2];
      xyDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[0] * radius[1];
      yzDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[1] * radius[2];
      zxDerivWinFuncs[i] /= spline->getSplineIntegral() * radius[2] * radius[0];
    }

    // norming by particle weight obtained from FEM mesh
    else if(winFuncNorming == 2 && windowFuncType != 1) {
      windowFuncs[i] /= ptcls[m].getWeight();
      xDerivWinFuncs[i] /= ptcls[m].getWeight();
      yDerivWinFuncs[i] /= ptcls[m].getWeight();
      zDerivWinFuncs[i] /= ptcls[m].getWeight();
      xxDerivWinFuncs[i] /= ptcls[m].getWeight();
      yyDerivWinFuncs[i] /= ptcls[m].getWeight();
      zzDerivWinFuncs[i] /= ptcls[m].getWeight();
      xyDerivWinFuncs[i] /= ptcls[m].getWeight();
      yzDerivWinFuncs[i] /= ptcls[m].getWeight();
      zxDerivWinFuncs[i] /= ptcls[m].getWeight();
    }

    // norming by particle weight obtained from FEM mesh
    else if(winFuncNorming == 2 && windowFuncType == 1) {
      windowFuncs[i] /= ptcls[m].getWeight();
      xDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[0];
      yDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[1];
      zDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[2];
      xxDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[0] * radius[0];
      yyDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[1] * radius[1];
      zzDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[2] * radius[2];
      xyDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[0] * radius[1];
      yzDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[1] * radius[2];
      zxDerivWinFuncs[i] /= ptcls[m].getWeight() * radius[2] * radius[0];
    }

    else if(windowFuncType == 1) {
      xDerivWinFuncs[i] /= radius[0];
      yDerivWinFuncs[i] /= radius[1];
      zDerivWinFuncs[i] /= radius[2];
      xxDerivWinFuncs[i] /= radius[0] * radius[0];
      yyDerivWinFuncs[i] /= radius[1] * radius[1];
      zzDerivWinFuncs[i] /= radius[2] * radius[2];
      xyDerivWinFuncs[i] /= radius[0] * radius[1];
      yzDerivWinFuncs[i] /= radius[1] * radius[2];
      zxDerivWinFuncs[i] /= radius[2] * radius[0];
    }

  }

}

