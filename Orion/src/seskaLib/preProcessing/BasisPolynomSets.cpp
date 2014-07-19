// Return a vector containing the polynom basis for a shapefuntion's set.

// Polynom-Type (Pascal) number of terms l = (p + 1)*(p + 2)*(p + 3)

#include "BasisPolynomSets.h"

dbVector PascalLinear::P(double& x,double& y,double& z){

  using namespace std;

  dbVector basis(4);

  basis[0] = 1;
  basis[1] = x;
  basis[2] = y;
  basis[3] = z;

  return (basis);
}

dbVector PascalLinear::dPx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(4);

  basis[0] = 0;
  basis[1] = 1.0;
  basis[2] = 0;
  basis[3] = 0;

  return (basis);
}

dbVector PascalLinear:: dPy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(4);

  basis[0] = 0;
  basis[1] = 0;
  basis[2] = 1.0;
  basis[3] = 0;

  return (basis);
}

dbVector PascalLinear:: dPz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(4);

  basis[0] = 0;
  basis[1] = 0;
  basis[2] = 0;
  basis[3] = 1.0;

  return (basis);
}

/***********************************************************************/
/***********************************************************************/

dbVector PascalQuadratic::P(double& x,double& y,double& z){

  using namespace std;

  dbVector basis(10);

  basis[0] = 1;
  basis[1] = x;
  basis[2] = y;
  basis[3] = z;
  basis[4] = x*y;
  basis[5] = y*z;
  basis[6] = z*x;
  basis[7] = pow(x,2);
  basis[8] = pow(y,2);
  basis[9] = pow(z,2);

  return (basis);
}

/***********************************************************************/
dbVector PascalQuadratic::dPx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 1;       // x
  basis[2] = 0;       // y
  basis[3] = 0;       // z
  basis[4] = y;       // xy
  basis[5] = 0;       // yz
  basis[6] = z;       // zx
  basis[7] = 2.0*x;   // x^2
  basis[8] = 0;       // y^2
  basis[9] = 0;       // z^2

  return (basis);
}

dbVector PascalQuadratic:: dPy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 1;       // y
  basis[3] = 0;       // z
  basis[4] = x;       // xy
  basis[5] = z;       // yz
  basis[6] = 0;       // zx
  basis[7] = 0;       // x^2
  basis[8] = 2.0*y;   // y^2
  basis[9] = 0;       // z^2

  return (basis);
}

dbVector PascalQuadratic:: dPz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 0;       // y
  basis[3] = 1;       // z
  basis[4] = 0;       // xy
  basis[5] = y;       // yz
  basis[6] = x;       // zx
  basis[7] = 0;       // x^2
  basis[8] = 0;       // y^2
  basis[9] = 2.0*z;   // z^2

  return (basis);
}

/***********************************************************************/
dbVector PascalQuadratic::dPxx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 0;       // y
  basis[3] = 0;       // z
  basis[4] = 0;       // xy
  basis[5] = 0;       // yz
  basis[6] = 0;       // zx
  basis[7] = 2;       // x^2
  basis[8] = 0;       // y^2
  basis[9] = 0;       // z^2

  return (basis);
}

dbVector PascalQuadratic:: dPyy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 0;       // y
  basis[3] = 0;       // z
  basis[4] = 0;       // xy
  basis[5] = 0;       // yz
  basis[6] = 0;       // zx
  basis[7] = 0;       // x^2
  basis[8] = 2;       // y^2
  basis[9] = 0;       // z^2

  return (basis);
}

dbVector PascalQuadratic:: dPzz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 0;       // y
  basis[3] = 0;       // z
  basis[4] = 0;       // xy
  basis[5] = 0;       // yz
  basis[6] = 0;       // zx
  basis[7] = 0;       // x^2
  basis[8] = 0;       // y^2
  basis[9] = 2;       // z^2

  return (basis);
}

dbVector PascalQuadratic::dPxy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 0;       // y
  basis[3] = 0;       // z
  basis[4] = 1;       // xy
  basis[5] = 0;       // yz
  basis[6] = 0;       // zx
  basis[7] = 0;       // x^2
  basis[8] = 0;       // y^2
  basis[9] = 0;       // z^2

  return (basis);
}

dbVector PascalQuadratic:: dPyz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 0;       // y
  basis[3] = 0;       // z
  basis[4] = 0;       // xy
  basis[5] = 1;       // yz
  basis[6] = 0;       // zx
  basis[7] = 0;       // x^2
  basis[8] = 0;       // y^2
  basis[9] = 0;       // z^2

  return (basis);
}

dbVector PascalQuadratic:: dPzx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(10);

  basis[0] = 0;       // 1
  basis[1] = 0;       // x
  basis[2] = 0;       // y
  basis[3] = 0;       // z
  basis[4] = 0;       // xy
  basis[5] = 0;       // yz
  basis[6] = 1;       // zx
  basis[7] = 0;       // x^2
  basis[8] = 0;       // y^2
  basis[9] = 0;       // z^2

  return (basis);
}

/***********************************************************************/
/***********************************************************************/

dbVector PascalCubic::P(double& x,double& y,double& z){

  using namespace std;

  dbVector basis(20);

  basis[0] = 1;            // 1
  basis[1] = x;            // x
  basis[2] = y;            // y
  basis[3] = z;            // z
  basis[4] = x*y;           // xy
  basis[5] = y*z;           // yz
  basis[6] = z*x;           // zx
  basis[7] = pow(x,2);     // x^2
  basis[8] = pow(y,2);     // y^2
  basis[9] = pow(z,2);     // z^2
  basis[10] = x*y*z;       // x*y*z
  basis[11] = pow(x,2)*y;  // x^2*y
  basis[12] = pow(x,2)*z;  // x^2*z
  basis[13] = pow(y,2)*x;  // y^2*x
  basis[14] = pow(y,2)*z;  // y^2*z 
  basis[15] = pow(z,2)*x;  // z^2*x 
  basis[16] = pow(z,2)*y;  // z^2*y
  basis[17] = pow(x,3);    // x^3
  basis[18] = pow(y,3);    // y^3
  basis[19] = pow(z,3);    // z^3

  return (basis);
}

/***********************************************************************/
dbVector PascalCubic::dPx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 1;            // x
  basis[2] = 0;            // y
  basis[3] = 0;            // z
  basis[4] = y;            // xy
  basis[5] = 0;            // yz
  basis[6] = z;            // zx
  basis[7] = 2.0*x;        // x^2
  basis[8] = 0;            // y^2
  basis[9] = 0;            // z^2
  basis[10] = y*z;         // x*y*z
  basis[11] = 2.0*x*y;     // x^2*y
  basis[12] = 2.0*x*z;     // x^2*z
  basis[13] = pow(y,2);    // y^2*x
  basis[14] = 0;           // y^2*z 
  basis[15] = pow(z,2);    // z^2*x 
  basis[16] = 0;           // z^2*y
  basis[17] = 3.0*pow(x,2);// x^3
  basis[18] = 0;           // y^3
  basis[19] = 0;           // z^3

  return (basis);
}

dbVector PascalCubic:: dPy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 1;            // y
  basis[3] = 0;            // z
  basis[4] = x;            // xy
  basis[5] = z;            // yz
  basis[6] = 0;            // zx
  basis[7] = 0;            // x^2
  basis[8] = 2.0*y;        // y^2
  basis[9] = 0;            // z^2
  basis[10] = x*z;         // x*y*z
  basis[11] = pow(x,2);    // x^2*y
  basis[12] = 0;           // x^2*z
  basis[13] = 2.0*y*x;     // y^2*x
  basis[14] = 2.0*y*z;     // y^2*z 
  basis[15] = 0;           // z^2*x 
  basis[16] = pow(z,2);    // z^2*y
  basis[17] = 0;           // x^3
  basis[18] = 3.0*pow(y,2);// y^3
  basis[19] = 0;           // z^3

  return (basis);
}

dbVector PascalCubic:: dPz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 0;            // y
  basis[3] = 1;            // z
  basis[4] = 0;            // xy
  basis[5] = y;            // yz
  basis[6] = x;            // zx
  basis[7] = 0;            // x^2
  basis[8] = 0;            // y^2
  basis[9] = 2.0*z;        // z^2
  basis[10] = x*y;         // x*y*z
  basis[11] = 0;           // x^2*y
  basis[12] = pow(x,2);    // x^2*z
  basis[13] = 0;           // y^2*x
  basis[14] = pow(y,2);    // y^2*z 
  basis[15] = 2.0*z*x;     // z^2*x 
  basis[16] = 2.0*z*y;     // z^2*y
  basis[17] = 0;           // x^3
  basis[18] = 0;           // y^3
  basis[19] = 3.0*pow(z,2);// z^3

  return (basis);
}

/***********************************************************************/
dbVector PascalCubic::dPxx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 0;            // y
  basis[3] = 0;            // z
  basis[4] = 0;            // xy
  basis[5] = 0;            // yz
  basis[6] = 0;            // zx
  basis[7] = 2;            // x^2
  basis[8] = 0;            // y^2
  basis[9] = 0;            // z^2
  basis[10] = 0;           // x*y*z
  basis[11] = 2.0*y;       // x^2*y
  basis[12] = 2.0*z;       // x^2*z
  basis[13] = 0;           // y^2*x
  basis[14] = 0;           // y^2*z 
  basis[15] = 0;           // z^2*x 
  basis[16] = 0;           // z^2*y
  basis[17] = 6.0*x;       // x^3
  basis[18] = 0;           // y^3
  basis[19] = 0;           // z^3

  return (basis);
}

dbVector PascalCubic:: dPyy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 0;            // y
  basis[3] = 0;            // z
  basis[4] = 0;            // xy
  basis[5] = 0;            // yz
  basis[6] = 0;            // zx
  basis[7] = 0;            // x^2
  basis[8] = 2;            // y^2
  basis[9] = 0;            // z^2  
  basis[10] = 0;           // x*y*z
  basis[11] = 0;           // x^2*y
  basis[12] = 0;           // x^2*z
  basis[13] = 2.0*x;       // y^2*x
  basis[14] = 2.0*z;       // y^2*z 
  basis[15] = 0;           // z^2*x 
  basis[16] = 0;           // z^2*y
  basis[17] = 0;           // x^3
  basis[18] = 6.0*y;       // y^3
  basis[19] = 0;           // z^3

  return (basis);
}

dbVector PascalCubic:: dPzz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 0;            // y
  basis[3] = 0;            // z
  basis[4] = 0;            // xy
  basis[5] = 0;            // yz
  basis[6] = 0;            // zx
  basis[7] = 0;            // x^2
  basis[8] = 0;            // y^2
  basis[9] = 2;            // z^2
  basis[10] = 0;           // x*y*z
  basis[11] = 0;           // x^2*y
  basis[12] = 0;           // x^2*z
  basis[13] = 0;           // y^2*x
  basis[14] = 0;           // y^2*z 
  basis[15] = 2.0*x;       // z^2*x 
  basis[16] = 2.0*y;       // z^2*y
  basis[17] = 0;           // x^3
  basis[18] = 0;           // y^3
  basis[19] = 6.0*z;       // z^3

  return (basis);
}

dbVector PascalCubic::dPxy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 0;            // y
  basis[3] = 0;            // z
  basis[4] = 1;            // xy
  basis[5] = 0;            // yz
  basis[6] = 0;            // zx
  basis[7] = 0;            // x^2
  basis[8] = 0;            // y^2
  basis[9] = 0;            // z^2
  basis[10] = z;           // x*y*z
  basis[11] = 2.0*x;       // x^2*y
  basis[12] = 0;           // x^2*z
  basis[13] = 2.0*y;       // y^2*x
  basis[14] = 0;           // y^2*z 
  basis[15] = 0;           // z^2*x 
  basis[16] = 0;           // z^2*y
  basis[17] = 0;           // x^3
  basis[18] = 0;           // y^3
  basis[19] = 0;           // z^3

  return (basis);
}

dbVector PascalCubic:: dPyz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 0;            // y
  basis[3] = 0;            // z
  basis[4] = 0;            // xy
  basis[5] = 1;            // yz
  basis[6] = 0;            // zx
  basis[7] = 0;            // x^2
  basis[8] = 0;            // y^2
  basis[9] = 0;            // z^2
  basis[10] = x;           // x*y*z
  basis[11] = 0;           // x^2*y
  basis[12] = 0;           // x^2*z
  basis[13] = 0;           // y^2*x
  basis[14] = 2.0*y;       // y^2*z 
  basis[15] = 0;           // z^2*x 
  basis[16] = 2.0*z;       // z^2*y
  basis[17] = 0;           // x^3
  basis[18] = 0;           // y^3
  basis[19] = 0;           // z^3

  return (basis);
}

dbVector PascalCubic:: dPzx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;            // 1
  basis[1] = 0;            // x
  basis[2] = 0;            // y
  basis[3] = 0;            // z
  basis[4] = 0;            // xy
  basis[5] = 0;            // yz
  basis[6] = 1;            // zx
  basis[7] = 0;            // x^2
  basis[8] = 0;            // y^2
  basis[9] = 0;            // z^2
  basis[10] = y;           // x*y*z
  basis[11] = 0;           // x^2*y
  basis[12] = 2.0*x;       // x^2*z
  basis[13] = 0;           // y^2*x
  basis[14] = 0;           // y^2*z 
  basis[15] = 2.0*z;       // z^2*x 
  basis[16] = 0;           // z^2*y
  basis[17] = 0;           // x^3
  basis[18] = 0;           // y^3
  basis[19] = 0;           // z^3


  return (basis);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

// Polynom-Type 2 (Lagrange) number of terms l = (p + 1)^3

dbVector LagrangeLinear::P(double& x,double& y,double& z){

  using namespace std;

  dbVector basis(8);

  basis[0] = 1;
  basis[1] = x;
  basis[2] = y;
  basis[3] = z;
  basis[4] = x*y;
  basis[5] = y*z;
  basis[6] = z*x;
  basis[7] = x*y*z;

  return (basis);
}

/***********************************************************************/
dbVector LagrangeLinear::dPx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = 0;
  basis[1] = 1.0;
  basis[2] = 0;
  basis[3] = 0;
  basis[4] = y;
  basis[5] = 0;
  basis[6] = z;
  basis[7] = y*z;

  return (basis);
}

dbVector LagrangeLinear:: dPy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = 0;
  basis[1] = 0;
  basis[2] = 1.0;
  basis[3] = 0;
  basis[4] = x;
  basis[5] = z;
  basis[6] = 0;
  basis[7] = x*z;

  return (basis);
}

dbVector LagrangeLinear:: dPz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = 0;
  basis[1] = 0;
  basis[2] = 0;
  basis[3] = 1.0;
  basis[4] = 0;
  basis[5] = y;
  basis[6] = x;
  basis[7] = x*y;

  return (basis);
}

/***********************************************************************/
dbVector LagrangeLinear::dPxy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = 0.0;
  basis[1] = 0.0;
  basis[2] = 0.0;
  basis[3] = 0.0;
  basis[4] = 1.0;
  basis[5] = 0.0;
  basis[6] = 0.0;
  basis[7] = z;

  return (basis);
}

dbVector LagrangeLinear:: dPyz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = 0.0;
  basis[1] = 0.0;
  basis[2] = 0.0;
  basis[3] = 0.0;
  basis[4] = 0.0;
  basis[5] = 1.0;
  basis[6] = 0.0;
  basis[7] = x;

  return (basis);
}

dbVector LagrangeLinear:: dPzx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = 0.0;
  basis[1] = 0.0;
  basis[2] = 0.0;
  basis[3] = 0.0;
  basis[4] = 0.0;
  basis[5] = 0.0;
  basis[6] = 1.0;
  basis[7] = y;

  return (basis);
}

/***********************************************************************/
/***********************************************************************/
dbVector LagrangeQuadratic::P(double& x,double& y,double& z){

  using namespace std;

  dbVector basis(27);

  basis[0] = 1;                           // 1
  basis[1] = x;                           // x
  basis[2] = y;                           // y
  basis[3] = z;                           // z
  basis[4] = x*y;                         // x*y
  basis[5] = y*z;                         // y*z
  basis[6] = z*x;                         // z*x
  basis[7] = x*y*z;                       // x*y*z
  basis[8] = pow(x,2);                    // x^2
  basis[9] = pow(y,2);                    // y^2
  basis[10] = pow(z,2);                   // z^2
  basis[11] = pow(x,2)*y;                 // x^2*y
  basis[12] = pow(x,2)*z;                 // x^2*z
  basis[13] = pow(y,2)*x;                 // y^2*x
  basis[14] = pow(y,2)*z;                 // y^2*z
  basis[15] = pow(z,2)*x;                 // z^2*x
  basis[16] = pow(z,2)*y;                 // z^2*y
  basis[17] = pow(x,2)*y*z;               // x^2*y*z
  basis[18] = pow(y,2)*x*z;               // y^2*x*z
  basis[19] = pow(z,2)*x*y;               // z^2*x*y
  basis[20] = pow(x,2)*pow(y,2);          // x^2*y^2
  basis[21] = pow(y,2)*pow(z,2);          // y^2*z^2
  basis[22] = pow(z,2)*pow(x,2);          // z^2*x^2
  basis[23] = pow(x,2)*pow(y,2)*z;        // x^2*y^2*z
  basis[24] = pow(y,2)*pow(z,2)*x;        // y^2*z^2*x
  basis[25] = pow(z,2)*pow(x,2)*y;        // z^2*x^2*y
  basis[26] = pow(x,2)*pow(y,2)*pow(z,2); // x^2*y^2*z^2

  return (basis);
}

/***********************************************************************/
dbVector LagrangeQuadratic::dPx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                       // 1
  basis[1] = 1;                       // x
  basis[2] = 0;                       // y
  basis[3] = 0;                       // z
  basis[4] = y;                       // x*y
  basis[5] = 0;                       // y*z
  basis[6] = z;                       // z*x
  basis[7] = y*z;                     // x*y*z
  basis[8] = 2.0*x;                   // x^2
  basis[9] = 0;                       // y^2
  basis[10] = 0;                      // z^2
  basis[11] = 2.0*x*y;                // x^2*y
  basis[12] = 2.0*x*z;                // x^2*z
  basis[13] = pow(y,2);               // y^2*x
  basis[14] = 0;                      // y^2*z
  basis[15] = pow(z,2);               // z^2*x
  basis[16] = 0;                      // z^2*y
  basis[17] = 2.0*x*y*z;              // x^2*y*z
  basis[18] = pow(y,2)*z;             // y^2*x*z
  basis[19] = pow(z,2)*y;             // z^2*x*y
  basis[20] = 2.0*x*pow(y,2);         // x^2*y^2
  basis[21] = 0;                      // y^2*z^2
  basis[22] = 2.0*pow(z,2)*x;         // z^2*x^2
  basis[23] = 2.0*x*pow(y,2)*z;       // x^2*y^2*z
  basis[24] = pow(y,2)*pow(z,2);      // y^2*z^2*x
  basis[25] = 2.0*pow(z,2)*x*y;       // z^2*x^2*y
  basis[26] = 2.0*x*pow(y,2)*pow(z,2);// x^2*y^2*z^2

  return (basis);
}

dbVector LagrangeQuadratic:: dPy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                       // 1
  basis[1] = 0;                       // x
  basis[2] = 1;                       // y
  basis[3] = 0;                       // z
  basis[4] = x;                       // x*y
  basis[5] = z;                       // y*z
  basis[6] = 0;                       // z*x
  basis[7] = x*z;                     // x*y*z
  basis[8] = 0;                       // x^2
  basis[9] = 2.0*y;                   // y^2
  basis[10] = 0;                      // z^2
  basis[11] = pow(x,2);               // x^2*y
  basis[12] = 0;                      // x^2*z
  basis[13] = 2.0*y*x;                // y^2*x
  basis[14] = 2.0*y*z;                // y^2*z
  basis[15] = 0;                      // z^2*x
  basis[16] = pow(z,2);               // z^2*y
  basis[17] = pow(x,2)*z;             // x^2*y*z
  basis[18] = 2.0*y*x*z;              // y^2*x*z
  basis[19] = pow(z,2)*x;             // z^2*x*y
  basis[20] = 2.0*pow(x,2)*y;         // x^2*y^2
  basis[21] = 2.0*y*pow(z,2);         // y^2*z^2
  basis[22] = 0;                      // z^2*x^2
  basis[23] = 2.0*pow(x,2)*y*z;       // x^2*y^2*z
  basis[24] = 2.0*y*pow(z,2)*x;       // y^2*z^2*x
  basis[25] = pow(z,2)*pow(x,2);      // z^2*x^2*y
  basis[26] = 2.0*pow(x,2)*y*pow(z,2);// x^2*y^2*z^2

  return (basis);
}

dbVector LagrangeQuadratic:: dPz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                       // 1
  basis[1] = 0;                       // x
  basis[2] = 0;                       // y
  basis[3] = 1;                       // z
  basis[4] = 0;                       // x*y
  basis[5] = y;                       // y*z
  basis[6] = x;                       // z*x
  basis[7] = x*y;                     // x*y*z
  basis[8] = 0;                       // x^2
  basis[9] = 0;                       // y^2
  basis[10] = 2.0*z;                  // z^2
  basis[11] = 0;                      // x^2*y
  basis[12] = pow(x,2);               // x^2*z
  basis[13] = 0;                      // y^2*x
  basis[14] = pow(y,2);               // y^2*z
  basis[15] = 2.0*z*x;                // z^2*x
  basis[16] = 2.0*z*y;                // z^2*y
  basis[17] = pow(x,2)*y;             // x^2*y*z
  basis[18] = pow(y,2)*x;             // y^2*x*z
  basis[19] = 2.0*z*x*y;              // z^2*x*y
  basis[20] = 0;                      // x^2*y^2
  basis[21] = 2.0*pow(y,2)*z;         // y^2*z^2
  basis[22] = 2.0*z*pow(x,2);         // z^2*x^2
  basis[23] = pow(x,2)*pow(y,2);     // x^2*y^2*z
  basis[24] = 2.0*pow(y,2)*z*x;       // y^2*z^2*x
  basis[25] = 2.0*z*pow(x,2)*y;       // z^2*x^2*y
  basis[26] = 2.0*pow(x,2)*pow(y,2)*z;// x^2*y^2*z^2

  return (basis);
}

/***********************************************************************/
dbVector LagrangeQuadratic::dPxx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                     // 1
  basis[1] = 0;                     // x
  basis[2] = 0;                     // y
  basis[3] = 0;                     // z
  basis[4] = 0;                     // x*y
  basis[5] = 0;                     // y*z
  basis[6] = 0;                     // z*x
  basis[7] = 0;                     // x*y*z
  basis[8] = 2;                     // x^2
  basis[9] = 0;                     // y^2
  basis[10] = 0;                    // z^2
  basis[11] = 2.0*y;                // x^2*y
  basis[12] = 2.0*z;                // x^2*z
  basis[13] = 0;                    // y^2*x
  basis[14] = 0;                    // y^2*z
  basis[15] = 0;                    // z^2*x
  basis[16] = 0;                    // z^2*y
  basis[17] = 2.0*y*z;              // x^2*y*z
  basis[18] = 0;                    // y^2*x*z
  basis[19] = 0;                    // z^2*x*y
  basis[20] = 2.0*pow(y,2);         // x^2*y^2
  basis[21] = 0;                    // y^2*z^2
  basis[22] = 2.0*pow(z,2);         // z^2*x^2
  basis[23] = 2.0*pow(y,2)*z;       // x^2*y^2*z
  basis[24] = 0;                    // y^2*z^2*x
  basis[25] = 2.0*pow(z,2)*y;       // z^2*x^2*y
  basis[26] = 2.0*pow(y,2)*pow(z,2);// x^2*y^2*z^2

  return (basis);
}

dbVector LagrangeQuadratic:: dPyy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                     // 1
  basis[1] = 0;                     // x
  basis[2] = 0;                     // y
  basis[3] = 0;                     // z
  basis[4] = 0;                     // x*y
  basis[5] = 0;                     // y*z
  basis[6] = 0;                     // z*x
  basis[7] = 0;                     // x*y*z
  basis[8] = 0;                      // x^2
  basis[9] = 2;                     // y^2
  basis[10] = 0;                    // z^2
  basis[11] = 0;                    // x^2*y
  basis[12] = 0;                    // x^2*z
  basis[13] = 2.0*x;                // y^2*x
  basis[14] = 2.0*z;                // y^2*z
  basis[15] = 0;                    // z^2*x
  basis[16] = 0;                    // z^2*y
  basis[17] = 0;                    // x^2*y*z
  basis[18] = 2.0*x*z;              // y^2*x*z
  basis[19] = 0;                    // z^2*x*y
  basis[20] = 2.0*pow(x,2);         // x^2*y^2
  basis[21] = 2.0*pow(z,2);         // y^2*z^2
  basis[22] = 0;                    // z^2*x^2
  basis[23] = 2.0*pow(x,2)*z;       // x^2*y^2*z
  basis[24] = 2.0*pow(z,2)*x;       // y^2*z^2*x
  basis[25] = 0;                    // z^2*x^2*y
  basis[26] = 2.0*pow(x,2)*pow(z,2);// x^2*y^2*z^2

  return (basis);
}

dbVector LagrangeQuadratic:: dPzz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                     // 1
  basis[1] = 0;                     // x
  basis[2] = 0;                     // y
  basis[3] = 0;                     // z
  basis[4] = 0;                     // x*y
  basis[5] = 0;                     // y*z
  basis[6] = 0;                     // z*x
  basis[7] = 0;                     // x*y*z
  basis[8] = 0;                     // x^2
  basis[9] = 0;                     // y^2
  basis[10] = 2;                    // z^2
  basis[11] = 0;                    // x^2*y
  basis[12] = 0;                    // x^2*z
  basis[13] = 0;                    // y^2*x
  basis[14] = 0;                    // y^2*z
  basis[15] = 2.0*x;                // z^2*x
  basis[16] = 2.0*y;                // z^2*y
  basis[17] = 0;                    // x^2*y*z
  basis[18] = 0;                    // y^2*x*z
  basis[19] = 2.0*x*y;              // z^2*x*y
  basis[20] = 0;                    // x^2*y^2
  basis[21] = 2.0*pow(y,2);         // y^2*z^2
  basis[22] = 2.0*pow(x,2);         // z^2*x^2
  basis[23] = 0;                    // x^2*y^2*z
  basis[24] = 2.0*pow(y,2)*x;       // y^2*z^2*x
  basis[25] = 2.0*pow(x,2)*y;       // z^2*x^2*y
  basis[26] = 2.0*pow(x,2)*pow(y,2);// x^2*y^2*z^2

  return (basis);
}

dbVector LagrangeQuadratic::dPxy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                // 1
  basis[1] = 0;                // x
  basis[2] = 0;                // y
  basis[3] = 0;                // z
  basis[4] = 1;                // x*y
  basis[5] = 0;                // y*z
  basis[6] = 0;                // z*x
  basis[7] = z;                // x*y*z
  basis[8] = 0;                // x^2
  basis[9] = 0;                // y^2
  basis[10] = 0;               // z^2
  basis[11] = 2.0*x;           // x^2*y
  basis[12] = 0;               // x^2*z
  basis[13] = 2.0*y;           // y^2*x
  basis[14] = 0;               // y^2*z
  basis[15] = 0;               // z^2*x
  basis[16] = 0;               // z^2*y
  basis[17] = 2.0*x*z;         // x^2*y*z
  basis[18] = 2.0*y*z;         // y^2*x*z
  basis[19] = pow(z,2);        // z^2*x*y
  basis[20] = 4.0*x*y;         // x^2*y^2
  basis[21] = 0;               // y^2*z^2
  basis[22] = 0;               // z^2*x^2
  basis[23] = 4.0*x*y*z;       // x^2*y^2*z
  basis[24] = 2.0*y*pow(z,2);  // y^2*z^2*x
  basis[25] = 2.0*pow(z,2)*x;  // z^2*x^2*y
  basis[26] = 4.0*x*y*pow(z,2);// x^2*y^2*z^2

  return (basis);
}

dbVector LagrangeQuadratic:: dPyz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                // 1
  basis[1] = 0;                // x
  basis[2] = 0;                // y
  basis[3] = 0;                // z
  basis[4] = 0;                // x*y
  basis[5] = 1;                // y*z
  basis[6] = 0;                // z*x
  basis[7] = x;                // x*y*z
  basis[8] = 0;                // x^2
  basis[9] = 0;                // y^2
  basis[10] = 0;               // z^2
  basis[11] = 0;               // x^2*y
  basis[12] = 0;               // x^2*z
  basis[13] = 0;               // y^2*x
  basis[14] = 2.0*y;           // y^2*z
  basis[15] = 0;               // z^2*x
  basis[16] = 2.0*z;           // z^2*y
  basis[17] = pow(x,2);        // x^2*y*z
  basis[18] = 2.0*y*x;         // y^2*x*z
  basis[19] = 2.0*z*x;         // z^2*x*y
  basis[20] = 0;               // x^2*y^2
  basis[21] = 4.0*y*z;         // y^2*z^2
  basis[22] = 0;               // z^2*x^2
  basis[23] = 2.0*pow(x,2)*y;  // x^2*y^2*z
  basis[24] = 4.0*y*z*x;       // y^2*z^2*x
  basis[25] = 2.0*z*pow(x,2);  // z^2*x^2*y
  basis[26] = 4.0*pow(x,2)*y*z;// x^2*y^2*z^2

  return (basis);
}

dbVector LagrangeQuadratic:: dPzx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(27);

  basis[0] = 0;                // 1
  basis[1] = 0;                // x
  basis[2] = 0;                // y
  basis[3] = 0;                // z
  basis[4] = 0;                // x*y
  basis[5] = 0;                // y*z
  basis[6] = 1;                // z*x
  basis[7] = y;                // x*y*z
  basis[8] = 0;                // x^2
  basis[9] = 0;                // y^2
  basis[10] = 0;               // z^2
  basis[11] = 0;               // x^2*y
  basis[12] = 2.0*x;           // x^2*z
  basis[13] = 0;               // y^2*x
  basis[14] = 0;               // y^2*z
  basis[15] = 2.0*z;           // z^2*x
  basis[16] = 0;               // z^2*y
  basis[17] = 2.0*x*y;         // x^2*y*z
  basis[18] = pow(y,2);        // y^2*x*z
  basis[19] = 2.0*z*y;         // z^2*x*y
  basis[20] = 0;               // x^2*y^2
  basis[21] = 0;               // y^2*z^2
  basis[22] = 4.0*z*x;         // z^2*x^2
  basis[23] = 2.0*x*pow(y,2);  // x^2*y^2*z
  basis[24] = 2.0*pow(y,2)*z;  // y^2*z^2*x
  basis[25] = 4.0*z*x*y;       // z^2*x^2*y
  basis[26] = 4.0*x*pow(y,2)*z;// x^2*y^2*z^2

  return (basis);
}


/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

// Polynom-Type 3 (Serendipity)

dbVector SerendipityQuadratic::P(double& x,double& y,double& z){

  using namespace std;

  dbVector basis(20);

  basis[0] = 1;                           // 1
  basis[1] = x;                           // x
  basis[2] = y;                           // y
  basis[3] = z;                           // z
  basis[4] = x*y;                         // x*y
  basis[5] = y*z;                         // y*z
  basis[6] = z*x;                         // z*x
  basis[7] = x*y*z;                       // x*y*z
  basis[8] = pow(x,2);                    // x^2
  basis[9] = pow(y,2);                    // y^2
  basis[10] = pow(z,2);                   // z^2
  basis[11] = pow(x,2)*y;                 // x^2*y
  basis[12] = pow(x,2)*z;                 // x^2*z
  basis[13] = pow(y,2)*x;                 // y^2*x
  basis[14] = pow(y,2)*z;                 // y^2*z
  basis[15] = pow(z,2)*x;                 // z^2*x
  basis[16] = pow(z,2)*y;                 // z^2*y
  basis[17] = pow(x,2)*y*z;               // x^2*y*z
  basis[18] = pow(y,2)*x*z;               // y^2*x*z
  basis[19] = pow(z,2)*x*y;               // z^2*x*y

  return (basis);
}

/***********************************************************************/
dbVector SerendipityQuadratic::dPx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                       // 1
  basis[1] = 1;                       // x
  basis[2] = 0;                       // y
  basis[3] = 0;                       // z
  basis[4] = y;                       // x*y
  basis[5] = 0;                       // y*z
  basis[6] = z;                       // z*x
  basis[7] = y*z;                     // x*y*z
  basis[8] = 2.0*x;                   // x^2
  basis[9] = 0;                       // y^2
  basis[10] = 0;                      // z^2
  basis[11] = 2.0*x*y;                // x^2*y
  basis[12] = 2.0*x*z;                // x^2*z
  basis[13] = pow(y,2);               // y^2*x
  basis[14] = 0;                      // y^2*z
  basis[15] = pow(z,2);               // z^2*x
  basis[16] = 0;                      // z^2*y
  basis[17] = 2.0*x*y*z;              // x^2*y*z
  basis[18] = pow(y,2)*z;             // y^2*x*z
  basis[19] = pow(z,2)*y;             // z^2*x*y

  return (basis);
}

dbVector SerendipityQuadratic:: dPy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                       // 1
  basis[1] = 0;                       // x
  basis[2] = 1;                       // y
  basis[3] = 0;                       // z
  basis[4] = x;                       // x*y
  basis[5] = z;                       // y*z
  basis[6] = 0;                       // z*x
  basis[7] = x*z;                     // x*y*z
  basis[8] = 0;                       // x^2
  basis[9] = 2.0*y;                   // y^2
  basis[10] = 0;                      // z^2
  basis[11] = pow(x,2);               // x^2*y
  basis[12] = 0;                      // x^2*z
  basis[13] = 2.0*y*x;                // y^2*x
  basis[14] = 2.0*y*z;                // y^2*z
  basis[15] = 0;                      // z^2*x
  basis[16] = pow(z,2);               // z^2*y
  basis[17] = pow(x,2)*z;             // x^2*y*z
  basis[18] = 2.0*y*x*z;              // y^2*x*z
  basis[19] = pow(z,2)*x;             // z^2*x*y

  return (basis);
}

dbVector SerendipityQuadratic:: dPz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                       // 1
  basis[1] = 0;                       // x
  basis[2] = 0;                       // y
  basis[3] = 1;                       // z
  basis[4] = 0;                       // x*y
  basis[5] = y;                       // y*z
  basis[6] = x;                       // z*x
  basis[7] = x*y;                     // x*y*z
  basis[8] = 0;                       // x^2
  basis[9] = 0;                       // y^2
  basis[10] = 2.0*z;                  // z^2
  basis[11] = 0;                      // x^2*y
  basis[12] = pow(x,2);               // x^2*z
  basis[13] = 0;                      // y^2*x
  basis[14] = pow(y,2);               // y^2*z
  basis[15] = 2.0*z*x;                // z^2*x
  basis[16] = 2.0*z*y;                // z^2*y
  basis[17] = pow(x,2)*y;             // x^2*y*z
  basis[18] = pow(y,2)*x;             // y^2*x*z
  basis[19] = 2.0*z*x*y;              // z^2*x*y

  return (basis);
}

/***********************************************************************/
dbVector SerendipityQuadratic::dPxx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                     // 1
  basis[1] = 0;                     // x
  basis[2] = 0;                     // y
  basis[3] = 0;                     // z
  basis[4] = 0;                     // x*y
  basis[5] = 0;                     // y*z
  basis[6] = 0;                     // z*x
  basis[7] = 0;                     // x*y*z
  basis[8] = 2;                     // x^2
  basis[9] = 0;                     // y^2
  basis[10] = 0;                    // z^2
  basis[11] = 2.0*y;                // x^2*y
  basis[12] = 2.0*z;                // x^2*z
  basis[13] = 0;                    // y^2*x
  basis[14] = 0;                    // y^2*z
  basis[15] = 0;                    // z^2*x
  basis[16] = 0;                    // z^2*y
  basis[17] = 2.0*y*z;              // x^2*y*z
  basis[18] = 0;                    // y^2*x*z
  basis[19] = 0;                    // z^2*x*y

  return (basis);
}

dbVector SerendipityQuadratic:: dPyy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                     // 1
  basis[1] = 0;                     // x
  basis[2] = 0;                     // y
  basis[3] = 0;                     // z
  basis[4] = 0;                     // x*y
  basis[5] = 0;                     // y*z
  basis[6] = 0;                     // z*x
  basis[7] = 0;                     // x*y*z
  basis[8] = 0;                      // x^2
  basis[9] = 2;                     // y^2
  basis[10] = 0;                    // z^2
  basis[11] = 0;                    // x^2*y
  basis[12] = 0;                    // x^2*z
  basis[13] = 2.0*x;                // y^2*x
  basis[14] = 2.0*z;                // y^2*z
  basis[15] = 0;                    // z^2*x
  basis[16] = 0;                    // z^2*y
  basis[17] = 0;                    // x^2*y*z
  basis[18] = 2.0*x*z;              // y^2*x*z
  basis[19] = 0;                    // z^2*x*y

  return (basis);
}

dbVector SerendipityQuadratic:: dPzz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                     // 1
  basis[1] = 0;                     // x
  basis[2] = 0;                     // y
  basis[3] = 0;                     // z
  basis[4] = 0;                     // x*y
  basis[5] = 0;                     // y*z
  basis[6] = 0;                     // z*x
  basis[7] = 0;                     // x*y*z
  basis[8] = 0;                     // x^2
  basis[9] = 0;                     // y^2
  basis[10] = 2;                    // z^2
  basis[11] = 0;                    // x^2*y
  basis[12] = 0;                    // x^2*z
  basis[13] = 0;                    // y^2*x
  basis[14] = 0;                    // y^2*z
  basis[15] = 2.0*x;                // z^2*x
  basis[16] = 2.0*y;                // z^2*y
  basis[17] = 0;                    // x^2*y*z
  basis[18] = 0;                    // y^2*x*z
  basis[19] = 2.0*x*y;              // z^2*x*y

  return (basis);
}

dbVector SerendipityQuadratic::dPxy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                // 1
  basis[1] = 0;                // x
  basis[2] = 0;                // y
  basis[3] = 0;                // z
  basis[4] = 1;                // x*y
  basis[5] = 0;                // y*z
  basis[6] = 0;                // z*x
  basis[7] = z;                // x*y*z
  basis[8] = 0;                // x^2
  basis[9] = 0;                // y^2
  basis[10] = 0;               // z^2
  basis[11] = 2.0*x;           // x^2*y
  basis[12] = 0;               // x^2*z
  basis[13] = 2.0*y;           // y^2*x
  basis[14] = 0;               // y^2*z
  basis[15] = 0;               // z^2*x
  basis[16] = 0;               // z^2*y
  basis[17] = 2.0*x*z;         // x^2*y*z
  basis[18] = 2.0*y*z;         // y^2*x*z
  basis[19] = pow(z,2);        // z^2*x*y

  return (basis);
}

dbVector SerendipityQuadratic:: dPyz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                // 1
  basis[1] = 0;                // x
  basis[2] = 0;                // y
  basis[3] = 0;                // z
  basis[4] = 0;                // x*y
  basis[5] = 1;                // y*z
  basis[6] = 0;                // z*x
  basis[7] = x;                // x*y*z
  basis[8] = 0;                // x^2
  basis[9] = 0;                // y^2
  basis[10] = 0;               // z^2
  basis[11] = 0;               // x^2*y
  basis[12] = 0;               // x^2*z
  basis[13] = 0;               // y^2*x
  basis[14] = 2.0*y;           // y^2*z
  basis[15] = 0;               // z^2*x
  basis[16] = 2.0*z;           // z^2*y
  basis[17] = pow(x,2);        // x^2*y*z
  basis[18] = 2.0*y*x;         // y^2*x*z
  basis[19] = 2.0*z*x;         // z^2*x*y

  return (basis);
}

dbVector SerendipityQuadratic:: dPzx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(20);

  basis[0] = 0;                // 1
  basis[1] = 0;                // x
  basis[2] = 0;                // y
  basis[3] = 0;                // z
  basis[4] = 0;                // x*y
  basis[5] = 0;                // y*z
  basis[6] = 1;                // z*x
  basis[7] = y;                // x*y*z
  basis[8] = 0;                // x^2
  basis[9] = 0;                // y^2
  basis[10] = 0;               // z^2
  basis[11] = 0;               // x^2*y
  basis[12] = 2.0*x;           // x^2*z
  basis[13] = 0;               // y^2*x
  basis[14] = 0;               // y^2*z
  basis[15] = 2.0*z;           // z^2*x
  basis[16] = 0;               // z^2*y
  basis[17] = 2.0*x*y;         // x^2*y*z
  basis[18] = pow(y,2);        // y^2*x*z
  basis[19] = 2.0*z*y;         // z^2*x*y

  return (basis);
}


/***********************************************************************/
/***********************************************************************/
dbVector BernsteinLinear::P(double& x,double& y,double& z){

  using namespace std;

  dbVector basis(8);

  basis[0] = x*y*z;
  basis[1] = x*y*(1.0-z);
  basis[2] = x*(1.0-y)*z;
  basis[3] = x*(1.0-y)*(1.0-z);
  basis[4] = (1.0-x)*y*z;
  basis[5] = (1.0-x)*y*(1.0-z);
  basis[6] = (1.0-x)*(1.0-y)*z;
  basis[7] = (1.0-x)*(1.0-y)*(1.0-z);

  return (basis);
}

dbVector BernsteinLinear::dPx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = y*z;
  basis[1] = y*(1.0-z);
  basis[2] = (1.0-y)*z;
  basis[3] = (1.0-y)*(1.0-z);
  basis[4] = (-1.0)*y*z;
  basis[5] = (-1.0)*y*(1.0-z);
  basis[6] = (-1.0)*(1.0-y)*z;
  basis[7] = (-1.0)*(1.0-y)*(1.0-z);

  return (basis);
}

dbVector BernsteinLinear:: dPy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = x*z;
  basis[1] = x*(1.0-z);
  basis[2] = (-1.0)*x*z;
  basis[3] = (-1.0)*x*(1.0-z);
  basis[4] = (1.0-x)*z;
  basis[5] = (1.0-x)*(1.0-z);
  basis[6] = (-1.0)*(1.0-x)*z;
  basis[7] = (-1.0)*(1.0-x)*(1.0-z);

  return (basis);
}

dbVector BernsteinLinear:: dPz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = x*y;
  basis[1] = (-1.0)*x*y;
  basis[2] = x*(1.0-y);
  basis[3] = (-1.0)*x*(1.0-y);
  basis[4] = (1.0-x)*y;
  basis[5] = (-1.0)*(1.0-x)*y;
  basis[6] = (1.0-x)*(1.0-y);
  basis[7] = (-1.0)*(1.0-x)*(1.0-y);

  return (basis);
}

dbVector BernsteinLinear:: dPxy(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = z;
  basis[1] = 1.0-z;
  basis[2] = (-1.0)*z;
  basis[3] = (-1.0)*(1.0-z);
  basis[4] = (-1.0)*z;
  basis[5] = (-1.0)*(1.0-z);
  basis[6] = z;
  basis[7] = 1.0-z;

  return (basis);
}

dbVector BernsteinLinear:: dPyz(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = x;
  basis[1] = (-1.0)*x;
  basis[2] = (-1.0)*x;
  basis[3] = x;
  basis[4] = 1.0-x;
  basis[5] = (-1.0)*(1.0-x);
  basis[6] = (-1.0)*(1.0-x);
  basis[7] = 1.0-x;

  return (basis);
}

dbVector BernsteinLinear:: dPzx(double& x,double& y,double& z) {

  using namespace std;

  dbVector basis(8);

  basis[0] = y;
  basis[1] = (-1.0)*y;
  basis[2] = 1.0-y;
  basis[3] = (-1.0)*(1.0-y);
  basis[4] = (-1.0)*y;
  basis[5] = y;
  basis[6] = (-1.0)*(1.0-y);
  basis[7] = 1.0-y;

  return (basis);
}
