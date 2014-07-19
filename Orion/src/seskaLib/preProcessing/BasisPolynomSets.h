/* Return a vector containing the polynom basis for a shapefuntions set.*/

#ifndef BasisPolynomSets_h_
#define BasisPolynomSets_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "mpi.h"
#include "petscsys.h"

// Polynom-Type 1 (Pascal)

namespace PascalLinear {

  dbVector P(double& x,double& y,double& z);
  dbVector dPx(double& x,double& y,double& z);
  dbVector dPy(double& x,double& y,double& z);
  dbVector dPz(double& x,double& y,double& z);
}

namespace PascalQuadratic {

  dbVector P(double& x,double& y,double& z);

  dbVector dPx(double& x,double& y,double& z);
  dbVector dPy(double& x,double& y,double& z);
  dbVector dPz(double& x,double& y,double& z);

  dbVector dPxx(double& x,double& y,double& z);
  dbVector dPyy(double& x,double& y,double& z);
  dbVector dPzz(double& x,double& y,double& z);

  dbVector dPxy(double& x,double& y,double& z);
  dbVector dPyz(double& x,double& y,double& z);
  dbVector dPzx(double& x,double& y,double& z);
}

namespace PascalCubic {

  dbVector P(double& x,double& y,double& z);

  dbVector dPx(double& x,double& y,double& z);
  dbVector dPy(double& x,double& y,double& z);
  dbVector dPz(double& x,double& y,double& z);

  dbVector dPxx(double& x,double& y,double& z);
  dbVector dPyy(double& x,double& y,double& z);
  dbVector dPzz(double& x,double& y,double& z);

  dbVector dPxy(double& x,double& y,double& z);
  dbVector dPyz(double& x,double& y,double& z);
  dbVector dPzx(double& x,double& y,double& z);
}

// Polynom-Type 2 (Lagrangian).

namespace LagrangeLinear {

  dbVector P(double& x,double& y,double& z);

  dbVector dPx(double& x,double& y,double& z);
  dbVector dPy(double& x,double& y,double& z);
  dbVector dPz(double& x,double& y,double& z);

  dbVector dPxy(double& x,double& y,double& z);
  dbVector dPyz(double& x,double& y,double& z);
  dbVector dPzx(double& x,double& y,double& z);
}

namespace LagrangeQuadratic {

  dbVector P(double& x,double& y,double& z);

  dbVector dPx(double& x,double& y,double& z);
  dbVector dPy(double& x,double& y,double& z);
  dbVector dPz(double& x,double& y,double& z);

  dbVector dPxx(double& x,double& y,double& z);
  dbVector dPyy(double& x,double& y,double& z);
  dbVector dPzz(double& x,double& y,double& z);

  dbVector dPxy(double& x,double& y,double& z);
  dbVector dPyz(double& x,double& y,double& z);
  dbVector dPzx(double& x,double& y,double& z);
}

// Polynom-Type 3 (Serendipity).

namespace SerendipityQuadratic {

  dbVector P(double& x,double& y,double& z);

  dbVector dPx(double& x,double& y,double& z);
  dbVector dPy(double& x,double& y,double& z);
  dbVector dPz(double& x,double& y,double& z);

  dbVector dPxx(double& x,double& y,double& z);
  dbVector dPyy(double& x,double& y,double& z);
  dbVector dPzz(double& x,double& y,double& z);

  dbVector dPxy(double& x,double& y,double& z);
  dbVector dPyz(double& x,double& y,double& z);
  dbVector dPzx(double& x,double& y,double& z);
}

// Polynom-Type 4 (Bernstein).

namespace BernsteinLinear {

  dbVector P(double& x,double& y,double& z);

  dbVector dPx(double& x,double& y,double& z);
  dbVector dPy(double& x,double& y,double& z);
  dbVector dPz(double& x,double& y,double& z);

  dbVector dPxy(double& x,double& y,double& z);
  dbVector dPyz(double& x,double& y,double& z);
  dbVector dPzx(double& x,double& y,double& z);
}

#endif
