// Calculate the basis polynom set for a certain gauss point - one for each
// particle.

#ifndef BasisPolynom_h_
#define BasisPolynom_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"

class BasisPolynom {

  public:

    BasisPolynom() {};
    ~BasisPolynom() {};

    dbMatrix& getBasis() { return basis; };

    dbMatrix& getXDerivBasis() { return xDerivBasis; };
    dbMatrix& getYDerivBasis() { return yDerivBasis; };
    dbMatrix& getZDerivBasis() { return zDerivBasis; };

    dbMatrix& getXXDerivBasis() { return xxDerivBasis; };
    dbMatrix& getYYDerivBasis() { return yyDerivBasis; };
    dbMatrix& getZZDerivBasis() { return zzDerivBasis; };

    dbMatrix& getXYDerivBasis() { return xyDerivBasis; };
    dbMatrix& getYZDerivBasis() { return yzDerivBasis; };
    dbMatrix& getZXDerivBasis() { return zxDerivBasis; };

  protected:
    dbMatrix basis;

    dbMatrix xDerivBasis;
    dbMatrix yDerivBasis;
    dbMatrix zDerivBasis;

    dbMatrix xxDerivBasis;
    dbMatrix yyDerivBasis;
    dbMatrix zzDerivBasis;

    dbMatrix xyDerivBasis;
    dbMatrix yzDerivBasis;
    dbMatrix zxDerivBasis;

 private:

    /* Calculate basis polynom set.*/
    virtual void calcPolynomSet() {};

    /* Calculation of first order derivations of the basis polynom */
    virtual void calcPoly1stDerivs() {};

    /* Calculation of second order derivations of the basis polynom */
    virtual void calcPoly2ndDerivs() {};

};

#endif
