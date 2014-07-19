/*
 * ParticleExt.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef PARTICLEEXT_H_
#define PARTICLEEXT_H_

#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "InputFileData.h"
#include "IntegrationPointX.h"
#include "mpi.h"
#include "Particle.h"


class ParticleExt: public Particle {

  public:

    ParticleExt(int usedDOF);
    ParticleExt(Particle ptcle);

    ~ParticleExt();

    void setMaterialID(int matID){this->getMaterialID() =matID;};

    void setStepDOFMat(dbMatrix mat){stepDegreesOfFreedomMat = mat;};
    dbMatrix& getStepDOFMat(){return stepDegreesOfFreedomMat;};

    void dummyFunction();


    private:

    dbVector degreesOfFreedom;
    dbMatrix stepDegreesOfFreedomMat;
    dbVector stepDegreesOfFreedom; // DOF obtained during one time step
    dbVector deltaDOF;

};
#endif /* PARTICLEEXT_H_ */
