// Stores all properties of a single particle.

#include "ParticleExt.h"

ParticleExt::ParticleExt(int usedDOF) : Particle(usedDOF){

  stepDegreesOfFreedomMat = dbMatrix(0);
  degreesOfFreedom = dbVector(usedDOF);
}

ParticleExt::ParticleExt(Particle ptcle) : Particle(ptcle){

  stepDegreesOfFreedomMat = dbMatrix(0);
  degreesOfFreedom = dbVector(this->getDOF().size());
}



ParticleExt::~ParticleExt() {}

void ParticleExt::dummyFunction(){
	this->getCoords();
}
