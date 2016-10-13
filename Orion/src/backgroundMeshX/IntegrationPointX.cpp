/*
 * IntegrationPointX.cpp
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

// Stores all properties of a single integration point.
#include "IntegrationPointX.h"


IntegrationPointX::IntegrationPointX() : globalID(0),materialID(0),lc0(0),
				       depolarisationTime(0),
                       activeTension(0),
				       activeTensionVariation(0),
				       intermediateJacobian(0) {

  coords = dbVector(3);

  firstDerivShapeFuncs = dbMatrix(3);
  secondDerivShapeFuncs = dbMatrix(6);

}

void IntegrationPointX::setCoords(double coord1,double coord2,double coord3) {

  coords[0] = coord1;
  coords[1] = coord2;
  coords[2] = coord3;
}

void IntegrationPointX::setWeight(int ID,double value) {

  if(ID >= weights.size())
    weights.resize(ID+1);

  weights[ID] = value;
}

double& IntegrationPointX::getIntWeight(int ID) {

  if(ID < weights.size())
    return weights[ID];

  else {
    weights.resize(ID+1);
    return weights[ID];
  }

}

void IntegrationPointX::setSupportPtcle(int particle) {

  supportingPtcls.push_back(particle);
}

void IntegrationPointX::setLocalMinPtcleDistance(double localMinPtcleDist) {

    localMinPtcleDistance = localMinPtcleDist;
}

void IntegrationPointX::checkSuppPtclsSize(int size) {

  if(supportingPtcls.size() < (unsigned int)size ||
     supportingPtcls.size() > (unsigned int)size)
    supportingPtcls.resize(size);
}

void IntegrationPointX::setGlobalID(int idx) {

  globalID = idx;
}

void IntegrationPointX::setShapeFuncsSize(int& size) {

  shapeFuncs = dbVector(size);
}

void IntegrationPointX::setFirstDerivShapesSize(int& size) {

  for(int i=0;i<3;i++)
    firstDerivShapeFuncs[i].resize(size);

}

void IntegrationPointX::setSecondDerivShapesSize(int& size) {

  for(int i=0;i<6;i++)
    secondDerivShapeFuncs[i].resize(size);

}

void IntegrationPointX::setShapeFuncs(int& size,dbVector& sFuncs) {

  if((unsigned int)size > shapeFuncs.size())
    shapeFuncs.resize(size);

  for(int i=0;i<size;i++)
    shapeFuncs[i] = sFuncs[i];
}

void IntegrationPointX::setFirstDerivShapes(int& size,
				     dbMatrix& firstDerivShapes) {

  if((unsigned int)size > firstDerivShapeFuncs[0].size()) {
    firstDerivShapeFuncs[0].resize(size);
    firstDerivShapeFuncs[1].resize(size);
    firstDerivShapeFuncs[2].resize(size);
  }

  for(int i=0;i<size;i++) {
    firstDerivShapeFuncs[0][i] = firstDerivShapes[0][i];
    firstDerivShapeFuncs[1][i] = firstDerivShapes[1][i];
    firstDerivShapeFuncs[2][i] = firstDerivShapes[2][i];
  }
}

void IntegrationPointX::setSecondDerivShapes(int& size,
				      dbMatrix& secondDerivShapes) {

  if((unsigned int)size > secondDerivShapeFuncs[0].size()) {
    secondDerivShapeFuncs[0].resize(size);
    secondDerivShapeFuncs[1].resize(size);
    secondDerivShapeFuncs[2].resize(size);
    secondDerivShapeFuncs[3].resize(size);
    secondDerivShapeFuncs[4].resize(size);
    secondDerivShapeFuncs[5].resize(size);
  }

  for(int i=0;i<size;i++) {
    secondDerivShapeFuncs[0][i] = secondDerivShapes[0][i];
    secondDerivShapeFuncs[1][i] = secondDerivShapes[1][i];
    secondDerivShapeFuncs[2][i] = secondDerivShapes[2][i];

    secondDerivShapeFuncs[3][i] = secondDerivShapes[3][i];
    secondDerivShapeFuncs[4][i] = secondDerivShapes[4][i];
    secondDerivShapeFuncs[5][i] = secondDerivShapes[5][i];
  }
}

void IntegrationPointX::setGlobalDOF(int idx,int value) {

  globalDOF[idx] = value;
}

void IntegrationPointX::setGlobalDOFSize(int dof) {

  using namespace std;

  globalDOF = intVector(supportingPtcls.size()*dof);
}

void IntegrationPointX::setSurfaceNormal(dbVector& normal) {

  using namespace std;

  surfaceNormal = normal;
}

/***********************************************************************/
// Neumann boundary conditions
double& IntegrationPointX::getSurfacePressure(int ID) {

  if(ID < surfacePressureLoads.size())
    return surfacePressureLoads[ID];

  else {
    surfacePressureLoads.resize(ID+1);
    return surfacePressureLoads[ID];
  }

}

blVector& IntegrationPointX::getBodyForceDOF(int ID) {

  if(ID < bodyForceDOF.size())
    return bodyForceDOF[ID];

  else {
    bodyForceDOF.resize(ID+1);
    return bodyForceDOF[ID];
  }

}

dbVector& IntegrationPointX::getBodyForce(int ID) {

  if(ID < bodyForceLoads.size())
    return bodyForceLoads[ID];

  else {
    bodyForceLoads.resize(ID+1);
   return bodyForceLoads[ID];
  }

}

blVector& IntegrationPointX::getTractionDOF(int ID) {

  if(ID < tractionDOF.size())
    return tractionDOF[ID];

  else {
    tractionDOF.resize(ID+1);
    return tractionDOF[ID];
  }

}

dbVector& IntegrationPointX::getTraction(int ID) {

  if(ID < tractionLoads.size())
    return tractionLoads[ID];

  else {
    tractionLoads.resize(ID+1);
    return tractionLoads[ID];
  }

}


blVector& IntegrationPointX::getLineForceDOF(int ID) {

  if(ID < lineForceDOF.size())
    return lineForceDOF[ID];

  else {
    lineForceDOF.resize(ID+1);
    return lineForceDOF[ID];
  }

}

dbVector& IntegrationPointX::getLineForce(int ID) {

  if(ID < lineForceLoads.size())
    return lineForceLoads[ID];

  else {
    lineForceLoads.resize(ID+1);
    return lineForceLoads[ID];
  }

}

blVector& IntegrationPointX::getBodyMomentDOF(int ID) {

  if(ID < bodyMomentDOF.size())
    return bodyMomentDOF[ID];

  else {
    bodyMomentDOF.resize(ID+1);
    return bodyMomentDOF[ID];
  }

}

dbVector& IntegrationPointX::getBodyMoment(int ID) {

  if(ID < bodyMomentLoads.size())
    return bodyMomentLoads[ID];

  else {
    bodyMomentLoads.resize(ID+1);
   return bodyMomentLoads[ID];
  }

}

blVector& IntegrationPointX::getSurfaceMomentDOF(int ID) {

  if(ID < surfaceMomentDOF.size())
    return surfaceMomentDOF[ID];

  else {
    surfaceMomentDOF.resize(ID+1);
    return surfaceMomentDOF[ID];
  }

}

dbVector& IntegrationPointX::getSurfaceMoment(int ID) {

  if(ID < surfaceMomentLoads.size())
    return surfaceMomentLoads[ID];

  else {
    surfaceMomentLoads.resize(ID+1);
    return surfaceMomentLoads[ID];
  }

}

blVector& IntegrationPointX::getSurfaceElectricChargeDOF(int ID) {

  if(ID < surfaceElectricChargeDOF.size())
    return surfaceElectricChargeDOF[ID];

  else {
    surfaceElectricChargeDOF.resize(ID+1);
    return surfaceElectricChargeDOF[ID];
  }

}

dbVector& IntegrationPointX::getSurfaceElectricCharge(int ID) {

  if(ID < surfaceElectricChargeLoads.size())
    return surfaceElectricChargeLoads[ID];

  else {
    surfaceElectricChargeLoads.resize(ID+1);
    return surfaceElectricChargeLoads[ID];
  }

}

blVector& IntegrationPointX::getBodyElectricChargeDOF(int ID) {

  if(ID < bodyElectricChargeDOF.size())
    return bodyElectricChargeDOF[ID];

  else {
    bodyElectricChargeDOF.resize(ID+1);
    return bodyElectricChargeDOF[ID];
  }

}

dbVector& IntegrationPointX::getBodyElectricCharge(int ID) {

  if(ID < bodyElectricChargeLoads.size())
    return bodyElectricChargeLoads[ID];

  else {
    bodyElectricChargeLoads.resize(ID+1);
    return bodyElectricChargeLoads[ID];
  }

}

/***********************************************************************/
// Dirichlet boundary conditions

// deformation boundary conditions
blVector& IntegrationPointX::getDeformationBoundDOF(int ID) {

  if(ID < deformationBoundDOF.size())
    return deformationBoundDOF[ID];

  else {
    deformationBoundDOF.resize(ID+1);
    return deformationBoundDOF[ID];
  }

}

dbVector& IntegrationPointX::getDeformationBoundConds(int ID) {

  if(ID < deformationBoundConds.size())
    return deformationBoundConds[ID];

  else {
    deformationBoundConds.resize(ID+1);
   return deformationBoundConds[ID];
  }

}

dbVector& IntegrationPointX::getDeltaDeformationBoundConds(int ID) {

  if(ID < deltaDeformationBoundConds.size())
    return deltaDeformationBoundConds[ID];

  else {
    deltaDeformationBoundConds.resize(ID+1);
    return deltaDeformationBoundConds[ID];
  }

}

dbVector& IntegrationPointX::getInitialDeformationBoundConds(int ID) {

  if(ID < initialDeformationBoundConds.size())
    return initialDeformationBoundConds[ID];

  else {
    initialDeformationBoundConds.resize(ID+1);
    return initialDeformationBoundConds[ID];
  }

}

// ----------------------------------------------------------------------
// electric boundary conditions

blVector& IntegrationPointX::getElectricBoundDOF(int ID) {

  if(ID < electricBoundDOF.size())
    return electricBoundDOF[ID];

  else {
    electricBoundDOF.resize(ID+1);
    return electricBoundDOF[ID];
  }

}

dbVector& IntegrationPointX::getElectricBoundConds(int ID) {

  if(ID < electricBoundConds.size())
    return electricBoundConds[ID];

  else {
    electricBoundConds.resize(ID+1);
   return electricBoundConds[ID];
  }

}

dbVector& IntegrationPointX::getDeltaElectricBoundConds(int ID) {

  if(ID < deltaElectricBoundConds.size())
    return deltaElectricBoundConds[ID];

  else {
    deltaElectricBoundConds.resize(ID+1);
    return deltaElectricBoundConds[ID];
  }

}

dbVector& IntegrationPointX::getInitialElectricBoundConds(int ID) {

  if(ID < initialElectricBoundConds.size())
    return initialElectricBoundConds[ID];

  else {
    initialElectricBoundConds.resize(ID+1);
    return initialElectricBoundConds[ID];
  }

}


// ----------------------------------------------------------------------
// depolarisation boundary conditions

blVector& IntegrationPointX::getDepolarisationBoundDOF(int ID) {

  if(ID < depolarisationBoundDOF.size())
    return depolarisationBoundDOF[ID];

  else {
    depolarisationBoundDOF.resize(ID+1);
    return depolarisationBoundDOF[ID];
  }

}

dbVector& IntegrationPointX::getDepolarisationBoundConds(int ID) {

  if(ID < depolarisationBoundConds.size())
    return depolarisationBoundConds[ID];

  else {
    depolarisationBoundConds.resize(ID+1);
   return depolarisationBoundConds[ID];
  }

}

dbVector& IntegrationPointX::getDeltaDepolarisationBoundConds(int ID) {

  if(ID < deltaDepolarisationBoundConds.size())
    return deltaDepolarisationBoundConds[ID];

  else {
    deltaDepolarisationBoundConds.resize(ID+1);
    return deltaDepolarisationBoundConds[ID];
  }

}

dbVector& IntegrationPointX::getInitialDepolarisationBoundConds(int ID) {

  if(ID < initialDepolarisationBoundConds.size())
    return initialDepolarisationBoundConds[ID];

  else {
    initialDepolarisationBoundConds.resize(ID+1);
    return initialDepolarisationBoundConds[ID];
  }

}

/***********************************************************************/
// Return the deformation gradient tensor
dbMatrix& IntegrationPointX::getDeformationGradient() {

  using namespace std;

  // Check first if the deformation gradient is already set.
  if(deformationGradient.size() > 0)
    return deformationGradient;

  else {
    dbMatrix delta = getKroneckerSymbol(3);
    deformationGradient = delta;
    return deformationGradient;
  }

}

/***********************************************************************/
// Return the intermediate deformation gradient which is a history
// variable needed when undeformed configuration is unknown; transfers
// undeformed to deformed configuration
dbMatrix& IntegrationPointX::getIntermediateDefGradient() {

  using namespace std;

  // Check first if the deformation gradient is already set.
  if(intermediateDefGradient.size() > 0)
    return intermediateDefGradient;

  else {
    dbMatrix delta = getKroneckerSymbol(3);
    intermediateDefGradient = delta;
    return intermediateDefGradient;
  }

}

double& IntegrationPointX::getIntermediateJacobian() {

  using namespace std;

  // Check first if the deformation gradient is already set.
  if(intermediateJacobian != 0)
    return intermediateJacobian;

  else {
    calcDetDouble(getIntermediateDefGradient(),intermediateJacobian);
    return intermediateJacobian;
  }

}

/***********************************************************************/
// Return the rotation tensor
dbMatrix& IntegrationPointX::getRotationTens() {

  using namespace std;

  // Check first if the rotation tensor is already set.
  if(rotationTensor.size() > 0)
    return rotationTensor;

  else {
    dbMatrix delta = getKroneckerSymbol(3);
    rotationTensor = delta;
    return rotationTensor;
  }

}




