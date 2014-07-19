#include "FEMElement.h"

FEMElement::FEMElement(int usedDOF) : elementOrder(0),elementType(0),
				      globalID(0),materialID(0) {

  numOfNodeDOF = usedDOF;

}

void FEMElement::setGlobalID(int globalIdx) {

  globalID = globalIdx;
}

dbVector& FEMElement::getSurfaceNormal(int ID) { 

  if(ID < surfaceNormals.size())
    return surfaceNormals[ID]; 

  else {
    surfaceNormals.resize(ID+1);
    return surfaceNormals[ID];  
  }

}

dbMatrix3& FEMElement::getJacobian(std::vector<Particle>& particles,
				   std::ofstream& logFile) {

  if(jacobian.size() == 0) {
    resizeArray(jacobian,VolumeGaussSet->coord.size());

    // loop over all element Gauss points
    for(int i=0;i<VolumeGaussSet->coord.size();i++)
      FEVolumeSet->jacobian(particles,nodes,VolumeGaussSet->coord[i],
			    jacobian[i],logFile);

  }

  return jacobian;
}

/***********************************************************************/
// loading boundary conditions

double& FEMElement::getSurfacePressure(int ID) { 

  if(ID < surfacePressureLoads.size())
    return surfacePressureLoads[ID]; 

  else {
    surfacePressureLoads.resize(ID+1);
    return surfacePressureLoads[ID]; 
  }

}

blVector& FEMElement::getBodyForceDOF(int ID) { 

  if(ID < bodyForceDOF.size())
    return bodyForceDOF[ID]; 

  else {
    bodyForceDOF.resize(ID+1);
    return bodyForceDOF[ID]; 
  }

}

dbVector& FEMElement::getBodyForce(int ID) { 

  if(ID < bodyForceLoads.size())
    return bodyForceLoads[ID]; 

  else {
    bodyForceLoads.resize(ID+1);
   return bodyForceLoads[ID]; 
  }

}

blVector& FEMElement::getTractionDOF(int ID) { 

  if(ID < tractionDOF.size())
    return tractionDOF[ID]; 

  else {
    tractionDOF.resize(ID+1);
    return tractionDOF[ID]; 
  }

}

dbVector& FEMElement::getTraction(int ID) { 

  if(ID < tractionLoads.size())
    return tractionLoads[ID]; 

  else {
    tractionLoads.resize(ID+1);
    return tractionLoads[ID]; 
  }

}


blVector& FEMElement::getLineForceDOF(int ID) { 

  if(ID < lineForceDOF.size())
    return lineForceDOF[ID]; 

  else {
    lineForceDOF.resize(ID+1);
    return lineForceDOF[ID]; 
  }

}

dbVector& FEMElement::getLineForce(int ID) { 

  if(ID < lineForceLoads.size())
    return lineForceLoads[ID]; 

  else {
    lineForceLoads.resize(ID+1);
    return lineForceLoads[ID]; 
  }

}

blVector& FEMElement::getBodyMomentDOF(int ID) { 

  if(ID < bodyMomentDOF.size())
    return bodyMomentDOF[ID]; 

  else {
    bodyMomentDOF.resize(ID+1);
    return bodyMomentDOF[ID]; 
  }

}

dbVector& FEMElement::getBodyMoment(int ID) { 

  if(ID < bodyMomentLoads.size())
    return bodyMomentLoads[ID]; 

  else {
    bodyMomentLoads.resize(ID+1);
   return bodyMomentLoads[ID]; 
  }

}

blVector& FEMElement::getSurfaceMomentDOF(int ID) { 

  if(ID < surfaceMomentDOF.size())
    return surfaceMomentDOF[ID]; 

  else {
    surfaceMomentDOF.resize(ID+1);
    return surfaceMomentDOF[ID]; 
  }

}

dbVector& FEMElement::getSurfaceMoment(int ID) { 

  if(ID < surfaceMomentLoads.size())
    return surfaceMomentLoads[ID]; 

  else {
    surfaceMomentLoads.resize(ID+1);
    return surfaceMomentLoads[ID]; 
  }

}

blVector& FEMElement::getSurfaceElectricChargeDOF(int ID) { 

  if(ID < surfaceElectricChargeDOF.size())
    return surfaceElectricChargeDOF[ID]; 

  else {
    surfaceElectricChargeDOF.resize(ID+1);
    return surfaceElectricChargeDOF[ID]; 
  }

}

dbVector& FEMElement::getSurfaceElectricCharge(int ID) { 

  if(ID < surfaceElectricChargeLoads.size())
    return surfaceElectricChargeLoads[ID]; 

  else {
    surfaceElectricChargeLoads.resize(ID+1);
    return surfaceElectricChargeLoads[ID]; 
  }

}

blVector& FEMElement::getBodyElectricChargeDOF(int ID) { 

  if(ID < bodyElectricChargeDOF.size())
    return bodyElectricChargeDOF[ID]; 

  else {
    bodyElectricChargeDOF.resize(ID+1);
    return bodyElectricChargeDOF[ID]; 
  }

}

dbVector& FEMElement::getBodyElectricCharge(int ID) { 

  if(ID < bodyElectricChargeLoads.size())
    return bodyElectricChargeLoads[ID]; 

  else {
    bodyElectricChargeLoads.resize(ID+1);
    return bodyElectricChargeLoads[ID]; 
  }

}

/***********************************************************************/
// deformation boundary conditions
blVector& FEMElement::getLineDeformationBoundDOF(int ID) { 

  if(ID < lineDefBoundDOF.size())
    return lineDefBoundDOF[ID]; 

  else {
    lineDefBoundDOF.resize(ID+1);
   return lineDefBoundDOF[ID]; 
  }

}

dbVector& FEMElement::getLineDeformationBoundConds(int ID) { 

  if(ID < lineDefBoundConds.size())
    return lineDefBoundConds[ID]; 

  else {
    lineDefBoundConds.resize(ID+1);
   return lineDefBoundConds[ID]; 
  }

}

blVector& FEMElement::getSurfaceDeformationBoundDOF(int ID) { 

  if(ID < surfaceDefBoundDOF.size())
    return surfaceDefBoundDOF[ID]; 

  else {
    surfaceDefBoundDOF.resize(ID+1);
   return surfaceDefBoundDOF[ID];  
  }

}

dbVector& FEMElement::getSurfaceDeformationBoundConds(int ID) { 

  if(ID < surfaceDefBoundConds.size())
    return surfaceDefBoundConds[ID]; 

  else {
    surfaceDefBoundConds.resize(ID+1);
    return surfaceDefBoundConds[ID]; 
  }

}

/***********************************************************************/
// electrical boundary conditions

blVector& FEMElement::getLineElectricBoundDOF(int ID) { 

  if(ID < lineElectricBoundDOF.size())
    return lineElectricBoundDOF[ID]; 

  else {
    lineElectricBoundDOF.resize(ID+1);
   return lineElectricBoundDOF[ID]; 
  }

}

dbVector& FEMElement::getLineElectricBoundConds(int ID) { 

  if(ID < lineElectricBoundConds.size())
    return lineElectricBoundConds[ID]; 

  else {
    lineElectricBoundConds.resize(ID+1);
   return lineElectricBoundConds[ID]; 
  }

}

blVector& FEMElement::getSurfaceElectricBoundDOF(int ID) { 

  if(ID < surfaceElectricBoundDOF.size())
    return surfaceElectricBoundDOF[ID]; 

  else {
    surfaceElectricBoundDOF.resize(ID+1);
   return surfaceElectricBoundDOF[ID];  
  }

}

dbVector& FEMElement::getSurfaceElectricBoundConds(int ID) { 

  if(ID < surfaceElectricBoundConds.size())
    return surfaceElectricBoundConds[ID]; 

  else {
    surfaceElectricBoundConds.resize(ID+1);
    return surfaceElectricBoundConds[ID]; 
  }

}

/***********************************************************************/
// depolarisation boundary conditions

blVector& FEMElement::getLineDepolarisationBoundDOF(int ID) { 

  if(ID < lineDepolarisationBoundDOF.size())
    return lineDepolarisationBoundDOF[ID]; 

  else {
    lineDepolarisationBoundDOF.resize(ID+1);
   return lineDepolarisationBoundDOF[ID]; 
  }

}

dbVector& FEMElement::getLineDepolarisationBoundConds(int ID) { 

  if(ID < lineDepolarisationBoundConds.size())
    return lineDepolarisationBoundConds[ID]; 

  else {
    lineDepolarisationBoundConds.resize(ID+1);
   return lineDepolarisationBoundConds[ID]; 
  }

}

blVector& FEMElement::getSurfaceDepolarisationBoundDOF(int ID) { 

  if(ID < surfaceDepolarisationBoundDOF.size())
    return surfaceDepolarisationBoundDOF[ID]; 

  else {
    surfaceDepolarisationBoundDOF.resize(ID+1);
   return surfaceDepolarisationBoundDOF[ID];  
  }

}

dbVector& FEMElement::getSurfaceDepolarisationBoundConds(int ID) { 

  if(ID < surfaceDepolarisationBoundConds.size())
    return surfaceDepolarisationBoundConds[ID]; 

  else {
    surfaceDepolarisationBoundConds.resize(ID+1);
    return surfaceDepolarisationBoundConds[ID]; 
  }

}

/***********************************************************************/
// micro boundary conditions

blVector& FEMElement::getLineMicroBoundDOF(int ID) { 

  if(ID < lineMicroBoundDOF.size())
    return lineMicroBoundDOF[ID]; 

  else {
    lineMicroBoundDOF.resize(ID+1);
   return lineMicroBoundDOF[ID]; 
  }

}

dbVector& FEMElement::getLineMicroBoundConds(int ID) { 

  if(ID < lineMicroBoundConds.size())
    return lineMicroBoundConds[ID]; 

  else {
    lineMicroBoundConds.resize(ID+1);
   return lineMicroBoundConds[ID]; 
  }

}

blVector& FEMElement::getSurfaceMicroBoundDOF(int ID) { 

  if(ID < surfaceMicroBoundDOF.size())
    return surfaceMicroBoundDOF[ID]; 

  else {
    surfaceMicroBoundDOF.resize(ID+1);
   return surfaceMicroBoundDOF[ID];  
  }

}

dbVector& FEMElement::getSurfaceMicroBoundConds(int ID) { 

  if(ID < surfaceMicroBoundConds.size())
    return surfaceMicroBoundConds[ID]; 

  else {
    surfaceMicroBoundConds.resize(ID+1);
    return surfaceMicroBoundConds[ID]; 
  }

}


