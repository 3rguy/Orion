// Stores all properties of a single particle.

#include "ParticleX.h"

ParticleX::ParticleX(int usedDOF) : ID(0),depolarisationTime(0),
                                  activeTension(0),activeTensionVariation(0),
				  intermediateJacobian(0) {
  coords = dbVector(3); 
  influenceRadii = dbVector(3);

  stepDegreesOfFreedomMat = dbMatrix(0);
  degreesOfFreedom = dbVector(usedDOF);
  suppNBoundPtcls = intMatrix(usedDOF);
  suppBoundPtcls = intMatrix(usedDOF);
  nBoundShapeFuncs = dbMatrix(usedDOF);
  boundShapeFuncs = dbMatrix(usedDOF);

  firstDerivShapeFuncs.resize(3);
  secondDerivShapeFuncs.resize(6);

  localMinPtcleDistance = 0;
  beta = 0;
  betaDerivs = dbVector(3);
}

ParticleX::~ParticleX() {}

void ParticleX::setID(int idx) {

  ID = idx;
}

void ParticleX::setWeight(int ID,double value) {

  if(ID >= weights.size())
    weights.resize(ID+1);

  weights[ID] = value;
}

double& ParticleX::getIntWeight(int ID) {

  if(ID < weights.size())
    return weights[ID];

  else {
    weights.resize(ID+1);
    return weights[ID];
  }

}

void ParticleX::setElems(int elem) {


  int pos = findIntVecPos(elem,0,elems.size(),elems);

  // Key was not found.
  if(pos == -1)
    elems.push_back(elem);

}

void ParticleX::setCoords(double coord1,double coord2,double coord3) {

  coords[0] = coord1;
  coords[1] = coord2;
  coords[2] = coord3;
}

void ParticleX::setRadii(double& rx,double& ry,double& rz) {

  influenceRadii[0] = rx;
  influenceRadii[1] = ry;
  influenceRadii[2] = rz;
}

void ParticleX::setDOF(int idx,double& value) {


  degreesOfFreedom[idx] = value;
}

void ParticleX::setDOFs(dbVector& dofs) {

  degreesOfFreedom = dofs;
}
/*
  void ParticleX::setOldDOF(int idx,double& value) {


  oldDegreesOfFreedom[idx] = value;
  }

  void ParticleX::setOldDOFs(dbVector& dofs) {

  oldDegreesOfFreedom = dofs;
  }
*/

void ParticleX::setLocalMinPtcleDistance(double localMinPtcleDist) {

  localMinPtcleDistance = localMinPtcleDist;
}

void ParticleX::setBeta(double betaVal) {

  beta = betaVal;
}


void ParticleX::setBetaDerivs(dbVector betaDer) {

  betaDerivs = betaDer;
}


void ParticleX::setShapeFuncsSize(int& size,int& order) {

  if(order == 0)
    shapeFuncs.resize(size);

  else if(order == 1) {
    shapeFuncs.resize(size);

    for(int i=0;i<firstDerivShapeFuncs.size();i++)
      firstDerivShapeFuncs[i].resize(size);

  }
  else if(order == 2) {
    shapeFuncs.resize(size);

    for(int i=0;i<firstDerivShapeFuncs.size();i++)
      firstDerivShapeFuncs[i].resize(size);

    for(int i=0;i<secondDerivShapeFuncs.size();i++)
      secondDerivShapeFuncs[i].resize(size);

  }
  else {
    std::cerr<<"Shape functions derivation order is not support in "
	     <<"ParticleX::setShapeFuncsSize!"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

void ParticleX::setShapeFuncs(int& size,dbVector& sFuncs) {

  if((unsigned int)size > shapeFuncs.size())
    shapeFuncs.resize(size);

  for(int i=0;i<size;i++)
    shapeFuncs[i] = sFuncs[i];
}

void ParticleX::setFirstDerivShapes(int& size,
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

void ParticleX::setSecondDerivShapes(int& size,
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

void ParticleX::setInflSpheres(int& size,intVector& iSpheres) {

  using namespace std;

  influencingSpheres.insert(influencingSpheres.begin(),iSpheres.begin(),
			    iSpheres.begin()+size);
}

/***********************************************************************/
void ParticleX::setNBShapeFuncs(int& DOF,int& supportSize,
                               dbVector& shapes) {

  if((unsigned int)supportSize > nBoundShapeFuncs[DOF].size())
    nBoundShapeFuncs[DOF].resize(supportSize);

  for(int i=0;i<supportSize;i++)
    nBoundShapeFuncs[DOF][i] = shapes[i];
}

void ParticleX::setBShapeFunc(int& DOF,int& idx,double& value) {

  if((unsigned int)idx < boundShapeFuncs[DOF].size())
    boundShapeFuncs[DOF][idx] = value;
  else
    boundShapeFuncs[DOF].push_back(value);

}

void ParticleX::setNBShapeFunc(int& DOF,int& idx,double& value) {

  if((unsigned int)idx < nBoundShapeFuncs[DOF].size())
    nBoundShapeFuncs[DOF][idx] = value;
  else
    nBoundShapeFuncs[DOF].push_back(value);

}

/***********************************************************************/
double& ParticleX::getSurfacePressure(int ID) {

  if(ID < surfacePressureLoads.size())
    return surfacePressureLoads[ID];

  else {
    surfacePressureLoads.resize(ID+1);
    return surfacePressureLoads[ID];
  }

}

blVector& ParticleX::getBodyForceDOF(int ID) {

  if(ID < bodyForceDOF.size())
    return bodyForceDOF[ID];

  else {
    bodyForceDOF.resize(ID+1);
    return bodyForceDOF[ID];
  }

}

dbVector& ParticleX::getBodyForce(int ID) {

  if(ID < bodyForceLoads.size())
    return bodyForceLoads[ID];

  else {
    bodyForceLoads.resize(ID+1);
    return bodyForceLoads[ID];
  }

}

blVector& ParticleX::getTractionDOF(int ID) {

  if(ID < tractionDOF.size())
    return tractionDOF[ID];

  else {
    tractionDOF.resize(ID+1);
    return tractionDOF[ID];
  }

}

dbVector& ParticleX::getTraction(int ID) {

  if(ID < tractionLoads.size())
    return tractionLoads[ID];

  else {
    tractionLoads.resize(ID+1);
    return tractionLoads[ID];
  }

}


blVector& ParticleX::getLineForceDOF(int ID) {

  if(ID < lineForceDOF.size())
    return lineForceDOF[ID];

  else {
    lineForceDOF.resize(ID+1);
    return lineForceDOF[ID];
  }

}

dbVector& ParticleX::getLineForce(int ID) {

  if(ID < lineForceLoads.size())
    return lineForceLoads[ID];

  else {
    lineForceLoads.resize(ID+1);
    return lineForceLoads[ID];
  }

}

blVector& ParticleX::getPointForceDOF(int ID) {

  if(ID < pointForceDOF.size())
    return pointForceDOF[ID];

  else {
    pointForceDOF.resize(ID+1);
    return pointForceDOF[ID];
  }

}

dbVector& ParticleX::getPointForce(int ID) {

  if(ID < pointForceLoads.size())
    return pointForceLoads[ID];

  else {
    pointForceLoads.resize(ID+1);
    return pointForceLoads[ID];
  }

}

dbVector& ParticleX::getSurfaceNormal(int ID) {

  if(ID < surfaceNormals.size())
    return surfaceNormals[ID];

  else {
    surfaceNormals.resize(ID+1);
    return surfaceNormals[ID];
  }

}

/***********************************************************************/
// Dirichlet boundary conditions

// deformation boundary conditions
blVector& ParticleX::getDeformationBoundDOF(int ID) {

  if(ID < deformationBoundDOF.size())
    return deformationBoundDOF[ID];

  else {
    deformationBoundDOF.resize(ID+1);
    return deformationBoundDOF[ID];
  }

}

dbVector& ParticleX::getDeformationBoundConds(int ID) {

  if(ID < deformationBoundConds.size())
    return deformationBoundConds[ID];

  else {
    deformationBoundConds.resize(ID+1);
    return deformationBoundConds[ID];
  }

}

dbVector& ParticleX::getDeltaDeformationBoundConds(int ID) {

  if(ID < deltaDeformationBoundConds.size())
    return deltaDeformationBoundConds[ID];

  else {
    deltaDeformationBoundConds.resize(ID+1);
    return deltaDeformationBoundConds[ID];
  }

}

dbVector& ParticleX::getInitialDeformationBoundConds(int ID) {

  if(ID < initialDeformationBoundConds.size())
    return initialDeformationBoundConds[ID];

  else {
    initialDeformationBoundConds.resize(ID+1);
    return initialDeformationBoundConds[ID];
  }

}

// ----------------------------------------------------------------------
// electric boundary conditions

blVector& ParticleX::getElectricBoundDOF(int ID) {

  if(ID < electricBoundDOF.size())
    return electricBoundDOF[ID];

  else {
    electricBoundDOF.resize(ID+1);
    return electricBoundDOF[ID];
  }

}

dbVector& ParticleX::getElectricBoundConds(int ID) {

  if(ID < electricBoundConds.size())
    return electricBoundConds[ID];

  else {
    electricBoundConds.resize(ID+1);
    return electricBoundConds[ID];
  }

}

dbVector& ParticleX::getDeltaElectricBoundConds(int ID) {

  if(ID < deltaElectricBoundConds.size())
    return deltaElectricBoundConds[ID];

  else {
    deltaElectricBoundConds.resize(ID+1);
    return deltaElectricBoundConds[ID];
  }

}

dbVector& ParticleX::getInitialElectricBoundConds(int ID) {

  if(ID < initialElectricBoundConds.size())
    return initialElectricBoundConds[ID];

  else {
    initialElectricBoundConds.resize(ID+1);
    return initialElectricBoundConds[ID];
  }

}


// ----------------------------------------------------------------------
// depolarisation boundary conditions

blVector& ParticleX::getDepolarisationBoundDOF(int ID) {

  if(ID < depolarisationBoundDOF.size())
    return depolarisationBoundDOF[ID];

  else {
    depolarisationBoundDOF.resize(ID+1);
    return depolarisationBoundDOF[ID];
  }

}

dbVector& ParticleX::getDepolarisationBoundConds(int ID) {

  if(ID < depolarisationBoundConds.size())
    return depolarisationBoundConds[ID];

  else {
    depolarisationBoundConds.resize(ID+1);
    return depolarisationBoundConds[ID];
  }

}

dbVector& ParticleX::getDeltaDepolarisationBoundConds(int ID) {

  if(ID < deltaDepolarisationBoundConds.size())
    return deltaDepolarisationBoundConds[ID];

  else {
    deltaDepolarisationBoundConds.resize(ID+1);
    return deltaDepolarisationBoundConds[ID];
  }

}

dbVector& ParticleX::getInitialDepolarisationBoundConds(int ID) {

  if(ID < initialDepolarisationBoundConds.size())
    return initialDepolarisationBoundConds[ID];

  else {
    initialDepolarisationBoundConds.resize(ID+1);
    return initialDepolarisationBoundConds[ID];
  }

}

// ----------------------------------------------------------------------
// micro boundary conditions

blVector& ParticleX::getMicroBoundDOF(int ID) {

  if(ID < microBoundDOF.size())
    return microBoundDOF[ID];

  else {
    microBoundDOF.resize(ID+1);
    return microBoundDOF[ID];
  }

}

dbVector& ParticleX::getMicroBoundConds(int ID) {

  if(ID < microBoundConds.size())
    return microBoundConds[ID];

  else {
    microBoundConds.resize(ID+1);
    return microBoundConds[ID];
  }

}

dbVector& ParticleX::getDeltaMicroBoundConds(int ID) {

  if(ID < deltaMicroBoundConds.size())
    return deltaMicroBoundConds[ID];

  else {
    deltaMicroBoundConds.resize(ID+1);
    return deltaMicroBoundConds[ID];
  }

}

dbVector& ParticleX::getInitialMicroBoundConds(int ID) {

  if(ID < initialMicroBoundConds.size())
    return initialMicroBoundConds[ID];

  else {
    initialMicroBoundConds.resize(ID+1);
    return initialMicroBoundConds[ID];
  }

}

/***********************************************************************/
// Return the deformation gradient tensor 
dbMatrix& ParticleX::getDeformationGradient() {

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
dbMatrix& ParticleX::getIntermediateDefGradient() {

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

double& ParticleX::getIntermediateJacobian() {

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
/***********************************************************************/
// Check if a point is supported by this particle.
bool ParticleX::querySupported(InputFileData* InputData,
                              dbVector& pointCoords,
                              std::map<std::string,double>& modelData,
                              std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  bool plusMinusDependent =
    (bool)InputData->getValue("plusMinusDirectionDependentRadius");
  int choice = (int)InputData->getValue("shapefunctionType");
  bool supported = true;
  double dist;
  int counter=0;



  // Maximum Entropy
  if(choice == 6){

    for (int i=0; i<usedDims; i++){
      dist = fabs(pointCoords[i]-coords[i]);
      if (dist < influenceRadii[i]){
	counter +=1;
      }
    }
    if (counter == 3){
      supported = true;
    }
    else {
      supported = false;
    }
  }

  // all other shape functions
  else {

    if(!plusMinusDependent) {

      // check if point is within the particle's support zone
      for(int i=0;i<usedDims;i++) {

	if(fabs(pointCoords[i] - coords[i]) >= influenceRadii[i])
	  supported = false;

#ifdef _supportDebugMode_
	logFile<<pointCoords[i]<<"-("<<coords[i]<<") = "
	       <<fabs(pointCoords[i] - coords[i])<<" <=> "
	       <<influenceRadii[i]<<"("<<supported<<")"<<endl;
#endif

      }

    }

    else {

      double dist;

      // choose the radii for the positive or negative coordinate direction
      // according to the position of the point within the support-zone

      for(int i=0;i<usedDims;i++) {

	dist = pointCoords[i] - coords[i];

	if(dist > 0 && fabs(dist) >= influenceRadii[i])
	  supported = false;

	else if(dist < 0 && fabs(dist) >= influenceRadii[i+usedDims])
	  supported = false;
      }
    }
  }

  return supported;
}

/***********************************************************************/
/***********************************************************************/
// Check if a point is supported by this particle.
bool ParticleX::querySupported(InputFileData* InputData,
                              dbVector& pointCoords,
                              double& ptcleDist,
                              std::map<std::string,double>& modelData,
                              std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  bool plusMinusDependent =
    (bool)InputData->getValue("plusMinusDirectionDependentRadius");
  int choice = (int)InputData->getValue("shapefunctionType");
  bool supported = true;
  double dist,dist_sqr=0;
  int counter=0;

    

  // Calculate the distance between two particles
  for (int j=0;j<usedDims;j++){
    dist_sqr += pow(pointCoords[j]-coords[j],2);
  }
  ptcleDist = sqrt(dist_sqr);

  // Maximum Entropy
  if(choice == 6){

    for (int i=0; i<usedDims; i++){
      dist = fabs(pointCoords[i]-coords[i]);
      if (dist < influenceRadii[i]){
	counter +=1;
      }
    }
    if (counter == 3){
      supported = true;
    }
    else {
      supported = false;
    }
  }

  // all other shape functions
  else {

    if(!plusMinusDependent) {

      // check if point is within the particle's support zone
      for(int i=0;i<usedDims;i++) {

	if(fabs(pointCoords[i] - coords[i]) >= influenceRadii[i])
	  supported = false;

#ifdef _supportDebugMode_
	logFile<<pointCoords[i]<<"-"<<coords[i]<<" <=> "<<influenceRadii[i]<<endl;
#endif

      }

    }

    else {

      double dist;

      // choose the radii for the positive or negative coordinate direction
      // according to the position of the point within the support-zone

      for(int i=0;i<usedDims;i++) {

	dist = pointCoords[i] - coords[i];

	if(dist > 0 && fabs(dist) >= influenceRadii[i])
	  supported = false;

	else if(dist < 0 && fabs(dist) >= influenceRadii[i+usedDims])
	  supported = false;
      }
    }
  }

  return supported;
}



/**********************************************************************/
/**********************************************************************/
// Clear all class arrays which are not locally needed.
void ParticleX::clearArrays() {

  using namespace std;

  resizeArray(coords,0);
  resizeArray(influenceRadii,0);
  resizeArray(degreesOfFreedom,0);
  resizeArray(supportingPtcls,0);
  resizeArray(influencingSpheres,0);

  resizeArray(shapeFuncs,0);
  resizeArray(firstDerivShapeFuncs,0);
  resizeArray(secondDerivShapeFuncs,0);

  resizeArray(elems,0);

  //   for(int j=0;firstDerivShapeFuncs.size();j++) {
  //     resizeArray(firstDerivShapeFuncs[j],0);

  //   for(int j=0;secondDerivShapeFuncs.size();j++)
  //     resizeArray(secondDerivShapeFuncs[j],0);

}

