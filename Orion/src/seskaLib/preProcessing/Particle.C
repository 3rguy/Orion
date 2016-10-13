// Stores all properties of a single particle.

#include "Particle.h"

Particle::Particle(int usedDOF) :
    ID(0), depolarisationTime(0), activeTension(0), activeTensionVariation(0),
    intermediateJacobian(0) {
  coords = dbVector(3);
  influenceRadii = dbVector(3);

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

Particle::~Particle() {
}

void
Particle::setID(int idx) {

  ID = idx;
}

void
Particle::setWeight(int ID,double value) {

  if(ID >= weights.size()) weights.resize(ID + 1);

  weights[ID] = value;
}

double&
Particle::getIntWeight(int ID) {

  if(ID < weights.size()) return weights[ID];

  else {
    weights.resize(ID + 1);
    return weights[ID];
  }

}

void
Particle::setElems(int elem) {

  int pos = findIntVecPos(elem,0,elems.size(),elems);

  // Key was not found.
  if(pos == -1) elems.push_back(elem);

}

void
Particle::setCoords(double coord1,double coord2,double coord3) {

  coords[0] = coord1;
  coords[1] = coord2;
  coords[2] = coord3;
}

void
Particle::setRadii(double& rx,double& ry,double& rz) {

  influenceRadii[0] = rx;
  influenceRadii[1] = ry;
  influenceRadii[2] = rz;
}

void
Particle::setDOF(int idx,double& value) {

  degreesOfFreedom[idx] = value;
}

void
Particle::setDOFs(dbVector& dofs) {

  degreesOfFreedom = dofs;
}

void
Particle::setOldDOF(int idx,double& value) {
  
  oldDegreesOfFreedom[idx] = value;
}

void
Particle::setOldDOFs(dbVector& dofs) {
  
  oldDegreesOfFreedom = dofs;
}

void
Particle::setLocalMinPtcleDistance(double localMinPtcleDist) {

  localMinPtcleDistance = localMinPtcleDist;
}

void
Particle::setBeta(double betaVal) {

  beta = betaVal;
}

void
Particle::setBetaDerivs(dbVector betaDer) {

  betaDerivs = betaDer;
}

void
Particle::setShapeFuncsSize(int& size,int& order) {

  if(order == 0) shapeFuncs.resize(size);

  else if(order == 1) {
    shapeFuncs.resize(size);

    for(int i = 0;i < firstDerivShapeFuncs.size();i++)
      firstDerivShapeFuncs[i].resize(size);

  }
  else if(order == 2) {
    shapeFuncs.resize(size);

    for(int i = 0;i < firstDerivShapeFuncs.size();i++)
      firstDerivShapeFuncs[i].resize(size);

    for(int i = 0;i < secondDerivShapeFuncs.size();i++)
      secondDerivShapeFuncs[i].resize(size);

  }
  else {
    std::cerr << "Shape functions derivation order is not support in "
        << "Particle::setShapeFuncsSize!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

void
Particle::setShapeFuncs(int& size,dbVector& sFuncs) {

  if((unsigned int) size > shapeFuncs.size()) shapeFuncs.resize(size);

  for(int i = 0;i < size;i++)
    shapeFuncs[i] = sFuncs[i];
}

void
Particle::setFirstDerivShapes(int& size,dbMatrix& firstDerivShapes) {

  if((unsigned int) size > firstDerivShapeFuncs[0].size()) {
    firstDerivShapeFuncs[0].resize(size);
    firstDerivShapeFuncs[1].resize(size);
    firstDerivShapeFuncs[2].resize(size);
  }

  for(int i = 0;i < size;i++) {
    firstDerivShapeFuncs[0][i] = firstDerivShapes[0][i];
    firstDerivShapeFuncs[1][i] = firstDerivShapes[1][i];
    firstDerivShapeFuncs[2][i] = firstDerivShapes[2][i];
  }
}

void
Particle::setSecondDerivShapes(int& size,dbMatrix& secondDerivShapes) {

  if((unsigned int) size > secondDerivShapeFuncs[0].size()) {
    secondDerivShapeFuncs[0].resize(size);
    secondDerivShapeFuncs[1].resize(size);
    secondDerivShapeFuncs[2].resize(size);
    secondDerivShapeFuncs[3].resize(size);
    secondDerivShapeFuncs[4].resize(size);
    secondDerivShapeFuncs[5].resize(size);
  }

  for(int i = 0;i < size;i++) {
    secondDerivShapeFuncs[0][i] = secondDerivShapes[0][i];
    secondDerivShapeFuncs[1][i] = secondDerivShapes[1][i];
    secondDerivShapeFuncs[2][i] = secondDerivShapes[2][i];

    secondDerivShapeFuncs[3][i] = secondDerivShapes[3][i];
    secondDerivShapeFuncs[4][i] = secondDerivShapes[4][i];
    secondDerivShapeFuncs[5][i] = secondDerivShapes[5][i];
  }
}

void
Particle::setInflSpheres(int& size,intVector& iSpheres) {

  using namespace std;

  influencingSpheres.insert(influencingSpheres.begin(),iSpheres.begin(),
                            iSpheres.begin() + size);
}

/***********************************************************************/
void
Particle::setNBShapeFuncs(int& DOF,int& supportSize,dbVector& shapes) {

  if((unsigned int) supportSize > nBoundShapeFuncs[DOF].size()) nBoundShapeFuncs[DOF].resize(
      supportSize);

  for(int i = 0;i < supportSize;i++)
    nBoundShapeFuncs[DOF][i] = shapes[i];
}

void
Particle::setBShapeFunc(int& DOF,int& idx,double& value) {

  if((unsigned int) idx < boundShapeFuncs[DOF].size()) boundShapeFuncs[DOF][idx] =
    value;
  else boundShapeFuncs[DOF].push_back(value);

}

void
Particle::setNBShapeFunc(int& DOF,int& idx,double& value) {

  if((unsigned int) idx < nBoundShapeFuncs[DOF].size()) nBoundShapeFuncs[DOF][idx] =
    value;
  else nBoundShapeFuncs[DOF].push_back(value);

}

/***********************************************************************/
double&
Particle::getSurfacePressure(int ID) {

  if(ID < surfacePressureLoads.size()) return surfacePressureLoads[ID];

  else {
    surfacePressureLoads.resize(ID + 1);
    return surfacePressureLoads[ID];
  }

}

blVector&
Particle::getBodyForceDOF(int ID) {

  if(ID < bodyForceDOF.size()) return bodyForceDOF[ID];

  else {
    bodyForceDOF.resize(ID + 1);
    return bodyForceDOF[ID];
  }

}

dbVector&
Particle::getBodyForce(int ID) {

  if(ID < bodyForceLoads.size()) return bodyForceLoads[ID];

  else {
    bodyForceLoads.resize(ID + 1);
    return bodyForceLoads[ID];
  }

}

blVector&
Particle::getTractionDOF(int ID) {

  if(ID < tractionDOF.size()) return tractionDOF[ID];

  else {
    tractionDOF.resize(ID + 1);
    return tractionDOF[ID];
  }

}

dbVector&
Particle::getTraction(int ID) {

  if(ID < tractionLoads.size()) return tractionLoads[ID];

  else {
    tractionLoads.resize(ID + 1);
    return tractionLoads[ID];
  }

}

blVector&
Particle::getLineForceDOF(int ID) {

  if(ID < lineForceDOF.size()) return lineForceDOF[ID];

  else {
    lineForceDOF.resize(ID + 1);
    return lineForceDOF[ID];
  }

}

dbVector&
Particle::getLineForce(int ID) {

  if(ID < lineForceLoads.size()) return lineForceLoads[ID];

  else {
    lineForceLoads.resize(ID + 1);
    return lineForceLoads[ID];
  }

}

blVector&
Particle::getPointForceDOF(int ID) {

  if(ID < pointForceDOF.size()) return pointForceDOF[ID];

  else {
    pointForceDOF.resize(ID + 1);
    return pointForceDOF[ID];
  }

}

dbVector&
Particle::getPointForce(int ID) {

  if(ID < pointForceLoads.size()) return pointForceLoads[ID];

  else {
    pointForceLoads.resize(ID + 1);
    return pointForceLoads[ID];
  }

}

dbVector&
Particle::getSurfaceNormal(int ID) {

  if(ID < surfaceNormals.size()) return surfaceNormals[ID];

  else {
    surfaceNormals.resize(ID + 1);
    return surfaceNormals[ID];
  }

}

/***********************************************************************/
// Dirichlet boundary conditions
// deformation boundary conditions
blVector&
Particle::getDeformationBoundDOF(int ID) {

  if(ID < deformationBoundDOF.size()) return deformationBoundDOF[ID];

  else {
    deformationBoundDOF.resize(ID + 1);
    return deformationBoundDOF[ID];
  }

}

dbVector&
Particle::getDeformationBoundConds(int ID) {

  if(ID < deformationBoundConds.size()) return deformationBoundConds[ID];

  else {
    deformationBoundConds.resize(ID + 1);
    return deformationBoundConds[ID];
  }

}

dbVector&
Particle::getDeltaDeformationBoundConds(int ID) {

  if(ID < deltaDeformationBoundConds.size()) return deltaDeformationBoundConds[ID];

  else {
    deltaDeformationBoundConds.resize(ID + 1);
    return deltaDeformationBoundConds[ID];
  }

}

dbVector&
Particle::getInitialDeformationBoundConds(int ID) {

  if(ID < initialDeformationBoundConds.size()) return initialDeformationBoundConds[ID];

  else {
    initialDeformationBoundConds.resize(ID + 1);
    return initialDeformationBoundConds[ID];
  }

}

// ----------------------------------------------------------------------
// electric boundary conditions

blVector&
Particle::getElectricBoundDOF(int ID) {

  if(ID < electricBoundDOF.size()) return electricBoundDOF[ID];

  else {
    electricBoundDOF.resize(ID + 1);
    return electricBoundDOF[ID];
  }

}

dbVector&
Particle::getElectricBoundConds(int ID) {

  if(ID < electricBoundConds.size()) return electricBoundConds[ID];

  else {
    electricBoundConds.resize(ID + 1);
    return electricBoundConds[ID];
  }

}

dbVector&
Particle::getDeltaElectricBoundConds(int ID) {

  if(ID < deltaElectricBoundConds.size()) return deltaElectricBoundConds[ID];

  else {
    deltaElectricBoundConds.resize(ID + 1);
    return deltaElectricBoundConds[ID];
  }

}

dbVector&
Particle::getInitialElectricBoundConds(int ID) {

  if(ID < initialElectricBoundConds.size()) return initialElectricBoundConds[ID];

  else {
    initialElectricBoundConds.resize(ID + 1);
    return initialElectricBoundConds[ID];
  }

}

// ----------------------------------------------------------------------
// TPM boundary conditions

blVector&
Particle::getTPMBoundDOF(int ID) {

  if(ID < TPMBoundDOF.size()) return TPMBoundDOF[ID];

  else {
    TPMBoundDOF.resize(ID + 1);
    return TPMBoundDOF[ID];
  }

}

dbVector&
Particle::getTPMBoundConds(int ID) {

  if(ID < TPMBoundConds.size()) return TPMBoundConds[ID];

  else {
    TPMBoundConds.resize(ID + 1);
    return TPMBoundConds[ID];
  }

}

dbVector&
Particle::getDeltaTPMBoundConds(int ID) {

  if(ID < deltaTPMBoundConds.size()) return deltaTPMBoundConds[ID];

  else {
    deltaTPMBoundConds.resize(ID + 1);
    return deltaTPMBoundConds[ID];
  }

}

dbVector&
Particle::getInitialTPMBoundConds(int ID) {

  if(ID < initialTPMBoundConds.size()) return initialTPMBoundConds[ID];

  else {
    initialTPMBoundConds.resize(ID + 1);
    return initialTPMBoundConds[ID];
  }

}

// ----------------------------------------------------------------------
// depolarisation boundary conditions

blVector&
Particle::getDepolarisationBoundDOF(int ID) {

  if(ID < depolarisationBoundDOF.size()) return depolarisationBoundDOF[ID];

  else {
    depolarisationBoundDOF.resize(ID + 1);
    return depolarisationBoundDOF[ID];
  }

}

dbVector&
Particle::getDepolarisationBoundConds(int ID) {

  if(ID < depolarisationBoundConds.size()) return depolarisationBoundConds[ID];

  else {
    depolarisationBoundConds.resize(ID + 1);
    return depolarisationBoundConds[ID];
  }

}

dbVector&
Particle::getDeltaDepolarisationBoundConds(int ID) {

  if(ID < deltaDepolarisationBoundConds.size()) return deltaDepolarisationBoundConds[ID];

  else {
    deltaDepolarisationBoundConds.resize(ID + 1);
    return deltaDepolarisationBoundConds[ID];
  }

}

dbVector&
Particle::getInitialDepolarisationBoundConds(int ID) {

  if(ID < initialDepolarisationBoundConds.size()) return initialDepolarisationBoundConds[ID];

  else {
    initialDepolarisationBoundConds.resize(ID + 1);
    return initialDepolarisationBoundConds[ID];
  }

}

// ----------------------------------------------------------------------
// micro boundary conditions

blVector&
Particle::getMicroBoundDOF(int ID) {

  if(ID < microBoundDOF.size()) return microBoundDOF[ID];

  else {
    microBoundDOF.resize(ID + 1);
    return microBoundDOF[ID];
  }

}

dbVector&
Particle::getMicroBoundConds(int ID) {

  if(ID < microBoundConds.size()) return microBoundConds[ID];

  else {
    microBoundConds.resize(ID + 1);
    return microBoundConds[ID];
  }

}

dbVector&
Particle::getDeltaMicroBoundConds(int ID) {

  if(ID < deltaMicroBoundConds.size()) return deltaMicroBoundConds[ID];

  else {
    deltaMicroBoundConds.resize(ID + 1);
    return deltaMicroBoundConds[ID];
  }

}

dbVector&
Particle::getInitialMicroBoundConds(int ID) {

  if(ID < initialMicroBoundConds.size()) return initialMicroBoundConds[ID];

  else {
    initialMicroBoundConds.resize(ID + 1);
    return initialMicroBoundConds[ID];
  }

}

/***********************************************************************/
// Return the deformation gradient tensor 
dbMatrix&
Particle::getDeformationGradient() {

  using namespace std;

  // Check first if the deformation gradient is already set.
  if(deformationGradient.size() > 0) return deformationGradient;

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
dbMatrix&
Particle::getIntermediateDefGradient() {

  using namespace std;

  // Check first if the deformation gradient is already set.
  if(intermediateDefGradient.size() > 0) return intermediateDefGradient;

  else {
    dbMatrix delta = getKroneckerSymbol(3);
    intermediateDefGradient = delta;
    return intermediateDefGradient;
  }

}

double&
Particle::getIntermediateJacobian() {

  using namespace std;

  // Check first if the deformation gradient is already set.
  if(intermediateJacobian != 0) return intermediateJacobian;

  else {
    calcDetDoubleDense(getIntermediateDefGradient(),intermediateJacobian);
    return intermediateJacobian;
  }

}

/***********************************************************************/
/***********************************************************************/
// Check if a point is supported by this particle.
bool
Particle::querySupported(InputFileData* InputData,dbVector& pointCoords,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int windowFuncShape = (int) InputData->getValue("windowFunctionShape");
  bool plusMinusDependent = (bool) InputData->getValue(
      "plusMinusDirectionDependentRadius");
  int choice = (int) InputData->getValue("shapefunctionType");
  bool supported = true;
  double dist;
  int counter = 0;

  // Maximum Entropy
  if(choice == 6) {

    for(int i = 0;i < usedDims;i++) {
      dist = fabs(pointCoords[i] - coords[i]);
      if(dist < influenceRadii[i]) {
        counter += 1;
      }
    }
    if(counter == 3) {
      supported = true;
    }
    else {
      supported = false;
    }
  }

  /*********************************************************************/
  // all other shape functions
  else {

    // prismatic window function shape
    if(windowFuncShape == 1) {

      if( !plusMinusDependent) {

        // check if point is within the particle's support zone
        for(int i = 0;i < usedDims;i++) {

          if(fabs(pointCoords[i] - coords[i]) >= influenceRadii[i]) supported =
            false;

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

        for(int i = 0;i < usedDims;i++) {

          dist = pointCoords[i] - coords[i];

          if(dist > 0 && fabs(dist) >= influenceRadii[i]) supported = false;

          else if(dist < 0 && fabs(dist) >= influenceRadii[i + usedDims]) supported =
            false;
        }

      }

    }
    // -----------------------------------------------------------------
    // spherical window function
    else {

      if(calcDistance(pointCoords,coords) >= influenceRadii[0]) supported =
        false;

    }

  }

  return supported;
}

/***********************************************************************/
/***********************************************************************/
// Check if a point is supported by this particle.
bool
Particle::querySupported(InputFileData* InputData,dbVector& pointCoords,
                         double& ptcleDist,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int windowFuncShape = (int) InputData->getValue("windowFunctionShape");
  bool plusMinusDependent = (bool) InputData->getValue(
      "plusMinusDirectionDependentRadius");
  int choice = (int) InputData->getValue("shapefunctionType");
  bool supported = true;
  double dist,dist_sqr = 0;
  int counter = 0;

  // Calculate the distance between two particles
  for(int j = 0;j < usedDims;j++) {
    dist_sqr += pow(pointCoords[j] - coords[j],2);
  }
  ptcleDist = sqrt(dist_sqr);

  // Maximum Entropy
  if(choice == 6) {

    for(int i = 0;i < usedDims;i++) {
      dist = fabs(pointCoords[i] - coords[i]);
      if(dist < influenceRadii[i]) {
        counter += 1;
      }
    }
    if(counter == 3) {
      supported = true;
    }
    else {
      supported = false;
    }
  }
  /*********************************************************************/
  // all other shape functions
  else {

    // prismatic window function shape
    if(windowFuncShape == 1) {

      if( !plusMinusDependent) {

        // check if point is within the particle's support zone
        for(int i = 0;i < usedDims;i++) {

          if(fabs(pointCoords[i] - coords[i]) >= influenceRadii[i]) supported =
            false;

#ifdef _supportDebugMode_
          logFile<<pointCoords[i]<<"-"<<coords[i]<<" <=> "<<influenceRadii[i]<<endl;
#endif

        }

      }

      else {

        double dist;

        // choose the radii for the positive or negative coordinate direction
        // according to the position of the point within the support-zone

        for(int i = 0;i < usedDims;i++) {

          dist = pointCoords[i] - coords[i];

          if(dist > 0 && fabs(dist) >= influenceRadii[i]) supported = false;

          else if(dist < 0 && fabs(dist) >= influenceRadii[i + usedDims]) supported =
            false;
        }
      }

    }
    // -----------------------------------------------------------------
    // spherical window function
    else {

      ptcleDist = calcDistance(pointCoords,coords);

      if(ptcleDist >= influenceRadii[0]) supported = false;

    }

  }

  return supported;
}

/**********************************************************************/
/**********************************************************************/
// Clear all class arrays which are not locally needed.
void
Particle::clearArrays() {

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

