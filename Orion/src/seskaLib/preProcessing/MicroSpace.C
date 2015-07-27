// Manipulate all properties of a microspace.

#include "MicroSpace.h"

// retrieve the micro-length scale parameter
dbMatrix& MicroSpace::getMicroLengths(InputFileData* InputData,
                                      std::ofstream& logFile) {
  
  using namespace std;

  double value1 = InputData->getMatValue(materialID,"microLength1");
  double value2 = InputData->getMatValue(materialID,"microLength2");
  double value3 = InputData->getMatValue(materialID,"microLength3");

  microLengths.resize(0);

  if(value1 != 0) {
    microLengths.resize(microLengths.size() + 1);
    microLengths[microLengths.size() - 1].resize(2);
    microLengths[microLengths.size() - 1][0] = 0;
    microLengths[microLengths.size() - 1][1] = value1;
  }
  if(value2 != 0) {
    microLengths.resize(microLengths.size() + 1);
    microLengths[microLengths.size() - 1].resize(2);
    microLengths[microLengths.size() - 1][0] = 1;
    microLengths[microLengths.size() - 1][1] = value2;
  }
  if(value3 != 0) {
    microLengths.resize(microLengths.size() + 1);
    microLengths[microLengths.size() - 1].resize(2);
    microLengths[microLengths.size() - 1][0] = 2;
    microLengths[microLengths.size() - 1][1] = value3;
  }

  return microLengths;
}

/************************************************************************/
/************************************************************************/
// retrieve the micro-rotation vector
dbVector& MicroSpace::getMicroRotationVec(InputFileData* InputData,
                                          std::ofstream& logFile) {

  using namespace std;

  microRotVec.resize(3);

  microRotVec[0] = InputData->getMatValue(materialID,"microRotation1");
  microRotVec[1] = InputData->getMatValue(materialID,"microRotation2");
  microRotVec[2] = InputData->getMatValue(materialID,"microRotation3");

  return microRotVec;

}

/************************************************************************/
/************************************************************************/
// retrieve the volume ratio of each microSpace
void MicroSpace::setMicroVolRatio(InputFileData* InputData,
                                  std::ofstream& logFile) {

  using namespace std;

  microVolRatio = InputData->getMatValue(materialID,"microVolumeRatio");

  if(microVolRatio == 0) logFile << "Warning: material ID " << materialID
      << " has a zero volume " << "fraction." << endl;

}

/************************************************************************/
/************************************************************************/
// set initial micro-directions (Anisotropy)
// where matID is the global (macro) material ID, i.e. of superimposed 
// layers
void MicroSpace::setInitialMicroDirections(
    InputFileData* InputData,dbVector& coords,dbMatrix& macroDirectors,
    int& macroMatID,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {
  
  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int anisotropyType = (int) InputData->getMatValue(materialID,
                                                    "anisotropyType");

  int type;
  dbVector dummyDirector;

  switch(anisotropyType) {

  // isotropy
  case 0:
    break;

    // simple transverse isotropy: i = e1 (Cartesian basis vector)
  case 1:

    allocateArray(initialMicroDirectors,1,usedDims);
    initialMicroDirectors[0][0] = 1.0;
    break;

    // simple transverse isotropy: i = e2 (Cartesian basis vector)
  case 2:

    allocateArray(initialMicroDirectors,1,usedDims);
    initialMicroDirectors[0][1] = 1.0;
    break;

    // simple transverse isotropy: i = e3 (Cartesian basis vector)
  case 3:

    allocateArray(initialMicroDirectors,1,usedDims);
    initialMicroDirectors[0][2] = 1.0;
    break;

    // simple transverse isotropy: i = Qe1 (rotated Cartesian basis vector)
  case 4:

    allocateArray(initialMicroDirectors,1,usedDims);
    allocateArray(dummyDirector,usedDims);
    dummyDirector[0] = 1.0;

    getMicroRotationVec(InputData,logFile);
    calcRotationMatrix(microRotVec,microMetrics,logFile);
    innerTensorProduct(microMetrics,dummyDirector,initialMicroDirectors[0],
                       false,logFile);
    break;

    // anisotropy within myocardial models (single preferred direction)
  case 5:

    allocateArray(initialMicroDirectors,1,usedDims);

    type = (int) InputData->getValue("cardiacMechanicsProblem");

    //colli rotating function 1994
    //in a block with thickness in "z" direction increasing from zero
    //thickness and alpha0 constants in input file
    if(type == 1) {

      double thickness = InputData->getValue("thickness");
      modelData["thickness"] = thickness;

      double alpha0 = InputData->getValue("alpha0");
      modelData["alpha0"] = alpha0;

      double pi = getPi();
      double x3 = coords[2];
      double gamma = (pi / 4)
        * ( -sin(4 * pi * pow((x3 / thickness - 0.5),3)) + 1) + alpha0;

      //fibre directions
      initialMicroDirectors[0][0] = cos(gamma);
      initialMicroDirectors[0][1] = sin(gamma);
      initialMicroDirectors[0][2] = 0;

    }
    //circumferential fibre direction in an annulus
    //"z" direction normal to annulus
    else if(type == 2) {

      double pi = getPi();
      double r = sqrt(pow(coords[0],2) + pow(coords[1],2));

      //fibre directions
      initialMicroDirectors[0][0] = -coords[1] / r;
      initialMicroDirectors[0][1] = coords[0] / r;
      initialMicroDirectors[0][2] = 0;

    }

    //ellipsoidal fibre directions Rijckens
    //sheet tangent is in longitudinal direction
    else if(type == 3 || type == 4) {

      double muEpi = InputData->getValue("muEpi");
      modelData["muEpi"] = muEpi;
      double muEndo = InputData->getValue("muEndo");
      modelData["muEndo"] = muEndo;
      modelData["prolateC"] = InputData->getValue("prolateC");

      double p_1 = InputData->getValue("p1");
      modelData["p1"] = p_1;
      double p_2 = InputData->getValue("p2");
      modelData["p2"] = p_2;
      double p_3 = InputData->getValue("p3");
      modelData["p3"] = p_3;

      double pi = getPi();

      dbVector proSpheroidCoords;
      cartToProSpheroid(coords,proSpheroidCoords,modelData,logFile);

      double v = (proSpheroidCoords[0] - (muEpi + muEndo) / 2) * 2
        / (muEpi - muEndo);
      double u = -(proSpheroidCoords[1] - (pi / 2)) * 2 / pi;

      double alpha_h = p_1 + p_2 * v;
      double alpha_t = p_3 * u * (1 - pow(v,2));

      dbMatrix unitVectors;
      getProSpheroidVecs(proSpheroidCoords,unitVectors,modelData);

      // fibre direction
      initialMicroDirectors[0][0] = tan(alpha_t) * unitVectors[0][0]
        + tan(alpha_h) * unitVectors[1][0] + unitVectors[2][0];
      initialMicroDirectors[0][1] = tan(alpha_t) * unitVectors[0][1]
        + tan(alpha_h) * unitVectors[1][1] + unitVectors[2][1];
      initialMicroDirectors[0][2] = tan(alpha_t) * unitVectors[0][2]
        + tan(alpha_h) * unitVectors[1][2] + unitVectors[2][2];

      // normalise
      double norm = computeNorm(initialMicroDirectors[0],2,logFile);
      initialMicroDirectors[0][0] /= norm;
      initialMicroDirectors[0][1] /= norm;
      initialMicroDirectors[0][2] /= norm;

    }

    break;

    // anisotropy within myocardial models (multiple preferred directions)
  case 6:

    allocateArray(initialMicroDirectors,3,usedDims);

    type = (int) InputData->getValue("cardiacMechanicsProblem");

    // parallel to Cartesian axes
    if(type == 0) {
      initialMicroDirectors = getKroneckerSymbol(usedDims);
    }

    //colli rotating function 1994
    //in a block with thickness in "z" direction increasing from zero
    //thickness and alpha0 constants in input file
    else if(type == 1) {

      double thickness = InputData->getValue("thickness");
      modelData["thickness"] = thickness;

      double alpha0 = InputData->getValue("alpha0");
      modelData["alpha0"] = alpha0;

      double pi = getPi();
      double x3 = coords[2];
      double gamma = (pi / 4)
        * ( -sin(4 * pi * pow((x3 / thickness - 0.5),3)) + 1) + alpha0;

      //fibre directions
      initialMicroDirectors[0][0] = cos(gamma);
      initialMicroDirectors[0][1] = sin(gamma);
      initialMicroDirectors[0][2] = 0;

      //sheet tangent directions
      initialMicroDirectors[1][0] = -sin(gamma);
      initialMicroDirectors[1][1] = cos(gamma);
      initialMicroDirectors[1][2] = 0;

      //sheet normal directions
      initialMicroDirectors[2][0] = 0;
      initialMicroDirectors[2][1] = 0;
      initialMicroDirectors[2][2] = 1;

    }
    //circumferential fibre direction in an annulus
    //"z" direction normal to annulus
    else if(type == 2) {

      double pi = getPi();
      double r = sqrt(pow(coords[0],2) + pow(coords[1],2));

      //fibre directions
      initialMicroDirectors[0][0] = -coords[1] / r;
      initialMicroDirectors[0][1] = coords[0] / r;
      initialMicroDirectors[0][2] = 0;

      //sheet tangent directions
      initialMicroDirectors[1][0] = coords[0] / r;
      initialMicroDirectors[1][1] = coords[1] / r;
      initialMicroDirectors[1][2] = 0;

      //sheet normal directions
      initialMicroDirectors[2][0] = 0;
      initialMicroDirectors[2][1] = 0;
      initialMicroDirectors[2][2] = 1;

    }

    //ellipsoidal fibre directions Rijckens
    //sheet tangent is in longitudinal direction
    else if(type == 3 || type == 4) {

      double muEpi = InputData->getValue("muEpi");
      modelData["muEpi"] = muEpi;
      double muEndo = InputData->getValue("muEndo");
      modelData["muEndo"] = muEndo;
      modelData["prolateC"] = InputData->getValue("prolateC");

      double p_1 = InputData->getValue("p1");
      modelData["p1"] = p_1;
      double p_2 = InputData->getValue("p2");
      modelData["p2"] = p_2;
      double p_3 = InputData->getValue("p3");
      modelData["p3"] = p_3;

      double pi = getPi();

      dbVector proSpheroidCoords;
      cartToProSpheroid(coords,proSpheroidCoords,modelData,logFile);

      double v = (proSpheroidCoords[0] - (muEpi + muEndo) / 2) * 2
        / (muEpi - muEndo);
      double u = -(proSpheroidCoords[1] - (pi / 2)) * 2 / pi;

      double alpha_h = p_1 + p_2 * v;
      double alpha_t = p_3 * u * (1 - pow(v,2));

      dbMatrix unitVectors;
      getProSpheroidVecs(proSpheroidCoords,unitVectors,modelData);

      //fibre directions
      initialMicroDirectors[0][0] = tan(alpha_t) * unitVectors[0][0]
        + tan(alpha_h) * unitVectors[1][0] + unitVectors[2][0];
      initialMicroDirectors[0][1] = tan(alpha_t) * unitVectors[0][1]
        + tan(alpha_h) * unitVectors[1][1] + unitVectors[2][1];
      initialMicroDirectors[0][2] = tan(alpha_t) * unitVectors[0][2]
        + tan(alpha_h) * unitVectors[1][2] + unitVectors[2][2];

      //normalise
      double norm = computeNorm(initialMicroDirectors[0],2,logFile);
      initialMicroDirectors[0][0] /= norm;
      initialMicroDirectors[0][1] /= norm;
      initialMicroDirectors[0][2] /= norm;
      //sheet normal directions

      dbVector sheetNormal;
      allocateArray(sheetNormal,usedDims);
      crossProduct(initialMicroDirectors[0],unitVectors[1],sheetNormal);
      initialMicroDirectors[2][0] = sheetNormal[0];
      initialMicroDirectors[2][1] = sheetNormal[1];
      initialMicroDirectors[2][2] = sheetNormal[2];

      //sheet tangent directions

      dbVector sheetTangent;
      allocateArray(sheetTangent,usedDims);
      crossProduct(initialMicroDirectors[2],initialMicroDirectors[0],
                   sheetTangent);
      initialMicroDirectors[1][0] = sheetTangent[0];
      initialMicroDirectors[1][1] = sheetTangent[1];
      initialMicroDirectors[1][2] = sheetTangent[2];

    }

    break;

    // micro-directors follow cylindrical macro tangent space
  case 7:
    break;

  case 8:
    // Testing Anisotropy type for LV; only fibre direction nb
  {
    allocateArray(initialMicroDirectors,3,usedDims);

    double muEpi = InputData->getValue("muEpi");
    modelData["muEpi"] = muEpi;
    double muEndo = InputData->getValue("muEndo");
    modelData["muEndo"] = muEndo;
    modelData["prolateC"] = InputData->getValue("prolateC");

    double p_1 = InputData->getValue("p1");
    modelData["p1"] = p_1;
    double p_2 = InputData->getValue("p2");
    modelData["p2"] = p_2;
    double p_3 = InputData->getValue("p3");
    modelData["p3"] = p_3;

    double pi = getPi();

    dbVector proSpheroidCoords;
    cart2ProSpheroid(coords,proSpheroidCoords,modelData,logFile);

    //double v = (proSpheroidCoords[0] -(muEpi+muEndo)/2)*2/(muEpi-muEndo);
    double v = 1.6160 * (proSpheroidCoords[0] - muEndo) / (muEpi - muEndo);
    double u = -(proSpheroidCoords[1] - (pi / 2)) * 2 / pi;

    double alpha_h = (110 * cos(v) - 56) * pi / 180;
    double alpha_t = p_3 * u * (1 - pow(v,2));

    dbMatrix unitVectors;
    getProSpheroidVecs(proSpheroidCoords,unitVectors,modelData);

    //fibre directions
    initialMicroDirectors[0][0] = tan(alpha_t) * unitVectors[0][0]
      + tan(alpha_h) * unitVectors[1][0] + unitVectors[2][0];

    initialMicroDirectors[0][1] = tan(alpha_t) * unitVectors[0][1]
      + tan(alpha_h) * unitVectors[1][1] + unitVectors[2][1];

    initialMicroDirectors[0][2] = tan(alpha_t) * unitVectors[0][2]
      + tan(alpha_h) * unitVectors[1][2] + unitVectors[2][2];

    //normalise
    double norm0 = computeNorm(initialMicroDirectors[0],2,logFile);
    initialMicroDirectors[0][0] /= norm0;
    initialMicroDirectors[0][1] /= norm0;
    initialMicroDirectors[0][2] /= norm0;

    //sheet directions

    dbVector sheet;
    allocateArray(sheet,usedDims);
    crossProduct(initialMicroDirectors[0],unitVectors[0],sheet);
    initialMicroDirectors[1][0] = sheet[0];
    initialMicroDirectors[1][1] = sheet[1];
    initialMicroDirectors[1][2] = sheet[2];

    double norm1 = computeNorm(initialMicroDirectors[1],2,logFile);
    initialMicroDirectors[1][0] /= norm1;
    initialMicroDirectors[1][1] /= norm1;
    initialMicroDirectors[1][2] /= norm1;

    //normal directions

    dbVector normal;
    allocateArray(normal,usedDims);
    crossProduct(initialMicroDirectors[1],initialMicroDirectors[0],normal);
    initialMicroDirectors[2][0] = normal[0];
    initialMicroDirectors[2][1] = normal[1];
    initialMicroDirectors[2][2] = normal[2];

    double norm2 = computeNorm(initialMicroDirectors[2],2,logFile);
    initialMicroDirectors[2][0] /= norm2;
    initialMicroDirectors[2][1] /= norm2;
    initialMicroDirectors[2][2] /= norm2;

    break;
  }

    // anisotropy directions read in from file at specified points and 
    // then approximated at particles and Gauss points using MLS
  case 9:

    initialMicroDirectors = macroDirectors;
    break;

  case 10:
    // mean fibre direction with slight variations
    // made specifically for a 3x3 cube to replicate the experiments
    // of Dokos et al 2002.
  {

    int meanDirection = (int) InputData->getMatValue(materialID,
                                                     "meanDirection");
    allocateArray(initialMicroDirectors,1,usedDims);

    switch(meanDirection) {

    case 1: {
      //fibre directions
      initialMicroDirectors[0][0] = 1.0;

      initialMicroDirectors[0][1] = 0.0;

      initialMicroDirectors[0][2] = tan(15) * (1 - coords[1] / 1.5);
      break;
    }

    case 2: {
      //fibre directions
      initialMicroDirectors[0][0] = tan(15) * (1 - coords[2] / 1.5);

      initialMicroDirectors[0][1] = 1.0;

      initialMicroDirectors[0][2] = 0.0;
      break;
    }

    case 3: {
      //fibre directions
      logFile << "In anisotropyType 9, meanDirection 3" << endl;
      initialMicroDirectors[0][0] = tan(15) * (1 - coords[1] / 1.5);

      initialMicroDirectors[0][1] = 0.0;

      initialMicroDirectors[0][2] = 1.0;
      break;
    }

    default:
      logFile
          << "In MicroSpace::setMicroDirections anisotropyType 9, meanDirection"
          << meanDirection << " is not supported!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;

    }

    //normalise
    double norm0 = computeNorm(initialMicroDirectors[0],2,logFile);
    initialMicroDirectors[0][0] /= norm0;
    initialMicroDirectors[0][1] /= norm0;
    initialMicroDirectors[0][2] /= norm0;

    break;
  }

  default:
    logFile << "In MicroSpace::setMicroDirections anisotropy type "
        << anisotropyType << " is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/************************************************************************/
/************************************************************************/
// set the structural tensor (Anisotropy)
void MicroSpace::setStructuralTensors(InputFileData* InputData,
                                      std::map<std::string,double>& calcData,
                                      std::map<std::string,double>& modelData,
                                      std::ofstream& logFile) {
  
  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];

  if(structuralTensors.size() != initialMicroDirectors.size())

  allocateArray(structuralTensors,initialMicroDirectors.size(),usedDims,
                usedDims);

  // loop over all directions and set the corresponding structural tensors
  for(int dir = 0;dir < structuralTensors.size();dir++)

    for(int i = 0;i < structuralTensors[dir].size();i++)

      for(int j = 0;j < structuralTensors[dir].size();j++)

        structuralTensors[dir][i][j] = initialMicroDirectors[dir][i]
          * initialMicroDirectors[dir][j];

}

/************************************************************************/
/************************************************************************/
// set the micrometrics
void MicroSpace::setMicroMetrics(InputFileData* InputData,dbVector& coords,
                                 std::map<std::string,double>& calcData,
                                 std::map<std::string,double>& modelData,
                                 std::ofstream& logFile) {
  
  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];

  if((int) InputData->getMatValue(materialID,"anisotropyType") == 1
    || (int) InputData->getMatValue(materialID,"anisotropyType") == 2
    || (int) InputData->getMatValue(materialID,"anisotropyType") == 3) {
    
    InputData->setValue("variableMicrospace",0.0);
    
  }
  // The micro-continuum is rotated by multiplying the (Cartesian) basis 
  // vectors with a constant rotation matrix which is function of a 
  // corresponding rotation vector (Rodriugues formula!).
  // 
  // micro-directors are rotated by a constant rotation matrix throughout 
  // the macrospace
  else if((int) InputData->getMatValue(materialID,"anisotropyType") == 4) {
    
    getMicroRotationVec(InputData,logFile);
    calcRotationMatrix(microRotVec,microMetrics,logFile);

    InputData->setValue("variableMicrospace",1.0);
    
#ifdef _constitutiveDebugMode_
    for(int k=0;k<microRotVec.size();k++)
    logFile<<"mRotVec["<<k<<"] = "<<microRotVec[k]<<" ";
    logFile<<endl;
    for(int k=0;k<microMetrics.size();k++)
    for(int l=0;l<microMetrics[k].size();l++)
    logFile<<"mMetrics["<<k<<"]["<<l<<"]="<<microMetrics[k][l]<<endl;
#endif
    
  }

  // micro-directors follow cylindrical macro tangent space
  //
  // microDirection1 = micro-director parallel to longitudinal axis
  // microDirection2 = micro-director parallel to radial direction
  // microDirection3 = micro-director parallel to circumferential direction
  
  else if((int) InputData->getMatValue(materialID,"anisotropyType") == 7) {
    
    // set micro rotation vector  
    dbVector& microRotVec = getMicroRotationVec(InputData,logFile);

    // set micro metrics and micro metrics derivatives
    setCylinderMetrics(microRotVec,coords,microMetrics,microMetricsDerivs,
                       logFile);

    // convert the given cylindrical axes IDs to Cartesian
    intVector idx(3);
    
    // cylinder axis = X1 => tangent = X3
    if(microRotVec[0] != 0 && microRotVec[1] == 0 && microRotVec[2] == 0) {
      
      idx[0] = 0;
      idx[1] = 1;
      idx[2] = 2;
      
    }
    
    // cylinder axis = X2 => tangent = X1
    else if(microRotVec[0] == 0 && microRotVec[1] != 0 && microRotVec[2] == 0) {
      
      idx[0] = 1;
      idx[1] = 2;
      idx[2] = 0;
      
    }
    
    // cylinder axis = X3 => tangent = X2
    else if(microRotVec[0] == 0 && microRotVec[1] == 0 && microRotVec[2] != 0) {
      
      idx[0] = 2;
      idx[1] = 0;
      idx[2] = 1;
      
    }
    
    else {

      cerr << "In MicroSpace::setMicroMetrics chosen cylinder axis "
          << "is not supported!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    dbMatrix oldMicroLengths = microLengths;
    
    for(int k = 0;k < microLengths.size();k++) {
      microLengths[k][0] = idx[(int) oldMicroLengths[k][0]];
      microLengths[k][1] = oldMicroLengths[k][1];
    }

#ifdef _constitutiveDebugMode_
    for(int k=0;k<microRotVec.size();k++)
    logFile<<"mRotVec["<<k<<"] = "<<microRotVec[k]<<" ";
    logFile<<endl;
    logFile<<"-------"<<endl;
    for(int k=0;k<microMetrics.size();k++)
    for(int l=0;l<microMetrics[k].size();l++)
    logFile<<"mMetrics["<<k<<"]["<<l<<"]="<<microMetrics[k][l]<<endl;
    logFile<<"-------"<<endl;
    for(int i=0;i<microMetricsDerivs.size();i++)
    for(int j=0;j<microMetricsDerivs[i].size();j++)
    for(int k=0;k<microMetricsDerivs[i][j].size();k++)
    logFile<<"mMetricsDeriv["<<i<<"]["<<j<<"]["<<k<<"]="
    <<microMetricsDerivs[i][j][k]<<endl;
    logFile<<"-------"<<endl;
    logFile<<"re-defined microLengths:"<<endl;
    for(int k=0;k<microLengths.size();k++)
    logFile<<"dof "<<microLengths[k][0]
    <<" value "<<microLengths[k][1]<<endl;
#endif
    
  }
  
}

/************************************************************************/
/************************************************************************/
// set the material properties
void MicroSpace::setMaterialProperties(InputFileData* InputData,
                                       std::map<std::string,double>& calcData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {
  
  using namespace std;
  
  switch(constitutiveLawID) {

  // Saint-Venant-Kirchhoff
  case 1:

    microDependentMaterial = false;
    break;

    // Saint-Venant-Kirchhoff
  case 2:

    microDependentMaterial = false;
    break;

    // Saint-Venant-Kirchhoff
  case 3:

    microDependentMaterial = false;
    break;

    // Arruda and Boyce hyperelasticity
  case 4:

    microDependentMaterial = true;
    break;
    
    // Neo-Hookean hyperelasticity
  case 5:

    microDependentMaterial = true;
    break;

    // right Cauchy-Green deformation tensor-based viscoplasticity 
  case 7:

    microDependentMaterial = true;
    break;
    
    // electro-active polymer
  case 8:

    if((int) modelData["mechanicalConstitutiveLaw"] == 1)

    microDependentMaterial = false;
    
    else

    microDependentMaterial = true;

    break;

    // nonlinear hypoelasticity based on power law
  case 9:

    microDependentMaterial = true;
    break;

    // Mises plasticity based on stretch tensor
  case 16:

    microDependentMaterial = true;
    break;
    
    // Linear cosseratFibre based
  case 17:

    microDependentMaterial = false;
    break;

    // Neohookean hyperelasticity based on stretch tensor
  case 18:

    microDependentMaterial = true;
    break;

    // Power law hypoelasticity based on stretch tensor
  case 19:

    microDependentMaterial = true;
    break;

    // Nonlinear cosseratFibre
  case 21:

    microDependentMaterial = false;
    break;

    // stretch cardiac mechanics
  case 22:

    microDependentMaterial = false;
    break;

    // invariant based cosseratfibre
  case 23:

    microDependentMaterial = false;
    break;

    // Cosserat just the fibre
  case 30:

    microDependentMaterial = false;
    break;

  default:
    logFile << "In MicroSpace::setMaterialProperties chosen\n"
        << "constitutive law '" << constitutiveLawID << "' is "
        << "not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

#ifdef _constitutiveDebugMode_
  logFile<<"microDependentMaterial = "<<microDependentMaterial<<endl;
#endif

}
