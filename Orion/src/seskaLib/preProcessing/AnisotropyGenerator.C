// Project the preferred a material direction onto two opposite planes
// of a structure, each with a user-specified normal vector 'z' and
// rotated by a used-specified angle
#include "AnisotropyGenerator.h"

AnisotropyGenerator::AnisotropyGenerator(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,
    std::map<std::string,double>& constitutiveData,std::ofstream& logFile) {

  using namespace std;

  map<string,vector<AnisotropyCondition> >& anisotropyConditions =
    InputData->getAnisotropyConditions();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef _anisotropyInterpolationDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"############ Compute Surface Anisotropy #############"<<endl;
  logFile<<"#####################################################"<<endl;
#endif

  string conditionName = "Surface-Material-Direction-Angle-Condition";
  pi = getPi();

  if(rank == 0) cout << "generating surface preferred material directions"
      << endl;
  logFile << "generating surface preferred material directions" << endl;


  initPreSetup(MeshlessData,InputData,calcData,modelData,logFile);


  /**************************************************************************/
  // loop over all anisotropy conditions
  for(int c = 0;c < anisotropyConditions[conditionName].size();c++) {
    AnisotropyCondition& condition = anisotropyConditions[conditionName][c];

#ifdef _anisotropyInterpolationDebugMode_
    logFile<<"**************************************************"<< endl;
    logFile<<condition.getType()<<" "<<condition.getID()<<endl;
#endif

    // Calculate all surface normals
    calcAllSurfaceNormals(MeshlessData,InputData,condition,calcData,modelData,
                          logFile);

    // Compute the averaged surface normal at each node
    calcNodalNormal(MeshlessData,InputData,condition,calcData,modelData,
                    logFile);

    // Calculate the circumferential direction at each node
    calcNodalCircumDirections(MeshlessData,InputData,condition,calcData,
                              modelData,logFile);

    // Calculate the outward pointing normal at each node
    calcOutwardNormal(MeshlessData,InputData,condition,calcData,modelData,
                      logFile);

    // Calculate the sheet normal at each node
    calcSheetNormal(MeshlessData,InputData,condition,calcData,modelData,
                    logFile);

    // Calculate the fibre direction projection vector at each node
    calcFibreDirectProjection(MeshlessData,InputData,condition,calcData,
                              modelData,logFile);

    // Calculate the fibre direction vector at each node
    calcFibreDirection(MeshlessData,InputData,condition,calcData,modelData,
                       logFile);

    // Calculate the fibre direction vector at each node
    calcOrthogonalFibreDirection(MeshlessData,InputData,condition,calcData,
                                 modelData,logFile);

    // store material directions for for SESKA
    storeMaterialDirections(MeshlessData,InputData,condition,
                            condPtcleRootLists[c],calcData,modelData,
                            constitutiveData,logFile);

    // Generate result for GiD
    if(rank == 0) saveResultsToFile_flexible(MeshlessData,InputData,condition,
                                             seskaCondPtcleIDLists[c],calcData,
                                             modelData,logFile);

  }

  femResFile.close();

  if(rank == 0) cout
      << "finished generating surface preferred material directions" << endl;
  logFile << "finished generating surface preferred material directions"
      << endl;

}

/******************************************************************************/
/******************************************************************************/
// pre-setup
void AnisotropyGenerator::initPreSetup(MeshlessApproximation* MeshlessData,
                                       InputFileData* InputData,
                                       std::map<std::string,double>& calcData,
                                       std::map<std::string,double>& modelData,
                                       std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  intVector& ptcleRootList = MeshlessData->getPtcleRootList();
  vector<Particle>& ptcls = MeshlessData->getParticlesVec();
  intVector& newIdx = MeshlessData->getNewIdx();
  std::map<std::string,std::vector<AnisotropyCondition> >& anisotropyConditions =
    InputData->getAnisotropyConditions();
  int nodesPerSurfaceElem =
    InputData->getBackGroundMeshInfo()["nodesPerSurfaceElement"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank == 0) system("cp fem.msh fem_fibre.msh");

  int ID,surfaceID,elemCount,ptcleCount;
  int currentElem,condPtcle,ptcle,ptcleNum;
  string conditionName = "Surface-Material-Direction-Angle-Condition";
  double value;

  int globalElemNum,localElemNum,number;

  resizeArray(seskaCondPtcleIDLists,anisotropyConditions[conditionName].size(),
              ptcls.size());
  resizeArray(condPtcleRootLists,anisotropyConditions[conditionName].size());

  // loop over all anisotropy conditions and store particle coordinates
  // anisotropy angles
  for(int c = 0;c < anisotropyConditions[conditionName].size();c++) {
    AnisotropyCondition& condition = anisotropyConditions[conditionName][c];

    vector<ConditionElement>& condElems = condition.getElements();
    vector<ConditionParticle>& condPtcls = condition.getParticles();
    dbVector& condValues = condition.getConditionValues();
    dbVector& rotAngles = condition.getRotationAngle();

    // -----------------------------------------------------------------
    // determine a local portion of surface elements
    globalElemNum = condElems.size();
    number = (int) ceil((double) globalElemNum / size);

    if(number * (rank + 1) <= globalElemNum) {
      localElemNum = number;
    }
    else if(number * rank >= globalElemNum) {
      localElemNum = 0;
      logFile << "In AnisotropyGenerator::initPreSetup " << localElemNum
          << " local volume elements!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {
      localElemNum = globalElemNum - number * rank;
    }

    // -----------------------------------------------------------------
    elemCount = 0;
    ptcleCount = 0;
    seskaCondPtcleIDLists[c].assign(seskaCondPtcleIDLists[c].size(), -1);
    vector<ConditionElement> localElems(localElemNum);
    condPtcls = vector<ConditionParticle>(globalElemNum * nodesPerSurfaceElem,
                                          ConditionParticle());

    intVector elemPtcls;
    dbMatrix coordMat;

    // loop over all elements of condition
    for(int i = 0;i < condElems.size();i++) {
      ConditionElement& globalElem = condElems[i];
      intVector& elemNodes = globalElem.getNodes();

      if(elemPtcls.size() < elemNodes.size()) {
        resizeArray(elemPtcls,elemNodes.size());
        resizeArray(coordMat,elemNodes.size());
      }
      clearArray(elemPtcls);
      clearArray(coordMat);

      // loop over all nodes of current element and store the nodes, coordinates and ID
      for(int j = 0;j < elemNodes.size();j++) {
        ptcle = newIdx[elemNodes[j] - 1]; // SESKA ordering

        // not yet set
        if(seskaCondPtcleIDLists[c][ptcle] == -1) {
          seskaCondPtcleIDLists[c][ptcle] = ptcleCount;
          condPtcls[ptcleCount].getID() = ptcle;
          condPtcls[ptcleCount].getCoords() = ptcls[ptcle].getCoords();
          condPtcls[ptcleCount].getMaterialDirectionAngle() =
            condValues[condElems[i].getSurfaceID()];
          condPtcls[ptcleCount].getMaterialRotationAngle() =
            rotAngles[condElems[i].getSurfaceID()];

          // broadcast coordinates
          if(condPtcls[ptcleCount].getCoords().size() == 0) resizeArray(
              condPtcls[ptcleCount].getCoords(),usedDims);
          MPI_Bcast( &condPtcls[ptcleCount].getCoords()[0],usedDims,MPI_DOUBLE,
                    ptcleRootList[ptcle],MPI_COMM_WORLD);
          ptcleCount++;
        }
        elemPtcls[j] = seskaCondPtcleIDLists[c][ptcle]; // condition ordering
        coordMat[j] = condPtcls[seskaCondPtcleIDLists[c][ptcle]].getCoords();
      }

      // store only the local elements
      if(i >= rank * number && i < rank * number + localElemNum) {
        localElems[elemCount].getNodalCoords() = coordMat;
        localElems[elemCount].getNodes() = elemPtcls;

        // transfer global element properties to local one
        localElems[elemCount].getID() = condElems[i].getID();
        localElems[elemCount].getSurfaceID() = condElems[i].getSurfaceID();
        localElems[elemCount].getSurfaceTangentDirection() =
          condElems[i].getSurfaceTangentDirection();
        localElems[elemCount].getConditionValues() =
          condElems[i].getConditionValues();

        elemCount++;
      }

    }

    condElems = localElems;
    resizeArray(condPtcls,ptcleCount);

    // set condition particle root list
    resizeArray(condPtcleRootLists[c],condPtcls.size());

    for(int i = 0;i < condElems.size();i++) {
      ConditionElement& elem = condElems[i];
      intVector& elemNodes = elem.getNodes();
      for(int j = 0;j < elemNodes.size();j++)
        condPtcleRootLists[c][elemNodes[j]] = rank;
    }
    for(int i = 0;i < condPtcleRootLists[c].size();i++) {
      int localValue = condPtcleRootLists[c][i];
      MPI_Allreduce( &localValue, &condPtcleRootLists[c][i],1,MPI_INT,MPI_MAX,
                    MPI_COMM_WORLD);
    }

#ifdef _anisotropyInterpolationDebugMode_
    logFile << "**************************************************" << endl;
    logFile << "***** local anisotropy surface elements **********" << endl;
    for(int i = 0;i < condElems.size();i++) {
      ConditionElement& elem = condElems[i];
      intVector& nodes = elem.getNodes();
      dbMatrix& nodalCoords = elem.getNodalCoords();
      logFile << "-------------------------------" << endl;
      logFile << i + 1 << ") ELEM " << elem.getID() << ": angle="
      <<elem.getConditionValues()[0] <<" nodalCoords:"<<endl;
      for(int j = 0;j < nodalCoords.size();j++) {
        logFile<< " ptcle " << nodes[j] << ": ";
        for(int k = 0;k < nodalCoords[j].size();k++) {
          logFile << nodalCoords[j][k] << " ";
        }
        logFile<<endl;
      }
      for(int j = 0;j < nodes.size();j++) {
        logFile << j + 1 << " ptcle " << nodes[j] << " ("
        << condPtcls[nodes[j]].getID() << ") coords=";
        dbVector& coords = condPtcls[nodes[j]].getCoords();
        for(int k = 0;k < coords.size();k++) {
          logFile << coords[k] << " ";
        }
        logFile << endl;
      }
    }
    logFile << "**************************************************" << endl;
    logFile << "***** local anisotropy surface particles *********" << endl;
    for(int i = 0;i < condPtcls.size();i++) {
      dbVector& coords = condPtcls[i].getCoords();
      logFile << i + 1 << ") PTCLE " << condPtcls[i].getID() << " root="
      <<condPtcleRootLists[c][i]<< ": coords=";
      for(int j = 0;j < coords.size();j++)
      logFile << coords[j] << " ";
      logFile << "angle="<<condPtcls[i].getMaterialDirectionAngle() <<endl;
    }
#endif

  }

  string fileName = "fem_fibre.res";
  femResFile.open(fileName.c_str());
  femResFile.precision(12);
  femResFile.setf(ios_base::scientific,ios_base::floatfield);

  string headerLine = "GiD Post Results File 1.0";
  femResFile << headerLine << endl;
}

/*****************************************************************************/
/*****************************************************************************/
// Calculate all surface normal
void AnisotropyGenerator::calcAllSurfaceNormals(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  vector<ConditionElement>& elems = condition.getElements();
  vector<ConditionParticle>& ptcls = condition.getParticles();

  // Compute the surface normal of all surfaces

  for(int i = 0;i < elems.size();i++) {

    elems[i].getSurfaceNormal() = calcSurfaceNormal(MeshlessData,InputData,
                                                    elems[i].getNodalCoords(),
                                                    calcData,modelData,logFile);

  }

#ifdef _anisotropyInterpolationDebugMode_
  logFile << "*****************************************************" << endl;
  logFile << "**************** all surface normals ****************" << endl;
  logFile << "*****************************************************" << endl;
  for(int i = 0;i < elems.size();i++) {
    ConditionElement& elem = elems[i];
    dbVector& normal = elems[i].getSurfaceNormal();
    intVector& nodes = elems[i].getNodes();
    dbMatrix& coordMat = elems[i].getNodalCoords();
    logFile<<"-------------------------------------------"<<endl;
    logFile << i + 1 << ") ELEMENT " << elem.getID() << ": surfaceID="
    << elem.getSurfaceID() << " normal=";
    for(int j = 0;j < normal.size();j++)
    logFile << normal[j] << " ";
    logFile << endl;
    for(int j=0;j<nodes.size();j++) {
      logFile<<"Node "<<nodes[j]<<" coords=";
      for(int k=0;k<coordMat[j].size();k++)
      logFile<<coordMat[j][k]<<" ";
      logFile<<endl;
    }
  }
  logFile << endl;
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the surface normal from a set of nodal coordinates
dbVector AnisotropyGenerator::calcSurfaceNormal(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    dbMatrix& surfacePointsCoords,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(surfacePointsCoords.size() < usedDims) {
    logFile << "In AnisotropyGenerator::calcSurfaceNormal: Too few coordinates"
        " to compute surface normal" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  dbVector vecOne(surfacePointsCoords[0].size());
  dbVector vecTwo(surfacePointsCoords[0].size());

  for(int i = 0;i < surfacePointsCoords[0].size();i++) {

    if(abs(surfacePointsCoords[1][i] - surfacePointsCoords[0][i]) < 1e-14) vecOne[i] =
      0;
    else vecOne[i] = surfacePointsCoords[1][i] - surfacePointsCoords[0][i];

    if(abs(surfacePointsCoords[2][i] - surfacePointsCoords[0][i]) < 1e-14) vecTwo[i] =
      0;
    else vecTwo[i] = surfacePointsCoords[2][i] - surfacePointsCoords[0][i];
  }

  dbVector surfaceNormalVec;
  crossProduct(vecOne,vecTwo,surfaceNormalVec);

  // Invert the direction of the surface normal
  for(int i = 0;i < surfaceNormalVec.size();i++) {
    surfaceNormalVec[i] = surfaceNormalVec[i] * -1;
  }

  normalizeVector(surfaceNormalVec);

  return surfaceNormalVec;
}

/*****************************************************************************/
/*****************************************************************************/
// Compute the averaged surface normal at each node
void AnisotropyGenerator::calcNodalNormal(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  vector<ConditionElement>& elems = condition.getElements();
  vector<ConditionParticle>& ptcls = condition.getParticles();

  for(int i = 0;i < elems.size();i++) {
    intVector& nodes = elems[i].getNodes();

    for(int j = 0;j < nodes.size();j++) {
      ptcls[nodes[j]].getElems().push_back(i);
      ptcls[nodes[j]].getSurfaceID() = elems[i].getSurfaceID();
      ptcls[nodes[j]].getSurfaceTangentDirection() =
        elems[i].getSurfaceTangentDirection();
    }

  }

  // Enforce epicardial nodes
  for(int i = 0;i < elems.size();i++) {

    if(elems[i].getSurfaceID() == 0) {
      intVector& nodes = elems[i].getNodes();

      for(int j = 0;j < nodes.size();j++) {
        ptcls[nodes[j]].getSurfaceID() = elems[i].getSurfaceID();
        ptcls[nodes[j]].getSurfaceTangentDirection() =
          elems[i].getSurfaceTangentDirection();
      }

    }

  }

  calcAveragedNodalNormals(MeshlessData,InputData,condition,calcData,modelData,
                           logFile);

#ifdef _anisotropyInterpolationDebugMode_
  logFile << "*****************************************************" << endl;
  logFile << "************ all surface particle normals ***********" << endl;
  logFile << "*****************************************************" << endl;
  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];
    intVector& ptcleElems = ptcls[i].getElems();
    logFile << i + 1 << ") PTCLE " << ptcle.getID() << ": surfaceID="
    << ptcle.getSurfaceID() << " elems=";
    for(int j = 0;j < ptcleElems.size();j++)
    logFile << ptcleElems[j] << " ";
    if(ptcle.getAllSurfaceNormals().size()==0) {
      logFile<<endl;
      continue;
    }
    dbVector& normal = ptcle.getAllSurfaceNormals()[0];
    logFile << "tangentDirection="
    << ptcle.getSurfaceTangentDirection() << " normal=";
    for(int j = 0;j < normal.size();j++)
    logFile << normal[j] << " ";
    logFile << endl;
  }
  logFile << endl;
#endif

}

/*****************************************************************************/
/*****************************************************************************/
// calculate the nodal normals by averaging the normal of its
// surrounding/neighbouring surfaces.
void AnisotropyGenerator::calcAveragedNodalNormals(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();

  for(int i = 0;i < ptcls.size();i++) {
    intVector surfElems = ptcls[i].getElems();

    if(ptcls[i].getElems().size() == 0) continue;

    eliminateSpecificSurface(InputData,condition,surfElems,int(1),logFile);

    dbVector avgNormals(usedDims,0);
    int numOfSurfaces = surfElems.size();

    // Calculate the average normal from the selected surfaces
    for(int j = 0;j < numOfSurfaces;j++) {
      dbVector& sNormal = elems[surfElems[j]].getSurfaceNormal();

      for(int k = 0;k < sNormal.size();k++) {
        avgNormals[k] += sNormal[k] / numOfSurfaces;
      }
    }

    normalizeVector(avgNormals);

    // store the averaged surface normal at each node
    ptcls[i].getAllSurfaceNormals().resize(1);
    ptcls[i].getAllSurfaceNormals()[0] = avgNormals;
  }

}

/*****************************************************************************/
/*****************************************************************************/
// calculate the nodal normals by averaging the normal of its
// surrounding/neighbouring surfaces.
void AnisotropyGenerator::calcNodalCircumDirections(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();
  dbVector& condNormal = condition.getSurfaceNormal();

  //double orientationIndexEpi = 1.0;
  //double orientationIndexEndo = -1.0;

  double value;
  dbVector rotVec;

  // For all surface nodes
  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];
    dbVector& n = ptcle.getAllSurfaceNormals()[0];

    if(ptcls[i].getElems().size() == 0) continue;

    dbVector c;
    crossProduct(condNormal,n,c);

    if(computeNorm(c,2,logFile) == 0) {
      logFile << "In AnisotropyGenerator::calcNodalCircumDirections specified "
          << "input parameter 'normalDirection' is not admissible" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // give all circumferential direction a counter-clockwise or orientation
    // clockwise orientation about condNormal
    crossProduct(n,c,rotVec);
    scalarProduct(condNormal,rotVec,value,logFile);

    if((ptcle.getSurfaceTangentDirection() == 1 && value < 0)
      || (ptcle.getSurfaceTangentDirection() == -1 && value > 0)) transform(
        c.begin(),c.end(),c.begin(),bind1st(std::multiplies<double>(), -1.0));

    // other counter-clockwise
    else if(ptcle.getSurfaceTangentDirection() != 1
      && ptcle.getSurfaceTangentDirection() != -1) {
      logFile << "In AnisotropyGenerator::calcNodalCircumDirections specified "
          << "input parameter 'tangentDirection' is not admissible" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    ptcle.getAllSurfaceNormals().resize(2);
    ptcle.getAllSurfaceNormals()[1] = c;
  }

#ifdef _anisotropyInterpolationDebugMode_
  logFile << "*****************************************************" << endl;
  logFile << "*********** all surface particle tangents ***********" << endl;
  logFile << "*****************************************************" << endl;
  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];
    if(ptcle.getAllSurfaceNormals().size()==0)
    continue;
    dbVector& normal = ptcle.getAllSurfaceNormals()[0];
    dbVector& tangent = ptcle.getAllSurfaceNormals()[1];
    logFile << i + 1 << ") PTCLE " << ptcle.getID() << ": surfaceID="
    << ptcle.getSurfaceID() << " tangentDirection="
    << ptcle.getSurfaceTangentDirection() << " normal=";
    for(int j = 0;j < normal.size();j++)
    logFile << normal[j] << " ";
    logFile << " tangent=";
    for(int j = 0;j < tangent.size();j++)
    logFile << tangent[j] << " ";
    logFile << endl;
  }
  logFile << endl;
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the outward pointing normal at each node
void AnisotropyGenerator::calcOutwardNormal(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();
  dbVector& condNormal = condition.getSurfaceNormal();

  // For all surface nodes
  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];

    if(ptcls[i].getElems().size() == 0) continue;

    dbVector n_cz;
    crossProduct(ptcle.getAllSurfaceNormals()[1],condNormal,n_cz);

    ptcle.getAllSurfaceNormals().resize(3);
    ptcle.getAllSurfaceNormals()[2] = n_cz;
  }

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the sheet normal at each node
void AnisotropyGenerator::calcSheetNormal(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();

//	logFile << "##### calcSheetNormal #####" << endl;

  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];

    if(ptcls[i].getElems().size() == 0) continue;

    dbVector& n = ptcle.getAllSurfaceNormals()[0];
    dbVector& n_cz = ptcle.getAllSurfaceNormals()[2];
    double val;

    scalarProduct(n_cz,n,val,logFile);

    double sign_val = copysign(1.0,val);

    dbVector s(n.size());
    std::transform(n.begin(),n.end(),s.begin(),
                   std::bind1st(std::multiplies<double>(),sign_val));

    normalizeVector(s);

    ptcle.getAllSurfaceNormals().resize(4);
    ptcle.getAllSurfaceNormals()[3] = s;
  }

}

/*****************************************************************************/
/*****************************************************************************/
// project the fibres
void AnisotropyGenerator::calcFibreDirectProjection(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();
  dbVector& condNormal = condition.getSurfaceNormal();

  // For all surface nodes
  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];

    if(ptcls[i].getElems().size() == 0) continue;

    dbVector& c = ptcle.getAllSurfaceNormals()[1];
    dbVector p(c.size());

    double cos_alpha = cos(ptcle.getMaterialDirectionAngle() * pi / 180.0);
    double sin_alpha = sin(ptcle.getMaterialDirectionAngle() * pi / 180.0);

    for(int j = 0;j < c.size();j++) {
      p[j] = (cos_alpha * c[j]) + (sin_alpha * condNormal[j]);
    }

    ptcle.getAllSurfaceNormals().resize(5);
    ptcle.getAllSurfaceNormals()[4] = p;
  }

}

/*****************************************************************************/
/*****************************************************************************/
// calculate the fibre direction at each node
void AnisotropyGenerator::calcFibreDirection(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();

  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];

    if(ptcls[i].getElems().size() == 0) continue;

    dbVector n_cz = ptcle.getAllSurfaceNormals()[2];
    dbVector s = ptcle.getAllSurfaceNormals()[3];
    dbVector p = ptcle.getAllSurfaceNormals()[4];

    double p_dot_s;
    scalarProduct(p,s,p_dot_s,logFile);

    double n_cz_dot_s;
    scalarProduct(n_cz,s,n_cz_dot_s,logFile);

    dbVector f(n_cz.size());
    for(int j = 0;j < f.size();j++) {
      f[j] = ( -1 * p_dot_s * n_cz[j]) + (n_cz_dot_s * p[j]);
    }

    normalizeVector(f);

    ptcle.getAllSurfaceNormals().resize(6);
    ptcle.getAllSurfaceNormals()[5] = f;
  }

#ifdef _anisotropyInterpolationDebugMode_
  logFile << "*****************************************************" << endl;
  logFile << "***** all surface particles fibre directions ********" << endl;
  logFile << "*****************************************************" << endl;
  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];
    if(ptcle.getAllSurfaceNormals().size()==0)
    continue;
    dbVector& normal = ptcle.getAllSurfaceNormals()[5];
    logFile << i + 1 << ") PTCLE " << ptcle.getID() << ": surfaceID="
    << ptcle.getSurfaceID() << " fibre=";
    for(int j = 0;j < normal.size();j++)
    logFile << normal[j] << " ";
    logFile << endl;
  }
  logFile << endl;
#endif

}

/*****************************************************************************/
/*****************************************************************************/

void AnisotropyGenerator::calcOrthogonalFibreDirection(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();

  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];

    if(ptcls[i].getElems().size() == 0) continue;

    dbVector s = ptcle.getAllSurfaceNormals()[3];
    dbVector f = ptcle.getAllSurfaceNormals()[5];

    dbVector m;
    crossProduct(s,f,m);

    normalizeVector(m);

    ptcle.getAllSurfaceNormals().resize(7);
    ptcle.getAllSurfaceNormals()[6] = m;
  }

}

/*****************************************************************************/
/*****************************************************************************/
// calculate the outward pointing normal at each node
void AnisotropyGenerator::storeMaterialDirections(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,intVector& ptcleRootList,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,
    std::map<std::string,double>& constitutiveData,
    std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int& conditionID = condition.getID();
  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();
//  double rotationAngle = condition.getRotationAngle() * pi / 180.0;
//  dbVector& rotationAxis = condition.getRotationAxis();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  dbVector rotationAxis;
  dbVector localVector;
  dbVector rotVec(usedDims);
  dbMatrix rotMat;

  // -------------------------------------------------------------------------
  // loop over all ptcls
  for(int i = 0;i < ptcls.size();i++) {

    ConditionParticle& ptcle = ptcls[i];
    dbMatrix3& matDirs = ptcle.getMaterialDirections();
    if(matDirs.size() <= conditionID)
      allocateArray(matDirs,conditionID+1);

    if(ptcle.getAllSurfaceNormals().size() < 7) resizeArray(
        ptcle.getAllSurfaceNormals(),7,usedDims);

    // broadcast local results
    for(int k = 0;k < ptcle.getAllSurfaceNormals().size();k++)
      MPI_Bcast( &ptcle.getAllSurfaceNormals()[k][0],usedDims,MPI_DOUBLE,
                ptcleRootList[i],MPI_COMM_WORLD);

    double rotationAngle = ptcle.getMaterialRotationAngle() * pi / 180.0;

    // ---------------
    // rigid-body rotation of directions
    if(rotationAngle != 0) {

      int rotationAxisType = condition.getRotationAxisType();

      switch(rotationAxisType) {
      case 0:
        // rotate about an arbitrary axis (user-specified)
        rotationAxis = condition.getRotationAxis();
        normalizeVector(rotationAxis);
        break;

      case 1:
        // rotate about the x-axis
        rotationAxis = dbVector(3,0);
        rotationAxis[0] = 1;
        rotationAxis[1] = 0;
        rotationAxis[2] = 0;
        break;

      case 2:
        // rotate about the y-axis
        rotationAxis = dbVector(3,0);
        rotationAxis[0] = 0;
        rotationAxis[1] = 1;
        rotationAxis[2] = 0;
        break;

      case 3:
        // rotate about the z-axis
        rotationAxis = dbVector(3,0);
        rotationAxis[0] = 0;
        rotationAxis[1] = 0;
        rotationAxis[2] = 1;
        break;

      case 4:
        // rotate about the fibre axis
        rotationAxis = ptcle.getAllSurfaceNormals()[5];
        break;

      case 5:
        // rotate about the sheet axis
        rotationAxis = ptcle.getAllSurfaceNormals()[3];
        break;

      case 6:
        // rotate about the sheet-normal axis
        rotationAxis = ptcle.getAllSurfaceNormals()[6];
        break;

      case 7:
        // rotate about the averaged surface normal axis
        rotationAxis = ptcle.getAllSurfaceNormals()[0];
        break;

      case 8:
        // rotate about the circumferential axis
        rotationAxis = ptcle.getAllSurfaceNormals()[1];
        break;

      default:
        cout
            << "In AnisotropyGenerator::storeMaterialDirections, rotationAxisType is invalid."
            << endl;
        logFile
            << "In AnisotropyGenerator::storeMaterialDirections, rotationAxisType is invalid."
            << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
        break;

      }

      dbVector& f = ptcle.getAllSurfaceNormals()[5];
      dbVector& s = ptcle.getAllSurfaceNormals()[3];
      dbVector& m = ptcle.getAllSurfaceNormals()[6];

      dbVector f_rot,s_rot,m_rot;

      rotateVecAboutAxis(f,rotationAxis,rotationAngle,f_rot,logFile);
      rotateVecAboutAxis(s,rotationAxis,rotationAngle,s_rot,logFile);
      rotateVecAboutAxis(m,rotationAxis,rotationAngle,m_rot,logFile);

      ptcle.getAllSurfaceNormals()[5] = f_rot;
      ptcle.getAllSurfaceNormals()[3] = s_rot;
      ptcle.getAllSurfaceNormals()[6] = m_rot;

    }

    // ---------------
    // cardiac tissue orthotropy
    if((bool) constitutiveData["orthotropicMaterial"]) {
      matDirs[conditionID].resize(3);
      matDirs[conditionID][0] = ptcle.getAllSurfaceNormals()[5]; // Fibre direction
      matDirs[conditionID][1] = ptcle.getAllSurfaceNormals()[3]; // Sheet direction
      matDirs[conditionID][2] = ptcle.getAllSurfaceNormals()[6]; // Sheet-normal direction
    }
    // transverse isotropy
    else {
      matDirs[conditionID].resize(1);
      matDirs[conditionID][0] = ptcle.getAllSurfaceNormals()[5];
    }
  }

#ifdef _anisotropyInterpolationDebugMode_
  logFile << "**************************************************" << endl;
  logFile << "********** stored surface anisotropy *************" << endl;
  int count = 0;
  for(int i = 0;i < ptcls.size();i++) {
    ConditionParticle& ptcle = ptcls[i];
    dbVector& coords = ptcle.getCoords();
    dbMatrix& dirs = ptcle.getMaterialDirections()[conditionID];
    logFile << count + 1 << ") PTCLE " << ptcle.getID() << ": " << endl;
    for(int j = 0;j < dirs.size();j++) {
      for(int k = 0;k < dirs[j].size();k++) {
        logFile << "V[" << j << "][" << k << "]=" << dirs[j][k] << endl;
      }
      logFile<<"norm="<<computeNorm(dirs[j],2)<<endl;
    }
    count++;
    dbVector nVec(3);
    nVec[0]=dirs[0][1]*dirs[1][2]-dirs[0][2]*dirs[1][1];
    nVec[1]=dirs[0][2]*dirs[1][0]-dirs[0][0]*dirs[1][2];
    nVec[2]=dirs[0][0]*dirs[1][1]-dirs[0][1]*dirs[1][0];
    logFile<<"n="<<nVec[0]<<" "<<nVec[1]<<" "<<nVec[2]<<endl;
    crossProduct(dirs[0],dirs[1],nVec);
    logFile<<"n="<<nVec[0]<<" "<<nVec[1]<<" "<<nVec[2]<<endl;
  }
#endif
}

/*****************************************************************************/
/*****************************************************************************/
// Rotate the vector v around vector k using a prescribed angle
void AnisotropyGenerator::rotateVecAboutAxis(dbVector& v, dbVector& k, double& angle,
                                             dbVector& rVec, ofstream& logFile){

  // Using Rodrigues' rotation formula
  // (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula)
  // v = vector to be rotated, k = axis
  // v_rot = v cos(angle) + (k x v)sin(angle) + k(k.v)(1-cos(angle))

  double cos_rAngle = cos(angle);
  double sin_rAngle = sin(angle);

  dbVector k_cross_v;
  crossProduct(k,v,k_cross_v);

  double k_dot_v;
  scalarProduct(k,v,k_dot_v,logFile);

  dbVector v_cos(v.size(),0);
  dbVector k_dot_s_sin(v.size(),0);
  dbVector k_k_dot_s_cos(v.size(),0);
  dbVector v_rot(v.size(),0);
  for(int i=0; i<v.size();i++){
    v_rot[i] = v[i]*cos_rAngle + k_cross_v[i]*sin_rAngle + k[i]*k_dot_v*(1-cos_rAngle);
  }

  rVec = v_rot;

}

/*****************************************************************************/
/*****************************************************************************/
// Outputting calculated matrix to fem.res file
void AnisotropyGenerator::saveResultsToFile_flexible(
    MeshlessApproximation* MeshlessData,InputFileData* InputData,
    AnisotropyCondition& condition,intVector& seskaCondPtcleID,
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  intVector& newIdx = MeshlessData->getNewIdx();
  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();

  intVector resultsToSave(7);
  resultsToSave[0] = 0;
  resultsToSave[1] = 1;
  resultsToSave[2] = 2;
  resultsToSave[3] = 3;
  resultsToSave[4] = 4;
  resultsToSave[5] = 5;
  resultsToSave[6] = 6;

  map<int,string> resultInfo;
  resultInfo[0] = "averaged_nodal_normal_n";
  resultInfo[1] = "circumferential_direction_c";
  resultInfo[2] = "outward_pointing_normal_n_cz";
  resultInfo[3] = "sheet_normal_s";
  resultInfo[4] = "fibre_direction_projection_p";
  resultInfo[5] = "fibre_direction_f";
  resultInfo[6] = "fibre_direction_m";

  // -------------------------------------------------------------------------

  vector<string> vectorResultName(resultsToSave.size());
  vector<dbMatrix> vectorResultMatrix(resultsToSave.size(),dbMatrix());
  vector<double> stepValue(resultsToSave.size(),condition.getID());

  string cID = convertIntToString(condition.getID());

  for(int i = 0;i < resultsToSave.size();i++) {
    vectorResultName[i] = resultInfo.find(resultsToSave[i])->second;

    dbMatrix resultMat(newIdx.size(),dbVector());
    dbVector dummyVec(usedDims,0);

    // loop over all particles (mesh file numbering)
    for(int j = 0;j < newIdx.size();j++) {

      if(seskaCondPtcleID[newIdx[j]] != -1) {
        ConditionParticle& condPtcle = ptcls[seskaCondPtcleID[newIdx[j]]];
        dbVector& resultVec = condPtcle.getAllSurfaceNormals()[resultsToSave[i]];
        resultMat[j] = resultVec;
      }
      else resultMat[j] = dummyVec;
    }

    vectorResultMatrix[i] = resultMat;
  }

  string scalarResultToSave = "fibre_angle";
  vector<string> scalarResultName(1);
  scalarResultName[0] = scalarResultToSave;

  vector<dbVector> scalarResultMatrix(scalarResultName.size());
  for(int i = 0;i < scalarResultMatrix.size();i++) {
    for(int j = 0;j < ptcls.size();j++) {
      scalarResultMatrix[i].push_back(ptcls[j].getMaterialDirectionAngle());
    }
  }

  saveResultsToFile_res_format_flexible_resultTypes(vectorResultName,
                                                    vectorResultMatrix,
                                                    scalarResultName,
                                                    scalarResultMatrix,
                                                    stepValue,InputData,
                                                    logFile);

}

/*****************************************************************************/
/*****************************************************************************/
// Outputting calculated matrix to fem.res file
void AnisotropyGenerator::saveResultsToFile_res_format_flexible_resultTypes(
    std::vector<std::string>& vectorResultName,
    std::vector<dbMatrix>& vectorResultMatrix,
    std::vector<std::string>& scalarResultName,
    std::vector<dbVector>& scalarResultVector,dbVector& stepValue,
    InputFileData* InputData,ofstream& logFile) {

  using namespace std;

  string my_location = "OnNodes";
  string analysis_name = " ";

  string resultType = "Vector";
  for(int i = 0;i < vectorResultName.size();i++) {

    // Setting up the Result line
    femResFile << "Result " << "\"" << vectorResultName[i] << "\" \""
        << analysis_name << "\" " << stepValue[i] << " " << resultType << " "
        << my_location << endl;

    // Setting up the result type header line
    int numDofs = vectorResultMatrix[i][0].size();
    string my_type_header = "ComponentNames \"DOF_1\",";
    if(numDofs > 1) {
      for(int i = 2;i < numDofs + 1;i++) {
        stringstream i_ss;
        i_ss << i;
        my_type_header += " \"DOF_" + i_ss.str() + "\",";
      }
    }
    femResFile << my_type_header << endl;

    // Outputting results
    femResFile << "Values" << endl;
    int nodeNum = 1;
    for(int r = 0;r < vectorResultMatrix[i].size();r++) {

      femResFile << nodeNum;

      for(int s = 0;s < numDofs;s++) {
        femResFile << " " << vectorResultMatrix[i][r][s];
      }

      femResFile << endl;
      nodeNum++;

    }
    femResFile << "End Values" << endl;
  }

  // --------------------------------------------------------------------
  // Output scalar values
  resultType = "Scalar";
  for(int i = 0;i < scalarResultName.size();i++) {

    // Setting up the Result line
    femResFile << "Result " << "\"" << scalarResultName[i] << "\" \""
        << analysis_name << "\" " << stepValue[i] << " " << resultType << " "
        << my_location << endl;

    // Setting up the result type header line
    string my_type_header = "ComponentNames \"DOF_1\",";
    femResFile << my_type_header << endl;

    // Outputting results
    femResFile << "Values" << endl;
    int nodeNum = 1;
    for(int r = 0;r < scalarResultVector[i].size();r++) {

      femResFile << nodeNum << " " << scalarResultVector[i][r] << endl;

      nodeNum++;

    }
    femResFile << "End Values" << endl;
  }

}

/*****************************************************************************/
/*****************************************************************************/
//	enforce epicardium properties on boundary nodes (boundary between
//	epicardium and endocardium mesh)
void AnisotropyGenerator::eliminateSpecificSurface(InputFileData* InputData,
                                                   AnisotropyCondition& condition,
                                                   intVector& surfElems,
                                                   int surfID,
                                                   std::ofstream& logFile) {

  using namespace std;

  vector<ConditionParticle>& ptcls = condition.getParticles();
  vector<ConditionElement>& elems = condition.getElements();

  bool notSame = false;

  for(int j = 1;j < surfElems.size();j++) {

    if(elems[surfElems[j - 1]].getSurfaceID()
      != elems[surfElems[j]].getSurfaceID()) {
      notSame = true;
    }
  }
  if(notSame == true) {

    intVector dummyVec;
    for(int k = 0;k < surfElems.size();k++) {
      if(elems[surfElems[k]].getSurfaceID() == surfID) {
        dummyVec.push_back(surfElems[k]);
      }
    }
    surfElems = dummyVec;
  }

}
