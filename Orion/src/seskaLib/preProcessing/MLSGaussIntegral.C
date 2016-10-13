#include "MLSGaussIntegral.h"

MLSGaussIntegral::MLSGaussIntegral(InputFileData* InputData,
				   std::ofstream& logFile) 
  : MLSShapeFuncSet(InputData,logFile),
    ParticleDistribution(InputData,logFile),
    BackgroundMesh(InputData,logFile){}


/**********************************************************************/
/**********************************************************************/
// Calculate for all Gauss points and all particles a vector containing 
// the calculated shape functions and their derivatives for all 
// supporting particles.
void MLSGaussIntegral::setAllShapeFuncs(InputFileData* InputData,
					std::map<std::string,double>& calcData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile,
					PetscViewer& viewerMPI,
					PetscViewer& viewerSEQ) {

  using namespace std;

  // Compute derivatives of the influence radius field.
//   if((int)InputData->getValue("shapefunctionType") == 6)
    
//     setInflRadiusDerivsOnGauss(InputData,calcData,modelData,
// 			       logFile,viewerMPI,viewerSEQ);

  // Calculate for all particles a vector containing the calculated
  // shape functions and their derivations for all supporting particles.
  setShapeFuncsOnPtcls(InputData,calcData,modelData,logFile,
		       viewerMPI,viewerSEQ);

  // Calculate for all inner and boundary gauss points a vector 
  // containing the calculated shape functions and their derivations
  // for all supporting particles.
  setShapeFuncsOnGauss(InputData,calcData,modelData,logFile,viewerMPI,
		       viewerSEQ);
  setShapeFuncsOnBGauss(InputData,calcData,modelData,logFile,viewerMPI,
			viewerSEQ);

}

// /**********************************************************************/
// /**********************************************************************/
// // Compute derivatives of the influence radius field.
// void MLSGaussIntegral::setInflRadiusDerivsOnGauss(InputFileData* InputData,
// 						  std::map<std::string,double>& calcData,
// 						  std::map<std::string,double>& modelData,
// 						  std::ofstream& logFile,
// 						  PetscViewer& viewerMPI,
// 						  PetscViewer& viewerSEQ) {

//   using namespace std;

//   int usedDims = (int)modelData["usedDimensions"];

// #ifdef _geometryDebugMode_
//   logFile<<"##########################################################"<<endl;
//   logFile<<"**** computation of influence radius derivatives *********"<<endl;
// #endif

//   int ptcle;

//   for(int i=0;i<localGaussPtsNum;i++) {
//     dbMatrix& radiusDerivs = gaussPoints[i].getInflRadiusDerivs();
//     dbMatrix& feShapeDerivs = gaussPoints[i].getFEShapeFuncDerivs(); 
//     intVector& elemInfo = gaussPoints[i].getElementInfo();
//     int& elemID = elemInfo[0];

//     intVector& nodes = nodesElements[elemID].getNodes();

//     allocateArray(radiusDerivs,usedDims,usedDims);

//     // Loop over all element's particles to approximate dr/dx_k.
//     for(int j=0;j<nodes.size();j++) {
//       ptcle = nodes[j]-1;

//       // loop over radii
//       for(int k=0;k<usedDims;k++)
	
// 	// loop over derivatives
// 	for(int l=0;l<usedDims;l++)

// 	  radiusDerivs[k][l] += 
// 	    feShapeDerivs[k][j]*particles[ptcle].getRadius(k);


//     }

//   }


// #ifdef _geometryDebugMode_
//   logFile<<"******************************************************"<<endl;
//   logFile<<"********** influence radius derivatives **************"<<endl;
//   for(int i=0;i<localGaussPtsNum;i++) {
//     logFile<<"GAUSSPOINT "<<i<<" ("<<gaussPoints[i].getGlobalID()<<")"
// 	   <<endl;
//     dbMatrix& radiusDerivs = gaussPoints[i].getInflRadiusDerivs();
//     for(int k=0;k<usedDims;k++)
//       for(int l=0;l<usedDims;l++)
// 	logFile<<"dN["<<k<<"]["<<l<<"]="<<radiusDerivs[k][l]<<endl;
//   }
//   logFile<<"******************************************************"<<endl;
//   int globalElemNum = elementRootList.size();
//   int localElemNum = elemGaussIdx.size();
//   std::map<std::string,double>& backGroundMeshInfo = 
//     InputData->getBackGroundMeshInfo();
//   int gaussPtsPerVolElem = 
//     backGroundMeshInfo["gaussPointsPerVolumeElement"];
//   int idx,node,sIdx1;
//   dbVector localRadii;
//   dbMatrix3 globalRadii;
//   for(int i=0;i<nodesElements.size();i++) {
//     FEMElement& elem = nodesElements[i];
//     intVector& intPts = nodesElements[i].getVolumeIntegrationPts();
//     intVector& nodes = elem.getNodes();
//     int nodesPerElem = nodes.size();
//     dbMatrix& sFuncs = elem.getVolumeShapeFuncOrds();
//     for(int j=0;j<intPts.size();j++) {
//       idx = intPts[j];
//       dbVector& radius = gaussPoints[idx].getInflRadiusField();      
//       radius.resize(usedDims);
//       for(int k=0;k<nodes.size();k++) {
// 	node = nodes[k]-1;
// 	radius[0] += 
// 	  sFuncs[j][k]*particles[node].getRadius(0);	
// 	radius[1] += 
// 	  sFuncs[j][k]*particles[node].getRadius(1);
// 	radius[2] += 
// 	  sFuncs[j][k]*particles[node].getRadius(2);
//       }
//       if(localRadii.size() == 0) 
// 	allocateArray(localRadii,
// 		      localElemNum*gaussPtsPerVolElem*radius.size());
//       for(int k=0;k<radius.size();k++)
// 	localRadii[sIdx1+k] = radius[k];
//       sIdx1 += radius.size();
//     } 
//   }
//   globalRadii.resize(1);
//   mergeLocalElemGaussVectors(InputData,localRadii,usedDims,
// 			     globalRadii,calcData,modelData,
// 			     logFile);
//   system("cp fem.res radius.res");
//   ofstream radiusMeshFile;
//   radiusMeshFile.open("radius.msh");
//   radiusMeshFile<<"MESH  dimension 3  ElemType Point Nnode 1"<<endl;
//   radiusMeshFile<<"Coordinates"<<endl;
//   for(int i=0;i<surfaceGaussPtsNum;i++) {
//     int& currentPoint = surfaceGaussPtsIdx[i][0];
//     dbVector& coords = boundGaussPoints[currentPoint].getCoords();
//     radiusMeshFile<<i+1<<" ";
//     for(int j=0;j<coords.size();j++)
//       radiusMeshFile<<coords[j]<<" ";
//     radiusMeshFile<<endl;
//   }
//   radiusMeshFile<<"end coordinates"<<endl;
//   radiusMeshFile<<endl;
//   radiusMeshFile<<"Elements"<<endl;
//   for(int i=0;i<surfaceGaussPtsNum;i++)
//     radiusMeshFile<<i+1<<" "<<i+1<<"  1"<<endl;
//   radiusMeshFile<<"end elements"<<endl;
//   std::ofstream radiusResFile;
//   radiusResFile.open("radius.res");
//   radiusResFile.precision(12);
//   radiusResFile.setf(ios_base::scientific,ios_base::floatfield);
//   radiusResFile<<"GiD Post Results File 1.0"<<endl;

//   GaussPointSet* GaussSet;
//   if((int)backGroundMeshInfo["elemType"] == 1) {      
//     switch((int)backGroundMeshInfo["gaussPointsPerVolumeElement"]) { 
//     case 1:
//       GaussSet = new GaussSetTetra1();
//       break;
//     case 4:
//       GaussSet = new GaussSetTetra4();
//       break;
//     case 5:
//       GaussSet = new GaussSetTetra5();
//       break;
//     default:
//       logFile<<"In function MLSGaussIntegral::setShapeFuncsOnGauss element\n "
// 	     <<(int)backGroundMeshInfo["gaussPointsPerVolumeElement"]
// 	     <<" Gauss points per tetrahedral element\n " 
// 	     <<" are not supported!"<<endl;
//       MPI_Abort(MPI_COMM_WORLD,1);
//     }
//     radiusResFile<<"GaussPoints \"Volume Gauss points\" ElemType Tetrahedra"<<endl;
    
//   }
//   else if((int)backGroundMeshInfo["elemType"] == 2) {
//     switch((int)backGroundMeshInfo["gaussPointsPerVolumeElement"]) {
//     case 1:
//       GaussSet = new GaussSetCube1();
//       break;
//     case 8:
//       GaussSet = new GaussSetCube8();
//       break;
//     case 27:
//       GaussSet = new GaussSetCube27();
//       break;
//     case 64:
//       GaussSet = new GaussSetCube64();
//       break;
//     case 125:
//       GaussSet = new GaussSetCube125();
//       break;
//     default:
//       logFile<<"In function MLSGaussIntegral::setShapeFuncsOnGauss element\n "
// 	     <<(int)backGroundMeshInfo["gaussPointsPerVolumeElement"]
// 	     <<" Gauss points per hexahedral element\n " 
// 	     <<" are not supported!"<<endl;
//       MPI_Abort(MPI_COMM_WORLD,1);
//     }
//     radiusResFile<<"GaussPoints \"Volume Gauss points\" ElemType Hexahedra"<<endl;
//   }
  
//   else {
//     logFile<<"In function MLSGaussIntegral::setShapeFuncsOnGauss element "
// 	   <<"type "<<(int)backGroundMeshInfo["elemType"]
// 	   <<"\n is not supported!"<<endl;
//     MPI_Abort(MPI_COMM_WORLD,1);
//   }
//   radiusResFile<<"Number Of Gauss Points: "
// 	       <<(int)backGroundMeshInfo["gaussPointsPerVolumeElement"]<<endl;
//   radiusResFile<<"Natural Coordinates: given"<<endl;
//   for(int i=0;i<(GaussSet->coord).size();i++) {
//     for(int j=0;j<(GaussSet->coord[i]).size();j++)
//       radiusResFile<<GaussSet->coord[i][j]<<" ";
//     radiusResFile<<endl;
//   }
//   radiusResFile<<"end gausspoints"<<endl;
//   delete GaussSet;

//   radiusResFile<<"Result \"radius field \" \" \" "
// 	       <<calcData["currentTime"]
// 	       <<"  Vector OnGaussPoints \"Volume Gauss points\""<<endl;
//   radiusResFile<<"ComponentNames ";
//   for(int i=0;i<3;i++)
//     radiusResFile<<"\"r_"<<i+1<<"\", ";
//   radiusResFile<<endl;
//   radiusResFile<<"Values"<<endl;
//   for(int i=0;i<globalRadii.size();i++) {
//     radiusResFile<<i+1<<" ";
//     for(int gPt=0;gPt<globalRadii[i].size();gPt++) {
//       for(int j=0;j<globalRadii[i][gPt].size();j++) {
// 	radiusResFile<<globalRadii[i][gPt][j]<<" ";
//       }
//       radiusResFile<<endl;
//       radiusResFile.flush();
//     }
//   }
//   radiusResFile<<"End Values"<<endl;
//   globalRadii.resize(0);
//   radiusResFile.flush();
// #endif

//   vector<FEMElement>(0,FEMElement(0)).swap(nodesElements);

// }

/************************************************************************/
/************************************************************************/
// Calculate for all gauss points a vector containing the calculated
// shape functions and their derivations for all supporting particles.
void MLSGaussIntegral::setShapeFuncsOnGauss(InputFileData* InputData,
					    std::map<std::string,double>& calcData,
					    std::map<std::string,double>& modelData,
					    std::ofstream& logFile,
					    PetscViewer& viewerMPI,
					    PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  unsigned int derivationOrder = 
    (unsigned int)modelData["shapesDerivationOrderOnIntPoints"];

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  dbMatrix secondDerivShapes(6);

  int supportSize;
  
  // Loop over all local gauss points to set their shape function 
  // vectors.
  int rank; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int basisTermNum;
  int maxSupportSize=0;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"******************* gausspoint shapes ********************"<<endl;
#endif


  for(int i=0;i<localGaussPtsNum;i++) {
    supportSize = gaussPoints[i].getSupportCounts();
    
    if(supportSize > maxSupportSize) {
      shapeFuncs.resize(supportSize);

      if(derivationOrder > 0)

	for(int j=0;j<3;j++)
	  firstDerivShapes[j].resize(supportSize);

      if(derivationOrder > 1)

	for(int j=0;j<6;j++)
	  secondDerivShapes[j].resize(supportSize);

    }
    
#ifdef _geometryDebugMode_
    logFile<<"Gauss point "<<gaussPoints[i].getGlobalID()<<endl;
#endif

    double& x = gaussPoints[i].getCoord(0);
    double& y = gaussPoints[i].getCoord(1);
    double& z = gaussPoints[i].getCoord(2);

    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,supportSize,gaussPoints[i].getSupportPtcls(),
		     particles,x,y,z,shapeFuncs,basisTermNum,modelData,
		     logFile,viewerSEQ);
    
      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);

    }
    else if(derivationOrder == 1) {
      
      calcShapeFuncs(InputData,supportSize,gaussPoints[i].getSupportPtcls(),
		     particles,x,y,z,shapeFuncs,
		     firstDerivShapes,modelData,
		     logFile,viewerSEQ);
    
      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      gaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else if(derivationOrder == 2) {

      calcShapeFuncs(InputData,supportSize,gaussPoints[i].getSupportPtcls(),
		     particles,x,y,z,shapeFuncs,firstDerivShapes,
		     secondDerivShapes,modelData,logFile,viewerSEQ);
    
      // Store the calculated shape function set.
      gaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      gaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);
      gaussPoints[i].setSecondDerivShapes(supportSize,secondDerivShapes);

    }
    else {
      logFile <<"In MLSGaussIntegral::setShapeFuncsOnGauss shape function "
	      <<"derivations higher than second order are not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

#ifdef _geometryDebugMode_
  double PUM;
  logFile<<"******************************************************"<<endl;
  logFile<<"******** calculated shape sets on gauss pts **********"<<endl;
  for(int i=0;i<localGaussPtsNum;i++) {
    PUM = 0;
    logFile<<"GAUSSPOINT "<<gaussPoints[i].getGlobalID()<<":: ";
    for(int j=0;j<gaussPoints[i].getSupportCounts();j++) {
      PUM += gaussPoints[i].getShapeFunc(j);
      logFile<<gaussPoints[i].getSupportPtcle(j)<<": "
	     <<gaussPoints[i].getShapeFunc(j)<<" | ";
    }
    if(fabs(1.0-PUM) > 0.000000001)
      logFile<<" PUM= "<<PUM<<endl;
    else
      logFile<<endl;
    if(derivationOrder > 0) {
      dbMatrix& firstDerivs = gaussPoints[i].getFirstDerivShapes();
      logFile<<"first order derivations"<<endl;
      for(int k=0;k<3;k++) {
	for(int j=0;j<firstDerivs[k].size();j++)
	  logFile<<gaussPoints[i].getSupportPtcle(j)<<": "
		 <<firstDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
    if(derivationOrder > 1) {
      logFile<<"second order derivations"<<endl;
      dbMatrix& secondDerivs = gaussPoints[i].getSecondDerivShapes();
      for(int k=0;k<6;k++) {
	for(int j=0;j<secondDerivs[k].size();j++)
	  logFile<<gaussPoints[i].getSupportPtcle(j)<<": "
		 <<secondDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
  }
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
}

/**********************************************************************/
/**********************************************************************/
// Calculate for all particles a vector containing the calculated
// shape functions and their derivations for all supported particles.
void MLSGaussIntegral::setShapeFuncsOnPtcls(InputFileData* InputData,
					    std::map<std::string,double>& calcData,
					    std::map<std::string,double>& modelData,
					    std::ofstream& logFile,
					    PetscViewer& viewerMPI,
					    PetscViewer& viewerSEQ) {

  using namespace std;


  int rank,size; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int derivationOrder = 
    (int)modelData["shapesDerivationOrderOnPtcls"];
  derivationOrder = 1;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"********************* particle shapes ********************"<<endl;
  logFile<<"**********************************************************"<<endl;
  logFile<<"derivation order="<<derivationOrder<<endl;
#endif

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  dbMatrix secondDerivShapes(6);
  int basisTermNum,supportSize;

  // Loop over processor's locally needed particles portion.
  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    supportSize = suppPtcls.size();
    
    if(shapeFuncs.size() < supportSize)
      shapeFuncs.resize(supportSize);
    
    if(derivationOrder > 0)

      for(int j=0;j<firstDerivShapes.size();j++)
	firstDerivShapes[j].resize(supportSize);
    
    if(derivationOrder > 1)
      
      for(int j=0;j<secondDerivShapes.size();j++)
	secondDerivShapes[j].resize(supportSize);
    
#ifdef _geometryDebugMode_
    logFile<<i<<".) PARTICLE "<<ptcle<<endl;
#endif

    double& x = particles[ptcle].getCoord(0);
    double& y = particles[ptcle].getCoord(1);
    double& z = particles[ptcle].getCoord(2);

    // compute the shape functions
    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,supportSize,suppPtcls,particles,x,y,z,
		     shapeFuncs,basisTermNum,modelData,logFile,viewerSEQ);

      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);

    }

    // compute the shape functions and their first derivatives
    else if(derivationOrder == 1) {


      calcShapeFuncs(InputData,supportSize,suppPtcls,particles,x,y,z,
		     shapeFuncs,firstDerivShapes,modelData,logFile,
		     viewerSEQ);
      
      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);
      particles[ptcle].setFirstDerivShapes(supportSize,firstDerivShapes);

    }

    // compute the shape functions, their first and second derivatives
    else if(derivationOrder == 2) {

      calcShapeFuncs(InputData,supportSize,suppPtcls,particles,x,y,z,
		     shapeFuncs,firstDerivShapes,secondDerivShapes,
		     modelData,logFile,viewerSEQ);
    
      // Store the calculated shape function set.
      particles[ptcle].setShapeFuncs(supportSize,shapeFuncs);
      particles[ptcle].setFirstDerivShapes(supportSize,firstDerivShapes);
      particles[ptcle].setSecondDerivShapes(supportSize,secondDerivShapes);

    }
    else {
      logFile <<"In MLSGaussIntegral::setShapeFuncsOnPtcls shape function\n "
	      <<"derivations higher than second order are not supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }


  }

#if defined _geometryDebugMode_ || defined _sphericalShapeFuncDebugMode_
  double PUM;
  int zeroCounter;
  logFile<<"******************************************************"<<endl;
  logFile<<"*********** calculated shape sets on ptcls ***********"<<endl;
  for(int i=0;i<localGlobalPtcls.size();i++) {
    int& ptcle = localGlobalPtcls[i];
    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    dbVector& suppShapes = particles[ptcle].getShapeFuncs();
    dbMatrix& firstDerivs = particles[i].getFirstDerivShapes();
    dbMatrix& secondDerivs = particles[i].getSecondDerivShapes();
    PUM = 0;
    logFile<<"Ptcle "<<ptcle<<"("<<oldIdx[ptcle]<<"):: ";
    for(int j=0;j<suppPtcls.size();j++) {
      PUM += suppShapes[j];
      logFile<<suppPtcls[j]<<": "
	     <<suppShapes[j]<<" | ";
    }
    if(fabs(1.0-PUM) > 1.0e-13)
      logFile<<" PUM= "<<PUM<<endl;
    if(derivationOrder > 0) {
      zeroCounter = 0;
      logFile<<"first order derivations"<<endl;
      for(int k=0;k<firstDerivs.size();k++) {
	for(int j=0;j<firstDerivs[k].size();j++)
	  logFile<<suppPtcls[j]<<": "<<firstDerivs[k][j]<<" | ";
	logFile<<endl;
      }
      logFile<<"zeroCounter = "<<zeroCounter<<" of "<<suppPtcls.size()*3<<endl;
    }
    if(derivationOrder > 1) {
      logFile<<"second order derivations"<<endl;
      for(int k=0;k<secondDerivs.size();k++) {
	for(int j=0;j<secondDerivs[k].size();j++)
	  logFile<<suppPtcls[j]<<": "<<secondDerivs[k][j]<<" | ";
	logFile<<endl;
      }
    }
  }
#endif
// #ifdef _geometryDebugMode_
//   if(size == 1) {
//     for(int i=0;i<particles.size();i++) {
//       intVector& sPtclsI = particles[i].getSupportPtcls();
//       dbVector& shapesI = particles[i].getShapeFuncs();
//       dbMatrix& firstI = particles[i].getFirstDerivShapes();
//       dbMatrix& secondI = particles[i].getSecondDerivShapes();
//       for(int j=0;j<sPtclsI.size();j++) {
// 	int& sPtcleJ = sPtclsI[j];
// 	intVector& sPtclsJ = particles[sPtcleJ].getSupportPtcls();
// 	dbVector& shapesJ = particles[sPtcleJ].getShapeFuncs();
// 	dbMatrix& firstJ = particles[sPtcleJ].getFirstDerivShapes();
// 	dbMatrix& secondJ = particles[sPtcleJ].getSecondDerivShapes();
// 	int posI = findIntVecPos(i,0,sPtclsJ.size(),sPtclsJ);
// 	if(posI > -1)
// 	  if(shapesI[sPtcleJ] != shapesJ[posI]) {
// 	    logFile<<"non-matching shapes of ptcles "
// 		   <<i<<"("<<shapesI[sPtcleJ]<<") - "
// 		   <<sPtcleJ<<"("<<shapesJ[posI]<<")"<<endl;
// 	    //MPI_Abort(MPI_COMM_WORLD,1);
// 	  }
// 	for(int k=0;k<firstI.size();k++)
// 	  if(firstI[k][sPtcleJ] != firstJ[k][posI]) {
// 	    logFile<<"non-matching shapes deriv of ptcles "
// 		   <<i<<"("<<shapesI[sPtcleJ]<<") - "
// 		   <<sPtcleJ<<"("<<shapesJ[posI]<<")"<<endl;
// 	    //MPI_Abort(MPI_COMM_WORLD,1);
// 	  }
//       }
//     }
//   }
// #endif

}

/**********************************************************************/
/**********************************************************************/
// Calculate for all boundary gauss points a vector containing the 
// calculated shape functions for all supported particles.
void MLSGaussIntegral::setShapeFuncsOnBGauss(InputFileData* InputData,
					     std::map<std::string,double>& calcData,
					     std::map<std::string,double>& modelData,
					     std::ofstream& logFile,
					     PetscViewer& viewerMPI,
					     PetscViewer& viewerSEQ) {

  using namespace std;

#ifdef _geometryDebugMode_
  logFile<<"##########################################################"<<endl;
  logFile<<"*************** boundary gausspoint shapes ***************"<<endl;
#endif

  unsigned int derivationOrder = (unsigned int)
    modelData["shapesDerivationOrderOnBoundIntPoints"];

  int basisTermNum,supportSize;

  dbVector shapeFuncs;
  dbMatrix firstDerivShapes(3);
  dbMatrix secondDerivShapes(6);

  int m=0;
  int localMaxSupportSize = 0;

  // Loop over all local boundary gauss points to set their shape 
  // function vectors.
  for(int i=0;i<localBGaussPtsNum;i++) {
    supportSize = boundGaussPoints[i].getSupportCounts();

    if(supportSize > localMaxSupportSize) {
      localMaxSupportSize = supportSize;
      shapeFuncs.resize(supportSize);

      if(derivationOrder > 0)

	for(int j=0;j<3;j++)
	  firstDerivShapes[j].resize(supportSize);

      if(derivationOrder > 1)

	for(int j=0;j<6;j++)
	  secondDerivShapes[j].resize(supportSize);

    }

#ifdef _geometryDebugMode_
    logFile<<"Boundary gauss point "<<i<<endl;
#endif

    double& x = boundGaussPoints[i].getCoord(0);
    double& y = boundGaussPoints[i].getCoord(1);
    double& z = boundGaussPoints[i].getCoord(2);

    if(derivationOrder == 0) {

      calcShapeFuncs(InputData,supportSize,boundGaussPoints[i].getSupportPtcls(),
		     particles,x,y,z,shapeFuncs,basisTermNum,modelData,logFile,
		     viewerSEQ);

      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
    }
    else if(derivationOrder == 1) {

      calcShapeFuncs(InputData,supportSize,
		     boundGaussPoints[i].getSupportPtcls(),particles,
		     x,y,z,shapeFuncs,firstDerivShapes,modelData,logFile,
		     viewerSEQ);
    
      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      boundGaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);

    }
    else if(derivationOrder == 2) {

      calcShapeFuncs(InputData,supportSize,
		     boundGaussPoints[i].getSupportPtcls(),particles,
		     x,y,z,shapeFuncs,firstDerivShapes,secondDerivShapes,
		     modelData,logFile,viewerSEQ);
    
      // Store the calculated shape function set.
      boundGaussPoints[i].setShapeFuncs(supportSize,shapeFuncs);
      boundGaussPoints[i].setFirstDerivShapes(supportSize,firstDerivShapes);
      boundGaussPoints[i].setSecondDerivShapes(supportSize,secondDerivShapes);

    }
    else {
      logFile <<"In MLSGaussIntegral::setShapeFuncsOnBGauss shape function "
	      <<"derivations higher than second order are not supported!"<<endl;

      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

#ifdef _geometryDebugMode_
  double PUM;
  logFile<<"######################################################"<<endl;
  logFile<<"** calculated shapes on boundary integration points **"<<endl;
  for(int i=0;i<localBGaussPtsNum;i++) {
    PUM = 0;
    logFile<<"Bound GAUSSPOINT "<<i<<": support = "
	   <<boundGaussPoints[i].getSupportCounts()<<" ";
    for(int j=0;j<boundGaussPoints[i].getSupportCounts();j++) {
      PUM += boundGaussPoints[i].getShapeFunc(j);
      logFile<<boundGaussPoints[i].getSupportPtcle(j)<<": "<<boundGaussPoints[i].getShapeFunc(j)<<" | ";
    }
    if(fabs(1.0-PUM) > 0.000000001)
      logFile<<"PUM= "<<PUM<<endl;
    else
      logFile<<endl;
    if(derivationOrder > 0) {
      dbMatrix& firstDerivs = boundGaussPoints[i].getFirstDerivShapes();
      logFile<<"first order derivations"<<endl;
      for(int k=0;k<3;k++) {
	for(int j=0;j<firstDerivs[k].size();j++)
	  logFile<<boundGaussPoints[i].getSupportPtcle(j)<<": "
		 <<firstDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
    if(derivationOrder > 1) {
      dbMatrix& secondDerivs = boundGaussPoints[i].getSecondDerivShapes();
      logFile<<"second order derivations"<<endl;
      for(int k=0;k<6;k++) {
	for(int j=0;j<secondDerivs[k].size();j++)
	  logFile<<boundGaussPoints[i].getSupportPtcle(j)<<": "
		 <<secondDerivs[k][j]<<" |";
	logFile<<endl;
      }
    }
  }
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
}

/************************************************************************/
/************************************************************************/
// Check the shape function and their derivatives on all Gauss points.
void MLSGaussIntegral::checkGaussShapes(InputFileData* InputData,
					std::map<std::string,double>& modelData,
					std::ofstream& logFile) {
  
  using namespace std;

  int volDerivOrder = (int)modelData["shapesDerivationOrderOnIntPoints"];
  int boundDerivOrder = 
    (int)modelData["shapesDerivationOrderOnBoundIntPoints"];

  double shape;
  dbVector dShape;
  dbVector d2Shape;
  double limit = DBL_EPSILON*3.0e03;

  logFile<<"######################################################"<<endl;
  logFile<<"************** shape function check ******************"<<endl;
  logFile<<"limit = "<<limit<<endl;

  // Loop over all local volume gauss points.
  for(int i=0;i<localGaussPtsNum;i++) {

    logFile<<"********* vol GAUSS "<<i<<":"<<endl;

    shape = 0;
    clearArray(dShape); 
    clearArray(d2Shape);

    // loop over all shape functions.
    for(int j=0;j<gaussPoints[i].getSupportCounts();j++)
      shape += gaussPoints[i].getShapeFunc(j);
    
    if(fabs(1.0-shape) > limit)
      logFile<<"shape = "<<shape<<endl;
   
    // first derivatives
    if(volDerivOrder > 0) {
      dbMatrix& firstDerivs = gaussPoints[i].getFirstDerivShapes();

      if(dShape.size() < firstDerivs.size())
	dShape.resize(firstDerivs.size());

      // Loop over all kind of derivatives (dx,dy,dz)
      for(int k=0;k<firstDerivs.size();k++)

	for(int j=0;j<firstDerivs[k].size();j++)
	  dShape[k] += firstDerivs[k][j];
     
      for(int k=0;k<firstDerivs.size();k++)

	if(fabs(dShape[k]) > limit)
	  logFile<<"d"<<k<<"shape = "<<dShape[k]<<endl;
 
    }
    // second order derivatives
    if(volDerivOrder > 1) {
      dbMatrix& secondDerivs = gaussPoints[i].getSecondDerivShapes();

      if(d2Shape.size() < secondDerivs.size())
	d2Shape.resize(secondDerivs.size());

      // Loop over all kind of derivatives
      for(int k=0;k<secondDerivs.size();k++)

	for(int j=0;j<secondDerivs[k].size();j++)
	  d2Shape[k] += secondDerivs[k][j];
       
      for(int k=0;k<secondDerivs.size();k++)

	if(fabs(d2Shape[k]) > limit)
	  logFile<<"dd"<<k<<"shape = "<<d2Shape[k]<<endl;
    }

  }

  /*********************************************************************/
  // boundary Gauss points.
  for(int i=0;i<localBGaussPtsNum;i++) {

    logFile<<"********* bound GAUSS "<<i<<":"<<endl;

    shape = 0;
    clearArray(dShape); 
    clearArray(d2Shape);

    // loop over all shape functions.
    for(int j=0;j<boundGaussPoints[i].getSupportCounts();j++)
      shape += boundGaussPoints[i].getShapeFunc(j);
    
    if(fabs(1.0-shape) > limit)
      logFile<<"shape = "<<shape<<endl;
   
    // first derivatives
    if(volDerivOrder > 0) {
      dbMatrix& firstDerivs = boundGaussPoints[i].getFirstDerivShapes();

      if(dShape.size() < firstDerivs.size())
	dShape.resize(firstDerivs.size());

      // Loop over all kind of derivatives (dx,dy,dz)
      for(int k=0;k<firstDerivs.size();k++)

	for(int j=0;j<firstDerivs[k].size();j++)
	  dShape[k] += firstDerivs[k][j];
     
      for(int k=0;k<firstDerivs.size();k++)

	if(fabs(dShape[k]) > limit)
	  logFile<<"d"<<k<<"shape = "<<dShape[k]<<endl;
 
    }
    // second order derivatives
    if(volDerivOrder > 1) {
      dbMatrix& secondDerivs = boundGaussPoints[i].getSecondDerivShapes();

      if(d2Shape.size() < secondDerivs.size())
	d2Shape.resize(secondDerivs.size());

      // Loop over all kind of derivatives
      for(int k=0;k<secondDerivs.size();k++)

	for(int j=0;j<secondDerivs[k].size();j++)
	  d2Shape[k] += secondDerivs[k][j];
       
      for(int k=0;k<secondDerivs.size();k++)

	if(fabs(d2Shape[k]) > limit)
	  logFile<<"dd"<<k<<"shape = "<<d2Shape[k]<<endl;
    }

  }

}

