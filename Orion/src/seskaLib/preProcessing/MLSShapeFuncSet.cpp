#include "MLSShapeFuncSet.h"



/************************************************************************/
/************************************************************************/
// Calculate at a certain point for all its supporting particles their 
// shape functions.
void MLSShapeFuncSet::calcShapeFuncs(InputFileData* InputData,
				     int& supportSize,
				     intVector& sPtcls,
				     std::vector<Particle>& particles,
				     double& x,double& y,double& z,
				     dbVector& shapeFuncs,
				     int& basisTermNum,
				     std::map<std::string,double>& modelData,
				     std::ofstream& logFile,
				     PetscViewer& viewerSEQ) {

  using namespace std;

  int shapeFuncType = (int)InputData->getValue("shapefunctionType");

  switch(shapeFuncType) {

    // Calculate a Shepard shapefunction set.
  case 1:
    ShepardShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
				 x,y,z,shapeFuncs,basisTermNum,modelData,
				 logFile,viewerSEQ);

    break;

  // Calculate a EFG shapefunction set.
  case 2:

    EFGShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			     x,y,z,shapeFuncs,basisTermNum,modelData,
			     logFile,viewerSEQ);
    
    break;

  // Calculate RKPM  shapefunction set.
  case 3:

    RKPMShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			      x,y,z,shapeFuncs,basisTermNum,modelData,
			      logFile,viewerSEQ);
    
    break;

    // Calculate a orthogonalized MLS shapefunction set.
  case 4:
    OrthoShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			       x,y,z,shapeFuncs,basisTermNum,modelData,
			       logFile,viewerSEQ);

    break;

    // Calculate a MLS shapefunction set, where the influence zones
    // are different in negative and positive coordinate direction.
  case 5:
    AsymShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			      x,y,z,shapeFuncs,basisTermNum,modelData,
			      logFile,viewerSEQ);

    break;

  default:
    logFile<<"In MLSShapeFuncSet::calcShapeFuncs chosen shape function\n" 
	   <<"type isn't supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/************************************************************************/
/************************************************************************/
// Calculate at a certain point for all its supporting particles their 
// shape functions and their first order derivations.
void MLSShapeFuncSet::calcShapeFuncs(InputFileData* InputData,
				     int& supportSize,
				     intVector& sPtcls,
				     std::vector<Particle>& particles,
				     double& x,double& y,double& z,
				     dbVector& shapeFuncs,
				     dbMatrix& firstDerivShapes,
				     std::map<std::string,double>& modelData,
				     std::ofstream& logFile,
				     PetscViewer& viewerSEQ) {

  using namespace std;

  int shapeFuncType = (int)InputData->getValue("shapefunctionType");

  switch(shapeFuncType) {

    // Calculate a Shepard shapefunction set.
  case 1:
    ShepardShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
				 x,y,z,shapeFuncs,firstDerivShapes,
				 modelData,logFile,viewerSEQ);

    break;

  // Calculate a EFG shapefunction set.
  case 2:

    EFGShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			     x,y,z,shapeFuncs,firstDerivShapes[0],
			     firstDerivShapes[1],firstDerivShapes[2],
			     modelData,logFile,viewerSEQ);

    break;

  // Calculate RKPM  shapefunction set.
  case 3:

    RKPMShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			      x,y,z,shapeFuncs,firstDerivShapes[0],
			      firstDerivShapes[1],firstDerivShapes[2],
			      modelData,logFile,viewerSEQ);

    break;

    // Calculate a orthogonalized MLS shapefunction set.
  case 4:
    OrthoShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			       x,y,z,shapeFuncs,firstDerivShapes[0],
			       firstDerivShapes[1],firstDerivShapes[2],
			       modelData,logFile,viewerSEQ);

    break;

    // Calculate a MLS shapefunction set, where the influence zones
    // are different in negative and positive coordinate direction.
  case 5:
    AsymShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			      x,y,z,shapeFuncs,firstDerivShapes[0],
			      firstDerivShapes[1],firstDerivShapes[2],
			      modelData,logFile,viewerSEQ);

    break;

  default:
    logFile<<"In MLSShapeFuncSet::calcShapeFuncs chosen shape function\n" 
	   <<"type isn't supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/**********************************************************************/
// Calculate at a certain point for all its supporting particles their 
// shape functions and their first and second order derivations.
void MLSShapeFuncSet::calcShapeFuncs(InputFileData* InputData,
				     int& supportSize,
				     intVector& sPtcls,
				     std::vector<Particle>& particles,
				     double& x,double& y,double& z,
				     dbVector& shapeFuncs,
				     dbMatrix& firstDerivShapes,
				     dbMatrix& secondDerivShapes,
				     std::map<std::string,double>& modelData,
				     std::ofstream& logFile,
				     PetscViewer& viewerSEQ) {

  using namespace std;

  int shapeFuncType = (int)InputData->getValue("shapefunctionType");
  dbVector radius(3);

  switch(shapeFuncType) {

    // Calculate a Shepard shapefunction set.
  case 1:
    ShepardShapeFunc::calcShapes(InputData,supportSize,sPtcls,
				 particles,x,y,z,shapeFuncs,
				 firstDerivShapes,secondDerivShapes,
				 modelData,logFile,viewerSEQ);

    break;

  // Calculate a EFG shapefunction set.
  case 2:

    EFGShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			     x,y,z,shapeFuncs,firstDerivShapes[0],
			     firstDerivShapes[1],firstDerivShapes[2],
			     secondDerivShapes[0],secondDerivShapes[1],
			     secondDerivShapes[2],secondDerivShapes[3],
			     secondDerivShapes[4],secondDerivShapes[5],
			     modelData,logFile,viewerSEQ);
    break;

  // Calculate RKPM  shapefunction set.
  case 3:

    RKPMShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			      x,y,z,shapeFuncs,firstDerivShapes[0],
			      firstDerivShapes[1],firstDerivShapes[2],
			      secondDerivShapes[0],secondDerivShapes[1],
			      secondDerivShapes[2],secondDerivShapes[3],
			      secondDerivShapes[4],secondDerivShapes[5],
			      modelData,logFile,viewerSEQ);
    break;

    // Calculate a orthogonalized MLS shapefunction set.
  case 4:
    OrthoShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			       x,y,z,shapeFuncs,firstDerivShapes[0],
			       firstDerivShapes[1],firstDerivShapes[2],
			       secondDerivShapes[0],secondDerivShapes[1],
			       secondDerivShapes[2],secondDerivShapes[3],
			       secondDerivShapes[4],secondDerivShapes[5],
			       modelData,logFile,viewerSEQ);

    break;

    // Calculate a MLS shapefunction set, where the influence zones
    // are different in negative and positive coordinate direction.
  case 5:
    AsymShapeFunc::calcShapes(InputData,supportSize,sPtcls,particles,
			      x,y,z,shapeFuncs,firstDerivShapes[0],
			      firstDerivShapes[1],firstDerivShapes[2],
			      secondDerivShapes[0],secondDerivShapes[1],
			      secondDerivShapes[2],secondDerivShapes[3],
			      secondDerivShapes[4],secondDerivShapes[5],
			      modelData,logFile,viewerSEQ);
    
    break;

  default:
    logFile<<"In MLSShapeFuncSet::calcShapeFuncs chosen shape function\n" 
	   <<"type isn't supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
