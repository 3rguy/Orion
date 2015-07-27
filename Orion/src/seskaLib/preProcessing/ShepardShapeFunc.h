// Calculate for a given point a Shepard shapefunction set.

#ifndef ShepardShapeFunc_h_
#define ShepardShapeFunc_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "defs.h"
#include "InputFileData.h"
#include "mpi.h"
#include "Particle.h"
#include "petscksp.h"
#include "ShapefunctionSet.h"
#include "PrismaticWindowFunctionSet.h"
#include "SphericalWindowFunctionSet.h"
#include "WindowFunctionSet.h"


class ShepardShapeFunc : public ShapefunctionSet {

  public:

    ShepardShapeFunc() {};

    ShepardShapeFunc(InputFileData* InputData,
		     int& supportSize,
		     intVector& sPtcls,
		     std::vector<Particle>& ptcls,
		     double& x,double& y,double& z,
		     unsigned int derivationOrder,
		     std::map<std::string,double>& modelData,  
		     std::ofstream& logFile,
		     PetscViewer& viewerSEQ);

    ~ShepardShapeFunc() {};

    // Calculate at a gauss point for all its supporting particles their 
    // shape functions.
    void calcShapes(InputFileData* InputData,
		    int& supportSize,
		    intVector& sPtcls,
		    std::vector<Particle>& ptcls,
		    double& x,double& y,double& z,
		    std::map<std::string,double>& modelData,
		    std::ofstream& logFile,
		    PetscViewer& viewerSEQ);

    void calcShapes(InputFileData* InputData,
		    int& supportSize,
		    intVector& sPtcls,
		    std::vector<Particle>& ptcls,
		    double& x,double& y,double& z,
		    dbVector& shapes,int& basisTermNum,
		    std::map<std::string,double>& modelData,
		    std::ofstream& logFile,
		    PetscViewer& viewerSEQ);
    
    // Calculate at a gauss point for all its supporting particles their 
    // shape functions and their first order derivations.
    void calcShape1stDerivs(InputFileData* InputData,
			    int& supportSize,
			    intVector& sPtcls,
			    std::vector<Particle>& ptcls,
			    double& x,double& y,double& z,
			    std::map<std::string,double>& modelData,
			    std::ofstream& logFile,
			    PetscViewer& viewerSEQ);

    void calcShapes(InputFileData* InputData,
		    int& supportSize,
		    intVector& sPtcls,
		    std::vector<Particle>& ptcls,
		    double& x,double& y,double& z,
		    dbVector& shapes,
		    dbMatrix& dShapes,
		    std::map<std::string,double>& modelData,
		    std::ofstream& logFile,
		    PetscViewer& viewerSEQ);
    
    // Calculate at a point for all its supporting particles their 
    // shape functions and their first and second order derivations.
    void calcShape2ndDerivs(InputFileData* InputData,
			    int& supportSize,
			    intVector& sPtcls,
			    std::vector<Particle>& ptcls,
			    double& x,double& y,double& z,
			    std::map<std::string,double>& modelData,
			    std::ofstream& logFile,
			    PetscViewer& viewerSEQ);

    void calcShapes(InputFileData* InputData,
		    int& supportSize,
		    intVector& sPtcls,
		    std::vector<Particle>& ptcls,
		    double& x,double& y,double& z,
		    dbVector& shapes,
		    dbMatrix& dShapes,
		    dbMatrix& d2Shapes,
		    std::map<std::string,double>& modelData,
		    std::ofstream& logFile,
		    PetscViewer& viewerSEQ);

};



#endif
